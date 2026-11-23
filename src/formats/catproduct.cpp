// CATIA V5 Product (.CATProduct): the assembly structure of the product object stream.
//
// Objects carrying `_PartNumber` are product references (the document's
// own product and its local sub-assemblies). Objects carrying `_InstanceName` are instances.
// Each instance's relation header names its reference and its parent reference:
//   [link blob] <n> <parent>   reference in another document (blob: 11 03 <document uuid> ...)
//   <reference> <parent>       reference in this document
// The instance placement is the `_Position` (3 axis columns + translation, mm) of the
// PRDBAGREP object listed in its `_BagRepsList`; an instance without one is unplaced (identity).
// Linked documents are found by the UUID of their CATProdCont row, among the files that
// DecodeOptions::list_resources offers.
#include <algorithm>
#include <cctype>
#include <cstring>
#include <map>
#include <functional>
#include <memory>
#include <optional>
#include <set>
#include <string>

#include "cadcore/error.hpp"
#include "core/parallel.hpp"
#include "formats/catia/document.hpp"
#include "formats/catia/osm.hpp"
#include "formats/catia/spec.hpp"
#include "formats/registry.hpp"

namespace cadcore::formats {
namespace {

using catia::Object;
using catia::Stream;
using catia::Uuid;

[[noreturn]] void bad(const std::string& what, std::size_t at = ByteView::npos) { fail(ErrorCode::Malformed, "CATProduct: " + what, at); }

int sniff(ByteView d, std::string_view ext) {
    const bool magic = d.size >= 8 && std::memcmp(d.data, "V5_CFV2", 8) == 0;
    if (magic && ext == "catproduct") return 2;
    return ext == "catproduct" ? 1 : 0;
}

std::string hex(const Uuid& u) {
    static const char* k = "0123456789abcdef";
    std::string s;
    for (auto b : u) { s += k[b >> 4]; s += k[b & 15]; }
    return s;
}

std::string lower_ext(std::string_view path) {
    const auto dot = path.rfind('.');
    std::string e(dot == std::string_view::npos ? "" : path.substr(dot + 1));
    for (auto& c : e) c = char(std::tolower(static_cast<unsigned char>(c)));
    return e;
}

// Shared across the documents of one decode: the UUID index and decoded documents.
struct Context {
    Context(const DecodeOptions& o, Scene& s) : options(o), scene(s) {}
    const DecodeOptions& options;
    Scene& scene;
    bool indexed = false;
    std::map<Uuid, std::string> by_uuid;
    std::map<Uuid, std::uint32_t> node_of;  // decoded document -> its node
    std::set<Uuid> active;                  // documents being decoded (cycle guard)
    // Linked CATParts: their nodes exist as soon as the structure names them; the parts are
    // decoded afterwards in parallel and moved into their nodes in this order.
    struct Part { std::uint32_t node; std::string path; Scene scene; };
    std::vector<std::unique_ptr<Part>> parts;

    void index() {
        if (indexed) return;
        indexed = true;
        if (!options.list_resources || !options.load_resource) return;
        for (const std::string& path : options.list_resources()) {
            const std::string e = lower_ext(path);
            if (e != "catpart" && e != "catproduct") continue;
            const auto bytes = options.load_resource(path);
            if (!bytes) continue;
            const auto u = catia::document_uuid(ByteView(bytes->data(), bytes->size()));
            if (u && !by_uuid.count(*u)) by_uuid.emplace(*u, path);
        }
    }
};

std::uint32_t add_node(Scene& s, std::string name, const Matrix4& t = kIdentity) {
    Node n;
    n.name = std::move(name);
    n.transform = t;
    s.nodes.push_back(std::move(n));
    return std::uint32_t(s.nodes.size() - 1);
}

std::uint32_t decode_product(ByteView data, const std::string& label, Context& ctx);

// Node for a linked document, decoding it on first use.
std::uint32_t document_node(const Uuid& uuid, const std::string& instance, Context& ctx) {
    if (auto it = ctx.node_of.find(uuid); it != ctx.node_of.end()) return it->second;
    if (ctx.options.locate_document && !ctx.by_uuid.count(uuid))
        if (auto path = ctx.options.locate_document(hex(uuid))) ctx.by_uuid.emplace(uuid, *path);
    if (!ctx.options.locate_document) ctx.index();
    const auto it = ctx.by_uuid.find(uuid);
    if (it == ctx.by_uuid.end()) {
        ctx.scene.warnings.push_back("instance '" + instance + "': linked document " + hex(uuid) + " not found");
        return add_node(ctx.scene, "missing document " + hex(uuid));
    }
    if (ctx.active.count(uuid)) bad("assembly contains itself through '" + it->second + "'");
    const std::string& path = it->second;
    if (lower_ext(path) != "catproduct") {  // a part: read and decoded later, in parallel
        const std::uint32_t node = add_node(ctx.scene, path);
        ctx.parts.push_back(std::make_unique<Context::Part>(Context::Part{node, path, {}}));
        ctx.node_of.emplace(uuid, node);
        return node;
    }
    const auto bytes = ctx.options.load_resource(path);
    if (!bytes) {
        ctx.scene.warnings.push_back("instance '" + instance + "': cannot read '" + path + "'");
        return add_node(ctx.scene, path);
    }
    const ByteView v(bytes->data(), bytes->size());
    std::uint32_t node;
    ctx.active.insert(uuid);
    try {
        node = decode_product(v, path, ctx);
    } catch (const DecodeError& e) {
        if (std::string(e.what()).rfind(path, 0) == 0) throw;
        fail(e.code(), path + ": " + e.what(), e.offset());
    }
    ctx.active.erase(uuid);
    ctx.node_of.emplace(uuid, node);
    return node;
}

Matrix4 placement(const std::vector<double>& m, std::size_t at) {
    if (m.size() != 12) bad("_Position with " + std::to_string(m.size()) + " values", at);
    Matrix4 t = kIdentity;
    for (int c = 0; c < 3; ++c)
        for (int r = 0; r < 3; ++r) t[std::size_t(c * 4 + r)] = m[std::size_t(c * 3 + r)];
    for (int r = 0; r < 3; ++r) t[std::size_t(12 + r)] = m[std::size_t(9 + r)];
    return t;
}

// The stream of the product container: the one whose objects carry `_PartNumber` (product
// references: ASMPRODUCT, or application product classes such as CATTubTubingLine).
Stream& product_stream(std::vector<Stream>& streams) {
    Stream* found = nullptr;
    for (Stream& s : streams) {
        if (std::find(s.names.begin(), s.names.end(), "_PartNumber") == s.names.end()) continue;
        catia::parse_records(s);
        if (std::none_of(s.objects.begin(), s.objects.end(), [&](const Object& o) { return s.find(o, "_PartNumber"); })) continue;
        if (found) bad("more than one product object stream");
        found = &s;
    }
    if (!found) bad("no product object stream");
    return *found;
}

// The OSM streams of the file, or its SPEC containers (R11-era files) when no OSM stream holds
// product references.
// A SPEC container that may hold product references and failed is fatal; any other failed
// container is reported in `warnings` (its content is lost, the product structure is not).
std::vector<Stream> object_streams(ByteView data, std::vector<std::string>& warnings) {
    auto streams = catia::parse_streams(data);
    const bool osm_products = std::any_of(streams.begin(), streams.end(), [](const Stream& s) {
        return std::find(s.names.begin(), s.names.end(), "_PartNumber") != s.names.end();
    });
    if (osm_products) return streams;
    catia::SpecStreams spec = catia::parse_spec_streams(data);
    for (const auto& f : spec.failed) {
        if (f.has_products) fail(f.code, "product container not decoded: " + f.message, f.error_offset);
        warnings.push_back("object container at Data byte " + std::to_string(f.offset) + " (" + std::to_string(f.size) +
                           " bytes) not decoded, its content is lost: " + f.message + " (at byte " + std::to_string(f.error_offset) + ")");
    }
    return std::move(spec.streams);
}

// SPEC containers name no product classes; a reference's application content shows as
// attributes outside the product set.
bool product_attribute(std::string_view a) {
    static const std::set<std::string_view> k = {"_ViewsList", "_Connectors", "_RepresentedBy", "_Component", "_PartNumber", "_BagRepsList",
                                                 "_PrdVersion", "_PubConnectors", "_AttrForLiterals", "_RefProduct", "_BOMRefLiterals",
                                                 "_Nomenclature", "_Revision", "_Definition", "_Source", "_DescriptionRef"};
    return k.count(a) != 0;
}

// The document linked by the first object in `o`'s _BagRepsList whose relation header starts with a link blob.
std::optional<Uuid> bag_link(const Stream& s, const Object& o) {
    for (const auto& rel : catia::relations(s, o)) {
        if (s.name(rel.attr) != "_BagRepsList") continue;
        for (const std::uint32_t b : rel.objects) {
            const Object* bag = s.by_id(b);
            if (!bag) continue;
            const auto tok = catia::relation_header(*bag);
            if (!tok.empty() && tok[0].is_blob) return catia::link_target(tok[0].blob, bag->offset);
        }
    }
    return std::nullopt;
}

std::uint32_t decode_product(ByteView data, const std::string& label, Context& ctx) {
    Scene& scene = ctx.scene;
    auto streams = object_streams(data, scene.warnings);
    const Stream& s = product_stream(streams);

    // References.
    std::map<std::uint32_t, std::uint32_t> ref_node;  // object id -> node
    for (const Object& o : s.objects) {
        const catia::Chunk* pn = s.find(o, "_PartNumber");
        if (!pn) continue;
        if (s.spec_version) {
            std::string extra;
            for (const catia::Chunk& c : o.chunks)
                if (!product_attribute(s.name(c.attr))) extra += (extra.empty() ? "" : ", ") + s.name(c.attr);
            if (!extra.empty())
                scene.warnings.push_back("reference '" + s.name(o.name) + "' has application attributes (" + extra + "); its own content is not decoded");
        } else if (s.name(o.cls) != "ASMPRODUCT") {
            const std::string pn_name = pn->value.kind == catia::Value::Kind::Name ? s.name(std::uint32_t(pn->value.i)) : s.name(o.name);
            scene.warnings.push_back("reference '" + pn_name + "' is a " + s.name(o.cls) + " (an application product class, e.g. a Human Builder manikin or a piping run); its own content is not decoded");
        }
        std::string name = pn->value.kind == catia::Value::Kind::Name ? s.name(std::uint32_t(pn->value.i)) : s.name(o.name);
        ref_node.emplace(o.id, add_node(scene, name));
    }
    // References defined in another document (e.g. schematic/catalog components): a relation
    // header "<blob> <index>" without a parent. Their content is the linked document.
    std::vector<std::pair<std::uint32_t, Uuid>> linked_refs;
    for (const Object& o : s.objects) {
        if (o.flag == 8 || !(s.find(o, "_InstanceName") || s.find(o, "_PartNumber"))) continue;
        const auto tok = catia::relation_header(o);
        if (tok.size() == 2 && tok[0].is_blob && tok[1].is_number()) {
            if (!ref_node.count(o.id)) ref_node.emplace(o.id, add_node(scene, s.name(o.name)));
            linked_refs.emplace_back(o.id, catia::link_target(tok[0].blob, o.offset));
        }
    }

    // Link indices: a relation header "<blob> <index> ..." defines link <index>; later headers
    // may name the same index without repeating the blob.
    std::map<std::uint64_t, Uuid> links;
    for (const Object& o : s.objects) {
        if (o.flag == 8) continue;
        const auto tok = catia::relation_header(o);
        if (tok.size() >= 2 && tok[0].is_blob && tok[1].is_number()) {
            const Uuid u = catia::link_target(tok[0].blob, o.offset);
            const auto [it, fresh] = links.emplace(tok[1].number, u);
            if (!fresh && it->second != u) bad("link index " + std::to_string(tok[1].number) + " names two documents", o.offset);
        }
    }

    for (const auto& [id, uuid] : linked_refs) {
        const std::uint32_t node = ref_node.at(id);
        const std::uint32_t doc = document_node(uuid, scene.nodes[node].name, ctx);
        scene.nodes[node].children.push_back(doc);
    }

    // Instances. In SPEC containers an instance may lack _InstanceName (version 9): it is then an
    // object a reference lists in _Component, named by its object name; the header holds the
    // parent alone and the startup links the instantiated document.
    std::set<std::uint32_t> components;
    if (s.spec_version)
        for (const auto& [id, node] : ref_node)
            for (const auto& rel : catia::relations(s, *s.by_id(id)))
                if (s.name(rel.attr) == "_Component") components.insert(rel.objects.begin(), rel.objects.end());
    std::set<std::uint32_t> used_refs;
    for (const Object& o : s.objects) {
        const catia::Chunk* in = s.find(o, "_InstanceName");
        const bool startup_instance = !in && components.count(o.id);
        if (!in && !startup_instance) continue;
        if (std::any_of(linked_refs.begin(), linked_refs.end(), [&](const auto& l) { return l.first == o.id; })) continue;
        if (in && in->value.kind != catia::Value::Kind::Name && in->value.kind != catia::Value::Kind::Empty)
            bad("_InstanceName is not a name", o.offset);
        const std::string name = !in ? s.name(o.name) : in->value.kind == catia::Value::Kind::Name ? s.name(std::uint32_t(in->value.i)) : std::string();
        const auto tok = catia::relation_header(o);
        std::uint32_t parent = 0, child_node = 0;
        if (startup_instance) {
            if (tok.size() != 1 || !tok[0].is_number() || !o.startup.has_blob)
                bad("instance '" + name + "' (SPEC v" + std::to_string(s.spec_version) + ") without a parent header and a linked startup", o.offset);
            parent = std::uint32_t(tok[0].number);
            child_node = document_node(catia::link_target(o.startup.blob, o.offset), name, ctx);
        } else if (tok.size() == 3 && tok[0].is_blob && tok[1].is_number() && tok[2].is_number()) {
            parent = std::uint32_t(tok[2].number);
            child_node = document_node(catia::link_target(tok[0].blob, o.offset), name, ctx);
        } else if (tok.size() == 2 && tok[0].is_number() && tok[1].is_number()) {
            if (const auto r = ref_node.find(std::uint32_t(tok[0].number)); r != ref_node.end()) {
                used_refs.insert(r->first);
                child_node = r->second;
            } else if (const auto u = bag_link(s, o)) {
                // Older layout: the instance's bag of representations carries the document link.
                child_node = document_node(*u, name, ctx);
            } else {
                bad("instance '" + name + "' refers to " + std::to_string(tok[0].number) + ", neither a product reference nor a linked bag", o.offset);
            }
            parent = std::uint32_t(tok[1].number);
        } else {
            std::string form;
            for (const auto& t : tok) form += t.is_blob ? " blob" : t.is_guid ? " guid" : t.is_index ? " index" : " " + std::to_string(t.number);
            bad("instance '" + name + "' (object " + std::to_string(o.id) + ", " + s.name(o.cls) + ") with an unknown relation header:" + form, o.offset);
        }
        const auto p = ref_node.find(parent);
        if (p == ref_node.end()) bad("instance '" + name + "' has parent object " + std::to_string(parent) + ", not a product reference", o.offset);

        // Placement from the instance's bag of representations.
        Matrix4 t = kIdentity;
        int placed = 0;
        for (const auto& rel : catia::relations(s, o)) {
            if (s.name(rel.attr) != "_BagRepsList") continue;
            for (const std::uint32_t b : rel.objects) {
                const Object* bag = s.by_id(b);
                if (!bag) bad("instance '" + name + "' lists missing bag object " + std::to_string(b), o.offset);
                if (const catia::Chunk* pos = s.find(*bag, "_Position")) {
                    if (pos->value.kind != catia::Value::Kind::Doubles) bad("_Position is not a number array", bag->offset);
                    t = placement(pos->value.d, bag->offset);
                    ++placed;
                }
            }
        }
        if (placed > 1) bad("instance '" + name + "' has " + std::to_string(placed) + " positions", o.offset);
        if (startup_instance) {
            // In the one version-9 sample the instance's _RepresentedBy object carries an identity
            // _Position and its bag none; whether that object places the instance is not established.
            for (const auto& rel : catia::relations(s, o)) {
                if (s.name(rel.attr) != "_RepresentedBy") continue;
                for (const std::uint32_t r : rel.objects)
                    if (const Object* rep = s.by_id(r))
                        if (const catia::Chunk* pos = s.find(*rep, "_Position"))
                            if (pos->value.kind != catia::Value::Kind::Doubles || placement(pos->value.d, rep->offset) != kIdentity)
                                fail(ErrorCode::Unsupported, "CATProduct: instance '" + name + "': a non-identity _Position on its _RepresentedBy object (SPEC v" +
                                                                 std::to_string(s.spec_version) + "); its meaning is not established", rep->offset);
            }
        }
        const std::uint32_t inst = add_node(scene, name, t);
        scene.nodes[inst].children.push_back(child_node);
        scene.nodes[p->second].children.push_back(inst);
    }

    // A reference may carry its own shape: _BagRepsList -> PRDBAGREP, _Reps -> PRDREP whose
    // _Shapes link names the shape document. Those links are indices without a document blob,
    // and no table in the file maps them to a document yet; such shapes are reported, not dropped
    // silently.
    for (const auto& [id, node] : ref_node) {
        const Object* ref = s.by_id(id);
        if (!ref || !scene.nodes[node].children.empty()) continue;
        std::size_t shapes = 0;
        for (const auto& rel : catia::relations(s, *ref)) {
            if (s.name(rel.attr) != "_BagRepsList") continue;
            for (const std::uint32_t b : rel.objects)
                if (const Object* bag = s.by_id(b))
                    for (const auto& r2 : catia::relations(s, *bag))
                        if (s.name(r2.attr) == "_Reps")
                            for (const std::uint32_t rp : r2.objects)
                                if (const Object* rep = s.by_id(rp))
                                    if (const catia::Chunk* sh = s.find(*rep, "_Shapes"); sh && !sh->value.links.empty()) ++shapes;
        }
        if (shapes)
            scene.warnings.push_back("reference '" + scene.nodes[node].name + "': " + std::to_string(shapes) +
                                     " shape representation link(s) without a document link (not resolved, geometry missing)");
    }

    // The document's product is its first reference; others nothing instantiates are kept
    // out of the tree and reported.
    if (ref_node.empty()) bad("no product reference in '" + label + "'");
    const std::uint32_t root_id = ref_node.begin()->first;
    if (used_refs.count(root_id)) bad("the first product reference of '" + label + "' is instantiated");
    std::size_t unused = 0;
    for (const auto& [id, node] : ref_node)
        if (id != root_id && !used_refs.count(id) && scene.nodes[node].children.empty()) ++unused;
    if (unused) scene.notes.push_back(label + ": " + std::to_string(unused) + " product references are not instantiated (not shown)");
    return ref_node.begin()->second;
}

Scene decode(ByteView d, const DecodeOptions& options) {
    if (d.size < 8 || std::memcmp(d.data, "V5_CFV2", 8) != 0) bad("missing V5_CFV2 magic");
    Scene scene;
    scene.format = "catproduct";
    scene.nodes.push_back(Node{});  // nodes[0]: the root, filled below
    Context ctx(options, scene);
    if (const auto u = catia::document_uuid(d)) ctx.active.insert(*u);
    const std::uint32_t product = decode_product(d, "input", ctx);
    scene.nodes[0].name = scene.nodes[product].name;
    scene.nodes[0].children.push_back(product);
    std::vector<std::function<void()>> tasks;
    for (auto& p : ctx.parts)
        tasks.push_back([&options, part = p.get()] {
            const auto bytes = options.load_resource(part->path);
            if (!bytes) {
                part->scene.warnings.push_back("cannot read the linked part");
                return;
            }
            try {
                part->scene = catpart().decode(ByteView(bytes->data(), bytes->size()), options);
            } catch (const DecodeError& e) {
                fail(e.code(), part->path + ": " + e.what(), e.offset());
            }
        });
    run_all(tasks, options.threads);
    for (auto& p : ctx.parts) {
        const auto base = std::uint32_t(scene.meshes.size());
        for (auto& m : p->scene.meshes) scene.meshes.push_back(std::move(m));
        for (std::uint32_t i = base; i < scene.meshes.size(); ++i) scene.nodes[p->node].meshes.push_back(i);
        for (const auto& w : p->scene.warnings) scene.warnings.push_back(p->path + ": " + w);
        for (const auto& n : p->scene.notes) scene.notes.push_back(p->path + ": " + n);
    }
    return scene;
}

}  // namespace

const Decoder& catproduct() {
    static const SimpleDecoder d({"catproduct", "CATIA V5 Product (assembly)", {"catproduct"}}, sniff, decode);
    return d;
}

}  // namespace cadcore::formats
