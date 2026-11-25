// 3D XML (.3dxml): zip archive with Manifest.xml naming the root product structure file.
//
// Product structure (root file):
//   <ProductStructure root="R">
//     <Reference3D id name/>                                  product or part
//     <Instance3D id name> IsAggregatedBy(parent) IsInstanceOf(reference) RelativeMatrix(12) </Instance3D>
//     <ReferenceRep id format="TESSELLATED" associatedFile="urn:3DXML:file.3DRep"/>
//     <InstanceRep id> IsAggregatedBy(reference3D) IsInstanceOf(referenceRep) </InstanceRep>
// RelativeMatrix: x axis, y axis, z axis (3 values each), translation.
// Representations (.3DRep): XML (format TESSELLATED; PolygonalRepType: Faces/Face
// triangles|strips|fans, VertexBuffer Positions "x y z,x y z", Edges/Polyline) or binary
// (format UVR): a V5_CFV3 store whose Skeleton stream holds the tessellation (catia/tessellation.hpp).
// Colours: Face or Rep SurfaceAttributes/Color, or a MaterialId into CATMaterialRef.3dxml whose
// rendering domain carries DiffuseColor / Transparency.
#include <charconv>
#include <cmath>
#include <exception>
#include <functional>
#include <map>
#include <memory>
#include <string>

#include "cadcore/error.hpp"
#include "codec/xml.hpp"
#include "formats/catia/tessellation.hpp"
#include "formats/catia/store.hpp"
#include "core/parallel.hpp"
#include "formats/registry.hpp"
#include "formats/zip.hpp"

namespace cadcore::formats {
namespace {

// A resource the host loaded, as a Buffer like zip entries.
std::optional<Buffer> loaded(std::optional<std::vector<std::uint8_t>> v) {
    if (!v) return std::nullopt;
    return Buffer(v->begin(), v->end());
}

[[noreturn]] void bad(const std::string& w) { fail(ErrorCode::Malformed, "3DXML: " + w); }

int sniff(ByteView d, std::string_view ext) {
    if (ext == "3dxml" && d.starts_with("PK\x03\x04")) return 2;
    return ext == "3dxml" ? 1 : 0;
}

bool space(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }

// Numbers separated by spaces and/or commas.
template <class T>
void numbers(std::string_view s, std::vector<T>& out, const char* what) {
    std::size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && (space(s[i]) || s[i] == ',')) ++i;
        if (i == s.size()) break;
        T v{};
        auto r = std::from_chars(s.data() + i, s.data() + s.size(), v);
        if (r.ec != std::errc()) bad(std::string("bad number in ") + what + " near '" + std::string(s.substr(i, 16)) + "'");
        out.push_back(v);
        i = std::size_t(r.ptr - s.data());
        if (i < s.size() && !space(s[i]) && s[i] != ',') bad(std::string("bad separator in ") + what);
    }
}

// Groups separated by commas (strips, fans), each a list of indices.
std::vector<std::vector<std::uint32_t>> groups(std::string_view s, const char* what) {
    std::vector<std::vector<std::uint32_t>> out;
    std::size_t start = 0;
    while (start <= s.size()) {
        std::size_t end = s.find(',', start);
        if (end == std::string_view::npos) end = s.size();
        std::vector<std::uint32_t> g;
        numbers(s.substr(start, end - start), g, what);
        if (!g.empty()) out.push_back(std::move(g));
        start = end + 1;
    }
    return out;
}

std::uint8_t unit_u8(double v) { return std::uint8_t(v <= 0 ? 0 : v >= 1 ? 255 : std::lround(v * 255)); }

std::string strip_urn(std::string_view s) {
    constexpr std::string_view k = "urn:3DXML:";
    if (s.substr(0, k.size()) == k) s.remove_prefix(k.size());
    return std::string(s);
}

std::uint32_t child_id(const xml::Element& e, std::string_view tag) {
    const xml::Element* c = e.child(tag);
    if (!c) bad(std::string(e.local()) + " without " + std::string(tag));
    std::vector<std::uint32_t> v;
    numbers(c->text, v, std::string(tag).c_str());
    if (v.size() != 1) bad(std::string(tag) + " is not one id");
    return v[0];
}

struct Decoder3dxml {
    const ZipArchive& zip;
    const DecodeOptions& opt;
    Scene& scene;
    std::map<std::string, std::unique_ptr<xml::Document>> docs;   // parsed archive files
    std::map<std::string, Buffer> bytes;
    std::map<std::string, Color> material_colors;                  // "file#id" -> colour

    const xml::Document* doc(const std::string& name) {
        if (auto it = docs.find(name); it != docs.end()) return it->second.get();
        auto data = zip.read(name);
        if (!data && opt.load_resource) data = loaded(opt.load_resource(name));
        if (!data) return nullptr;
        auto& b = bytes[name] = std::move(*data);
        auto d = std::make_unique<xml::Document>(ByteView(b.data(), b.size()));
        return (docs[name] = std::move(d)).get();
    }

    std::optional<Color> material(std::string_view ref) {
        const std::string key = strip_urn(ref);
        if (auto it = material_colors.find(key); it != material_colors.end()) return it->second;
        const auto hash = key.find('#');
        if (hash == std::string::npos) return std::nullopt;
        const std::string file = key.substr(0, hash);
        const std::string id = key.substr(hash + 1);
        const xml::Document* d = doc(file);
        if (!d) return std::nullopt;
        // CATMatReference id -> MaterialDomainInstance (IsAggregatedBy id) -> MaterialDomain file.
        std::string domain;
        const xml::Element* lib = d->root().child("CATMaterialRef");
        if (!lib) return std::nullopt;
        for (const xml::Element* e : lib->children) {
            if (e->local() != "MaterialDomainInstance") continue;
            const xml::Element* agg = e->child("IsAggregatedBy");
            if (!agg || std::string(agg->text) != id) continue;
            const xml::Element* inst = e->child("IsInstanceOf");
            if (!inst) continue;
            for (const xml::Element* m : lib->children)
                if (m->local() == "MaterialDomain" && m->attr("id") == inst->text) domain = strip_urn(m->attr("associatedFile"));
        }
        if (domain.empty()) return std::nullopt;
        const xml::Document* r = doc(domain);
        if (!r) return std::nullopt;
        std::vector<double> diffuse;
        double transparency = 0;
        std::vector<const xml::Element*> stack{&r->root()};
        while (!stack.empty()) {
            const xml::Element* e = stack.back();
            stack.pop_back();
            for (const xml::Element* c : e->children) stack.push_back(c);
            if (e->local() != "Attr") continue;
            if (e->attr("Name") == "DiffuseColor") {
                std::string_view v = e->attr("Value");
                if (!v.empty() && v.front() == '[') v.remove_prefix(1);
                if (!v.empty() && v.back() == ']') v.remove_suffix(1);
                numbers(v, diffuse, "DiffuseColor");
            } else if (e->attr("Name") == "TextureImage" && !e->attr("Value").empty()) {
                scene.warnings.push_back("material " + key + ": texture image not applied (diffuse colour used)");
            } else if (e->attr("Name") == "Transparency") {
                std::vector<double> t;
                numbers(e->attr("Value"), t, "Transparency");
                if (t.size() == 1) transparency = t[0];
            }
        }
        if (diffuse.size() != 3) return std::nullopt;
        const Color c{unit_u8(diffuse[0]), unit_u8(diffuse[1]), unit_u8(diffuse[2]), unit_u8(1 - transparency)};
        material_colors.emplace(key, c);
        return c;
    }

    std::optional<Color> surface_color(const xml::Element* attrs, const std::string& where) {
        if (!attrs) return std::nullopt;
        if (const xml::Element* c = attrs->child("Color")) {
            std::vector<double> v;
            for (const char* k : {"red", "green", "blue", "alpha"}) numbers(c->attr(k), v, "Color");
            if (v.size() != 4) bad("Color needs red, green, blue, alpha (" + where + ")");
            return Color{unit_u8(v[0]), unit_u8(v[1]), unit_u8(v[2]), unit_u8(v[3])};
        }
        if (const xml::Element* ma = attrs->child("MaterialApplication"))
            if (const xml::Element* mid = ma->child("MaterialId")) {
                if (auto c = material(mid->attr("id"))) return c;
                scene.warnings.push_back("material " + std::string(mid->attr("id")) + " could not be resolved (" + where + ")");
            }
        return std::nullopt;
    }

    // One PolygonalRepType -> one mesh (+ polylines for its edges).
    void polygonal(const xml::Element& rep, const std::string& where, std::vector<std::uint32_t>& meshes, std::vector<std::uint32_t>& lines,
                   Color inherited) {
        const xml::Element* vb = rep.child("VertexBuffer");
        if (!vb) bad("PolygonalRep without VertexBuffer (" + where + ")");
        Mesh m;
        m.name = where;
        if (const xml::Element* p = vb->child("Positions")) numbers(p->text, m.positions, "Positions");
        if (m.positions.size() % 3) bad("Positions not a multiple of 3 (" + where + ")");
        if (const xml::Element* n = vb->child("Normals")) {
            numbers(n->text, m.normals, "Normals");
            if (m.normals.size() != m.positions.size()) bad("Normals count differs from Positions (" + where + ")");
        }
        const std::size_t nv = m.vertex_count();
        Color rep_color = surface_color(rep.child("SurfaceAttributes"), where).value_or(inherited);
        m.color = rep_color;
        auto check = [&](std::uint32_t i) { if (i >= nv) bad("vertex index out of range (" + where + ")"); };
        if (const xml::Element* faces = rep.child("Faces")) {
            for (const xml::Element* f : faces->children) {
                if (f->local() != "Face") continue;
                const auto first = std::uint32_t(m.triangle_count());
                std::vector<std::uint32_t> t;
                numbers(f->attr("triangles"), t, "triangles");
                if (t.size() % 3) bad("triangles not a multiple of 3 (" + where + ")");
                for (auto i : t) check(i);
                m.indices.insert(m.indices.end(), t.begin(), t.end());
                for (const auto& s : groups(f->attr("strips"), "strips")) {
                    for (auto i : s) check(i);
                    for (std::size_t k = 0; k + 2 < s.size(); ++k)
                        if (k % 2 == 0) m.indices.insert(m.indices.end(), {s[k], s[k + 1], s[k + 2]});
                        else m.indices.insert(m.indices.end(), {s[k + 1], s[k], s[k + 2]});
                }
                for (const auto& s : groups(f->attr("fans"), "fans")) {
                    for (auto i : s) check(i);
                    for (std::size_t k = 1; k + 1 < s.size(); ++k) m.indices.insert(m.indices.end(), {s[0], s[k], s[k + 1]});
                }
                const Color fc = surface_color(f->child("SurfaceAttributes"), where).value_or(rep_color);
                const auto count = std::uint32_t(m.triangle_count()) - first;
                if (count) m.groups.push_back({first, count, fc});
            }
        }
        if (const xml::Element* edges = rep.child("Edges")) {
            for (const xml::Element* pl : edges->children) {
                if (pl->local() != "Polyline") continue;
                Polyline line;
                numbers(pl->attr("vertices"), line.points, "Polyline");
                if (line.points.size() % 3) bad("Polyline not a multiple of 3 (" + where + ")");
                lines.push_back(std::uint32_t(scene.polylines.size()));
                scene.polylines.push_back(std::move(line));
            }
        }
        if (!m.indices.empty()) {
            meshes.push_back(std::uint32_t(scene.meshes.size()));
            scene.meshes.push_back(std::move(m));
        }
    }

    void rep_tree(const xml::Element& e, const std::string& where, std::vector<std::uint32_t>& meshes, std::vector<std::uint32_t>& lines, Color inherited) {
        for (const xml::Element* c : e.children) {
            if (c->local() != "Rep" && c->local() != "Root") continue;
            const std::string_view type = c->attr("xsi:type");
            if (type == "PolygonalRepType") {
                polygonal(*c, where, meshes, lines, inherited);
            } else if (type == "BagRepType" || c->local() == "Root") {
                const Color bag = surface_color(c->child("SurfaceAttributes"), where).value_or(inherited);
                rep_tree(*c, where, meshes, lines, bag);
            } else {
                scene.warnings.push_back("representation type '" + std::string(type) + "' not decoded (" + where + ")");
            }
        }
    }

    // Binary (UVR) representations decoded ahead in parallel, by file; an error is kept and raised
    // when its representation is reached, so warnings and errors keep their order.
    struct Uvr { bool missing = false; std::vector<catia::TessBody> bodies; std::exception_ptr error; };
    std::map<std::string, Uvr> uvr;

    void prefetch_uvr(const std::map<std::uint32_t, const xml::Element*>& reps) {
        for (const auto& [id, rr] : reps)
            if (rr->attr("format") == "UVR") uvr.emplace(strip_urn(rr->attr("associatedFile")), Uvr{});
        std::vector<std::function<void()>> tasks;
        for (auto& [file, u] : uvr)
            tasks.push_back([this, &file = file, &u = u] {
                try {
                    auto data = zip.read(file);
                    if (!data && opt.load_resource) data = loaded(opt.load_resource(file));
                    if (!data) { u.missing = true; return; }
                    const ByteView v(data->data(), data->size());
                    if (!v.starts_with(std::string_view("V5_CFV3\0", 8))) bad("UVR representation '" + file + "' is not a V5_CFV3 store");
                    const catia::Store store = catia::read_store(v);
                    const catia::StreamInfo* sk = store.find("Skeleton");
                    if (!sk) bad("UVR representation '" + file + "' has no Skeleton stream");
                    const auto stream = catia::read_stream(v, store, *sk);
                    u.bodies = catia::decode_tessellation(ByteView(stream.data(), stream.size()), catia::Generation::Cfv3);
                } catch (...) {
                    u.error = std::current_exception();
                }
            });
        run_all(tasks, opt.threads);
    }

    // Node holding the geometry of one ReferenceRep.
    std::uint32_t reference_rep(const xml::Element& rr) {
        const std::string file = strip_urn(rr.attr("associatedFile"));
        Node n;
        n.name = std::string(rr.attr("name"));
        const auto node = std::uint32_t(scene.nodes.size());
        scene.nodes.push_back(n);
        const std::string_view format = rr.attr("format");
        if (format != "TESSELLATED" && format != "UVR") {
            scene.warnings.push_back("representation '" + n.name + "' has format '" + std::string(format) + "' (not decoded)");
            return node;
        }
        if (format == "UVR") {
            Uvr& u = uvr.at(file);
            if (u.error) std::rethrow_exception(u.error);
            if (u.missing) {
                scene.warnings.push_back("representation file '" + file + "' is missing");
                return node;
            }
            if (u.bodies.empty()) scene.warnings.push_back("representation '" + file + "' holds no tessellation");
            catia::add_to_scene(std::move(u.bodies), file, scene, node);
            return node;
        }
        auto data = zip.read(file);
        if (!data && opt.load_resource) data = loaded(opt.load_resource(file));
        if (!data) {
            scene.warnings.push_back("representation file '" + file + "' is missing");
            return node;
        }
        bytes[file] = std::move(*data);
        xml::Document d(ByteView(bytes[file].data(), bytes[file].size()));
        std::vector<std::uint32_t> meshes, lines;
        rep_tree(d.root(), file, meshes, lines, Color{});
        scene.nodes[node].meshes = std::move(meshes);
        scene.nodes[node].polylines = std::move(lines);
        return node;
    }

    void run() {
        auto manifest = zip.read("Manifest.xml");
        if (!manifest) bad("no Manifest.xml");
        xml::Document md(ByteView(manifest->data(), manifest->size()));
        const xml::Element* root_el = md.root().child("Root");
        if (!root_el) bad("Manifest.xml without Root");
        const std::string root_file(root_el->text);
        const xml::Document* d = doc(root_file);
        if (!d) bad("root file '" + root_file + "' is missing");
        const xml::Element* ps = d->root().child("ProductStructure");
        if (!ps) bad("no ProductStructure in '" + root_file + "'");

        std::map<std::uint32_t, std::uint32_t> ref_node;     // Reference3D id -> node
        std::map<std::uint32_t, const xml::Element*> reps;    // ReferenceRep id -> element
        auto id_of = [](const xml::Element& e) {
            std::vector<std::uint32_t> v;
            numbers(e.attr("id"), v, "id");
            if (v.size() != 1) bad(std::string(e.local()) + " without a numeric id");
            return v[0];
        };
        for (const xml::Element* e : ps->children) {
            if (e->local() == "Reference3D") {
                Node n;
                n.name = std::string(e->attr("name"));
                ref_node[id_of(*e)] = std::uint32_t(scene.nodes.size());
                scene.nodes.push_back(n);
            } else if (e->local() == "ReferenceRep") {
                reps[id_of(*e)] = e;
            }
        }
        prefetch_uvr(reps);
        std::map<std::uint32_t, std::uint32_t> rep_node;
        for (const xml::Element* e : ps->children) {
            const std::string_view tag = e->local();
            if (tag == "Instance3D") {
                const auto parent = ref_node.find(child_id(*e, "IsAggregatedBy"));
                const auto ref = ref_node.find(child_id(*e, "IsInstanceOf"));
                if (parent == ref_node.end() || ref == ref_node.end()) bad("Instance3D " + std::string(e->attr("id")) + " links an unknown Reference3D");
                Node n;
                n.name = std::string(e->attr("name"));
                if (const xml::Element* m = e->child("RelativeMatrix")) {
                    std::vector<double> v;
                    numbers(m->text, v, "RelativeMatrix");
                    if (v.size() != 12) bad("RelativeMatrix with " + std::to_string(v.size()) + " values");
                    n.transform = {v[0], v[1], v[2], 0, v[3], v[4], v[5], 0, v[6], v[7], v[8], 0, v[9], v[10], v[11], 1};
                }
                n.children.push_back(ref->second);
                const auto node = std::uint32_t(scene.nodes.size());
                scene.nodes.push_back(n);
                scene.nodes[parent->second].children.push_back(node);
            } else if (tag == "InstanceRep") {
                const auto parent = ref_node.find(child_id(*e, "IsAggregatedBy"));
                const std::uint32_t rid = child_id(*e, "IsInstanceOf");
                const auto rep = reps.find(rid);
                if (parent == ref_node.end() || rep == reps.end()) bad("InstanceRep " + std::string(e->attr("id")) + " links an unknown reference");
                auto it = rep_node.find(rid);
                if (it == rep_node.end()) it = rep_node.emplace(rid, reference_rep(*rep->second)).first;
                scene.nodes[parent->second].children.push_back(it->second);
            }
        }
        std::vector<std::uint32_t> root_id;
        numbers(ps->attr("root"), root_id, "root");
        if (root_id.size() != 1 || !ref_node.count(root_id[0])) bad("ProductStructure root is not a Reference3D");
        scene.nodes[0].name = scene.nodes[ref_node[root_id[0]]].name;
        scene.nodes[0].children.push_back(ref_node[root_id[0]]);
    }
};

Scene decode(ByteView d, const DecodeOptions& opt) {
    ZipArchive zip(d, opt.threads);
    Scene scene;
    scene.format = "3dxml";
    scene.nodes.push_back(Node{});  // root, filled by run()
    Decoder3dxml dec{zip, opt, scene, {}, {}, {}, {}};
    dec.run();
    return scene;
}

}  // namespace

const Decoder& threedxml() {
    static const SimpleDecoder d({"3dxml", "3D XML (CATIA V5/V6 export)", {"3dxml"}}, sniff, decode);
    return d;
}

}  // namespace cadcore::formats
