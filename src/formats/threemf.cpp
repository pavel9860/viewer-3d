// 3MF (3D Manufacturing Format) — core spec + materials (basematerials, colorgroup)
// + production extension (p:path model parts).
//
// The package is an OPC zip. The root model is the target of the
// .../3dmanufacturing/2013/01/3dmodel relationship in _rels/.rels
// (fallback /3D/3dmodel.model). Objects with <mesh> become Meshes; triangle
// colours come from pid/p1 (or the object's pid/pindex) into basematerials or
// colorgroup resources. <components> and <build><item> become Nodes with the
// 3x4 row-vector transform "m00 m01 m02 m10 .. m32". Units map to millimetres.
// Numbers and indices are parsed strictly: an invalid value is a Malformed error.
#include <charconv>
#include <map>
#include <memory>
#include <string>

#include "codec/xml.hpp"
#include "core/parallel.hpp"
#include "core/text.hpp"
#include "formats/registry.hpp"
#include "formats/zip.hpp"

namespace cadcore::formats {
namespace {

constexpr std::string_view kModelRel = "http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel";

int sniff(ByteView d, std::string_view ext) {
    if (ext == "3mf" && d.starts_with("PK\x03\x04")) return 2;
    return ext == "3mf" ? 1 : 0;
}

std::string norm_path(std::string_view p) {
    while (!p.empty() && p[0] == '/') p.remove_prefix(1);
    return std::string(p);
}

using text::trim;

// The first n pieces, concatenated.
template <class T> std::vector<T> join(std::vector<Slot<std::vector<T>>>& piece, std::size_t n) {
    if (n == 1) return std::move(piece[0].value);
    std::size_t total = 0;
    for (std::size_t k = 0; k < n; ++k) total += piece[k].value.size();
    std::vector<T> out;
    out.reserve(total);
    for (std::size_t k = 0; k < n; ++k) out.insert(out.end(), piece[k].value.begin(), piece[k].value.end());
    return out;
}

float to_float(std::string_view s, const char* what) {
    s = trim(s);
    if (!s.empty() && s[0] == '+') s.remove_prefix(1);
    float v = 0;
    if (text::parse_real(s.data(), s.data() + s.size(), v) != s.data() + s.size()) fail(ErrorCode::Malformed, std::string("3MF: bad number for ") + what + ": '" + std::string(s) + "'");
    return v;
}

long long to_int(std::string_view s, const char* what) {
    s = trim(s);
    long long v = 0;
    auto r = std::from_chars(s.data(), s.data() + s.size(), v);
    if (r.ec != std::errc() || r.ptr != s.data() + s.size()) fail(ErrorCode::Malformed, std::string("3MF: bad integer for ") + what + ": '" + std::string(s) + "'");
    return v;
}

bool parse_color(std::string_view s, Color& out) {
    if ((s.size() != 7 && s.size() != 9) || s[0] != '#') return false;
    auto h = [](char c) { return c >= '0' && c <= '9' ? c - '0' : c >= 'a' && c <= 'f' ? c - 'a' + 10 : c >= 'A' && c <= 'F' ? c - 'A' + 10 : -1; };
    int v[4] = {0, 0, 0, 255};
    for (std::size_t k = 0; k * 2 + 1 < s.size(); ++k) {
        const int a = h(s[1 + 2 * k]), b = h(s[2 + 2 * k]);
        if (a < 0 || b < 0) return false;
        v[k] = a * 16 + b;
    }
    out = {std::uint8_t(v[0]), std::uint8_t(v[1]), std::uint8_t(v[2]), std::uint8_t(v[3])};
    return true;
}

Matrix4 parse_transform(std::string_view s) {
    Matrix4 m = kIdentity;
    s = trim(s);
    if (s.empty()) return m;
    double v[12];
    std::size_t p = 0;
    for (double& x : v) {
        while (p < s.size() && s[p] == ' ') ++p;
        std::size_t e = p;
        while (e < s.size() && s[e] != ' ') ++e;
        x = to_float(s.substr(p, e - p), "transform");
        p = e;
    }
    // Row-vector convention: p' = p * M  ->  column-major 4x4 for p' = M p.
    return {v[0], v[1], v[2], 0, v[3], v[4], v[5], 0, v[6], v[7], v[8], 0, v[9], v[10], v[11], 1};
}

double unit_scale(std::string_view u) {
    if (u.empty() || u == "millimeter") return 1;
    if (u == "micron") return 0.001;
    if (u == "centimeter") return 10;
    if (u == "inch") return 25.4;
    if (u == "foot") return 304.8;
    if (u == "meter") return 1000;
    fail(ErrorCode::Malformed, "unknown 3MF unit '" + std::string(u) + "'");
}

struct ModelFile {
    std::unique_ptr<xml::Document> doc;
    const xml::Element* model = nullptr;
    std::map<std::string, const xml::Element*> objects;          // id -> <object>
    std::map<std::string, std::vector<Color>> property_groups;   // id -> colours
};

class Decoder3mf {
public:
    Decoder3mf(ByteView d, Scene& s, unsigned threads) : zip_(d, threads), scene_(s), threads_(threads) {}

    void run() {
        std::string root = "3D/3dmodel.model";
        if (auto rels = zip_.read("_rels/.rels")) {
            xml::Document doc(ByteView(rels->data(), rels->size()));
            for (const xml::Element* r : doc.root().children)
                if (r->attr("Type") == kModelRel) root = norm_path(r->attr("Target"));
        }
        ModelFile& mf = model(root);
        scene_.unit_to_mm = unit_scale(mf.model->attr("unit"));
        scene_.nodes.emplace_back().name = "build";
        const xml::Element* build = mf.model->child("build");
        if (!build) fail(ErrorCode::Malformed, "3MF model without <build>");
        for (const xml::Element* item : build->children) {
            if (item->local() != "item") continue;
            const std::string path = item->has_attr("path") ? norm_path(item->attr("path")) : root;
            add_object(path, std::string(item->attr("objectid")), parse_transform(item->attr("transform")), 0, 0);
        }
    }

private:
    ModelFile& model(const std::string& path) {
        auto it = models_.find(path);
        if (it != models_.end()) return *it->second;
        auto bytes = zip_.read(path).value_or(Buffer{});
        if (bytes.empty()) fail(ErrorCode::Malformed, "3MF part '" + path + "' missing from package");
        auto mf = std::make_unique<ModelFile>();
        // vertex and triangle lists are parsed later, in parallel (for_each_element)
        mf->doc = std::make_unique<xml::Document>(std::move(bytes), std::vector<std::string_view>{"vertices", "triangles", "basematerials", "colorgroup"});
        mf->model = &mf->doc->root();
        if (mf->model->local() != "model") fail(ErrorCode::Malformed, "3MF part '" + path + "' root is not <model>");
        if (const xml::Element* res = mf->model->child("resources")) {
            for (const xml::Element* r : res->children) {
                const std::string_view ln = r->local();
                const std::string id(r->attr("id"));
                if (ln == "object") {
                    mf->objects[id] = r;
                } else if (ln == "basematerials" || ln == "colorgroup") {
                    std::vector<Slot<std::vector<Color>>> piece(kMaxThreads);
                    std::size_t n = 0;
                    const std::string_view attr = ln == "basematerials" ? "displaycolor" : "color";
                    mf->doc->for_each_element(*r, threads_, n, [&](std::size_t k, const xml::Element& e) {
                        Color c{200, 200, 200, 255};
                        parse_color(e.attr(attr), c);
                        piece[k].value.push_back(c);
                    });
                    mf->property_groups[id] = join(piece, n);
                }
            }
        }
        auto& ref = *mf;
        models_[path] = std::move(mf);
        return ref;
    }

    std::uint32_t mesh_for(const std::string& path, const std::string& id, ModelFile& mf, const xml::Element& obj) {
        const std::string key = path + "#" + id;
        if (auto it = mesh_ids_.find(key); it != mesh_ids_.end()) return it->second;
        const xml::Element* mesh_node = obj.child("mesh");
        Mesh m;
        m.name = std::string(obj.attr("name"));
        if (m.name.empty()) m.name = "object " + id;
        if (const xml::Element* verts = mesh_node->child("vertices")) {
            std::vector<Slot<std::vector<float>>> piece(kMaxThreads);
            std::size_t n = 0;
            mf.doc->for_each_element(*verts, threads_, n, [&](std::size_t k, const xml::Element& v) {
                if (v.local() != "vertex") return;
                piece[k].value.insert(piece[k].value.end(), {to_float(v.attr("x"), "x"), to_float(v.attr("y"), "y"), to_float(v.attr("z"), "z")});
            });
            m.positions = join(piece, n);
        }
        const auto nv = static_cast<long long>(m.vertex_count());
        const std::string obj_pid(obj.attr("pid"));
        const long long obj_pindex = obj.has_attr("pindex") ? to_int(obj.attr("pindex"), "pindex") : 0;
        auto lookup = [&](const std::string& pid, long long index, Color& out) {
            auto g = mf.property_groups.find(pid);
            if (g == mf.property_groups.end() || index < 0 || index >= static_cast<long long>(g->second.size())) return false;
            out = g->second[std::size_t(index)];
            return true;
        };
        Color obj_color{200, 200, 200, 255};
        lookup(obj_pid, obj_pindex, obj_color);
        m.color = obj_color;
        if (const xml::Element* tris = mesh_node->child("triangles")) {
            std::vector<Slot<std::vector<std::uint32_t>>> piece(kMaxThreads);
            std::vector<Slot<std::vector<Color>>> colors(kMaxThreads);
            std::size_t n = 0;
            mf.doc->for_each_element(*tris, threads_, n, [&](std::size_t k, const xml::Element& t) {
                if (t.local() != "triangle") return;
                const long long a = to_int(t.attr("v1"), "v1"), b = to_int(t.attr("v2"), "v2"), c = to_int(t.attr("v3"), "v3");
                if (a < 0 || b < 0 || c < 0 || a >= nv || b >= nv || c >= nv) fail(ErrorCode::Malformed, "3MF triangle index out of range in object " + id);
                piece[k].value.insert(piece[k].value.end(), {std::uint32_t(a), std::uint32_t(b), std::uint32_t(c)});
                Color col = obj_color;
                if (t.has_attr("pid") || t.has_attr("p1"))
                    lookup(t.has_attr("pid") ? std::string(t.attr("pid")) : obj_pid, t.has_attr("p1") ? to_int(t.attr("p1"), "p1") : obj_pindex, col);
                colors[k].value.push_back(col);
            });
            m.indices = join(piece, n);
            std::uint32_t tri = 0;
            for (std::size_t k = 0; k < n; ++k)
                for (const Color& col : colors[k].value) {
                    if (!m.groups.empty() && m.groups.back().color == col) m.groups.back().count++;
                    else m.groups.push_back({tri, 1, col});
                    ++tri;
                }
        }
        if (m.groups.size() == 1) m.color = m.groups[0].color;
        if (m.groups.size() <= 1) m.groups.clear();
        const auto idx = std::uint32_t(scene_.meshes.size());
        scene_.meshes.push_back(std::move(m));
        mesh_ids_[key] = idx;
        return idx;
    }

    void add_object(const std::string& path, const std::string& id, const Matrix4& t, std::uint32_t parent, int depth) {
        if (depth > 64) fail(ErrorCode::Malformed, "3MF component nesting deeper than 64 (cycle?)");
        ModelFile& mf = model(path);
        auto it = mf.objects.find(id);
        if (it == mf.objects.end()) fail(ErrorCode::Malformed, "3MF object " + id + " not found in " + path);
        const xml::Element& obj = *it->second;
        Node n;
        n.name = std::string(obj.attr("name"));
        n.transform = t;
        const auto nid = std::uint32_t(scene_.nodes.size());
        scene_.nodes.push_back(std::move(n));
        scene_.nodes[parent].children.push_back(nid);
        if (obj.child("mesh")) scene_.nodes[nid].meshes.push_back(mesh_for(path, id, mf, obj));
        if (const xml::Element* comps = obj.child("components")) {
            for (const xml::Element* c : comps->children) {
                if (c->local() != "component") continue;
                const std::string cpath = c->has_attr("path") ? norm_path(c->attr("path")) : path;
                add_object(cpath, std::string(c->attr("objectid")), parse_transform(c->attr("transform")), nid, depth + 1);
            }
        }
    }

    ZipArchive zip_;
    Scene& scene_;
    unsigned threads_;
    std::map<std::string, std::unique_ptr<ModelFile>> models_;
    std::map<std::string, std::uint32_t> mesh_ids_;
};

Scene decode(ByteView d, const DecodeOptions& options) {
    Scene scene;
    scene.format = "3mf";
    Decoder3mf(d, scene, options.threads).run();
    return scene;
}

}  // namespace

const Decoder& threemf() {
    static const SimpleDecoder d({"3mf", "3D Manufacturing Format", {"3mf"}}, sniff, decode);
    return d;
}

}  // namespace cadcore::formats
