// Wavefront OBJ (+ MTL colours).
//
// Supported statements: v (x y z [w] or x y z r g b), f (v, v/t, v//n, v/t/n,
// negative indices, polygons fan-triangulated), l (polylines), o / g (new mesh),
// usemtl (face group colour), mtllib (Kd, d, Tr via the resource loader),
// backslash line continuation. vt / vn / s and unknown statements are ignored.
// Each object becomes one Mesh whose vertices are the positions it references,
// renumbered locally.
#include <string>
#include <unordered_map>

#include "core/text.hpp"
#include "formats/registry.hpp"

namespace cadcore::formats {
namespace {

int sniff(ByteView d, std::string_view ext) {
    if (ext == "obj") return 1;
    return 0;  // OBJ has no signature; extension only
}

struct Material {
    Color color{200, 200, 200, 255};
};

void parse_mtl(std::string_view src, std::unordered_map<std::string, Material>& out) {
    text::Cursor c(src);
    c.set_hash_comments(true);
    Material* cur = nullptr;
    while (!c.eof()) {
        std::string_view kw = c.token();
        auto clamp = [](float v) { return std::uint8_t(v <= 0 ? 0 : v >= 1 ? 255 : v * 255.0f + 0.5f); };
        if (kw == "newmtl") {
            cur = &out[std::string(c.rest_of_line())];
            continue;
        }
        if (cur && kw == "Kd") {
            float r, g, b;
            if (c.try_number(r) && c.try_number(g) && c.try_number(b)) {
                cur->color.r = clamp(r); cur->color.g = clamp(g); cur->color.b = clamp(b);
            }
        } else if (cur && (kw == "d" || kw == "Tr")) {
            float v;
            if (c.try_number(v)) cur->color.a = clamp(kw == "d" ? v : 1.0f - v);
        }
        c.next_line();
    }
}

struct Builder {
    Scene& scene;
    const std::vector<float>& pos;      // global positions
    const std::vector<std::uint8_t>& vcol;  // global vertex colours (may be empty)
    std::vector<std::uint32_t> remap;   // global index -> local index + 1 (0 = unmapped) for the current mesh
    std::vector<std::uint32_t> touched;
    Mesh* mesh = nullptr;
    Color group_color{200, 200, 200, 255};

    void begin(std::string name) {
        for (auto g : touched) remap[g] = 0;
        touched.clear();
        scene.meshes.emplace_back();
        mesh = &scene.meshes.back();
        mesh->name = std::move(name);
        start_group();
    }
    void start_group() {
        if (!mesh->groups.empty() && mesh->groups.back().count == 0) {
            mesh->groups.back().color = group_color;
            return;
        }
        mesh->groups.push_back({std::uint32_t(mesh->triangle_count()), 0, group_color});
    }
    std::uint32_t local(std::uint32_t g) {
        if (remap.size() < pos.size() / 3) remap.resize(pos.size() / 3, 0);
        if (remap[g] == 0) {
            mesh->positions.insert(mesh->positions.end(), &pos[g * 3], &pos[g * 3] + 3);
            if (!vcol.empty()) mesh->colors.insert(mesh->colors.end(), &vcol[g * 4], &vcol[g * 4] + 4);
            remap[g] = std::uint32_t(mesh->vertex_count());
            touched.push_back(g);
        }
        return remap[g] - 1;
    }
};

Scene decode(ByteView d, const DecodeOptions& opt) {
    Scene scene;
    scene.format = "obj";
    // Line continuations are rare; normalise them once instead of handling them in the tokenizer.
    std::string utf8, joined;
    std::string_view src = text::as_utf8(d.text(), utf8);
    if (src.find("\\\n") != std::string_view::npos || src.find("\\\r\n") != std::string_view::npos) {
        joined.reserve(src.size());
        for (std::size_t i = 0; i < src.size(); ++i) {
            if (src[i] == '\\' && i + 1 < src.size() && (src[i + 1] == '\n' || (src[i + 1] == '\r' && i + 2 < src.size() && src[i + 2] == '\n'))) {
                joined.push_back(' ');
                i += src[i + 1] == '\r' ? 2 : 1;
            } else {
                joined.push_back(src[i]);
            }
        }
        src = joined;
    }

    std::vector<float> pos;
    std::vector<std::uint8_t> vcol;
    std::unordered_map<std::string, Material> materials;
    Builder b{scene, pos, vcol, {}, {}};
    std::vector<std::uint32_t> poly;
    std::string pending_name = "default";
    bool have_mesh = false;
    std::size_t texcoords = 0, normals = 0;

    text::Cursor c(src);
    c.set_hash_comments(true);
    std::vector<std::uint32_t> point_refs;
    while (!c.eof()) {
        std::string_view kw = c.token();
        if (kw.empty()) {
            c.next_line();
            continue;
        }
        if (kw == "v") {
            float xyz[3];
            for (float& v : xyz) v = c.number<float>("vertex");
            pos.insert(pos.end(), xyz, xyz + 3);
            float extra[4];
            int ne = 0;
            while (ne < 4 && c.try_number(extra[ne])) ++ne;
            if (ne >= 3) {  // x y z r g b extension
                if (vcol.empty()) vcol.assign((pos.size() / 3 - 1) * 4, 255);
                for (int k = 0; k < 3; ++k) vcol.push_back(std::uint8_t(extra[k] <= 0 ? 0 : extra[k] >= 1 ? 255 : extra[k] * 255 + 0.5f));
                vcol.push_back(255);
            } else if (!vcol.empty()) {
                vcol.insert(vcol.end(), {255, 255, 255, 255});
            }
        } else if (kw == "vt") {
            ++texcoords;
        } else if (kw == "vn") {
            ++normals;
        } else if (kw == "f" || kw == "l") {
            if (!have_mesh) {
                b.begin(pending_name);
                have_mesh = true;
            }
            poly.clear();
            const long long nverts = static_cast<long long>(pos.size() / 3);
            while (!c.at_eol()) {
                long long i = c.integer("face index");
                if (i < 0) i += nverts + 1;
                if (i < 1 || i > nverts) fail(ErrorCode::Malformed, "face index " + std::to_string(i) + " out of range", c.offset());
                poly.push_back(std::uint32_t(i - 1));
                while (c.peek() == '/') {  // skip /vt and /vn parts
                    c.consume('/');
                    long long ignored;
                    c.try_int(ignored);
                }
            }
            if (kw == "f") {
                if (poly.size() < 3) {
                    scene.warnings.push_back("face with fewer than 3 vertices skipped");
                } else {
                    std::uint32_t a = b.local(poly[0]);
                    for (std::size_t k = 1; k + 1 < poly.size(); ++k) {
                        b.mesh->indices.push_back(a);
                        b.mesh->indices.push_back(b.local(poly[k]));
                        b.mesh->indices.push_back(b.local(poly[k + 1]));
                        b.mesh->groups.back().count++;
                    }
                }
            } else if (poly.size() >= 2) {
                Polyline pl;
                for (auto g : poly) pl.points.insert(pl.points.end(), &pos[g * 3], &pos[g * 3] + 3);
                scene.polylines.push_back(std::move(pl));
            }
        } else if (kw == "p") {
            const long long nverts = static_cast<long long>(pos.size() / 3);
            while (!c.at_eol()) {
                long long i = c.integer("point index");
                if (i < 0) i += nverts + 1;
                if (i < 1 || i > nverts) fail(ErrorCode::Malformed, "point index out of range", c.offset());
                point_refs.push_back(std::uint32_t(i - 1));
            }
        } else if (kw == "o" || kw == "g") {
            pending_name = std::string(c.rest_of_line());
            if (have_mesh && b.mesh->triangle_count() == 0) {
                b.mesh->name = pending_name;
            } else {
                have_mesh = false;
            }
            continue;
        } else if (kw == "usemtl") {
            auto it = materials.find(std::string(c.rest_of_line()));
            b.group_color = it != materials.end() ? it->second.color : Color{200, 200, 200, 255};
            if (have_mesh) b.start_group();
            continue;
        } else if (kw == "mtllib") {
            std::string_view file = c.rest_of_line();
            if (!opt.load_resource) {
                scene.warnings.push_back("material library '" + std::string(file) + "' not loaded (no resource loader)");
            } else if (auto bytes = opt.load_resource(file)) {
                parse_mtl({reinterpret_cast<const char*>(bytes->data()), bytes->size()}, materials);
            } else if (auto dot = file.rfind('.'); dot != std::string_view::npos && file.substr(dot) != ".mtl" &&
                       (bytes = opt.load_resource(std::string(file.substr(0, dot)) + ".mtl"))) {
                // Exporters sometimes reference "x.mat" while writing "x.mtl"; assimp applies the same fallback.
                scene.notes.push_back("material library '" + std::string(file) + "' resolved as '" + std::string(file.substr(0, dot)) + ".mtl'");
                parse_mtl({reinterpret_cast<const char*>(bytes->data()), bytes->size()}, materials);
            } else {
                scene.warnings.push_back("material library '" + std::string(file) + "' not found");
            }
            continue;
        }
        c.next_line();
    }
    // Drop meshes that ended up empty (e.g. only "o" statements) and empty groups.
    std::vector<Mesh> kept;
    for (auto& m : scene.meshes) {
        std::vector<FaceGroup> g;
        for (auto& fg : m.groups)
            if (fg.count) g.push_back(fg);
        m.groups = std::move(g);
        if (m.groups.size() == 1) {
            m.color = m.groups[0].color;
            m.groups.clear();
        }
        if (!m.indices.empty()) kept.push_back(std::move(m));
    }
    scene.meshes = std::move(kept);
    if (scene.meshes.empty() && scene.polylines.empty()) {
        // Point cloud: the 'p' statements, or every vertex when there are none.
        if (pos.empty()) fail(ErrorCode::Malformed, "OBJ without vertices");
        Mesh pts;
        pts.name = "points";
        pts.topology = Topology::Points;
        if (point_refs.empty()) {
            pts.positions = pos;
            pts.colors = vcol;
        } else {
            for (auto g : point_refs) {
                pts.positions.insert(pts.positions.end(), &pos[g * 3], &pos[g * 3] + 3);
                if (!vcol.empty()) pts.colors.insert(pts.colors.end(), &vcol[g * 4], &vcol[g * 4] + 4);
            }
        }
        scene.meshes.push_back(std::move(pts));
    }
    return scene;
}

}  // namespace

const Decoder& obj() {
    static const SimpleDecoder d({"obj", "Wavefront OBJ", {"obj"}}, sniff, decode);
    return d;
}

}  // namespace cadcore::formats
