// Writers for the canonical models, one per format variant.
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>

#include "testgen/testgen.hpp"
#include "testgen/zipwriter.hpp"

namespace cadcore::testgen {
namespace {

using Bytes = std::vector<std::uint8_t>;
using V3 = std::array<float, 3>;

struct WorldTri {
    V3 p[3];
    Color c;
};

std::vector<Instance> instances_of(const Model& m) {
    if (!m.instances.empty()) return m.instances;
    std::vector<Instance> out;
    for (std::uint32_t i = 0; i < m.bodies.size(); ++i) out.push_back({i, kIdentity});
    return out;
}

V3 apply(const Matrix4& t, const V3& p) {
    V3 r;
    for (int k = 0; k < 3; ++k) r[k] = float(t[k] * p[0] + t[4 + k] * p[1] + t[8 + k] * p[2] + t[12 + k]);
    return r;
}

std::vector<WorldTri> world_triangles(const Model& m) {
    std::vector<WorldTri> out;
    for (const Instance& in : instances_of(m)) {
        const Body& b = m.bodies[in.body];
        for (const Face& f : b.faces)
            for (auto t : triangulate(f)) out.push_back({{apply(in.transform, b.points[t[0]]), apply(in.transform, b.points[t[1]]), apply(in.transform, b.points[t[2]])}, f.color});
    }
    return out;
}

void put_u32le(Bytes& b, std::uint32_t v) { for (int i = 0; i < 4; ++i) b.push_back(std::uint8_t(v >> (8 * i))); }
void put_u16le(Bytes& b, std::uint16_t v) { b.push_back(std::uint8_t(v)); b.push_back(std::uint8_t(v >> 8)); }
void put_u32be(Bytes& b, std::uint32_t v) { for (int i = 3; i >= 0; --i) b.push_back(std::uint8_t(v >> (8 * i))); }
void put_f32le(Bytes& b, float f) { std::uint32_t u; std::memcpy(&u, &f, 4); put_u32le(b, u); }
void put_f32be(Bytes& b, float f) { std::uint32_t u; std::memcpy(&u, &f, 4); put_u32be(b, u); }
void put_str(Bytes& b, const std::string& s) { b.insert(b.end(), s.begin(), s.end()); }
std::string num(double v) { char buf[32]; std::snprintf(buf, sizeof buf, "%.9g", v); return buf; }

}  // namespace

Files write_stl_binary(const Model& m) {
    auto tris = world_triangles(m);
    Bytes b(80, 0);
    std::string h = "cadcore testgen " + m.name;
    std::memcpy(b.data(), h.data(), std::min<std::size_t>(h.size(), 80));
    put_u32le(b, std::uint32_t(tris.size()));
    for (const auto& t : tris) {
        for (int k = 0; k < 3; ++k) put_f32le(b, 0);
        for (const auto& p : t.p)
            for (float v : p) put_f32le(b, v);
        put_u16le(b, 0);
    }
    return {{m.name + ".stl", b}};
}

Files write_stl_ascii(const Model& m) {
    std::string s;
    for (const Instance& in : instances_of(m)) {
        const Body& body = m.bodies[in.body];
        s += "solid " + body.name + "\n";
        for (const Face& f : body.faces)
            for (auto t : triangulate(f)) {
                s += "  facet normal 0 0 0\n    outer loop\n";
                for (auto i : t) {
                    V3 p = apply(in.transform, body.points[i]);
                    s += "      vertex " + num(p[0]) + " " + num(p[1]) + " " + num(p[2]) + "\n";
                }
                s += "    endloop\n  endfacet\n";
            }
        s += "endsolid " + body.name + "\n";
    }
    return {{m.name + ".stl", Bytes(s.begin(), s.end())}};
}

static std::string color_name(Color c) {
    char buf[16];
    std::snprintf(buf, sizeof buf, "c%02x%02x%02x%02x", c.r, c.g, c.b, c.a);
    return buf;
}

Files write_obj(const Model& m) {
    std::string obj = "# cadcore testgen\nmtllib " + m.name + ".mtl\n";
    std::map<std::string, Color> mats;
    std::size_t base = 1;
    int k = 0;
    for (const Instance& in : instances_of(m)) {
        const Body& body = m.bodies[in.body];
        obj += "o " + body.name + "_" + std::to_string(k++) + "\n";
        for (const auto& p0 : body.points) {
            V3 p = apply(in.transform, p0);
            obj += "v " + num(p[0]) + " " + num(p[1]) + " " + num(p[2]) + "\n";
        }
        for (const Face& f : body.faces) {
            std::string mn = color_name(f.color);
            mats[mn] = f.color;
            obj += "usemtl " + mn + "\n";
            for (auto t : triangulate(f))
                obj += "f " + std::to_string(base + t[0]) + " " + std::to_string(base + t[1]) + " " + std::to_string(base + t[2]) + "\n";
        }
        base += body.points.size();
    }
    std::string mtl;
    for (auto& [n, c] : mats)
        mtl += "newmtl " + n + "\nKd " + num(c.r / 255.0) + " " + num(c.g / 255.0) + " " + num(c.b / 255.0) + "\nd " + num(c.a / 255.0) + "\n";
    return {{m.name + ".obj", Bytes(obj.begin(), obj.end())}, {m.name + ".mtl", Bytes(mtl.begin(), mtl.end())}};
}

Files write_ply(const Model& m, int enc) {
    auto tris = world_triangles(m);
    std::string h = std::string("ply\nformat ") + (enc == 0 ? "ascii" : enc == 1 ? "binary_little_endian" : "binary_big_endian") +
                    " 1.0\ncomment cadcore testgen\nelement vertex " + std::to_string(tris.size() * 3) +
                    "\nproperty float x\nproperty float y\nproperty float z\nelement face " + std::to_string(tris.size()) +
                    "\nproperty list uchar int vertex_indices\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha\nend_header\n";
    Bytes b(h.begin(), h.end());
    if (enc == 0) {
        std::string s;
        for (const auto& t : tris)
            for (const auto& p : t.p) s += num(p[0]) + " " + num(p[1]) + " " + num(p[2]) + "\n";
        for (std::size_t i = 0; i < tris.size(); ++i) {
            const Color c = tris[i].c;
            s += "3 " + std::to_string(i * 3) + " " + std::to_string(i * 3 + 1) + " " + std::to_string(i * 3 + 2) + " " + std::to_string(c.r) + " " +
                 std::to_string(c.g) + " " + std::to_string(c.b) + " " + std::to_string(c.a) + "\n";
        }
        put_str(b, s);
    } else {
        auto f = [&](float v) { enc == 1 ? put_f32le(b, v) : put_f32be(b, v); };
        auto u = [&](std::uint32_t v) { enc == 1 ? put_u32le(b, v) : put_u32be(b, v); };
        for (const auto& t : tris)
            for (const auto& p : t.p) { f(p[0]); f(p[1]); f(p[2]); }
        for (std::size_t i = 0; i < tris.size(); ++i) {
            b.push_back(3);
            u(std::uint32_t(i * 3)); u(std::uint32_t(i * 3 + 1)); u(std::uint32_t(i * 3 + 2));
            const Color c = tris[i].c;
            b.insert(b.end(), {c.r, c.g, c.b, c.a});
        }
    }
    return {{m.name + ".ply", b}};
}

Files write_off(const Model& m, bool binary) {
    auto tris = world_triangles(m);
    if (!binary) {
        std::string s = "OFF\n# cadcore testgen\n" + std::to_string(tris.size() * 3) + " " + std::to_string(tris.size()) + " 0\n";
        for (const auto& t : tris)
            for (const auto& p : t.p) s += num(p[0]) + " " + num(p[1]) + " " + num(p[2]) + "\n";
        for (std::size_t i = 0; i < tris.size(); ++i) {
            const Color c = tris[i].c;
            s += "3 " + std::to_string(i * 3) + " " + std::to_string(i * 3 + 1) + " " + std::to_string(i * 3 + 2) + " " + std::to_string(c.r) + " " +
                 std::to_string(c.g) + " " + std::to_string(c.b) + " " + std::to_string(c.a) + "\n";
        }
        return {{m.name + ".off", Bytes(s.begin(), s.end())}};
    }
    Bytes b;
    put_str(b, "OFF BINARY\n");
    put_u32be(b, std::uint32_t(tris.size() * 3));
    put_u32be(b, std::uint32_t(tris.size()));
    put_u32be(b, 0);
    for (const auto& t : tris)
        for (const auto& p : t.p) { put_f32be(b, p[0]); put_f32be(b, p[1]); put_f32be(b, p[2]); }
    for (std::size_t i = 0; i < tris.size(); ++i) {
        put_u32be(b, 3);
        for (int k = 0; k < 3; ++k) put_u32be(b, std::uint32_t(i * 3 + k));
        put_u32be(b, 4);
        const Color c = tris[i].c;
        for (auto v : {c.r, c.g, c.b, c.a}) put_f32be(b, v / 255.0f);
    }
    return {{m.name + ".off", b}};
}

Files write_gltf(const Model& m, bool glb) {
    // One mesh per body, one primitive per face (material = face colour); a node per instance.
    Bytes bin;
    std::string accessors, views, meshes, materials, nodes;
    std::map<std::string, int> mat_index;
    int nacc = 0;
    auto add_view = [&](std::size_t offset, std::size_t length, int target) {
        int id = int(std::count(views.begin(), views.end(), '{'));
        if (!views.empty()) views += ",";
        views += "{\"buffer\":0,\"byteOffset\":" + std::to_string(offset) + ",\"byteLength\":" + std::to_string(length) + ",\"target\":" + std::to_string(target) + "}";
        return id;
    };
    for (std::size_t bi = 0; bi < m.bodies.size(); ++bi) {
        const Body& body = m.bodies[bi];
        float mn[3] = {1e30f, 1e30f, 1e30f}, mx[3] = {-1e30f, -1e30f, -1e30f};
        std::size_t off = bin.size();
        for (const auto& p : body.points)
            for (int k = 0; k < 3; ++k) { put_f32le(bin, p[k]); mn[k] = std::min(mn[k], p[k]); mx[k] = std::max(mx[k], p[k]); }
        int pv = add_view(off, bin.size() - off, 34962);
        int pa = nacc++;
        if (!accessors.empty()) accessors += ",";
        accessors += "{\"bufferView\":" + std::to_string(pv) + ",\"componentType\":5126,\"count\":" + std::to_string(body.points.size()) +
                     ",\"type\":\"VEC3\",\"min\":[" + num(mn[0]) + "," + num(mn[1]) + "," + num(mn[2]) + "],\"max\":[" + num(mx[0]) + "," + num(mx[1]) + "," + num(mx[2]) + "]}";
        std::string prims;
        for (const Face& f : body.faces) {
            auto tris = triangulate(f);
            std::size_t io = bin.size();
            for (auto t : tris)
                for (auto i : t) put_u32le(bin, i);
            int iv = add_view(io, bin.size() - io, 34963);
            int ia = nacc++;
            accessors += ",{\"bufferView\":" + std::to_string(iv) + ",\"componentType\":5125,\"count\":" + std::to_string(tris.size() * 3) + ",\"type\":\"SCALAR\"}";
            std::string cn = color_name(f.color);
            if (!mat_index.count(cn)) {
                int id = int(mat_index.size());
                mat_index[cn] = id;
                if (!materials.empty()) materials += ",";
                materials += "{\"name\":\"" + cn + "\",\"pbrMetallicRoughness\":{\"baseColorFactor\":[" + num(f.color.r / 255.0) + "," + num(f.color.g / 255.0) + "," +
                             num(f.color.b / 255.0) + "," + num(f.color.a / 255.0) + "]}}";
            }
            if (!prims.empty()) prims += ",";
            prims += "{\"attributes\":{\"POSITION\":" + std::to_string(pa) + "},\"indices\":" + std::to_string(ia) + ",\"material\":" + std::to_string(mat_index[cn]) + "}";
        }
        if (!meshes.empty()) meshes += ",";
        meshes += "{\"name\":\"" + body.name + "\",\"primitives\":[" + prims + "]}";
    }
    std::string node_list;
    auto inst = instances_of(m);
    for (std::size_t i = 0; i < inst.size(); ++i) {
        if (!nodes.empty()) { nodes += ","; node_list += ","; }
        std::string mat;
        for (int k = 0; k < 16; ++k) mat += (k ? "," : "") + num(inst[i].transform[k]);
        nodes += "{\"name\":\"instance" + std::to_string(i) + "\",\"mesh\":" + std::to_string(inst[i].body) + ",\"matrix\":[" + mat + "]}";
        node_list += std::to_string(i);
    }
    const std::string bin_name = m.name + ".bin";
    std::string json = "{\"asset\":{\"version\":\"2.0\",\"generator\":\"cadcore testgen\"},\"scene\":0,\"scenes\":[{\"nodes\":[" + node_list + "]}],\"nodes\":[" + nodes +
                       "],\"meshes\":[" + meshes + "],\"materials\":[" + materials + "],\"accessors\":[" + accessors + "],\"bufferViews\":[" + views +
                       "],\"buffers\":[{" + (glb ? std::string() : "\"uri\":\"" + bin_name + "\",") + "\"byteLength\":" + std::to_string(bin.size()) + "}]}";
    if (!glb) return {{m.name + ".gltf", Bytes(json.begin(), json.end())}, {bin_name, bin}};
    while (json.size() % 4) json += ' ';
    while (bin.size() % 4) bin.push_back(0);
    Bytes b;
    put_str(b, "glTF");
    put_u32le(b, 2);
    put_u32le(b, std::uint32_t(12 + 8 + json.size() + 8 + bin.size()));
    put_u32le(b, std::uint32_t(json.size()));
    put_str(b, "JSON");
    put_str(b, json);
    put_u32le(b, std::uint32_t(bin.size()));
    put_str(b, std::string("BIN\0", 4));
    b.insert(b.end(), bin.begin(), bin.end());
    return {{m.name + ".glb", b}};
}

Files write_3mf(const Model& m) {
    std::vector<Color> palette;
    auto color_id = [&](Color c) {
        for (std::size_t i = 0; i < palette.size(); ++i)
            if (palette[i] == c) return i;
        palette.push_back(c);
        return palette.size() - 1;
    };
    std::string objects;
    for (std::size_t bi = 0; bi < m.bodies.size(); ++bi) {
        const Body& body = m.bodies[bi];
        objects += "<object id=\"" + std::to_string(bi + 2) + "\" type=\"model\" name=\"" + body.name + "\"><mesh><vertices>";
        for (const auto& p : body.points) objects += "<vertex x=\"" + num(p[0]) + "\" y=\"" + num(p[1]) + "\" z=\"" + num(p[2]) + "\"/>";
        objects += "</vertices><triangles>";
        for (const Face& f : body.faces) {
            auto ci = color_id(f.color);
            for (auto t : triangulate(f))
                objects += "<triangle v1=\"" + std::to_string(t[0]) + "\" v2=\"" + std::to_string(t[1]) + "\" v3=\"" + std::to_string(t[2]) + "\" pid=\"1\" p1=\"" + std::to_string(ci) + "\"/>";
        }
        objects += "</triangles></mesh></object>";
    }
    std::string base = "<basematerials id=\"1\">";
    for (auto c : palette) {
        char buf[16];
        std::snprintf(buf, sizeof buf, "#%02X%02X%02X%02X", c.r, c.g, c.b, c.a);
        base += std::string("<base name=\"") + buf + "\" displaycolor=\"" + buf + "\"/>";
    }
    base += "</basematerials>";
    std::string items;
    for (const Instance& in : instances_of(m)) {
        const auto& t = in.transform;
        std::string tr = num(t[0]) + " " + num(t[1]) + " " + num(t[2]) + " " + num(t[4]) + " " + num(t[5]) + " " + num(t[6]) + " " + num(t[8]) + " " + num(t[9]) + " " +
                         num(t[10]) + " " + num(t[12]) + " " + num(t[13]) + " " + num(t[14]);
        items += "<item objectid=\"" + std::to_string(in.body + 2) + "\" transform=\"" + tr + "\"/>";
    }
    const std::string model = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<model unit=\"millimeter\" xml:lang=\"en-US\" xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/2015/02\"><resources>" +
                              base + objects + "</resources><build>" + items + "</build></model>";
    const std::string types = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\"><Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/><Default Extension=\"model\" ContentType=\"application/vnd.ms-package.3dmanufacturing-3dmodel+xml\"/></Types>";
    const std::string rels = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\"><Relationship Target=\"/3D/3dmodel.model\" Id=\"rel0\" Type=\"http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel\"/></Relationships>";
    // [Content_Types] stored, the rest deflated: both zip methods are exercised.
    const Bytes out = write_zip({{"[Content_Types].xml", types, false}, {"_rels/.rels", rels, true}, {"3D/3dmodel.model", model, true}});
    return {{m.name + ".3mf", out}};
}

// --- synthetic CATPart --------------------------------------------------------------------------

static void put_count(Bytes& b, std::uint32_t v) {
    if (v < 255) { b.push_back(std::uint8_t(v)); return; }
    b.push_back(0xff);
    put_u32le(b, v);
}

static void write_facet(Bytes& b, const Face& f, const std::vector<std::uint32_t>& remap, int width, bool u16_counts) {
    std::vector<std::uint32_t> tri_refs;
    std::vector<std::vector<std::uint32_t>> strips, fans;
    for (const Prim& p : f.prims) {
        std::vector<std::uint32_t> r;
        for (auto i : p.refs) r.push_back(remap[i]);
        if (p.kind == PrimKind::Triangles) tri_refs.insert(tri_refs.end(), r.begin(), r.end());
        else (p.kind == PrimKind::Strip ? strips : fans).push_back(r);
    }
    std::uint8_t marker = 0x40 | (tri_refs.empty() ? 0 : 1) | (strips.empty() ? 0 : 2) | (fans.empty() ? 0 : 4) | (f.has_normal ? 8 : 0);
    b.push_back(0x01);
    b.push_back(marker);
    if (!tri_refs.empty()) put_count(b, std::uint32_t(tri_refs.size() / 3));
    if (!strips.empty()) put_count(b, std::uint32_t(strips.size()));
    if (!fans.empty()) put_count(b, std::uint32_t(fans.size()));
    b.push_back(0xff);
    std::uint32_t n = std::uint32_t(tri_refs.size());
    for (auto& s : strips) n += std::uint32_t(s.size());
    for (auto& s : fans) n += std::uint32_t(s.size());
    put_u32le(b, n);
    if (f.has_normal)
        for (float v : f.normal) put_f32le(b, v);
    auto loop_count = [&](std::size_t v) {
        if (u16_counts) {
            if (v > 0xFFFF) throw std::runtime_error("u16 loop count cannot hold " + std::to_string(v));
            b.push_back(std::uint8_t(v >> 8));
            b.push_back(std::uint8_t(v));
        }
        else put_count(b, std::uint32_t(v));
    };
    for (auto& s : strips) loop_count(s.size());
    for (auto& s : fans) loop_count(s.size());
    auto put_ref = [&](std::uint32_t v) {
        for (int k = width - 1; k >= 0; --k) b.push_back(std::uint8_t(v >> (8 * k)));
    };
    for (auto v : tri_refs) put_ref(v);
    for (auto& s : strips) for (auto v : s) put_ref(v);
    for (auto& s : fans) for (auto v : s) put_ref(v);
}

static void write_attr(Bytes& b, const Face& f, bool solid) {
    std::uint8_t flags = std::uint8_t((solid ? 0x30 : 0x20) | (f.color.a < 255 ? 0x80 : 0));
    b.insert(b.end(), {flags, 0x04, 0x04, 0xff, f.color.a, f.color.b, f.color.g, f.color.r});
}

static Bytes catpart_body(const Body& body, const CatPartOptions& o) {
    // Elements in lexicographic order, like CATIA's point lists.
    std::vector<std::uint32_t> order(body.points.size());
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](auto a, auto b) { return body.points[a] < body.points[b]; });
    std::vector<std::uint32_t> remap(order.size());
    for (std::uint32_t i = 0; i < order.size(); ++i) remap[order[i]] = i;
    const auto ec = std::uint32_t(order.size());
    Bytes b = {0x04, 0x04, 0xfd, 0xff, 0, 0, 0, 1, 0, 0, 0, 0, 0xff, 0xff, 0x02, 0x00, 0x01};
    put_count(b, ec);
    put_count(b, ec);
    const bool tex = o.encoding == 2;
    b.insert(b.end(), {std::uint8_t(tex ? 1 : 0), 0, 0});
    if (tex) {
        b.push_back(2);
        for (std::uint32_t i = 0; i < ec; ++i) { put_f32le(b, float(i)); put_f32le(b, 0.5f); }
    }
    if (o.encoding == 1) {
        b.push_back(1);
        for (auto i : order)
            for (float v : body.points[i]) put_f32le(b, v);
    } else {
        b.push_back(0);
        std::vector<std::uint8_t> table((ec + 3) / 4, 0);
        std::vector<float> words;
        V3 prev{0, 0, 0};
        for (std::uint32_t k = 0; k < ec; ++k) {
            const V3& p = body.points[order[k]];
            int op;
            if (p == prev) op = 1;
            else if (p[0] == prev[0] && p[1] == prev[1]) { op = 2; words.push_back(p[2]); }
            else if (p[0] == prev[0]) { op = 3; words.insert(words.end(), {p[1], p[2]}); }
            else { op = 0; words.insert(words.end(), {p[0], p[1], p[2]}); }
            table[k >> 2] |= std::uint8_t(op << ((k & 3) * 2));
            prev = p;
        }
        b.insert(b.end(), table.begin(), table.end());
        put_count(b, std::uint32_t(words.size()));
        for (float v : words) put_f32le(b, v);
    }
    b.insert(b.end(), 16, 0);  // stands in for the per-element normal block
    const int width = ec > 0x10000 ? 3 : (o.force_u16 || o.u16_counts || ec > 0x100) ? 2 : 1;
    put_f32le(b, 0.05f);       // sag
    if (o.lod) {
        b.push_back(1);
        Face lod;
        Prim all{PrimKind::Triangles, {}};
        for (const Face& f : body.faces)
            for (auto t : triangulate(f)) all.refs.insert(all.refs.end(), t.begin(), t.end());
        lod.prims.push_back(all);
        write_facet(b, lod, remap, width, o.u16_counts);
        put_f32le(b, 3.0f);
        write_attr(b, body.faces[0], body.solid);
    } else {
        b.push_back(0);
    }
    b.insert(b.end(), {0x04, 0x02, 0x0a});
    put_count(b, std::uint32_t(body.faces.size()));
    for (const Face& f : body.faces) write_facet(b, f, remap, width, o.u16_counts);
    for (const Face& f : body.faces) write_attr(b, f, body.solid);
    b.insert(b.end(), {0x01, 0x01, 0x00});  // empty edge list
    return b;
}

Files write_catpart(const Model& m, const CatPartOptions& o) {
    Bytes body;
    for (const Body& bd : m.bodies) {
        Bytes x = catpart_body(bd, o);
        body.insert(body.end(), x.begin(), x.end());
    }
    const std::string catalog = std::string("CATIA_V5 CB0001\0", 16) + std::string(48, '\0');
    Bytes store;
    put_str(store, std::string("V5_CFV2\0", 8));
    put_u32be(store, std::uint32_t(32 + body.size()));
    put_u32be(store, std::uint32_t(catalog.size()));
    store.insert(store.end(), 8, 0xff);
    store.insert(store.end(), 8, 0);
    store.insert(store.end(), body.begin(), body.end());
    put_str(store, catalog);
    Bytes f;
    put_str(f, std::string("V5_CFV2\0", 8));
    const std::string filler = "cadcore testgen synthetic CATPart: " + m.name;
    put_u32be(f, std::uint32_t(32 + filler.size() + store.size()));
    put_u32be(f, std::uint32_t(catalog.size()));
    f.insert(f.end(), 8, 0xff);
    f.insert(f.end(), 8, 0);
    put_str(f, filler);
    f.insert(f.end(), store.begin(), store.end());
    put_str(f, catalog);
    return {{m.name + ".CATPart", f}};
}

// DXF: one 3DFACE per world triangle, true colour (420) and transparency (440) per face.
// format 0 = ASCII, 1 = binary with u16 group codes (R13+), 2 = binary with one-byte codes (R12).
Files write_dxf(const Model& m, int format) {
    auto tris = world_triangles(m);
    Bytes b;
    std::string a;
    if (format) put_str(b, std::string("AutoCAD Binary DXF\r\n\x1a\0", 22));
    auto code = [&](int c) {
        if (format == 0) { a += std::to_string(c) + "\n"; return; }
        if (format == 1) { put_u16le(b, std::uint16_t(c)); return; }
        if (c < 255) b.push_back(std::uint8_t(c));
        else { b.push_back(255); put_u16le(b, std::uint16_t(c)); }
    };
    auto str = [&](int c, const std::string& v) {
        code(c);
        if (format == 0) a += v + "\n";
        else { put_str(b, v); b.push_back(0); }
    };
    auto dbl = [&](int c, double v) {
        code(c);
        if (format == 0) a += num(v) + "\n";
        else { std::uint64_t u; std::memcpy(&u, &v, 8); for (int i = 0; i < 8; ++i) b.push_back(std::uint8_t(u >> (8 * i))); }
    };
    auto i32 = [&](int c, std::int32_t v) {
        code(c);
        if (format == 0) a += std::to_string(v) + "\n";
        else put_u32le(b, std::uint32_t(v));
    };
    str(0, "SECTION"); str(2, "HEADER");
    str(9, "$ACADVER"); str(1, format == 2 ? "AC1009" : "AC1018");
    str(0, "ENDSEC");
    str(0, "SECTION"); str(2, "ENTITIES");
    for (const auto& t : tris) {
        str(0, "3DFACE");
        str(8, "0");
        for (int k = 0; k < 4; ++k) {
            const auto& p = t.p[k < 3 ? k : 2];
            dbl(10 + k, p[0]); dbl(20 + k, p[1]); dbl(30 + k, p[2]);
        }
        i32(420, std::int32_t(std::uint32_t(t.c.r) << 16 | std::uint32_t(t.c.g) << 8 | t.c.b));
        if (t.c.a != 255) i32(440, std::int32_t(0x02000000u | t.c.a));
    }
    str(0, "ENDSEC");
    str(0, "EOF");
    if (format == 0) return {{m.name + ".dxf", Bytes(a.begin(), a.end())}};
    return {{m.name + ".dxf", b}};
}

const std::vector<Variant>& variants() {
    static const std::vector<Variant> v = {
        {"stl-binary", write_stl_binary},
        {"stl-ascii", write_stl_ascii},
        {"obj", write_obj},
        {"ply-ascii", [](const Model& m) { return write_ply(m, 0); }},
        {"ply-le", [](const Model& m) { return write_ply(m, 1); }},
        {"ply-be", [](const Model& m) { return write_ply(m, 2); }},
        {"off-ascii", [](const Model& m) { return write_off(m, false); }},
        {"off-binary", [](const Model& m) { return write_off(m, true); }},
        {"gltf", [](const Model& m) { return write_gltf(m, false); }},
        {"glb", [](const Model& m) { return write_gltf(m, true); }},
        {"3mf", write_3mf},
        {"dxf-ascii", [](const Model& m) { return write_dxf(m, 0); }},
        {"dxf-binary", [](const Model& m) { return write_dxf(m, 1); }},
        {"dxf-binary-r12", [](const Model& m) { return write_dxf(m, 2); }},
        {"catpart-table", [](const Model& m) { return write_catpart(m, {0, false, false, false}); }},
        {"catpart-raw", [](const Model& m) { return write_catpart(m, {1, false, false, false}); }},
        {"catpart-tex", [](const Model& m) { return write_catpart(m, {2, false, false, false}); }},
        {"catpart-u16-lod", [](const Model& m) { return write_catpart(m, {0, true, true, false}); }},
        {"catpart-r17-counts", [](const Model& m) { return write_catpart(m, {0, true, false, true}); }},
    };
    return v;
}

}  // namespace cadcore::testgen
