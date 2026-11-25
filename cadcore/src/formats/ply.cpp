// PLY — ascii, binary_little_endian, binary_big_endian.
//
// Header: "ply", "format <kind> 1.0", then elements with scalar or list
// properties of type (u)int8/16/32, float32/64 (and their alias names).
// Used: vertex.{x,y,z,red,green,blue,alpha} (u8 or float colours),
// face.{vertex_indices|vertex_index} (polygons, fan-triangulated),
// face.{red,green,blue,alpha}, edge.{vertex1,vertex2}. Every other element and
// property is skipped by size, so unknown content never breaks decoding.
#include <array>
#include <cstring>
#include <string>
#include <vector>

#include <algorithm>
#include <functional>

#include "core/parallel.hpp"
#include "core/text.hpp"
#include "formats/registry.hpp"

namespace cadcore::formats {
namespace {

enum class T : std::uint8_t { I8, U8, I16, U16, I32, U32, F32, F64 };

bool type_of(std::string_view s, T& t) {
    static const std::pair<std::string_view, T> names[] = {
        {"char", T::I8}, {"int8", T::I8}, {"uchar", T::U8}, {"uint8", T::U8}, {"short", T::I16}, {"int16", T::I16},
        {"ushort", T::U16}, {"uint16", T::U16}, {"int", T::I32}, {"int32", T::I32}, {"uint", T::U32}, {"uint32", T::U32},
        {"float", T::F32}, {"float32", T::F32}, {"double", T::F64}, {"float64", T::F64}};
    for (auto& [n, v] : names)
        if (n == s) { t = v; return true; }
    return false;
}
std::size_t size_of(T t) { return t == T::I8 || t == T::U8 ? 1 : t == T::I16 || t == T::U16 ? 2 : t == T::F64 ? 8 : 4; }

struct Property {
    std::string name;
    bool list = false;
    T count_type = T::U8;
    T type = T::F32;
};
struct Element {
    std::string name;
    std::uint64_t count = 0;
    std::vector<Property> props;
};
enum class Enc { Ascii, LE, BE };

int sniff(ByteView d, std::string_view ext) {
    if (d.starts_with("ply\n") || d.starts_with("ply\r\n")) return 2;
    return ext == "ply" ? 1 : 0;
}

// Reads one value of type t, as double, in the file's encoding.
class ValueReader {
public:
    ValueReader(ByteView d, std::size_t pos, Enc e) : r_(d, pos), c_(d.text().substr(pos), pos), enc_(e) {}
    double read(T t) {
        if (enc_ == Enc::Ascii) {
            c_.skip_ws_and_newlines();
            while (text::iequals_prefix(c_.token_peek(), "comment")) {
                c_.next_line();
                c_.skip_ws_and_newlines();
            }
            return c_.number<double>("PLY value");
        }
        // one bounds check per value, then a plain load (byte-swapped for big-endian files)
        const std::size_t n = size_of(t);
        std::uint8_t b[8];
        std::memcpy(b, r_.bytes(n), n);
        if (enc_ == Enc::BE) std::reverse(b, b + n);  // little-endian hosts only (x86, wasm)
        switch (t) {
            case T::I8: return double(std::int8_t(b[0]));
            case T::U8: return double(b[0]);
            case T::I16: { std::int16_t v; std::memcpy(&v, b, 2); return double(v); }
            case T::U16: { std::uint16_t v; std::memcpy(&v, b, 2); return double(v); }
            case T::I32: { std::int32_t v; std::memcpy(&v, b, 4); return double(v); }
            case T::U32: { std::uint32_t v; std::memcpy(&v, b, 4); return double(v); }
            case T::F32: { float v; std::memcpy(&v, b, 4); return double(v); }
            case T::F64: { double v; std::memcpy(&v, b, 8); return v; }
        }
        return 0;
    }
    void skip(T t, std::uint64_t n) {
        if (enc_ == Enc::Ascii) {
            for (std::uint64_t i = 0; i < n; ++i) read(t);
        } else {
            r_.skip(std::size_t(size_of(t) * n));
        }
    }
    void end_of_record() {
        if (enc_ == Enc::Ascii) c_.next_line();
    }
    std::size_t offset() const { return enc_ == Enc::Ascii ? c_.offset() : r_.pos(); }

private:
    Reader r_;
    text::Cursor c_;
    Enc enc_;
};

std::uint8_t to_u8(double v, T t) {
    if (t == T::F32 || t == T::F64) return std::uint8_t(v <= 0 ? 0 : v >= 1 ? 255 : v * 255 + 0.5);
    return std::uint8_t(v < 0 ? 0 : v > 255 ? 255 : v);
}

struct Layout {
    bool is_vertex = false, is_face = false, is_edge = false, vcolor = false, fcolor = false;
    int ix = -1, iy = -1, iz = -1, ir = -1, ig = -1, ib = -1, ia = -1, iidx = -1, ie1 = -1, ie2 = -1;
};

Layout layout(const Element& e) {
    Layout L;
    L.is_vertex = e.name == "vertex";
    L.is_face = e.name == "face";
    L.is_edge = e.name == "edge";
    for (int k = 0; k < int(e.props.size()); ++k) {
        const auto& n = e.props[std::size_t(k)].name;
        if (n == "x") L.ix = k; else if (n == "y") L.iy = k; else if (n == "z") L.iz = k;
        else if (n == "red" || n == "r" || n == "diffuse_red") L.ir = k;
        else if (n == "green" || n == "g" || n == "diffuse_green") L.ig = k;
        else if (n == "blue" || n == "b" || n == "diffuse_blue") L.ib = k;
        else if (n == "alpha") L.ia = k;
        else if (n == "vertex_indices" || n == "vertex_index") L.iidx = k;
        else if (n == "vertex1") L.ie1 = k; else if (n == "vertex2") L.ie2 = k;
    }
    if (L.is_vertex && (L.ix < 0 || L.iy < 0 || L.iz < 0)) fail(ErrorCode::Malformed, "vertex element without x/y/z");
    L.vcolor = L.is_vertex && L.ir >= 0 && L.ig >= 0 && L.ib >= 0;
    L.fcolor = L.is_face && L.ir >= 0 && L.ig >= 0 && L.ib >= 0;
    return L;
}

struct FaceRun { std::uint32_t triangles; Color color; };

// What a range of records contributes, in record order.
struct Out {
    // In place: when set, positions / colours / indices go straight into the mesh arrays here.
    float* pos_at = nullptr;
    std::uint8_t* col_at = nullptr;
    std::uint32_t* idx_at = nullptr;
    std::vector<float> pos;
    std::vector<std::uint8_t> col;
    std::vector<std::uint32_t> idx;
    std::vector<FaceRun> faces;  // faces with a colour
    std::vector<float> edges;    // vertex index pairs
};

void parse_records(const Element& e, const Layout& L, ValueReader& vr, std::uint64_t count, Out& o) {
    std::vector<std::uint32_t> poly;
    if (L.is_vertex && !o.pos_at) o.pos.reserve(std::size_t(count) * 3);
    if (L.is_face && !o.idx_at) o.idx.reserve(std::size_t(count) * 3);
    std::size_t tris = 0;  // triangles written in place
    for (std::uint64_t i = 0; i < count; ++i) {
        double vals[16] = {};
        const std::size_t before = o.idx_at ? tris * 3 : o.idx.size();
        double e1 = 0, e2 = 0;
        for (int k = 0; k < int(e.props.size()); ++k) {
            const Property& p = e.props[std::size_t(k)];
            if (p.list) {
                const double nd = vr.read(p.count_type);
                if (nd < 0 || nd > 1e7) fail(ErrorCode::Malformed, "implausible list length", vr.offset());
                const auto n = std::uint64_t(nd);
                if (L.is_face && k == L.iidx) {
                    poly.resize(n);
                    for (std::uint64_t j = 0; j < n; ++j) {
                        const double v = vr.read(p.type);
                        if (v < 0) fail(ErrorCode::Malformed, "negative face index", vr.offset());
                        poly[j] = std::uint32_t(v);
                    }
                    for (std::uint64_t j = 1; j + 1 < n; ++j) {
                        if (o.idx_at) {
                            std::uint32_t* t = o.idx_at + 3 * tris++;
                            t[0] = poly[0]; t[1] = poly[j]; t[2] = poly[j + 1];
                        } else {
                            o.idx.insert(o.idx.end(), {poly[0], poly[j], poly[j + 1]});
                        }
                    }
                } else {
                    vr.skip(p.type, n);
                }
            } else {
                const double v = vr.read(p.type);
                if (k < 16) vals[k] = v;
                if (k == L.ie1) e1 = v;
                if (k == L.ie2) e2 = v;
            }
        }
        vr.end_of_record();
        if (L.is_vertex) {
            const float xyz[3] = {float(vals[L.ix]), float(vals[L.iy]), float(vals[L.iz])};
            if (o.pos_at) std::memcpy(o.pos_at + 3 * i, xyz, 12);
            else o.pos.insert(o.pos.end(), xyz, xyz + 3);
            if (L.vcolor && o.col_at) {
                std::uint8_t* c = o.col_at + 4 * i;
                c[0] = to_u8(vals[L.ir], e.props[std::size_t(L.ir)].type);
                c[1] = to_u8(vals[L.ig], e.props[std::size_t(L.ig)].type);
                c[2] = to_u8(vals[L.ib], e.props[std::size_t(L.ib)].type);
                c[3] = L.ia >= 0 ? to_u8(vals[L.ia], e.props[std::size_t(L.ia)].type) : std::uint8_t(255);
            } else if (L.vcolor)
                o.col.insert(o.col.end(), {to_u8(vals[L.ir], e.props[std::size_t(L.ir)].type), to_u8(vals[L.ig], e.props[std::size_t(L.ig)].type),
                                           to_u8(vals[L.ib], e.props[std::size_t(L.ib)].type),
                                           L.ia >= 0 ? to_u8(vals[L.ia], e.props[std::size_t(L.ia)].type) : std::uint8_t(255)});
        } else if (L.fcolor) {
            const Color fc{to_u8(vals[L.ir], e.props[std::size_t(L.ir)].type), to_u8(vals[L.ig], e.props[std::size_t(L.ig)].type),
                           to_u8(vals[L.ib], e.props[std::size_t(L.ib)].type), L.ia >= 0 ? to_u8(vals[L.ia], e.props[std::size_t(L.ia)].type) : std::uint8_t(255)};
            o.faces.push_back({std::uint32_t(((o.idx_at ? tris * 3 : o.idx.size()) - before) / 3), fc});
        } else if (L.is_edge && L.ie1 >= 0 && L.ie2 >= 0) {
            o.edges.insert(o.edges.end(), {float(e1), float(e2)});
        }
    }
}


Scene decode(ByteView d, const DecodeOptions& options) {
    const std::size_t hdr_end_kw = d.find("end_header", 0, std::min<std::size_t>(d.size, 1 << 20));
    if (!d.starts_with("ply") || hdr_end_kw == ByteView::npos) fail(ErrorCode::Malformed, "missing PLY header");
    std::size_t body = d.find("\n", hdr_end_kw);
    if (body == ByteView::npos) fail(ErrorCode::Truncated, "header without body");
    ++body;

    Enc enc = Enc::Ascii;
    std::vector<Element> elements;
    text::Cursor h(d.text().substr(0, hdr_end_kw));
    h.next_line();
    while (!h.eof()) {
        std::string_view kw = h.token();
        if (kw == "format") {
            std::string_view f = h.token();
            enc = f == "ascii" ? Enc::Ascii : f == "binary_little_endian" ? Enc::LE : f == "binary_big_endian" ? Enc::BE
                                                                                         : (fail(ErrorCode::Unsupported, "PLY format " + std::string(f)), Enc::Ascii);
        } else if (kw == "element") {
            Element e;
            e.name = std::string(h.token());
            e.count = std::uint64_t(h.integer("element count"));
            elements.push_back(std::move(e));
        } else if (kw == "property") {
            if (elements.empty()) fail(ErrorCode::Malformed, "property before element", h.offset());
            Property p;
            std::string_view t = h.token();
            if (t == "list") {
                p.list = true;
                if (!type_of(h.token(), p.count_type)) fail(ErrorCode::Malformed, "bad list count type", h.offset());
                t = h.token();
            }
            if (!type_of(t, p.type)) fail(ErrorCode::Malformed, "bad property type '" + std::string(t) + "'", h.offset());
            p.name = std::string(h.token());
            elements.back().props.push_back(std::move(p));
        }
        h.next_line();
    }

    Scene scene;
    scene.format = enc == Enc::Ascii ? "ply-ascii" : enc == Enc::LE ? "ply-binary-le" : "ply-binary-be";
    Mesh mesh;
    std::vector<float> all_pos;
    std::vector<FaceRun> faces;
    const unsigned threads = options.threads;
    std::size_t pos = body;          // start of the current element's records
    std::vector<std::size_t> lines;  // ASCII: start of every body line (built on first use)
    for (const Element& e : elements) {
        const Layout L = layout(e);
        if (L.is_vertex) {
            if (e.count > d.size) fail(ErrorCode::Malformed, "vertex count exceeds file size");
            mesh.positions.reserve(std::size_t(e.count) * 3);
        }
        const std::uint64_t n = e.count;
        const std::size_t parts = n < 8192 ? 1 : thread_count(threads);
        std::vector<std::size_t> cut;   // record index at each part boundary
        std::vector<std::size_t> at;    // byte offset of each boundary
        std::vector<std::size_t> tri_at; // binary faces: triangles before each boundary
        for (std::size_t k = 0; k <= parts; ++k) cut.push_back(std::size_t(n * k / parts));
        std::vector<Slot<Out>> out(parts);
        std::size_t end = 0;
        bool done = false;
        if (parts > 1 && enc != Enc::Ascii) {
            // record offsets at the boundaries: fixed size, or a pass over the list counts only
            bool fixed = true;
            std::size_t rec = 0;
            for (const auto& p : e.props) { fixed &= !p.list; rec += size_of(p.type); }
            if (fixed) {
                if (std::uint64_t(rec) * n > d.size - pos) fail(ErrorCode::Truncated, "unexpected end of data", d.size);
                for (auto c : cut) at.push_back(pos + c * rec);
            } else {
                // record r: the fixed-size properties plus, per list, its count and count * item size
                const bool le = enc == Enc::LE;
                auto count_at = [&](std::size_t o, T t) -> std::uint64_t {
                    const std::size_t sz = size_of(t);
                    if (o + sz > d.size) fail(ErrorCode::Truncated, "unexpected end of data", o);
                    std::uint64_t v = 0;
                    for (std::size_t i = 0; i < sz; ++i) v |= std::uint64_t(d[o + (le ? i : sz - 1 - i)]) << (8 * i);
                    double x;  // as ValueReader reads it: signed and float counts included
                    switch (t) {
                        case T::I8: x = double(std::int8_t(v)); break;
                        case T::I16: x = double(std::int16_t(v)); break;
                        case T::I32: x = double(std::int32_t(v)); break;
                        case T::F32: { std::uint32_t u = std::uint32_t(v); float f; std::memcpy(&f, &u, 4); x = f; break; }
                        case T::F64: std::memcpy(&x, &v, 8); break;
                        default: x = double(v); break;
                    }
                    if (!(x >= 0 && x <= 1e7)) fail(ErrorCode::Malformed, "implausible list length", o + sz);
                    return std::uint64_t(x);
                };
                std::size_t o = pos, r = 0, t = 0;
                for (std::size_t k = 0; k <= parts; ++k) {
                    const std::size_t stop = k < parts ? cut[k] : std::size_t(n);
                    for (; r < stop; ++r)
                        for (std::size_t q = 0; q < e.props.size(); ++q) {
                            const auto& p = e.props[q];
                            if (!p.list) { o += size_of(p.type); continue; }
                            const std::uint64_t c = count_at(o, p.count_type);
                            if (L.is_face && int(q) == L.iidx && c >= 3) t += std::size_t(c) - 2;
                            o += size_of(p.count_type) + std::size_t(c) * size_of(p.type);
                        }
                    tri_at.push_back(t);
                    if (o > d.size) fail(ErrorCode::Truncated, "unexpected end of data", d.size);
                    at.push_back(o);
                }
            }
            // write in place: vertex records are one position (and colour) each; face triangles start
            // where the skim pass counted them
            const std::size_t pos0 = mesh.positions.size(), col0 = mesh.colors.size(), idx0 = mesh.indices.size();
            if (L.is_vertex) {
                mesh.positions.resize(pos0 + std::size_t(n) * 3);
                if (L.vcolor) mesh.colors.resize(col0 + std::size_t(n) * 4);
            }
            if (L.is_face && !tri_at.empty()) mesh.indices.resize(idx0 + tri_at.back() * 3);
            std::vector<std::function<void()>> tasks;
            for (std::size_t k = 0; k < parts; ++k)
                tasks.push_back([&, k] {
                    Out& o = out[k].value;
                    if (L.is_vertex) {
                        o.pos_at = mesh.positions.data() + pos0 + cut[k] * 3;
                        if (L.vcolor) o.col_at = mesh.colors.data() + col0 + cut[k] * 4;
                    }
                    if (L.is_face && !tri_at.empty()) o.idx_at = mesh.indices.data() + idx0 + tri_at[k] * 3;
                    ValueReader vr(d, at[k], enc);
                    parse_records(e, L, vr, cut[k + 1] - cut[k], o);
                });
            run_all(tasks, threads);
            end = at.back();
            done = true;
        } else if (parts > 1) {
            // ASCII: one record per line (as the format specifies) lets every part start at its line;
            // each part must end exactly at the next part's line, else (records across lines,
            // comments or blank lines in the body) the element is parsed sequentially
            if (lines.empty()) lines = text::line_starts(d.text(), body, threads);
            const auto first = std::size_t(std::lower_bound(lines.begin(), lines.end(), pos) - lines.begin());
            if (first < lines.size() && lines[first] == pos && first + n <= lines.size()) {
                auto line_at = [&](std::size_t r) { return first + r < lines.size() ? lines[first + r] : d.size; };
                std::vector<Slot<bool>> ok(parts);
                std::vector<std::function<void()>> tasks;
                for (std::size_t k = 0; k < parts; ++k)
                    tasks.push_back([&, k] {
                        try {
                            ValueReader vr(d, line_at(cut[k]), enc);
                            parse_records(e, L, vr, cut[k + 1] - cut[k], out[k].value);
                            ok[k].value = vr.offset() == line_at(cut[k + 1]);
                        } catch (const DecodeError&) {
                            ok[k].value = false;  // the sequential parse decides
                        }
                    });
                run_all(tasks, threads);
                done = std::all_of(ok.begin(), ok.end(), [](const Slot<bool>& x) { return x.value; });
                end = line_at(n);
            }
        }
        if (!done) {  // sequential
            out.assign(1, {});
            ValueReader vr(d, pos, enc);
            parse_records(e, L, vr, n, out[0].value);
            end = vr.offset();
        }
        pos = end;
        const std::size_t np = out.size();
        append_all(mesh.positions, np, [&](std::size_t k) -> auto& { return out[k].value.pos; }, threads);
        append_all(mesh.colors, np, [&](std::size_t k) -> auto& { return out[k].value.col; }, threads);
        append_all(mesh.indices, np, [&](std::size_t k) -> auto& { return out[k].value.idx; }, threads);
        append_all(faces, np, [&](std::size_t k) -> auto& { return out[k].value.faces; }, threads);
        append_all(all_pos, np, [&](std::size_t k) -> auto& { return out[k].value.edges; }, threads);
    }
    // face colour runs, in face order
    std::uint32_t tri = 0;
    for (const FaceRun& f : faces) {
        if (!mesh.groups.empty() && mesh.groups.back().color == f.color && mesh.groups.back().first + mesh.groups.back().count == tri)
            mesh.groups.back().count += f.triangles;
        else
            mesh.groups.push_back({tri, f.triangles, f.color});
        tri += f.triangles;
    }
    const std::size_t nv = mesh.vertex_count();
    for (auto i : mesh.indices)
        if (i >= nv) fail(ErrorCode::Malformed, "face index " + std::to_string(i) + " >= vertex count " + std::to_string(nv));
    for (std::size_t k = 0; k + 1 < all_pos.size(); k += 2) {
        auto a = std::size_t(all_pos[k]), b = std::size_t(all_pos[k + 1]);
        if (a >= nv || b >= nv) fail(ErrorCode::Malformed, "edge index out of range");
        Polyline pl;
        pl.points = {mesh.positions[a * 3], mesh.positions[a * 3 + 1], mesh.positions[a * 3 + 2],
                     mesh.positions[b * 3], mesh.positions[b * 3 + 1], mesh.positions[b * 3 + 2]};
        scene.polylines.push_back(std::move(pl));
    }
    if (mesh.groups.size() == 1) {
        mesh.color = mesh.groups[0].color;
        mesh.groups.clear();
    }
    if (mesh.indices.empty() && nv > 0) mesh.topology = Topology::Points;
    scene.meshes.push_back(std::move(mesh));
    return scene;
}

}  // namespace

const Decoder& ply() {
    static const SimpleDecoder d({"ply", "Polygon File Format (Stanford)", {"ply"}}, sniff, decode);
    return d;
}

}  // namespace cadcore::formats
