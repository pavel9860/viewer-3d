#include "codec/draco/attributes.hpp"

#include "codec/draco/kd_tree.hpp"

#include <algorithm>
#include <climits>
#include <cstring>
#include <cmath>
#include <memory>

namespace cadcore::codec::draco {
namespace {

// ---- enums ---------------------------------------------------------------------------------------------

constexpr int MESH_VERTEX_ATTRIBUTE = 0;
constexpr int TRAVERSAL_DEPTH_FIRST = 0, TRAVERSAL_PREDICTION_DEGREE = 1;
constexpr int SEQ_GENERIC = 0, SEQ_INTEGER = 1, SEQ_QUANTIZATION = 2, SEQ_NORMALS = 3;
constexpr int PREDICTION_NONE = -2, PREDICTION_DIFFERENCE = 0, PARALLELOGRAM = 1, MULTI_PARALLELOGRAM = 2,
              TEX_COORDS_DEPRECATED = 3, CONSTRAINED_MULTI_PARALLELOGRAM = 4, TEX_COORDS_PORTABLE = 5,
              GEOMETRIC_NORMAL = 6, NUM_PREDICTION_SCHEMES = 7;
constexpr int TRANSFORM_NONE = -1, TRANSFORM_WRAP = 1, TRANSFORM_OCTAHEDRON = 2, TRANSFORM_OCTAHEDRON_CANONICALIZED = 3,
              NUM_TRANSFORMS = 4;
constexpr int ATT_POSITION = 0;
constexpr int kNumNamedAttributeTypes = 9;  // POSITION .. WEIGHTS

int data_type_length(int dt) {
    static constexpr int kLen[12] = {0, 1, 1, 2, 2, 4, 4, 8, 8, 4, 8, 1};
    return dt >= 1 && dt <= 11 ? kLen[dt] : 0;
}

std::int32_t add_u(std::int32_t a, std::int32_t b) { return std::int32_t(std::uint32_t(a) + std::uint32_t(b)); }
std::int32_t sub_u(std::int32_t a, std::int32_t b) { return std::int32_t(std::uint32_t(a) - std::uint32_t(b)); }

int most_significant_bit(std::uint32_t n) {
    int msb = -1;
    while (n) { ++msb; n >>= 1; }
    return msb;
}

// ---- octahedral normals -----------------------------------------------------------------------------------

struct Octahedron {
    std::int32_t bits = 0, max_quantized = 0, max_value = 0, center = 0;
    float dequant_scale = 0;

    void set_bits(std::int32_t q) {
        if (q < 2 || q > 30) bad("normal quantization bits out of range");
        bits = q;
        max_quantized = std::int32_t((1u << q) - 1);
        max_value = max_quantized - 1;
        dequant_scale = 2.f / float(max_value);
        center = max_value / 2;
    }
    void set_max_quantized(std::int32_t mq) {
        if (mq % 2 == 0) bad("even maximum quantized normal value");
        set_bits(most_significant_bit(std::uint32_t(mq)) + 1);
    }
    bool in_diamond(std::int32_t s, std::int32_t t) const { return std::uint32_t(std::abs(s)) + std::uint32_t(std::abs(t)) <= std::uint32_t(center); }
    void invert_diamond(std::int32_t& s, std::int32_t& t) const {
        std::int32_t ss, st;
        if (s >= 0 && t >= 0) { ss = 1; st = 1; }
        else if (s <= 0 && t <= 0) { ss = -1; st = -1; }
        else { ss = s > 0 ? 1 : -1; st = t > 0 ? 1 : -1; }
        const std::uint32_t cs = std::uint32_t(ss * center), ct = std::uint32_t(st * center);
        std::uint32_t us = std::uint32_t(s), ut = std::uint32_t(t);
        us = us + us - cs;
        ut = ut + ut - ct;
        if (ss * st >= 0) { const std::uint32_t tmp = us; us = 0u - ut; ut = 0u - tmp; }
        else std::swap(us, ut);
        us += cs;
        ut += ct;
        s = std::int32_t(us) / 2;
        t = std::int32_t(ut) / 2;
    }
    std::int32_t mod_max(std::int32_t x) const {
        if (x > center) return x - max_quantized;
        if (x < -center) return x + max_quantized;
        return x;
    }
    void canonicalize_coords(std::int32_t& s, std::int32_t& t) const {
        if ((s == 0 && t == 0) || (s == 0 && t == max_value) || (s == max_value && t == 0)) { s = max_value; t = max_value; }
        else if (s == 0 && t > center) t = center - (t - center);
        else if (s == max_value && t < center) t = center + (center - t);
        else if (t == max_value && s < center) s = center + (center - s);
        else if (t == 0 && s > center) s = center - (s - center);
    }
    void vector_to_coords(const std::int32_t* v, std::int32_t& s, std::int32_t& t) const {
        if (v[0] >= 0) {
            s = v[1] + center;
            t = v[2] + center;
        } else {
            s = v[1] < 0 ? std::abs(v[2]) : max_value - std::abs(v[2]);
            t = v[2] < 0 ? std::abs(v[1]) : max_value - std::abs(v[1]);
        }
        canonicalize_coords(s, t);
    }
    void canonicalize_vector(std::int32_t* v) const {
        const std::int64_t abs_sum = std::int64_t(std::abs(std::int64_t(v[0]))) + std::abs(std::int64_t(v[1])) + std::abs(std::int64_t(v[2]));
        if (abs_sum == 0) {
            v[0] = center;
        } else {
            v[0] = std::int32_t(std::int64_t(v[0]) * center / abs_sum);
            v[1] = std::int32_t(std::int64_t(v[1]) * center / abs_sum);
            const std::int32_t r = center - std::abs(v[0]) - std::abs(v[1]);
            v[2] = v[2] >= 0 ? r : -r;
        }
    }
    void to_unit_vector(std::int32_t qs, std::int32_t qt, float* out) const {
        float y = float(qs) * dequant_scale - 1.f, z = float(qt) * dequant_scale - 1.f;
        const float x = 1.f - std::abs(y) - std::abs(z);
        float x_offset = -x;
        x_offset = x_offset < 0 ? 0 : x_offset;
        y += y < 0 ? x_offset : -x_offset;
        z += z < 0 ? x_offset : -x_offset;
        const float n2 = x * x + y * y + z * z;
        if (n2 < 1e-6f) {
            out[0] = out[1] = out[2] = 0;
        } else {
            const float d = 1.0f / std::sqrt(n2);
            out[0] = x * d;
            out[1] = y * d;
            out[2] = z * d;
        }
    }
};

// ---- prediction transforms ----------------------------------------------------------------------------------

struct Transform {
    int type = TRANSFORM_NONE;
    std::int32_t min_value = 0, max_value = 0, max_dif = 0;  // wrap
    Octahedron oct;                                          // octahedron transforms

    bool corrections_positive() const { return type == TRANSFORM_OCTAHEDRON || type == TRANSFORM_OCTAHEDRON_CANONICALIZED; }

    void decode(Stream& s) {
        if (type == TRANSFORM_WRAP) {
            min_value = s.i32();
            max_value = s.i32();
            if (min_value > max_value) bad("wrap transform minimum exceeds maximum");
            const std::int64_t dif = std::int64_t(max_value) - min_value;
            if (dif >= INT32_MAX) bad("wrap transform range too large");
            max_dif = std::int32_t(1 + dif);
        } else if (type == TRANSFORM_OCTAHEDRON_CANONICALIZED) {
            const std::int32_t mq = s.i32();
            s.i32();  // center value (derived from mq)
            oct.set_max_quantized(mq);
        } else if (type == TRANSFORM_OCTAHEDRON) {
            const std::int32_t mq = s.i32();
            if (s.ver < version(2, 2)) s.i32();
            oct.set_max_quantized(mq);
        }
    }

    void original(const std::int32_t* pred, const std::int32_t* corr, std::int32_t* out, int nc) const {
        if (type == TRANSFORM_WRAP) {
            for (int i = 0; i < nc; ++i) {
                const std::int32_t p = std::clamp(pred[i], min_value, max_value);
                std::int64_t v = std::int64_t(p) + corr[i];
                if (v > max_value) v -= max_dif;
                else if (v < min_value) v += max_dif;
                out[i] = std::int32_t(v);
            }
        } else if (type == TRANSFORM_OCTAHEDRON_CANONICALIZED) {
            std::int32_t px = sub_u(pred[0], oct.center), py = sub_u(pred[1], oct.center);
            const bool in_diamond = oct.in_diamond(px, py);
            if (!in_diamond) oct.invert_diamond(px, py);
            const bool bottom_left = (px == 0 && py == 0) || (px < 0 && py <= 0);
            int rot;
            if (px == 0) rot = py == 0 ? 0 : (py > 0 ? 3 : 1);
            else if (px > 0) rot = py >= 0 ? 2 : 1;
            else rot = py <= 0 ? 0 : 3;
            auto rotate = [](std::int32_t& x, std::int32_t& y, int r) {
                const std::int32_t a = x, b = y;
                switch (r) {
                    case 1: x = b; y = -a; break;
                    case 2: x = -a; y = -b; break;
                    case 3: x = -b; y = a; break;
                    default: break;
                }
            };
            if (!bottom_left) rotate(px, py, rot);
            std::int32_t ox = oct.mod_max(add_u(px, corr[0])), oy = oct.mod_max(add_u(py, corr[1]));
            if (!bottom_left) rotate(ox, oy, (4 - rot) % 4);
            if (!in_diamond) oct.invert_diamond(ox, oy);
            out[0] = add_u(ox, oct.center);
            out[1] = add_u(oy, oct.center);
        } else if (type == TRANSFORM_OCTAHEDRON) {
            std::int32_t px = sub_u(pred[0], oct.center), py = sub_u(pred[1], oct.center);
            const bool in_diamond = oct.in_diamond(px, py);
            if (!in_diamond) oct.invert_diamond(px, py);
            std::int32_t ox = oct.mod_max(add_u(px, corr[0])), oy = oct.mod_max(add_u(py, corr[1]));
            if (!in_diamond) oct.invert_diamond(ox, oy);
            out[0] = add_u(ox, oct.center);
            out[1] = add_u(oy, oct.center);
        } else {
            bad("prediction without a transform");
        }
    }
};

// ---- per-attribute and per-decoder state ----------------------------------------------------------------------

struct Attr {
    int att_type = 0, data_type = 0, num_components = 0;
    bool normalized = false;
    std::uint32_t unique_id = 0;
    int seq_type = 0;
    std::vector<std::uint32_t> point_to_entry;  // per point
    bool has_portable = false;
    int portable_components = 0;
    std::vector<std::int32_t> portable;  // integer (portable) values per entry
    bool has_values = false;
    std::vector<double> values;          // final values per entry (num_components each, 3 for normals)
    // transforms
    std::vector<float> min_values;
    float range = 0;
    int quantization_bits = 0;
};

struct Dec {
    int data_id = -1, decoder_type = MESH_VERTEX_ATTRIBUTE, traversal = TRAVERSAL_DEPTH_FIRST;
    std::vector<int> attrs;            // attribute ids
    std::vector<std::uint32_t> point_ids;  // sequence: entry -> point
    IVec data_to_corner, vertex_to_data;
    const AttributeCornerTable* att_table = nullptr;  // nullptr: position table
};

class AttributeDecoder {
public:
    AttributeDecoder(Stream& s, Connectivity& c, bool is_mesh, bool kd_tree) : s_(s), c_(c), mesh_(is_mesh), kd_(kd_tree) {}
    std::vector<Attribute> run();

private:
    // sequencing
    void generate_sequence(Dec& d);
    template <class T> void traverse(Dec& d, const T& table);
    template <class T> void depth_first(Dec& d, const T& t, std::int32_t corner, std::vector<bool>& face_vis, std::vector<bool>& vert_vis);
    template <class T> void prediction_degree(Dec& d, const T& t, std::int32_t corner, std::vector<bool>& face_vis, std::vector<bool>& vert_vis,
                                              std::vector<int>& degree);
    // values
    void decode_kd_tree(Dec& d);
    void decode_portable(Dec& d, Attr& a);
    void decode_integer(Dec& d, Attr& a, int method, int transform_type);
    void decode_transform_data(Attr& a);
    void transform_to_original(Attr& a);
    // prediction
    const Attr& parent_position(const Attr& self) const;
    template <class T>
    void predict(const Dec& d, const T& table, const Attr& a, int method, Transform& tr, std::vector<std::int32_t>& data, int nc);

    Stream& s_;
    Connectivity& c_;
    bool mesh_, kd_;
    std::vector<Attr> attrs_;
    std::vector<Dec> decs_;
    // prediction-scheme side data of the attribute being decoded
    std::vector<std::vector<bool>> crease_;
    std::vector<bool> orientations_;
    BitDecoder flips_;
    bool one_triangle_ = false;
};

// ---- traversal ---------------------------------------------------------------------------------------------------

template <class T>
void AttributeDecoder::depth_first(Dec& d, const T& t, std::int32_t corner, std::vector<bool>& face_vis, std::vector<bool>& vert_vis) {
    auto face_visited = [&](std::int32_t c) { return c == kInvalid || face_vis[std::size_t(c / 3)]; };
    auto visit = [&](std::int32_t v, std::int32_t c) {
        if (v < 0) bad("traversal reached an unmapped corner");
        vert_vis[std::size_t(v)] = true;
        d.point_ids.push_back(c_.faces[std::size_t(c)]);
        d.data_to_corner.push_back(c);
        d.vertex_to_data[std::size_t(v)] = std::int32_t(d.data_to_corner.size()) - 1;
    };
    if (face_visited(corner)) return;
    std::vector<std::int32_t> stack{corner};
    const std::int32_t nv = t.Vertex(t.Next(corner)), pv = t.Vertex(t.Previous(corner));
    if (nv < 0 || pv < 0) bad("traversal reached an unmapped corner");
    if (!vert_vis[std::size_t(nv)]) visit(nv, t.Next(corner));
    if (!vert_vis[std::size_t(pv)]) visit(pv, t.Previous(corner));
    while (!stack.empty()) {
        corner = stack.back();
        if (corner == kInvalid || face_vis[std::size_t(corner / 3)]) { stack.pop_back(); continue; }
        while (true) {
            face_vis[std::size_t(corner / 3)] = true;
            const std::int32_t v = t.Vertex(corner);
            if (v < 0) bad("traversal reached an unmapped corner");
            if (!vert_vis[std::size_t(v)]) {
                const bool on_boundary = t.IsOnBoundary(v);
                visit(v, corner);
                if (!on_boundary) {
                    corner = t.GetRightCorner(corner);
                    if (corner == kInvalid) bad("interior vertex without a right face");
                    continue;
                }
            }
            const std::int32_t rc = t.GetRightCorner(corner), lc = t.GetLeftCorner(corner);
            const bool rv = face_visited(rc), lv = face_visited(lc);
            if (rv) {
                if (lv) { stack.pop_back(); break; }
                corner = lc;
            } else if (lv) {
                corner = rc;
            } else {
                stack.back() = lc;
                stack.push_back(rc);
                break;
            }
        }
    }
}

template <class T>
void AttributeDecoder::prediction_degree(Dec& d, const T& t, std::int32_t corner, std::vector<bool>& face_vis, std::vector<bool>& vert_vis,
                                         std::vector<int>& degree) {
    constexpr int kMaxPriority = 3;
    std::vector<std::int32_t> stacks[kMaxPriority];
    int best = 0;
    auto face_visited = [&](std::int32_t c) { return c == kInvalid || face_vis[std::size_t(c / 3)]; };
    auto visit = [&](std::int32_t v, std::int32_t c) {
        if (v < 0) bad("traversal reached an unmapped corner");
        vert_vis[std::size_t(v)] = true;
        d.point_ids.push_back(c_.faces[std::size_t(c)]);
        d.data_to_corner.push_back(c);
        d.vertex_to_data[std::size_t(v)] = std::int32_t(d.data_to_corner.size()) - 1;
    };
    auto priority = [&](std::int32_t c) {
        const std::int32_t tip = t.Vertex(c);
        int p = 0;
        if (!vert_vis[std::size_t(tip)]) p = ++degree[std::size_t(tip)] > 1 ? 1 : 2;
        return std::min(p, kMaxPriority - 1);
    };
    auto pop = [&]() -> std::int32_t {
        for (int i = best; i < kMaxPriority; ++i)
            if (!stacks[i].empty()) {
                const std::int32_t r = stacks[i].back();
                stacks[i].pop_back();
                best = i;
                return r;
            }
        return kInvalid;
    };
    auto push = [&](std::int32_t c, int p) {
        stacks[p].push_back(c);
        if (p < best) best = p;
    };
    stacks[0].push_back(corner);
    const std::int32_t nv = t.Vertex(t.Next(corner)), pv = t.Vertex(t.Previous(corner)), tip = t.Vertex(corner);
    if (nv < 0 || pv < 0 || tip < 0) bad("traversal reached an unmapped corner");
    if (!vert_vis[std::size_t(nv)]) visit(nv, t.Next(corner));
    if (!vert_vis[std::size_t(pv)]) visit(pv, t.Previous(corner));
    if (!vert_vis[std::size_t(tip)]) visit(tip, corner);
    while ((corner = pop()) != kInvalid) {
        if (face_vis[std::size_t(corner / 3)]) continue;
        while (true) {
            face_vis[std::size_t(corner / 3)] = true;
            const std::int32_t v = t.Vertex(corner);
            if (v < 0) bad("traversal reached an unmapped corner");
            if (!vert_vis[std::size_t(v)]) visit(v, corner);
            const std::int32_t rc = t.GetRightCorner(corner), lc = t.GetLeftCorner(corner);
            const bool rv = face_visited(rc), lv = face_visited(lc);
            if (!lv) {
                const int p = priority(lc);
                if (rv && p <= best) { corner = lc; continue; }
                push(lc, p);
            }
            if (!rv) {
                const int p = priority(rc);
                if (p <= best) { corner = rc; continue; }
                push(rc, p);
            }
            break;
        }
    }
}

template <class T>
void AttributeDecoder::traverse(Dec& d, const T& t) {
    const std::size_t nverts = std::size_t(t.num_vertices());
    std::vector<bool> face_vis(std::size_t(t.num_faces()), false), vert_vis(nverts, false);
    std::size_t map_size = nverts;
    if (d.data_id >= 0) map_size = std::max(map_size, std::size_t(c_.table.num_vertices()));
    d.vertex_to_data.assign(map_size, -1);
    if (d.traversal == TRAVERSAL_DEPTH_FIRST) {
        for (std::int32_t f = 0; f < t.num_faces(); ++f) depth_first(d, t, 3 * f, face_vis, vert_vis);
    } else {
        std::vector<int> degree(nverts, 0);
        for (std::int32_t f = 0; f < t.num_faces(); ++f) prediction_degree(d, t, 3 * f, face_vis, vert_vis, degree);
    }
    // Point -> attribute entry for every attribute of this decoder.
    std::vector<std::uint32_t> p2e(c_.num_points, 0xffffffffu);
    for (std::int32_t c = 0; c < t.num_faces() * 3; ++c) {
        const std::uint32_t point = c_.faces[std::size_t(c)];
        const std::int32_t v = t.Vertex(c);
        if (v < 0 || d.vertex_to_data[std::size_t(v)] < 0) bad("corner vertex never traversed");
        const std::uint32_t entry = std::uint32_t(d.vertex_to_data[std::size_t(v)]);
        if (point >= c_.num_points || entry >= c_.num_points) bad("attribute mapping out of range");
        p2e[point] = entry;
    }
    for (const int id : d.attrs) attrs_[std::size_t(id)].point_to_entry = p2e;
}

void AttributeDecoder::generate_sequence(Dec& d) {
    if (!c_.edgebreaker) {  // linear sequence, identity mapping
        d.point_ids.resize(c_.num_points);
        std::vector<std::uint32_t> id(c_.num_points);
        for (std::uint32_t i = 0; i < c_.num_points; ++i) d.point_ids[i] = id[i] = i;
        for (const int a : d.attrs) attrs_[std::size_t(a)].point_to_entry = id;
        return;
    }
    if (d.att_table) traverse(d, *d.att_table);
    else traverse(d, c_.table);
}

// ---- prediction -----------------------------------------------------------------------------------------------

const Attr& AttributeDecoder::parent_position(const Attr& self) const {
    for (const Attr& a : attrs_) {
        if (a.att_type != ATT_POSITION) continue;
        if (&a == &self) bad("an attribute cannot predict from itself");
        if (a.num_components != 3) bad("parent position attribute must have 3 components");
        if (s_.ver < version(2, 0)) {
            if (!a.has_values) bad("parent position attribute is not decoded yet");
        } else if (!a.has_portable) {
            bad("parent position attribute has no portable values");
        }
        return a;
    }
    bad("prediction needs a POSITION attribute");
}

template <class T>
void AttributeDecoder::predict(const Dec& d, const T& t, const Attr& a, int method, Transform& tr, std::vector<std::int32_t>& data, int nc) {
    const IVec& d2c = d.data_to_corner;
    const IVec& v2d = d.vertex_to_data;
    const std::size_t n = d2c.size();
    if (n * std::size_t(nc) > data.size()) bad("prediction map larger than the value count");
    std::vector<std::int32_t> corr = data;  // corrections in, originals out (out may alias pred source)
    std::int32_t* out = data.data();
    std::vector<std::int32_t> pred(std::size_t(nc), 0);
    auto vtod = [&](std::int32_t v) -> std::int32_t {
        if (v < 0 || std::size_t(v) >= v2d.size()) bad("prediction vertex out of range");
        return v2d[std::size_t(v)];
    };
    auto parallelogram = [&](std::int32_t entry, std::int32_t ci, std::int32_t* res) {
        const std::int32_t oci = t.Opposite(ci);
        if (oci == kInvalid) return false;
        const std::int32_t vo = vtod(t.Vertex(oci)), vn = vtod(t.Vertex(t.Next(oci))), vp = vtod(t.Vertex(t.Previous(oci)));
        if (vo < 0 || vn < 0 || vp < 0) bad("parallelogram vertex without an entry");
        if (!(vo < entry && vn < entry && vp < entry)) return false;
        for (int c = 0; c < nc; ++c) {
            const std::int64_t r = std::int64_t(out[std::size_t(vn) * nc + c]) + out[std::size_t(vp) * nc + c] - out[std::size_t(vo) * nc + c];
            res[c] = std::int32_t(std::uint32_t(std::uint64_t(r)));
        }
        return true;
    };
    auto from_previous = [&](std::size_t p) { tr.original(out + (p - 1) * nc, corr.data() + p * nc, out + p * nc, nc); };

    if (method == PARALLELOGRAM || method == MULTI_PARALLELOGRAM || method == CONSTRAINED_MULTI_PARALLELOGRAM) {
        std::fill(pred.begin(), pred.end(), 0);
        if (n) tr.original(pred.data(), corr.data(), out, nc);
        if (method == PARALLELOGRAM) {
            for (std::size_t p = 1; p < n; ++p) {
                if (!parallelogram(std::int32_t(p), d2c[p], pred.data())) from_previous(p);
                else tr.original(pred.data(), corr.data() + p * nc, out + p * nc, nc);
            }
        } else if (method == MULTI_PARALLELOGRAM) {
            std::vector<std::int32_t> one(static_cast<std::size_t>(nc), 0);
            for (std::size_t p = 1; p < n; ++p) {
                const std::int32_t start = d2c[p];
                std::int32_t c = start;
                int count = 0;
                std::fill(pred.begin(), pred.end(), 0);
                std::size_t guard = 0;
                while (c != kInvalid) {
                    if (++guard > std::size_t(t.num_faces()) * 3 + 1) bad("corner cycle in multi-parallelogram");
                    if (parallelogram(std::int32_t(p), c, one.data())) {
                        for (int k = 0; k < nc; ++k) pred[std::size_t(k)] = add_u(pred[std::size_t(k)], one[std::size_t(k)]);
                        ++count;
                    }
                    c = t.SwingRight(c);
                    if (c == start) c = kInvalid;
                }
                if (count == 0) from_previous(p);
                else {
                    for (auto& v : pred) v /= count;
                    tr.original(pred.data(), corr.data() + p * nc, out + p * nc, nc);
                }
            }
        } else {
            constexpr int kMax = 4;
            std::vector<std::vector<std::int32_t>> pv(kMax, std::vector<std::int32_t>(std::size_t(nc), 0));
            std::vector<std::size_t> crease_pos(kMax, 0);
            for (std::size_t p = 1; p < n; ++p) {
                const std::int32_t start = d2c[p];
                std::int32_t c = start;
                int count = 0;
                bool first_pass = true;
                std::size_t guard = 0;
                while (c != kInvalid) {
                    if (++guard > std::size_t(t.num_faces()) * 6 + 2) bad("corner cycle in constrained multi-parallelogram");
                    if (parallelogram(std::int32_t(p), c, pv[std::size_t(count)].data())) {
                        if (++count == kMax) break;
                    }
                    c = first_pass ? t.SwingLeft(c) : t.SwingRight(c);
                    if (c == start) break;
                    if (c == kInvalid && first_pass) {
                        first_pass = false;
                        c = t.SwingRight(start);
                    }
                }
                int used = 0;
                std::fill(pred.begin(), pred.end(), 0);
                for (int i = 0; i < count; ++i) {
                    const int context = count - 1;
                    auto& flags = crease_[std::size_t(context)];
                    if (crease_pos[std::size_t(context)] >= flags.size()) bad("crease flags exhausted");
                    if (!flags[crease_pos[std::size_t(context)]++]) {
                        ++used;
                        for (int k = 0; k < nc; ++k) pred[std::size_t(k)] = add_u(pred[std::size_t(k)], pv[std::size_t(i)][std::size_t(k)]);
                    }
                }
                if (used == 0) from_previous(p);
                else {
                    for (auto& v : pred) v /= used;
                    tr.original(pred.data(), corr.data() + p * nc, out + p * nc, nc);
                }
            }
        }
        return;
    }

    if (method == TEX_COORDS_PORTABLE || method == TEX_COORDS_DEPRECATED) {
        if (nc != 2) bad("texture-coordinate prediction needs 2 components");
        const Attr& pos = parent_position(a);
        const bool portable = method == TEX_COORDS_PORTABLE || s_.ver >= version(2, 0);
        auto pos_i = [&](std::int32_t entry, std::int64_t* o) {
            const std::uint32_t point = d.point_ids.at(std::size_t(entry));
            const std::uint32_t e = pos.point_to_entry.at(point);
            for (int k = 0; k < 3; ++k) o[k] = pos.portable.at(std::size_t(e) * 3 + std::size_t(k));
        };
        auto pos_f = [&](std::int32_t entry, float* o) {
            const std::uint32_t point = d.point_ids.at(std::size_t(entry));
            const std::uint32_t e = pos.point_to_entry.at(point);
            for (int k = 0; k < 3; ++k) o[k] = portable ? float(pos.portable.at(std::size_t(e) * 3 + std::size_t(k))) : float(pos.values.at(std::size_t(e) * 3 + std::size_t(k)));
        };
        std::vector<bool>& orient = orientations_;
        for (std::size_t p = 0; p < n; ++p) {
            const std::int32_t data_id = std::int32_t(p);
            const std::int32_t ci = d2c[p];
            const std::int32_t next_id = vtod(t.Vertex(t.Next(ci))), prev_id = vtod(t.Vertex(t.Previous(ci)));
            bool done = false;
            if (prev_id < data_id && next_id < data_id) {
                if (method == TEX_COORDS_PORTABLE) {
                    const std::int64_t n_uv[2] = {out[next_id * 2], out[next_id * 2 + 1]};
                    const std::int64_t p_uv[2] = {out[prev_id * 2], out[prev_id * 2 + 1]};
                    if (p_uv[0] == n_uv[0] && p_uv[1] == n_uv[1]) {
                        pred[0] = std::int32_t(p_uv[0]);
                        pred[1] = std::int32_t(p_uv[1]);
                        done = true;
                    } else {
                        std::int64_t tip[3], nx[3], pr[3];
                        pos_i(data_id, tip);
                        pos_i(next_id, nx);
                        pos_i(prev_id, pr);
                        const std::int64_t pn[3] = {pr[0] - nx[0], pr[1] - nx[1], pr[2] - nx[2]};
                        const std::uint64_t pn2 = std::uint64_t(pn[0] * pn[0] + pn[1] * pn[1] + pn[2] * pn[2]);
                        if (pn2 != 0) {
                            const std::int64_t cn[3] = {tip[0] - nx[0], tip[1] - nx[1], tip[2] - nx[2]};
                            const std::int64_t cn_dot_pn = pn[0] * cn[0] + pn[1] * cn[1] + pn[2] * cn[2];
                            const std::int64_t pn_uv[2] = {p_uv[0] - n_uv[0], p_uv[1] - n_uv[1]};
                            const std::int64_t n_abs = std::max(std::abs(n_uv[0]), std::abs(n_uv[1]));
                            if (std::uint64_t(n_abs) > std::uint64_t(INT64_MAX) / pn2) bad("texture-coordinate prediction overflow");
                            const std::int64_t pn_uv_abs = std::max(std::abs(pn_uv[0]), std::abs(pn_uv[1]));
                            if (pn_uv_abs != 0 && std::abs(cn_dot_pn) > INT64_MAX / pn_uv_abs) bad("texture-coordinate prediction overflow");
                            const std::int64_t x_uv[2] = {std::int64_t(std::uint64_t(n_uv[0]) * pn2 + std::uint64_t(pn_uv[0]) * std::uint64_t(cn_dot_pn)),
                                                          std::int64_t(std::uint64_t(n_uv[1]) * pn2 + std::uint64_t(pn_uv[1]) * std::uint64_t(cn_dot_pn))};
                            const std::int64_t pn_abs = std::max({std::abs(pn[0]), std::abs(pn[1]), std::abs(pn[2])});
                            if (pn_abs != 0 && std::abs(cn_dot_pn) > INT64_MAX / pn_abs) bad("texture-coordinate prediction overflow");
                            std::int64_t x_pos[3];
                            for (int k = 0; k < 3; ++k) x_pos[k] = nx[k] + (cn_dot_pn * pn[k]) / std::int64_t(pn2);
                            const std::int64_t cx[3] = {tip[0] - x_pos[0], tip[1] - x_pos[1], tip[2] - x_pos[2]};
                            const std::uint64_t cx2 = std::uint64_t(cx[0] * cx[0] + cx[1] * cx[1] + cx[2] * cx[2]);
                            // IntSqrt
                            std::uint64_t number = cx2 * pn2, root = 0;
                            if (number != 0) {
                                std::uint64_t act = number;
                                root = 1;
                                while (act >= 2) { root *= 2; act /= 4; }
                                do { root = (root + number / root) / 2; } while (root * root > number);
                            }
                            const std::int64_t cx_uv[2] = {std::int64_t(std::uint64_t(pn_uv[1]) * root), std::int64_t(std::uint64_t(-pn_uv[0]) * root)};
                            if (orient.empty()) bad("texture-coordinate orientations exhausted");
                            const bool o = orient.back();
                            orient.pop_back();
                            for (int k = 0; k < 2; ++k) {
                                const std::int64_t num = o ? std::int64_t(std::uint64_t(x_uv[k]) + std::uint64_t(cx_uv[k]))
                                                           : std::int64_t(std::uint64_t(x_uv[k]) - std::uint64_t(cx_uv[k]));
                                pred[std::size_t(k)] = std::int32_t(num / std::int64_t(pn2));
                            }
                            done = true;
                        }
                    }
                } else {  // deprecated float predictor
                    const float n_uv[2] = {float(out[next_id * 2]), float(out[next_id * 2 + 1])};
                    const float p_uv[2] = {float(out[prev_id * 2]), float(out[prev_id * 2 + 1])};
                    auto to_int = [](double v) { return std::isnan(v) || v > INT_MAX || v < INT_MIN ? INT_MIN : int(v); };
                    if (p_uv[0] == n_uv[0] && p_uv[1] == n_uv[1]) {
                        pred[0] = to_int(p_uv[0]);
                        pred[1] = to_int(p_uv[1]);
                    } else {
                        float tip[3], nx[3], pr[3];
                        pos_f(data_id, tip);
                        pos_f(next_id, nx);
                        pos_f(prev_id, pr);
                        const float pn[3] = {pr[0] - nx[0], pr[1] - nx[1], pr[2] - nx[2]};
                        const float cn[3] = {tip[0] - nx[0], tip[1] - nx[1], tip[2] - nx[2]};
                        const float pn2 = pn[0] * pn[0] + pn[1] * pn[1] + pn[2] * pn[2];
                        float sv = 0, tv = 0;
                        if (s_.ver < version(1, 2) || pn2 > 0) {
                            sv = (pn[0] * cn[0] + pn[1] * cn[1] + pn[2] * cn[2]) / pn2;
                            const float r[3] = {cn[0] - pn[0] * sv, cn[1] - pn[1] * sv, cn[2] - pn[2] * sv};
                            tv = std::sqrt((r[0] * r[0] + r[1] * r[1] + r[2] * r[2]) / pn2);
                        }
                        const float pn_uv[2] = {p_uv[0] - n_uv[0], p_uv[1] - n_uv[1]};
                        const float pnus = pn_uv[0] * sv + n_uv[0], pnut = pn_uv[0] * tv, pnvs = pn_uv[1] * sv + n_uv[1], pnvt = pn_uv[1] * tv;
                        if (orient.empty()) bad("texture-coordinate orientations exhausted");
                        const bool o = orient.back();
                        orient.pop_back();
                        const float u = o ? pnus - pnvt : pnus + pnvt, v = o ? pnvs + pnut : pnvs - pnut;
                        pred[0] = to_int(std::floor(double(u) + 0.5));
                        pred[1] = to_int(std::floor(double(v) + 0.5));
                    }
                    done = true;
                }
            }
            if (!done) {
                std::int64_t off = 0;
                if (prev_id < data_id) off = std::int64_t(prev_id) * 2;
                if (next_id < data_id) {
                    off = std::int64_t(next_id) * 2;
                } else if (data_id > 0) {
                    off = std::int64_t(data_id - 1) * 2;
                } else {
                    pred[0] = pred[1] = 0;
                    off = -1;
                }
                if (off >= 0) {
                    pred[0] = out[off];
                    pred[1] = out[off + 1];
                }
            }
            tr.original(pred.data(), corr.data() + p * 2, out + p * 2, 2);
        }
        return;
    }

    if (method == GEOMETRIC_NORMAL) {
        if (nc != 2) bad("normal prediction needs 2 components");
        if (!tr.corrections_positive()) bad("geometric normal prediction without an octahedron transform");
        const Attr& pos = parent_position(a);
        const Octahedron& oct = tr.oct;
        auto pos_corner = [&](std::int32_t c, std::int64_t* o) {
            const std::int32_t entry = vtod(t.Vertex(c));
            if (entry < 0) bad("normal prediction vertex has no entry");
            const std::uint32_t point = d.point_ids.at(std::size_t(entry));
            const std::uint32_t e = pos.point_to_entry.at(point);
            for (int k = 0; k < 3; ++k)
                o[k] = s_.ver < version(2, 0) ? std::int64_t(pos.values.at(std::size_t(e) * 3 + std::size_t(k)))
                                              : std::int64_t(pos.portable.at(std::size_t(e) * 3 + std::size_t(k)));
        };
        for (std::size_t p = 0; p < n; ++p) {
            const std::int32_t ci = d2c[p];
            std::int64_t cent[3];
            pos_corner(ci, cent);
            std::uint64_t normal[3] = {0, 0, 0};
            for_each_vertex_corner(t, ci, [&](std::int32_t c) {
                const std::int32_t cn = one_triangle_ ? t.Next(ci) : t.Next(c);
                const std::int32_t cp = one_triangle_ ? t.Previous(ci) : t.Previous(c);
                std::int64_t pn[3], pp[3];
                pos_corner(cn, pn);
                pos_corner(cp, pp);
                const std::int64_t dn[3] = {pn[0] - cent[0], pn[1] - cent[1], pn[2] - cent[2]};
                const std::int64_t dp[3] = {pp[0] - cent[0], pp[1] - cent[1], pp[2] - cent[2]};
                const std::int64_t cr[3] = {dn[1] * dp[2] - dn[2] * dp[1], dn[2] * dp[0] - dn[0] * dp[2], dn[0] * dp[1] - dn[1] * dp[0]};
                for (int k = 0; k < 3; ++k) normal[k] += std::uint64_t(cr[k]);
            });
            std::int64_t nrm[3] = {std::int64_t(normal[0]), std::int64_t(normal[1]), std::int64_t(normal[2])};
            constexpr std::int64_t upper = 1 << 29;
            if (one_triangle_) {
                const std::int32_t abs_sum = std::int32_t(std::abs(nrm[0]) + std::abs(nrm[1]) + std::abs(nrm[2]));
                if (abs_sum > upper) { const std::int64_t q = abs_sum / upper; for (auto& x : nrm) x /= q; }
            } else {
                const std::int64_t abs_sum = std::abs(nrm[0]) + std::abs(nrm[1]) + std::abs(nrm[2]);
                if (abs_sum > upper) { const std::int64_t q = abs_sum / upper; for (auto& x : nrm) x /= q; }
            }
            std::int32_t v[3] = {std::int32_t(nrm[0]), std::int32_t(nrm[1]), std::int32_t(nrm[2])};
            oct.canonicalize_vector(v);
            if (flips_.next()) for (auto& x : v) x = -x;
            std::int32_t po[2];
            oct.vector_to_coords(v, po[0], po[1]);
            tr.original(po, corr.data() + p * 2, out + p * 2, 2);
        }
        return;
    }
    bad("unhandled prediction method " + std::to_string(method));
}

// ---- value decoding -------------------------------------------------------------------------------------------

void AttributeDecoder::decode_integer(Dec& d, Attr& a, int method, int transform_type) {
    const int nc = a.seq_type == SEQ_NORMALS ? 2 : a.num_components;
    const std::size_t num_entries = d.point_ids.size();
    const std::size_t num_values = num_entries * std::size_t(nc);
    if (num_values > (std::size_t(1) << 31)) bad("too many attribute values");

    // Prediction scheme selection ("CreateIntPredictionScheme" + "CreatePredictionSchemeForDecoder").
    bool scheme = false;
    int kind = PREDICTION_DIFFERENCE;
    Transform tr;
    if (method != PREDICTION_NONE) {
        const bool normal = a.seq_type == SEQ_NORMALS;
        const bool ok_transform = normal ? (transform_type == TRANSFORM_OCTAHEDRON || transform_type == TRANSFORM_OCTAHEDRON_CANONICALIZED)
                                         : transform_type == TRANSFORM_WRAP;
        if (ok_transform) {
            scheme = true;
            tr.type = transform_type;
            const bool mesh_method = method == PARALLELOGRAM || method == MULTI_PARALLELOGRAM || method == CONSTRAINED_MULTI_PARALLELOGRAM ||
                                     method == TEX_COORDS_PORTABLE || method == GEOMETRIC_NORMAL || method == TEX_COORDS_DEPRECATED;
            if (mesh_ && c_.edgebreaker && mesh_method) {
                // Octahedron transforms support only the geometric-normal mesh scheme.
                if (!normal || method == GEOMETRIC_NORMAL) kind = method;
            }
            if ((kind == TEX_COORDS_PORTABLE || kind == TEX_COORDS_DEPRECATED || kind == GEOMETRIC_NORMAL)) parent_position(a);
        }
    }

    std::vector<std::int32_t> data(num_values, 0);
    const int compressed = s_.u8();
    if (compressed > 0) {
        decode_symbols(s_, std::uint32_t(num_values), nc, reinterpret_cast<std::uint32_t*>(data.data()));
    } else {
        const int num_bytes = s_.u8();
        if (num_bytes < 1 || num_bytes > 4) bad("invalid raw integer size");
        if (std::uint64_t(num_bytes) * num_values > s_.remaining()) bad("raw integer values exceed the stream");
        for (auto& v : data) {
            std::uint32_t x = 0;
            const std::uint8_t* b = s_.bytes(std::size_t(num_bytes));
            for (int k = 0; k < num_bytes; ++k) x |= std::uint32_t(b[k]) << (8 * k);
            v = std::int32_t(x);
        }
    }
    if (num_values > 0 && !(scheme && tr.corrections_positive()))
        for (auto& v : data) {
            const std::uint32_t u = std::uint32_t(v);
            v = (u & 1) ? -std::int32_t(u >> 1) - 1 : std::int32_t(u >> 1);
        }
    if (scheme) {
        // DecodePredictionData
        crease_.assign(4, {});
        orientations_.clear();
        flips_ = BitDecoder();
        one_triangle_ = false;
        if (kind == CONSTRAINED_MULTI_PARALLELOGRAM) {
            if (s_.ver < version(2, 2) && s_.u8() != 0) unsupported("constrained multi-parallelogram mode");
            for (int i = 0; i < 4; ++i) {
                const std::uint32_t num_flags = s_.varint32();
                if (num_flags > std::uint32_t(c_.table.num_corners())) bad("too many crease flags");
                if (num_flags > 0) {
                    BitDecoder bd;
                    bd.start(s_);
                    crease_[std::size_t(i)].resize(num_flags);
                    for (std::uint32_t j = 0; j < num_flags; ++j) crease_[std::size_t(i)][j] = bd.next();
                }
            }
            tr.decode(s_);
        } else if (kind == TEX_COORDS_PORTABLE) {
            const std::int32_t num = s_.i32();
            if (num < 0 || std::size_t(num) > num_values) bad("invalid texture-coordinate orientation count");
            orientations_.resize(std::size_t(num));
            bool last = true;
            BitDecoder bd;
            bd.start(s_);
            for (std::int32_t i = 0; i < num; ++i) {
                if (!bd.next()) last = !last;
                orientations_[std::size_t(i)] = last;
            }
            tr.decode(s_);
        } else if (kind == TEX_COORDS_DEPRECATED) {
            const std::uint32_t num = s_.ver < version(2, 2) ? s_.u32() : s_.varint32();
            if (num == 0 || num > std::uint32_t(c_.table.num_corners())) bad("invalid texture-coordinate orientation count");
            orientations_.resize(num);
            bool last = true;
            BitDecoder bd;
            bd.start(s_);
            for (std::uint32_t i = 0; i < num; ++i) {
                if (!bd.next()) last = !last;
                orientations_[i] = last;
            }
            tr.decode(s_);
        } else if (kind == GEOMETRIC_NORMAL) {
            tr.decode(s_);
            if (s_.ver < version(2, 2)) {
                const int mode = s_.u8();
                if (mode > 1) bad("invalid normal prediction mode");
                one_triangle_ = mode == 0;
            }
            flips_.start(s_);
        } else {
            tr.decode(s_);
        }
        if (num_values > 0) {
            if (kind == PREDICTION_DIFFERENCE) {
                std::vector<std::int32_t> corr = data, zero(std::size_t(nc), 0);
                tr.original(zero.data(), corr.data(), data.data(), nc);
                for (std::size_t i = std::size_t(nc); i < num_values; i += std::size_t(nc))
                    tr.original(data.data() + i - std::size_t(nc), corr.data() + i, data.data() + i, nc);
            } else if (d.att_table) {
                predict(d, *d.att_table, a, kind, tr, data, nc);
            } else {
                predict(d, c_.table, a, kind, tr, data, nc);
            }
        }
    }
    a.portable = std::move(data);
    a.portable_components = nc;
    a.has_portable = true;
}

void AttributeDecoder::decode_transform_data(Attr& a) {
    if (a.seq_type == SEQ_QUANTIZATION) {
        a.min_values.resize(std::size_t(a.num_components));
        for (auto& v : a.min_values) v = s_.f32();
        a.range = s_.f32();
        a.quantization_bits = s_.u8();
        if (a.quantization_bits < 1 || a.quantization_bits > 30) bad("invalid quantization bits");
    } else if (a.seq_type == SEQ_NORMALS) {
        a.quantization_bits = s_.u8();
    }
}

void AttributeDecoder::transform_to_original(Attr& a) {
    const std::size_t n = a.portable.size() / std::size_t(std::max(1, a.portable_components));
    a.values.clear();
    if (a.seq_type == SEQ_INTEGER) {
        a.values.resize(a.portable.size());
        for (std::size_t i = 0; i < a.portable.size(); ++i) {
            const std::int32_t v = a.portable[i];
            switch (a.data_type) {
                case 1: a.values[i] = std::int8_t(v); break;
                case 2: a.values[i] = std::uint8_t(v); break;
                case 3: a.values[i] = std::int16_t(v); break;
                case 4: a.values[i] = std::uint16_t(v); break;
                case 5: a.values[i] = v; break;
                case 6: a.values[i] = std::uint32_t(v); break;
                default: bad("integer attribute with a non-integer data type");
            }
        }
    } else if (a.seq_type == SEQ_QUANTIZATION) {
        if (a.data_type != 9) bad("quantized attribute must be float32");
        const std::int32_t max_q = std::int32_t((1u << std::uint32_t(a.quantization_bits)) - 1);
        const float delta = a.range / float(max_q);
        a.values.resize(a.portable.size());
        const std::size_t nc = std::size_t(a.num_components);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t c = 0; c < nc; ++c) a.values[i * nc + c] = float(a.portable[i * nc + c]) * delta + a.min_values[c];
    } else if (a.seq_type == SEQ_NORMALS) {
        if (a.data_type != 9 || a.num_components != 3) bad("normal attribute must be 3 x float32");
        Octahedron oct;
        oct.set_bits(a.quantization_bits);
        a.values.resize(n * 3);
        for (std::size_t i = 0; i < n; ++i) {
            float v[3];
            oct.to_unit_vector(a.portable[i * 2], a.portable[i * 2 + 1], v);
            for (int k = 0; k < 3; ++k) a.values[i * 3 + std::size_t(k)] = v[k];
        }
    }
    a.has_values = true;
}

void AttributeDecoder::decode_portable(Dec& d, Attr& a) {
    if (a.seq_type == SEQ_GENERIC) {
        const int len = data_type_length(a.data_type);
        if (len == 0) bad("generic attribute with an invalid data type");
        const std::size_t n = d.point_ids.size(), nc = std::size_t(a.num_components);
        if (std::uint64_t(n) * nc * std::uint64_t(len) > s_.remaining()) bad("generic attribute values exceed the stream");
        a.values.resize(n * nc);
        for (std::size_t i = 0; i < n * nc; ++i) {
            const std::uint8_t* b = s_.bytes(std::size_t(len));
            double v = 0;
            switch (a.data_type) {
                case 1: v = std::int8_t(b[0]); break;
                case 2: case 11: v = b[0]; break;
                case 3: { std::int16_t x; std::memcpy(&x, b, 2); v = x; break; }
                case 4: { std::uint16_t x; std::memcpy(&x, b, 2); v = x; break; }
                case 5: { std::int32_t x; std::memcpy(&x, b, 4); v = x; break; }
                case 6: { std::uint32_t x; std::memcpy(&x, b, 4); v = x; break; }
                case 7: { std::int64_t x; std::memcpy(&x, b, 8); v = double(x); break; }
                case 8: { std::uint64_t x; std::memcpy(&x, b, 8); v = double(x); break; }
                case 9: { float x; std::memcpy(&x, b, 4); v = x; break; }
                case 10: { double x; std::memcpy(&x, b, 8); v = x; break; }
                default: break;
            }
            a.values[i] = v;
        }
        a.has_values = true;
        return;
    }
    const int method = s_.i8();
    if (method < PREDICTION_NONE || method >= NUM_PREDICTION_SCHEMES) bad("invalid prediction scheme " + std::to_string(method));
    int transform_type = TRANSFORM_NONE;
    if (method != PREDICTION_NONE) {
        transform_type = s_.i8();
        if (transform_type < TRANSFORM_NONE || transform_type >= NUM_TRANSFORMS) bad("invalid prediction transform");
    }
    if (s_.ver < version(2, 0)) decode_transform_data(a);  // legacy: parameters precede the values
    decode_integer(d, a, method, transform_type);
    if (s_.ver < version(2, 0)) transform_to_original(a);
}

// ---- kd-tree point clouds ------------------------------------------------------------------------------------

void AttributeDecoder::decode_kd_tree(Dec& d) {
    if (s_.ver < version(2, 3)) unsupported("kd-tree point clouds older than bitstream 2.3");
    const std::uint32_t n = c_.num_points;
    const int level = s_.u8();
    std::uint32_t dims = 0;
    for (const int id : d.attrs) {
        const Attr& a = attrs_[std::size_t(id)];
        const int dt = a.data_type;
        if (!(dt == 1 || dt == 2 || dt == 3 || dt == 4 || dt == 5 || dt == 6 || dt == 9)) bad("kd-tree attribute with an unsupported data type");
        dims += std::uint32_t(a.num_components);
    }
    std::vector<std::uint32_t> pts;
    decode_kd_tree_points(s_, level, dims, n, pts);
    // Transform parameters: float attributes (min, range, bits), then signed minimums.
    std::size_t signed_components = 0;
    for (const int id : d.attrs) {
        Attr& a = attrs_[std::size_t(id)];
        if (a.data_type == 9) {
            a.min_values.resize(std::size_t(a.num_components));
            for (auto& v : a.min_values) v = s_.f32();
            a.range = s_.f32();
            a.quantization_bits = s_.u8();
            if (a.quantization_bits > 31) bad("invalid kd-tree quantization bits");
        } else if (a.data_type == 1 || a.data_type == 3 || a.data_type == 5) {
            signed_components += std::size_t(a.num_components);
        }
    }
    std::vector<std::int32_t> mins(signed_components);
    for (auto& m : mins) {
        const std::uint32_t u = s_.varint32();
        m = (u & 1) ? -std::int32_t(u >> 1) - 1 : std::int32_t(u >> 1);
    }
    std::uint32_t offset = 0;
    std::size_t signed_pos = 0;
    std::vector<std::uint32_t> identity(n);
    for (std::uint32_t i = 0; i < n; ++i) identity[i] = i;
    for (const int id : d.attrs) {
        Attr& a = attrs_[std::size_t(id)];
        const std::size_t nc = std::size_t(a.num_components);
        a.values.resize(std::size_t(n) * nc);
        a.point_to_entry = identity;
        float delta = 0;
        if (a.data_type == 9) {
            const std::int32_t max_q = std::int32_t((1u << std::uint32_t(a.quantization_bits)) - 1);
            if (max_q <= 0) bad("invalid kd-tree quantization");
            delta = a.range / float(max_q);
        }
        for (std::uint32_t p = 0; p < n; ++p)
            for (std::size_t c = 0; c < nc; ++c) {
                const std::uint32_t u = pts[std::size_t(p) * dims + offset + c];
                double v = 0;
                switch (a.data_type) {
                    case 2: v = std::uint8_t(u); break;
                    case 4: v = std::uint16_t(u); break;
                    case 6: v = u; break;
                    case 1: case 3: case 5: {
                        const std::uint32_t raw = a.data_type == 1 ? std::uint8_t(u) : a.data_type == 3 ? std::uint16_t(u) : u;
                        if (raw > std::uint32_t(INT32_MAX)) bad("kd-tree signed value out of range");
                        const std::int32_t sv = std::int32_t(raw) + mins[signed_pos + c];
                        v = a.data_type == 1 ? double(std::int8_t(sv)) : a.data_type == 3 ? double(std::int16_t(sv)) : double(sv);
                        break;
                    }
                    case 9: v = float(std::int32_t(u)) * delta + a.min_values[c]; break;
                    default: break;
                }
                a.values[std::size_t(p) * nc + c] = v;
            }
        if (a.data_type == 1 || a.data_type == 3 || a.data_type == 5) signed_pos += nc;
        a.has_values = true;
        offset += std::uint32_t(nc);
    }
}

// ---- driver -----------------------------------------------------------------------------------------------------

std::vector<Attribute> AttributeDecoder::run() {
    const int num_decoders = s_.u8();
    decs_.resize(std::size_t(num_decoders));
    for (int i = 0; i < num_decoders; ++i) {
        Dec& d = decs_[std::size_t(i)];
        if (c_.edgebreaker) {
            d.data_id = s_.i8();
            d.decoder_type = s_.u8();
            if (d.data_id >= 0) {
                if (std::size_t(d.data_id) >= c_.attribute_data.size()) bad("attribute decoder refers to missing attribute data");
                AttributeData& ad = c_.attribute_data[std::size_t(d.data_id)];
                if (ad.decoder_id >= 0) bad("attribute data used by two decoders");
                ad.decoder_id = i;
            } else {
                if (c_.pos_data_decoder_id >= 0) bad("two decoders use the position data");
                c_.pos_data_decoder_id = i;
            }
            if (s_.ver >= version(1, 2)) {
                d.traversal = s_.u8();
                if (d.traversal > TRAVERSAL_PREDICTION_DEGREE) bad("unknown traversal method");
            }
            if (d.decoder_type == MESH_VERTEX_ATTRIBUTE) {
                if (d.data_id >= 0) c_.attribute_data[std::size_t(d.data_id)].connectivity_used = false;
            } else {
                if (d.traversal != TRAVERSAL_DEPTH_FIRST) bad("corner attributes support only depth-first traversal");
                if (d.data_id < 0) bad("corner attribute decoder without attribute data");
                d.att_table = &c_.attribute_data[std::size_t(d.data_id)].table;
            }
        }
    }
    for (Dec& d : decs_) {
        const std::uint32_t n = s_.ver < version(2, 0) ? s_.u32() : s_.varint32();
        if (n == 0) bad("attribute decoder without attributes");
        if (n > 5 * s_.remaining()) bad("attribute count exceeds the stream");
        for (std::uint32_t j = 0; j < n; ++j) {
            Attr a;
            a.att_type = s_.u8();
            a.data_type = s_.u8();
            a.num_components = s_.u8();
            a.normalized = s_.u8() > 0;
            if (a.att_type >= kNumNamedAttributeTypes) bad("unknown attribute type " + std::to_string(a.att_type));
            if (a.data_type < 1 || a.data_type > 11) bad("invalid attribute data type");
            if (a.num_components == 0) bad("attribute without components");
            a.unique_id = s_.ver < version(1, 3) ? s_.u16() : s_.varint32();
            d.attrs.push_back(int(attrs_.size()));
            attrs_.push_back(std::move(a));
        }
        for (const int id : d.attrs) {
            if (kd_) break;  // kd-tree decoders have no per-attribute decoder types
            Attr& a = attrs_[std::size_t(id)];
            a.seq_type = s_.u8();
            if (a.seq_type > SEQ_NORMALS) bad("unknown sequential attribute decoder " + std::to_string(a.seq_type));
            if (a.seq_type == SEQ_NORMALS && (a.num_components != 3 || a.data_type != 9)) bad("normal decoder needs 3 x float32");
        }
    }
    for (Dec& d : decs_) {
        if (kd_) {
            decode_kd_tree(d);
            continue;
        }
        generate_sequence(d);
        for (const int id : d.attrs) decode_portable(d, attrs_[std::size_t(id)]);
        if (s_.ver >= version(2, 0)) {
            for (const int id : d.attrs) decode_transform_data(attrs_[std::size_t(id)]);
            for (const int id : d.attrs)
                if (attrs_[std::size_t(id)].seq_type != SEQ_GENERIC) transform_to_original(attrs_[std::size_t(id)]);
        }
    }

    std::vector<Attribute> out;
    for (const Attr& a : attrs_) {
        Attribute r;
        r.unique_id = a.unique_id;
        r.type = a.att_type;
        r.components = a.num_components;
        const std::size_t nc = std::size_t(a.num_components);
        const std::size_t nentries = a.values.size() / nc;
        r.values.resize(std::size_t(c_.num_points) * nc);
        double scale = 1;
        if (a.normalized && a.seq_type != SEQ_QUANTIZATION && a.seq_type != SEQ_NORMALS && a.data_type < 9) {
            static constexpr double kMax[9] = {1, 127, 255, 32767, 65535, 2147483647.0, 4294967295.0, 9.2233720368547758e18, 1.8446744073709552e19};
            scale = 1.0 / kMax[a.data_type];
        }
        for (std::uint32_t p = 0; p < c_.num_points; ++p) {
            const std::uint32_t e = a.point_to_entry.at(p);
            if (e >= nentries) bad("point without an attribute value");
            for (std::size_t k = 0; k < nc; ++k) {
                double v = a.values[std::size_t(e) * nc + k] * scale;
                if (scale != 1 && v < -1) v = -1;  // signed normalized: max(c / max, -1)
                r.values[std::size_t(p) * nc + k] = float(v);
            }
        }
        out.push_back(std::move(r));
    }
    return out;
}

}  // namespace

std::vector<Attribute> decode_attributes(Stream& s, Connectivity& conn, bool is_mesh, bool kd_tree) {
    AttributeDecoder d(s, conn, is_mesh, kd_tree);
    return d.run();
}

}  // namespace cadcore::codec::draco
