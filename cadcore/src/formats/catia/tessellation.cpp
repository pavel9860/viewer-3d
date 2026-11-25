#include "formats/catia/tessellation.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstring>
#include <string>

#include "cadcore/error.hpp"

namespace cadcore::formats::catia {
namespace {

constexpr std::uint8_t kNode[] = {0x04, 0x04, 0xfd, 0xff};
constexpr std::size_t kTokenAt = 33;  // list token offset inside a node header
constexpr std::uint32_t kMaxCount = 50'000'000;

[[noreturn]] void bad(const std::string& what, std::size_t at) { fail(ErrorCode::Malformed, "CATIA tessellation: " + what, at); }

// Byte reader with the CFV3 stream-wide bit reader (refill words are taken inline).
class Reader {
public:
    Reader(ByteView s, Generation g) : s_(s), g_(g) {}

    std::size_t pos() const { return p_; }
    void seek(std::size_t p) { p_ = p; }
    Generation generation() const { return g_; }

    void need(std::size_t n) const {
        if (n > s_.size - p_) fail(ErrorCode::Truncated, "CATIA tessellation: stream ends inside a structure", p_);
    }
    std::uint8_t u8() { need(1); return s_.data[p_++]; }
    std::uint32_t u32() {
        need(4);
        std::uint32_t v;
        std::memcpy(&v, s_.data + p_, 4);  // little-endian hosts only (x86, wasm)
        p_ += 4;
        return v;
    }
    std::uint32_t count() {  // P1g: one byte, or ff + u32le
        const std::uint8_t v = u8();
        return v == 0xff ? u32() : v;
    }
    void floats(float* out, std::size_t n) {
        need(4 * n);
        if (n == 0) return;  // `out` may be the null data() of an empty vector
        std::memcpy(out, s_.data + p_, 4 * n);
        p_ += 4 * n;
    }
    void skip(std::size_t n) { need(n); p_ += n; }
    const std::uint8_t* here() const { return s_.data + p_; }

    std::uint32_t bits(int n) {
        std::uint32_t v = 0;
        for (int k = 0; k < n;) {
            if (left_ == 0) { word_ = u32(); left_ = 32; }
            const int take = std::min(n - k, left_);
            const std::uint32_t mask = take == 32 ? 0xffffffffu : ((1u << take) - 1);
            v |= (word_ & mask) << k;
            word_ = take == 32 ? 0 : word_ >> take;
            left_ -= take;
            k += take;
        }
        return v;
    }

private:
    ByteView s_;
    Generation g_;
    std::size_t p_ = 0;
    std::uint32_t word_ = 0;
    int left_ = 0;
};

// How one list stores its lengths and references.
struct Fields {
    Generation g;
    int width;  // CFV2: bytes per reference; CFV3: bits per field

    static Fields for_count(Generation g, std::uint32_t count) {
        if (g == Generation::Cfv3) return {g, std::max(1, int(std::bit_width(count)))};
        return {g, count < 0x100 ? 1 : count < 0x10000 ? 2 : count < 0x1000000 ? 3 : 4};
    }
    std::uint32_t length(Reader& r) const { return g == Generation::Cfv3 ? r.bits(width) : r.count(); }
    std::uint32_t ref(Reader& r) const {
        if (g == Generation::Cfv3) return r.bits(width);
        std::uint32_t v = 0;
        for (int k = 0; k < width; ++k) v = (v << 8) | r.u8();
        return v;
    }
};

void expect(Reader& r, std::uint8_t want, const char* what) {
    const std::size_t at = r.pos();
    if (r.u8() != want) bad(std::string("expected ") + what, at);
}

std::vector<float> points(Reader& r, std::uint32_t& count, std::uint32_t& normal_count) {
    const std::size_t at = r.pos();
    count = r.count();
    normal_count = r.count();
    if (count == 0 || count > kMaxCount || normal_count > count) bad("point list counts out of range", at);
    const std::uint8_t has_tex = r.u8();
    if (has_tex > 1 || r.u8() != 0 || r.u8() != 0) bad("point list flags", at);
    if (has_tex) {  // texture coordinates: no texture image is stored with them
        const std::uint8_t dims = r.u8();
        if (dims < 1 || dims > 4) bad("texture coordinate dimension", r.pos() - 1);
        r.skip(std::size_t(4) * dims * count);
    }
    const std::size_t eat = r.pos();
    const std::uint32_t encoding = r.generation() == Generation::Cfv3 ? r.bits(1) : r.u8();
    std::vector<float> xyz(std::size_t(count) * 3);
    if (encoding == 1) {
        r.floats(xyz.data(), xyz.size());
    } else if (encoding == 0) {
        const std::size_t tl = (std::size_t(count) + 3) / 4;
        r.need(tl);
        const std::uint8_t* table = r.here();
        r.skip(tl);
        static constexpr int kWords[4] = {3, 0, 1, 2};
        std::uint64_t need = 0;
        for (std::uint32_t i = 0; i < count; ++i) need += kWords[(table[i >> 2] >> ((i & 3) * 2)) & 3];
        const std::size_t wat = r.pos();
        if (r.count() != need) bad("point table word count does not match its opcodes", wat);
        std::vector<float> w(need);
        r.floats(w.data(), w.size());
        float x = 0, y = 0, z = 0;
        std::size_t q = 0;
        for (std::uint32_t i = 0; i < count; ++i) {
            switch ((table[i >> 2] >> ((i & 3) * 2)) & 3) {
                case 0: x = w[q]; y = w[q + 1]; z = w[q + 2]; q += 3; break;
                case 3: y = w[q]; z = w[q + 1]; q += 2; break;
                case 2: z = w[q]; q += 1; break;
                default: break;  // 1: same position as the previous point
            }
            xyz[i * 3] = x; xyz[i * 3 + 1] = y; xyz[i * 3 + 2] = z;
        }
    } else {
        bad("point encoding " + std::to_string(encoding), eat);
    }
    for (float v : xyz)
        if (!std::isfinite(v)) bad("non-finite coordinate", at);
    return xyz;
}

void skip_normals(Reader& r, std::uint32_t normal_count) {
    if (normal_count == 0) return;
    const std::uint32_t m = r.count();
    r.skip(std::size_t(m) * 2);
    r.skip((std::size_t(normal_count) + 1) / 2);
}

// Facet record -> triangles (strips alternate orientation, fans pivot on their first point).
std::vector<std::uint32_t> record(Reader& r, const Fields& f, std::uint32_t count) {
    const std::size_t at = r.pos();
    expect(r, 0x01, "facet record 01");
    const std::uint8_t mk = r.u8();
    if ((mk & 0xf0) != 0x40) bad("facet record marker", at + 1);
    std::uint32_t c[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k)
        if (mk & (1 << k)) c[k] = r.count();
    expect(r, 0xff, "ff before the reference count");
    const std::uint32_t n = r.u32();
    if (n > kMaxCount || c[0] > n || c[1] > n || c[2] > n) bad("facet record counts out of range", at);
    if (mk & 0x08) r.skip(12);  // planar face normal
    std::vector<std::uint32_t> strips(c[1]), fans(c[2]);
    std::uint64_t sum = std::uint64_t(3) * c[0];
    for (auto& v : strips) sum += v = f.length(r);
    for (auto& v : fans) sum += v = f.length(r);
    if (sum != n) bad("strip and fan lengths do not add up to the reference count", at);
    std::vector<std::uint32_t> refs(n);
    for (auto& v : refs)
        if ((v = f.ref(r)) >= count) bad("reference beyond the point list", at);
    std::vector<std::uint32_t> out;
    out.reserve(std::size_t(n) * 3);
    auto emit = [&](std::uint32_t a, std::uint32_t b, std::uint32_t d) {
        if (a != b && b != d && a != d) out.insert(out.end(), {a, b, d});
    };
    std::size_t k = 0;
    for (std::uint32_t t = 0; t < c[0]; ++t, k += 3) emit(refs[k], refs[k + 1], refs[k + 2]);
    for (auto len : strips) {
        for (std::uint32_t t = 0; t + 2 < len; ++t)
            t % 2 == 0 ? emit(refs[k + t], refs[k + t + 1], refs[k + t + 2]) : emit(refs[k + t + 1], refs[k + t], refs[k + t + 2]);
        k += len;
    }
    for (auto len : fans) {
        for (std::uint32_t t = 1; t + 1 < len; ++t) emit(refs[k], refs[k + t], refs[k + t + 1]);
        k += len;
    }
    return out;
}

struct Attribute {
    std::uint8_t flags;
    Color color;
};

Attribute attribute(Reader& r, bool edge) {
    const std::size_t at = r.pos();
    r.need(8);
    const std::uint8_t* a = r.here();
    const bool ok = a[2] == 0x04 && a[3] == 0xff && (a[1] == 0x04 || (edge && a[1] == 0x24));
    if (!ok) bad(edge ? "edge attribute" : "face attribute", at);
    r.skip(8);
    return {a[0], Color{a[7], a[6], a[5], a[4]}};
}

// count = 0: a polyline node, whose face list may only hold explicit polylines.
void face_list(Reader& r, TessBody& body, std::uint32_t count) {
    const Fields f = Fields::for_count(r.generation(), count);
    float sag;
    r.floats(&sag, 1);
    const std::uint8_t nlod = r.u8();
    for (std::uint8_t k = 0; k < nlod; ++k) {  // coarse whole-body meshes: consumed, not shown
        record(r, f, count);
        r.skip(4);
        attribute(r, false);
    }
    const std::uint8_t groups = r.u8();
    for (std::uint8_t g = 0; g < groups; ++g) {
        const std::size_t at = r.pos();
        const std::uint8_t mode = r.u8(), type = r.u8();
        const std::uint32_t items = r.count();
        if (mode != 1 && mode != 2) bad("group attribute mode " + std::to_string(mode), at);
        if (items > kMaxCount) bad("group item count out of range", at);
        if (type == 0x0a) {
            const std::size_t first = body.faces.size();
            for (std::uint32_t i = 0; i < items; ++i) body.faces.push_back({record(r, f, count), 0, {}});
            const std::size_t attrs = mode == 2 ? items : 1;
            for (std::size_t i = 0; i < attrs; ++i) {
                const Attribute a = attribute(r, false);
                const std::size_t lo = mode == 2 ? first + i : first, hi = mode == 2 ? lo + 1 : body.faces.size();
                for (std::size_t j = lo; j < hi; ++j) body.faces[j].flags = a.flags, body.faces[j].color = a.color;
            }
        } else if (type <= 2) {
            if (mode != 1) bad("polyline group with per-item attributes", at);
            TessLines lines;
            for (std::uint32_t i = 0; i < items; ++i) {
                const std::size_t iat = r.pos();
                const std::uint8_t kind = r.u8();
                const std::uint32_t n = r.count();
                if (n > kMaxCount) bad("polyline length out of range", iat);
                std::vector<float> pts(std::size_t(n) * 3);
                if (kind == 0x02) {  // indexed
                    for (std::uint32_t k = 0; k < n; ++k) {
                        const std::uint32_t v = f.ref(r);
                        if (v >= count) bad("polyline reference beyond the point list", iat);
                        std::copy_n(body.positions.data() + 3 * std::size_t(v), 3, pts.data() + 3 * std::size_t(k));
                    }
                } else if (kind == 0x03) {  // explicit points
                    expect(r, 0x01, "01 after an explicit polyline length");
                    r.floats(pts.data(), pts.size());
                    for (float v : pts)
                        if (!std::isfinite(v)) bad("non-finite polyline coordinate", iat);
                } else {
                    bad("polyline item kind " + std::to_string(kind), iat);
                }
                lines.polylines.push_back(std::move(pts));
            }
            lines.color = attribute(r, true).color;
            body.lines.push_back(std::move(lines));
        } else {
            bad("group type " + std::to_string(type), at);
        }
    }
    if (count == 0 && !body.faces.empty()) bad("facets in a node without points", body.offset);
}

}  // namespace

std::vector<TessBody> decode_tessellation(ByteView s, Generation generation) {
    Reader r(s, generation);
    std::vector<TessBody> bodies;
    const std::string_view node(reinterpret_cast<const char*>(kNode), 4);
    std::size_t pos = 0, n;
    while ((n = s.find(node, pos)) != ByteView::npos) {
        if (n + kTokenAt + 6 > s.size) bad("node header ends with the stream", n);
        const std::uint8_t* t = s.data + n + kTokenAt;
        // Only the list token ff ff {00|02} 00 {00|01} marks geometry; root, group and symbol
        // nodes carry other values there, and their children follow their header.
        const bool token = t[0] == 0xff && t[1] == 0xff && (t[2] == 0 || t[2] == 2) && t[3] == 0 && t[4] <= 1;
        const bool geometry = token && t[4] == 1, polyline = token && t[4] == 0 && t[5] != 0xff;
        if (!geometry && !polyline) {
            pos = n + 4;
            continue;
        }
        // Geometry and polyline nodes are framed: ff u32le length, 7 bytes before the tag.
        if (n < 7 || s[n - 7] != 0xff) bad("geometry node without framing", n);
        std::uint32_t framed;
        std::memcpy(&framed, s.data + n - 6, 4);
        const std::uint64_t end = n - 7 + std::uint64_t(framed);
        if (end > s.size) bad("framed node beyond the stream", n);
        r.seek(n + kTokenAt + 5);
        TessBody b;
        b.offset = n;
        if (geometry) {
            std::uint32_t count, normal_count;
            b.positions = points(r, count, normal_count);
            skip_normals(r, normal_count);
            face_list(r, b, count);
        } else {
            face_list(r, b, 0);
        }
        if (r.pos() != end) bad("node does not end at its framed length", r.pos());
        bodies.push_back(std::move(b));
        pos = std::size_t(end);
    }
    return bodies;
}

void add_to_scene(std::vector<TessBody>&& bodies, const std::string& name, Scene& scene, std::uint32_t node) {
    for (auto& b : bodies) {
        for (const auto& lines : b.lines)
            for (auto& pts : lines.polylines) {
                Polyline pl;
                pl.color = lines.color;
                pl.points = std::move(pts);
                scene.nodes[node].polylines.push_back(std::uint32_t(scene.polylines.size()));
                scene.polylines.push_back(std::move(pl));
            }
        if (b.faces.empty()) continue;
        Mesh m;
        m.name = name + " @" + std::to_string(b.offset);
        bool any_solid = false, any_surface = false;
        for (const auto& f : b.faces) {
            const auto first = std::uint32_t(m.triangle_count());
            m.indices.insert(m.indices.end(), f.triangles.begin(), f.triangles.end());
            const auto count = std::uint32_t(m.triangle_count()) - first;
            if (!count) continue;
            m.groups.push_back({first, count, f.color});
            ((f.flags & 0x10) ? any_solid : any_surface) = true;
        }
        if (m.indices.empty()) continue;
        m.kind = any_solid && !any_surface ? BodyKind::Solid : any_surface && !any_solid ? BodyKind::Surface : BodyKind::Unknown;
        m.color = m.groups[0].color;
        m.positions = std::move(b.positions);
        scene.nodes[node].meshes.push_back(std::uint32_t(scene.meshes.size()));
        scene.meshes.push_back(std::move(m));
    }
}

}  // namespace cadcore::formats::catia
