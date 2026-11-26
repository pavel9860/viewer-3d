#include "formats/catia/cgm_compact.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <optional>
#include <string>

#include "cadcore/error.hpp"
#include "formats/catia/cgm_compact_io.hpp"

namespace cadcore::formats::catia::cgm {
namespace {
using namespace compact_io;

int pm(int bit) { return bit ? 1 : -1; }

double load_f64_le(const std::uint8_t* p) { double d; std::memcpy(&d, p, 8); return d; }

// Knot vector: int degree, int n, u8 z (4: x = 0 omitted; 0: int x follows), f64 knots[n], u8 1
// (standard multiplicities: ends degree + 1, interior degree - 2). Written in the varint dialect's
// form: i32 degree, 0, x, n, 0, knots, i32 multiplicities. Returns the number of poles.
// Values a 2012+ record may take from others: surface parameter boxes, p-curve and edge-curve ranges.
struct Ctx {
    std::unordered_map<std::uint32_t, int> edge_uses;  // edge curve -> edges referencing it
    std::unordered_map<std::uint32_t, std::uint16_t> type;
    std::unordered_map<std::uint32_t, std::array<double, 4>> box;   // u0 u1 v0 v1
    std::unordered_map<std::uint32_t, std::array<double, 2>> range;  // p-curves, edge curves
    const std::array<double, 4>& surface_box(std::uint32_t id, std::size_t at) const {
        auto it = box.find(id);
        if (it == box.end()) {
            auto t = type.find(id);
            throw DecodeError(ErrorCode::Unsupported, "parameter box of surface class " + std::to_string(t == type.end() ? 0 : t->second) + " not established", at);
        }
        return it->second;
    }
    const std::array<double, 2>& curve_range(std::uint32_t id, std::size_t at) const {
        auto it = range.find(id);
        if (it == range.end()) {
            auto t = type.find(id);
            throw DecodeError(ErrorCode::Unsupported, "range of curve class " + std::to_string(t == type.end() ? 0 : t->second) + " not available", at);
        }
        return it->second;
    }
};

std::size_t knot_vector(In& in, Out& o, int& n_out) {
    const int deg = in.i32(), n = in.i32();
    const std::uint8_t z = in.u8();  // bit 2: x = 0 (omitted); bit 3: multiplicity byte omitted
    if (z & ~0x0cu) bad("knot-vector form byte", in.at);
    const int x = z & 4 ? 0 : in.i32();
    if (n < 2 || n > 100000 || deg < 1 || deg > 30) bad("knot vector size", in.at);
    o.i32(deg); o.i32(0); o.i32(x); o.i32(n); o.i32(0);
    for (int i = 0; i < n; ++i) o.f64(in.f64());
    // multiplicities: standard (ends degree + 1, interior degree - 2) when bit 3 is set or a byte 1 follows;
    // otherwise listed: raw-integer streams after a byte 0, compact streams directly
    std::vector<int> mult;
    bool standard = (z & 8) != 0;
    if (!standard) {
        in.need(1);
        if (in.b[in.k] == 1) { ++in.k; standard = true; }
        else {
            if (!in.ci && in.u8() != 0) bad("knot multiplicity form", in.at);
            for (int i = 0; i < n; ++i) mult.push_back(in.i32());
        }
    }
    if (standard) {
        if (n > 2 && deg < 3) throw DecodeError(ErrorCode::Unsupported, "interior knot multiplicity of degree " + std::to_string(deg) + " not established", in.at);
        for (int i = 0; i < n; ++i) mult.push_back(i == 0 || i == n - 1 ? deg + 1 : deg - 2);
    }
    std::size_t sum = 0;
    for (const int m : mult) {
        if (m < 1 || m > deg + 1) bad("knot multiplicity", in.at);
        o.i32(m);
        sum += std::size_t(m);
    }
    if (sum < std::size_t(deg) + 2) bad("knot multiplicities", in.at);
    n_out = n;
    return sum - std::size_t(deg) - 1;
}

// 216 fillet: int n, knot vector (n knots, degree 5), sections P[n][10] D1[n][10] D2[n][10] (f64), then a
// tail whose first four f64 are the parameter box u0 u1 v0 v1 (the rest is kept as written; not decoded).
void fillet(In& in, Out& o) {
    const int n = in.i32();
    o.i32(n);
    int nk = 0;
    knot_vector(in, o, nk);
    if (nk != n) bad("fillet section count differs from its knots", in.at);
    for (int i = 0; i < 30 * n; ++i) o.f64(in.f64());
    for (int i = 0; i < 4; ++i) o.f64(in.f64());
    o.raw(in.raw(in.left()));
}

// 211 revolution: o X Y Z (12 f64), box (4 f64), int, int, k, 3 f64, int, f64 (cgm_model.hpp).
void revolution(In& in, Out& o) {
    for (int i = 0; i < 16; ++i) o.f64(in.f64());
    o.i32(in.i32()); o.i32(in.i32());
    for (int i = 0; i < 4; ++i) o.f64(in.f64());
    o.i32(in.i32());
    o.f64(in.f64());
}

// Edge curves (191, 192, 193). A map-1 entry may omit its identity range; the record length decides
// (exactly one assignment must consume the body).
void edge_curve(In& in, Out& o, std::size_t n_refs) {
    (void)n_refs;
    const std::uint8_t h = in.u8();
    const std::uint32_t n = h & 8 ? 1 : in.count();
    if (n == 0 || n > 16) bad("edge-curve map count", in.at);
    const std::size_t k0 = in.k;
    std::optional<std::pair<std::size_t, std::vector<std::uint8_t>>> sol;
    int sols = 0;
    for (std::uint32_t opt = 0; opt < (1u << n); ++opt) {
        try {
            in.k = k0;
            Out t;
            t.count(n);
            bool tail = false;
            for (std::uint32_t j = 0; j < n; ++j) {
                const int m = h & 4 ? 0 : in.i16();
                const bool last = j + 1 == n;
                if (m != 1 && (opt >> j & 1)) throw 0;
                const double a0 = in.f64(), a1 = in.f64();
                if (m == 0) {
                    const double tol = last && (h & 1) ? 0.0 : in.f64();
                    t.i16(m); t.f64(a0); t.f64(a1); t.f64(tol);
                    continue;
                }
                const double tol = in.f64();
                t.i16(m); t.f64(a0); t.f64(a1); t.f64(tol);
                if (m == 1) {
                    if (opt >> j & 1) { t.f64(a0); t.f64(a1); t.f64(tol); continue; }
                    const double b0 = in.f64(), b1 = in.f64();
                    const double tb = last && (h & 1) && !tail ? 0.0 : in.f64();
                    t.f64(b0); t.f64(b1); t.f64(tb);
                    continue;
                }
                tail = true;
                t.f64(in.f64()); t.f64(in.f64());
                if (m == 3) { t.f64(in.f64()); t.f64(in.f64()); }
                else if (m == 5) {
                    const std::uint32_t N = in.count();
                    if (N > 100000) throw 0;
                    t.count(N);
                    for (std::uint32_t i = 0; i < 2 * N; ++i) t.f64(in.f64());
                } else if (m != 2) throw 0;
            }
            if (tail) { t.f64(in.f64()); t.f64(in.f64()); t.f64(h & 1 ? 0.0 : in.f64()); }
            if (in.left() == 8) in.raw(8);  // two i32 (raw-integer streams; not decoded)
            if (in.left() != 0) throw 0;
            ++sols;
            sol.emplace(in.k, std::move(t.p));
        } catch (int) {
        } catch (const DecodeError&) {
        }
    }
    if (sols != 1) bad(sols ? "edge-curve layout ambiguous" : "edge-curve layout", in.at);
    in.k = sol->first;
    o.p = std::move(sol->second);
    for (int i = 0; i < 4; ++i) o.u8(0);
}

// 2012+ edge curves: header h, byte x (0x80 | map count << 2 | bit 0: closing t0 = first t0 | bit 1:
// closing t1 = first t1; one curve with map 0 or 1: the first range is the p-curve's range), maps.
void edge_curve_derived(In& in, Out& o, const std::vector<std::uint32_t>& refs, Ctx& cx, std::uint32_t self) {
    const std::uint8_t h = in.u8(), x = in.u8();  // x = 0x80 | maps << 2 | closing derivations
    const std::uint32_t n = x >> 2 & 0x1f;
    if (!(x & 0x80) || n == 0 || ((h & 8) && n != 1)) bad("edge-curve derivation byte " + std::to_string(x), in.at);
    if (n == 0 || n > 16) bad("edge-curve map count", in.at);
    o.count(n);
    bool tail = false;
    double f0 = 0, f1 = 0;
    for (std::uint32_t j = 0; j < n; ++j) {
        const int m = h & 4 ? 0 : in.i16();
        const bool last = j + 1 == n;
        double a0, a1;
        if (refs.size() == 1 && m == 0 && (x & 3)) {
            const auto& pr = cx.curve_range(refs.at(0), in.at);
            a0 = x & 1 ? pr[0] : in.f64();
            a1 = x & 2 ? pr[1] : in.f64();
        } else {
            a0 = in.f64();
            a1 = in.f64();
        }
        if (j == 0) { f0 = a0; f1 = a1; }
        if (m == 0) {
            const double tol = last && (h & 1) ? 0.0 : in.f64();
            o.i16(m); o.f64(a0); o.f64(a1); o.f64(tol);
            continue;
        }
        const double tol = in.f64();
        o.i16(m); o.f64(a0); o.f64(a1); o.f64(tol);
        if (m == 1) {  // identity: x bits 0 / 1 take its range from the first
            const double b0 = x & 1 ? a0 : in.f64(), b1 = x & 2 ? a1 : in.f64();
            o.f64(b0); o.f64(b1); o.f64(last && (h & 1) && !tail ? 0.0 : in.f64());
            continue;
        }
        tail = true;
        o.f64(in.f64()); o.f64(in.f64());
        if (m == 3) { o.f64(in.f64()); o.f64(in.f64()); }
        else if (m == 5) {
            const std::uint32_t N = in.count();
            if (N > 100000) bad("edge-curve table size", in.at);
            o.count(N);
            for (std::uint32_t i = 0; i < 2 * N; ++i) o.f64(in.f64());
        } else if (m != 2) bad("edge-curve map " + std::to_string(m), in.at);
    }
    if (tail) { o.f64(x & 1 ? f0 : in.f64()); o.f64(x & 2 ? f1 : in.f64()); o.f64(h & 1 ? 0.0 : in.f64()); }
    for (int i = 0; i < 4; ++i) o.u8(0);
    // points take start / end from this range; when it is the infinite default (|t| >= 1e6), from
    // the first curve's own range
    auto it = cx.range.find(refs.at(0));
    const bool infinite = std::abs(f0) >= 1e6 || std::abs(f1) >= 1e6;
    // an edge curve of one edge whose range is wider than its first curve's (a full circle period
    // around a half): the edge spans the first curve's range
    bool wider = false;
    if (it != cx.range.end()) {
        auto u = cx.edge_uses.find(self);
        wider = u != cx.edge_uses.end() && u->second == 1 && std::abs(f1 - f0) > std::abs(it->second[1] - it->second[0]) + 1e-9;
    }
    cx.range[self] = (infinite || wider) && it != cx.range.end() ? it->second : std::array<double, 2>{f0, f1};
}

}  // namespace

std::optional<CompactLayout> compact_layout(const Stream& s) {
    if (s.dialect != Dialect::Compact || s.features.size() < 2) return std::nullopt;
    // 2009Fall streams (levels 83..85): feature byte 0 is 0b (compact integers) or 01 (raw)
    if (s.level >= 83 && s.level <= 85) {
        CompactLayout L;
        L.compact_ints = (s.features[0] & 8) != 0;
        L.fall2009 = true;
        return L;
    }
    const std::uint8_t f0 = s.features[0], f1 = s.features[1];
    // 2011Spring streams: feature bytes b5|bf 75|77 (bit 3 of byte 0: compact integers; bit 1 of
    // byte 1: edges without their edge-curve reference).
    if ((f0 == 0xb5 || f0 == 0xbf) && (f1 == 0x75 || f1 == 0x77)) return CompactLayout{(f0 & 8) != 0, (f1 & 2) == 0, false};
    // Streams with a leading format byte, then b5 75: 06 (2011Fall) and 08 (2012+, container GUIDs follow).
    // Both write compact integers and keep the edges' edge-curve reference; 08 also omits derived values.
    if ((f0 == 0x06 || f0 == 0x08) && s.features.size() >= 3 && f1 == 0xb5 && s.features[2] == 0x75)
        return CompactLayout{true, true, f0 == 0x08};
    return std::nullopt;
}

namespace {
Canonical canonical_in(const Record& r, const CompactLayout& L, Ctx* cx) {
    In in{r.body, r.id, L.compact_ints, r.offset};
    Out o;
    Canonical c;
    c.type = r.type;
    switch (r.type) {
        case 13: case 513: c.refs = in.refs(); break;
        case 14:
            if (L.derived) {
                // refs, code (cgm_compact.hpp); canonical: refs[0] (an edge curve or p-curve) and its parameter
                const auto refs = in.refs();
                if (refs.empty()) bad("point without curves", r.offset);
                const std::uint8_t cd = in.u8();
                // codes per ref: refs[0] (and refs[1] when bit 7) bits 2-3; then refs[1] bits 4-5 (no
                // bit 7) or refs[2] bits 0-1 (bit 7): 2 start, 1 end of the curve's range, 0 stored
                std::vector<std::pair<std::uint32_t, double>> qs;
                auto at_end = [&](std::uint32_t ref, int code) {
                    const auto& rg = cx->curve_range(ref, r.offset);
                    return code == 2 ? rg[0] : rg[1];
                };
                const int c0 = cd >> 2 & 3;
                if (c0 == 3) bad("point parameter code 3", r.offset);
                const double q0 = c0 == 0 ? in.f64() : at_end(refs[0], c0);
                qs.push_back({refs[0], q0});
                std::size_t next = 1;
                if (cd & 0x80) { if (refs.size() > 1) { qs.push_back({refs[1], c0 == 0 ? q0 : at_end(refs[1], c0)}); } next = 2; }
                const int c1 = cd & 0x80 ? (cd & 3) : (cd >> 4 & 3);
                if (next < refs.size() && (c1 == 1 || c1 == 2) && cx->range.count(refs[next])) qs.push_back({refs[next], at_end(refs[next], c1)});
                in.k = in.b.size;  // stored parameters of the other refs (not decoded)
                o.count(std::uint32_t(qs.size()));
                for (const auto& [ref, q] : qs) { c.refs.push_back(ref); o.f64(q); o.i32(1); }
                break;
            }
            c.refs = in.refs();
            o.count(std::uint32_t(c.refs.size()));
            for (std::size_t i = 0; i < c.refs.size(); ++i) { o.f64(in.f64()); o.i32(in.i32()); }
            break;
        case 510: c.refs = in.refs(1); o.u8(in.u8()); break;
        case 511:
            c.refs = in.refs(L.edge_curve_ref ? 5 : 4);
            o.u8(in.u8());
            break;
        case 512: c.refs = in.refs(); o.i16(in.i16()); break;
        case 520: {
            c.refs = in.refs();
            const std::size_t n = c.refs.size();
            // 2012+: some shells carry a 0 byte before the bits (the record length decides)
            if (L.derived && in.left() == (4 + n + 7) / 8 + 2 && in.b[in.k] == 0) in.u8();
            const auto bt = in.bits(4 + n);
            const int code = bt[0] | bt[1] << 1 | bt[2] << 2 | bt[3] << 3;
            const int kind = code == 1 ? 1 : code == 3 ? 3 : code == 8 ? -1 : 0;
            if (!kind) bad("shell kind code " + std::to_string(code), r.offset);
            o.count(std::uint32_t(n));
            o.i16(kind);
            for (std::size_t j = 0; j < n; ++j) o.i16(pm(bt[4 + j]));
            o.i16(in.i16());
            break;
        }
        case 521: {
            c.refs = in.refs();
            const std::uint32_t n = in.count();
            if (2 * std::size_t(n) + 1 != c.refs.size()) bad("loop count differs from its references", r.offset);
            const std::uint8_t h = in.u8();
            int a, b;
            switch (h) {
                case 0x41: a = 1; b = 1; break;
                case 0x50: a = 1; b = -1; break;
                case 0xc1: a = -1; b = 1; break;
                case 0xd0: a = -1; b = -1; break;
                default: bad("loop header byte", r.offset);
            }
            const auto bt = in.bits(3 * std::size_t(n));
            const int cc = in.i16();
            o.count(n);
            o.i16(a);
            o.i16(b);
            o.i16(-1);  // the first coedge's x (always -1)
            for (std::size_t j = 0; j < bt.size(); ++j) o.i16(pm(bt[j]));
            o.i16(cc);
            if (cc) {
                o.raw(in.raw(32));
                for (int i = 0; i < 4; ++i) o.i32(in.i32());
                o.u8(in.u8());
                o.raw(in.raw(24));
            }
            break;
        }
        case 524: c.refs = in.refs(1); o.i16(in.i16()); break;
        case 525: c.refs = in.refs(); o.i32(in.i32()); o.i16(in.i16()); break;
        case 500: case 502: c.refs = in.refs(); in.k = in.b.size; break;  // payload not used by the decoder
        case 191: case 192: case 193:
            c.refs = in.refs();
            if (L.derived) edge_curve_derived(in, o, c.refs, *cx, r.id);
            else edge_curve(in, o, c.refs.size());
            break;
        case 150: {
            c.refs = in.refs(1);
            if (L.derived) {
                const std::uint8_t f = in.u8();
                auto box = [&]() -> const std::array<double, 4>& { return cx->surface_box(c.refs[0], r.offset); };
                double o0, o1, d0, d1, t0, t1;
                if (f & 3) {
                    const bool vert = f & 1;
                    if (f & 0x20) throw DecodeError(ErrorCode::Unsupported, "line form " + std::to_string(f) + " not decoded", r.offset);
                    const double cc = f & 16 ? box()[vert ? 0 : 2] : in.f64();
                    t0 = f & 4 ? box()[vert ? 2 : 0] : in.f64();
                    t1 = f & 8 ? (f & 4 ? box()[vert ? 3 : 1] : -t0) : in.f64();
                    o0 = vert ? cc : 0; o1 = vert ? 0 : cc; d0 = vert ? 0 : 1; d1 = vert ? 1 : 0;
                } else if (f & 0xc0) {
                    const bool vert = f & 64;
                    // bit 2: t = minus the box range along the line; bits 4 / 5: o on the box start / end
                    // of the other parameter (o along = 0); else stored
                    if (f & 4) { const auto& b = box(); t0 = -b[vert ? 3 : 1]; t1 = -b[vert ? 2 : 0]; }
                    else { t0 = in.f64(); t1 = f & 8 ? -t0 : in.f64(); }
                    if (f & 0x30) {
                        const double cc = box()[(vert ? 0 : 2) + (f & 0x20 ? 1 : 0)];
                        o0 = vert ? cc : 0; o1 = vert ? 0 : cc;
                    } else { o0 = in.f64(); o1 = in.f64(); }
                    d0 = vert ? 0 : -1; d1 = vert ? -1 : 0;
                } else {
                    if (f & 0x30) throw DecodeError(ErrorCode::Unsupported, "line form " + std::to_string(f) + " not decoded", r.offset);
                    // bit 2: no range stored (the infinite default)
                    if (f & 4) { t0 = -1e7; t1 = 1e7; } else { t0 = in.f64(); t1 = f & 8 ? -t0 : in.f64(); }
                    o0 = in.f64(); o1 = in.f64();
                    d0 = in.f64(); d1 = in.f64();
                }
                for (double v : {o0, o1, d0, d1, t0, t1}) o.f64(v);
                cx->range[r.id] = {t0, t1};
                break;
            }
            const int form = in.i16();
            double v[6];
            if (form == 0) { for (double& x : v) x = in.f64(); }
            else if (form == 1) { v[0] = in.f64(); v[1] = 0; v[2] = 0; v[3] = 1; v[4] = in.f64(); v[5] = in.f64(); }
            else if (form == 2) { v[0] = 0; v[1] = in.f64(); v[2] = 1; v[3] = 0; v[4] = in.f64(); v[5] = in.f64(); }
            else bad("line form " + std::to_string(form), r.offset);
            for (double x : v) o.f64(x);
            break;
        }
        case 151: {
            c.refs = in.refs(1);
            const double c0 = in.f64(), c1 = in.f64();
            // r t0 t1, flag byte (1: scale 1), [scale], offset; some objects of raw-integer streams
            // carry two i32 (1, 1) after the centre
            auto fits = [&](std::size_t k) {
                return k + 24 < in.b.size && in.b[k + 24] <= 1 && in.b.size - k == 25 + (in.b[k + 24] == 0 ? 8u : 0u) + 8;
            };
            // the two ints are raw i32 (raw-integer streams) or compact integers (05 05)
            const std::size_t iw = L.compact_ints ? 2 : 8;
            bool plain = fits(in.k), ints = fits(in.k + iw);
            const bool ones = L.compact_ints ? in.left() >= 2 && in.b[in.k] == 5 && in.b[in.k + 1] == 5
                                             : in.left() >= 8 && load_u32(in.b.data + in.k, false) == 1 && load_u32(in.b.data + in.k + 4, false) == 1;
            if (plain && ints) plain = !ones;
            if (plain == ints) bad("circle layout", r.offset);
            if (ints) in.raw(iw);
            o.f64(c0); o.f64(c1); o.i32(1); o.i32(1);
            const double rad = in.f64(), t0 = in.f64(), t1 = in.f64();
            const std::uint8_t f = in.u8();
            const double s = f == 1 ? 1.0 : in.f64();
            o.f64(rad); o.f64(t0); o.f64(t1); o.f64(s); o.f64(in.f64());
            if (cx) cx->range[r.id] = {t0, t1};
            break;
        }
        case 160: {
            c.refs = in.refs(1);
            int n = 0;
            knot_vector(in, o, n);
            if (in.i16() != n) bad("Hermite value count", r.offset);
            const int x = in.i16();
            if (load_u32(o.p.data() + 8, false) != std::uint32_t(x)) bad("Hermite x field", r.offset);
            o.i32(n);
            for (int i = 0; i < 4 * n; ++i) o.f64(in.f64());
            const int has2 = in.i16();
            o.i16(has2);
            if (has2) for (int i = 0; i < 2 * n; ++i) o.f64(in.f64());
            const double t0 = in.f64(), t1 = in.f64();
            o.f64(t0); o.f64(t1);
            if (cx) cx->range[r.id] = {t0, t1};
            if (in.i16() != -2) bad("Hermite closing value", r.offset);
            break;
        }
        case 231: {
            int nu_k = 0, nv_k = 0;
            const std::size_t nu = knot_vector(in, o, nu_k), nv = knot_vector(in, o, nv_k);
            const int rat = in.i16();
            o.i16(rat);
            const std::size_t m = nu * nv * (rat ? 4 : 3);
            if (m > (in.left() / 8)) bad("NURBS poles past the record", r.offset);
            for (std::size_t i = 0; i < m; ++i) o.f64(in.f64());
            // trailing block (not decoded; 152 or 154 bytes in 2011, compressed in 2012+); it ends with
            // two ints -2
            const ByteView blk = in.raw(in.left());
            if (L.derived && (blk.size < 2 || blk[blk.size - 1] != 0x07 || blk[blk.size - 2] != 0x07)) bad("NURBS closing block", r.offset);
            o.raw(blk);
            if (cx) {
                // parameter box: the knot ranges
                auto knot = [&](std::size_t kv_at, int which) {
                    const std::size_t n = load_u32(o.p.data() + kv_at + 12, false);
                    return load_f64_le(o.p.data() + kv_at + 20 + 8 * (which ? n - 1 : 0));
                };
                const std::size_t v_at = 20 + std::size_t(nu_k) * 12;
                cx->box[r.id] = {knot(0, 0), knot(0, 1), knot(v_at, 0), knot(v_at, 1)};
            }
            break;
        }
        case 200: {
            if (L.derived) {
                const std::int64_t v = in.sint();
                if (v < 0 || v >> 13) throw DecodeError(ErrorCode::Unsupported, "plane header form not decoded", r.offset);
                double w[9];
                auto vec = [&](double* x, int zero, bool unit) {
                    int free = 0;
                    for (int j = 0; j < 3; ++j) free += !(zero >> j & 1);
                    if (unit && free != 1) throw DecodeError(ErrorCode::Unsupported, "plane header form not decoded", r.offset);
                    for (int j = 0; j < 3; ++j) x[j] = zero >> j & 1 ? 0.0 : unit ? 1.0 : in.f64();
                };
                vec(w, int(v & 7), false);
                vec(w + 3, int(v >> 4 & 7), v >> 3 & 1);
                vec(w + 6, int(v >> 8 & 7), v >> 7 & 1);
                for (double x : w) o.f64(x);
                if (v >> 11 & 1) { o.f64(1); o.f64(1); } else { o.f64(in.f64()); o.f64(in.f64()); }
                std::array<double, 4> bx{-1e7, 1e7, -1e7, 1e7};
                if (!(v >> 12 & 1)) for (double& x : bx) x = in.f64();
                for (double x : bx) o.f64(x);
                auto nrm = [](const double* x) { return x[0] * x[0] + x[1] * x[1] + x[2] * x[2]; };
                if (std::abs(nrm(w + 3) - 1) > 1e-9 || std::abs(nrm(w + 6) - 1) > 1e-9 || std::abs(w[3] * w[6] + w[4] * w[7] + w[5] * w[8]) > 1e-9)
                    throw DecodeError(ErrorCode::Unsupported, "plane header form not decoded (frame not orthonormal)", r.offset);
                cx->box[r.id] = bx;
                break;
            }
            for (int i = 0; i < 3; ++i) o.f64(in.f64());
            const std::uint8_t m = in.u8();
            for (int j = 0; j < 6; ++j) o.f64(m >> j & 1 ? 0.0 : in.f64());
            if (m & 64) { o.f64(1); o.f64(1); } else { o.f64(in.f64()); o.f64(in.f64()); }
            // default box: the stream's infinity (1e6 or 1e7 by stream); the decoder does not use it
            if (m & 128) { o.f64(-1e7); o.f64(1e7); o.f64(-1e7); o.f64(1e7); } else for (int i = 0; i < 4; ++i) o.f64(in.f64());
            break;
        }
        case 201: {
            for (int i = 0; i < 3; ++i) o.f64(in.f64());
            const std::uint8_t m = in.u8();
            if (m >> 6) bad("cylinder mask", r.offset);
            for (int j = 0; j < 6; ++j) o.f64(m >> j & 1 ? 0.0 : in.f64());
            double v[5];
            for (double& x : v) { x = in.f64(); o.f64(x); }
            const std::uint8_t f = in.u8();
            if (f >> 3) bad("cylinder flags", r.offset);
            o.f64(f & 1 ? 1.0 : in.f64());
            o.f64(f & 2 ? 1.0 : in.f64());
            o.f64(f & 4 ? v[1] : in.f64());
            if (cx) cx->box[r.id] = {v[1], v[2], v[3], v[4]};
            break;
        }
        case 216: fillet(in, o); break;
        case 10: case 100: case 202: case 203: case 204: case 243:
            o.raw(in.raw(in.left()));
            if (cx) {
                auto at = [&](int i) { return load_f64_le(o.p.data() + 8 * i); };
                // boxes (u0 u1 v0 v1) in the parameters the p-curves use: cone as stored; sphere and
                // torus store angles, their parameters are arc lengths (u = R a, v = r b)
                if (r.type == 202 && o.p.size() == 184) cx->box[r.id] = {at(14), at(15), at(16), at(17)};
                if (r.type == 203 && o.p.size() >= 17 * 8) cx->box[r.id] = {at(13) * at(12), at(14) * at(12), at(15) * at(12), at(16) * at(12)};
                if (r.type == 204 && o.p.size() >= 18 * 8) cx->box[r.id] = {at(14) * at(12), at(15) * at(12), at(16) * at(13), at(17) * at(13)};
            }
            break;
        // the decoder uses the references and the leading doubles of these (cgm_model.hpp)
        case 211: c.refs = in.refs(1); revolution(in, o); break;
        case 12: case 155: case 210: case 212: c.refs = in.refs(1); o.raw(in.raw(in.left())); break;
        case 214: c.refs = in.refs(2); o.raw(in.raw(in.left())); break;
        case 241: case 242: case 244: case 251: c.refs = in.refs(); o.raw(in.raw(in.left())); break;
        default:
            throw DecodeError(ErrorCode::Unsupported, "compact class " + std::to_string(r.type) + " (code " + std::to_string(r.code) + ") not decoded", r.offset);
    }
    if (in.left()) bad(std::to_string(in.left()) + " bytes past the class grammar", r.offset);
    c.payload = std::move(o.p);
    return c;
}
}  // namespace

Canonical canonical(const Record& r, const CompactLayout& L) {
    if (L.fall2009) return canonical_2009(r, L.compact_ints);
    if (L.derived) throw DecodeError(ErrorCode::Unsupported, "record with derived values needs its stream (canonical_records)", r.offset);
    return canonical_in(r, L, nullptr);
}

std::vector<CanonicalRecord> canonical_records(const Stream& s, const CompactLayout& L) {
    std::vector<CanonicalRecord> out(s.records.size());
    if (L.fall2009) {
        for (std::size_t i = 0; i < s.records.size(); ++i) {
            const Record& r = s.records[i];
            try {
                out[i].c = canonical_2009(r, L.compact_ints);
            } catch (const DecodeError& e) {
                out[i].error = "class " + std::to_string(r.type) + " record " + std::to_string(r.id) + ": " + e.what();
            }
        }
        return out;
    }
    Ctx cx;
    for (const Record& r : s.records) cx.type[r.id] = r.type;
    if (L.derived)
        for (const Record& r : s.records) {
            if (r.type != 511 || r.body.size == 0) continue;
            try { In in{r.body, r.id, true, r.offset}; ++cx.edge_uses[in.ref()]; } catch (const DecodeError&) {}
        }
    auto pass = [](std::uint16_t t) {
        if ((t >= 200 && t < 300) || t == 10 || t == 12 || t == 100) return 0;  // surfaces, points in space
        if (t >= 150 && t < 191) return 1;                                     // p-curves
        if (t >= 191 && t <= 193) return 2;                                    // edge curves
        if (t == 14 || t == 13) return 3;                                      // points on curves
        return 4;
    };
    for (int p = 0; p < 5; ++p) {
        for (std::size_t i = 0; i < s.records.size(); ++i) {
            const Record& r = s.records[i];
            if (pass(r.type) != p) continue;
            try {
                out[i].c = canonical_in(r, L, L.derived ? &cx : nullptr);
            } catch (const DecodeError& e) {
                out[i].error = "class " + std::to_string(r.type) + " record " + std::to_string(r.id) + ": " + e.what();
            }
        }
        // carriers (212, 241, 242, 244) use their first reference's parameters, offsets (214) their
        // base's (refs[1])
        for (bool more = p == 0; more;) {
            more = false;
            for (std::size_t i = 0; i < s.records.size(); ++i) {
                const std::uint16_t t = s.records[i].type;
                const std::size_t j = t == 214 ? 1 : 0;
                if ((t != 212 && t != 214 && t != 241 && t != 242 && t != 244) || !out[i].error.empty() || out[i].c.refs.size() <= j || cx.box.count(s.records[i].id)) continue;
                if (auto it = cx.box.find(out[i].c.refs[j]); it != cx.box.end()) { cx.box[s.records[i].id] = it->second; more = true; }
            }
        }
        // extrusions (210): u along the profile (its range), v from the payload (d[3] v0 v1); lines on
        // them decode once the profile's range is known, so the p-curve pass runs again for them
        if (p == 1 && L.derived) {
            bool added = false;
            for (std::size_t i = 0; i < s.records.size(); ++i) {
                const Record& r = s.records[i];
                if (r.type != 210 || !out[i].error.empty() || out[i].c.refs.empty() || out[i].c.payload.size() < 40 || cx.box.count(r.id)) continue;
                std::uint32_t prof = out[i].c.refs[0];
                if (auto t = cx.type.find(prof); t != cx.type.end() && t->second >= 191 && t->second <= 193) {
                    // an edge curve: its first curve's range
                    for (std::size_t m = 0; m < s.records.size(); ++m)
                        if (s.records[m].id == prof && out[m].error.empty() && !out[m].c.refs.empty()) prof = out[m].c.refs[0];
                }
                auto rg = cx.range.find(prof);
                if (rg == cx.range.end()) continue;
                const std::uint8_t* q = out[i].c.payload.data();
                cx.box[r.id] = {rg->second[0], rg->second[1], load_f64_le(q + 24), load_f64_le(q + 32)};
                added = true;
            }
            if (added)
                for (std::size_t i = 0; i < s.records.size(); ++i) {
                    const Record& r = s.records[i];
                    if (pass(r.type) != 1 || out[i].error.empty()) continue;
                    try {
                        out[i].c = canonical_in(r, L, &cx);
                        out[i].error.clear();
                    } catch (const DecodeError& e) {
                        out[i].error = "class " + std::to_string(r.type) + " record " + std::to_string(r.id) + ": " + e.what();
                    }
                }
        }
    }
    return out;
}

}  // namespace cadcore::formats::catia::cgm
