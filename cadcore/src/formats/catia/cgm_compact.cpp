#include "formats/catia/cgm_compact.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <string>
#include <unordered_map>

#include "cadcore/error.hpp"
#include "formats/catia/cgm_compact_io.hpp"

namespace cadcore::formats::catia::cgm {
namespace {
using namespace compact_io;

using Box = std::array<double, 4>;    // u0 u1 v0 v1
using Range = std::array<double, 2>;  // t0 t1

int pm(int bit) { return bit ? 1 : -1; }
bool is_edge_curve(std::uint16_t t) { return t >= 191 && t <= 193; }

// Values 2012+ records take from others: surface parameter boxes, curve ranges.
struct Ctx {
    std::unordered_map<std::uint32_t, int> edge_uses;  // edge curve -> edges referencing it
    std::unordered_map<std::uint32_t, std::uint16_t> type;
    std::unordered_map<std::uint32_t, Box> box;
    std::unordered_map<std::uint32_t, Range> range;
    std::unordered_map<std::uint32_t, std::uint32_t> first;  // edge curve -> its first curve
    std::uint16_t curve_type(std::uint32_t id) const {  // through edge curves
        for (int d = 0; d < 8; ++d) {
            auto f = first.find(id);
            if (f == first.end()) break;
            id = f->second;
        }
        auto t = type.find(id);
        return t == type.end() ? 0 : t->second;
    }

    template <class T>
    const T& get(const std::unordered_map<std::uint32_t, T>& m, std::uint32_t id, const char* what, const char* missing, std::size_t at) const {
        if (auto it = m.find(id); it != m.end()) return it->second;
        auto t = type.find(id);
        throw DecodeError(ErrorCode::Unsupported, std::string(what) + std::to_string(t == type.end() ? 0 : t->second) + missing, at);
    }
    const Box& surface_box(std::uint32_t id, std::size_t at) const { return get(box, id, "parameter box of surface class ", " not established", at); }
    const Range& curve_range(std::uint32_t id, std::size_t at) const { return get(range, id, "range of curve class ", " not available", at); }
};

struct Knots { int n = 0; std::size_t poles = 0; Range range{}; };

// int degree, int n, u8 z (bit 2: x = 0 omitted, else int x; bit 3: standard multiplicities), f64
// knots[n], multiplicities: standard (ends degree + 1, interior degree - 2) when bit 3 or a byte 1,
// else listed (after a byte 0 in raw-integer streams). Canonical: i32 degree 0 x n 0, knots, i32 mults.
Knots knot_vector(In& in, Out& o) {
    const int deg = in.i32(), n = in.i32();
    const std::uint8_t z = in.u8();
    if (z & ~0x0cu) bad("knot-vector form byte", in.at);
    const int x = z & 4 ? 0 : in.i32();
    if (n < 2 || n > 100000 || deg < 1 || deg > 30) bad("knot vector size", in.at);
    o.i32(deg); o.i32(0); o.i32(x); o.i32(n); o.i32(0);
    Knots kv;
    kv.n = n;
    for (int i = 0; i < n; ++i) {
        const double t = in.f64();
        o.f64(t);
        if (i == 0) kv.range[0] = t;
        kv.range[1] = t;
    }
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
    kv.poles = sum - std::size_t(deg) - 1;
    return kv;
}

// 216 (m = 10), 243 (m = 6): int n, knot vector, P D1 D2 [n][m], box, rest kept. Returns the box.
Box sections(In& in, Out& o, int m) {
    const int n = in.i32();
    o.i32(n);
    if (knot_vector(in, o).n != n) bad("section count differs from its knots", in.at);
    f64s(in, o, 3 * std::size_t(m) * std::size_t(n));
    Box b;
    for (double& x : b) { x = in.f64(); o.f64(x); }
    o.raw(in.raw(in.left()));
    return b;
}

// Poles (dim f64 each, [weights]) after a knot vector; `rational` read for 3D curves only.
void poles(In& in, Out& o, std::size_t np, int dim, bool rational, std::size_t at, const char* what) {
    const std::size_t m = np * std::size_t(dim + (rational ? 1 : 0));
    if (m > in.left() / 8) bad(std::string(what) + " poles past the record", at);
    f64s(in, o, m);
}

// 155 item: u8 0x41, int code, [1, 4: varint inputs], values (1: 3, 2: 4, 3: 0, 4: 3, 7: 5), children
// (3: 2, 4: one per input). Canonical: i16 code, [i16 inputs, i16 0], values, children.
void expr_item(In& in, Out& o, int depth) {
    if (depth > 32) bad("expression too deep", in.at);
    if (in.u8() != 0x41) bad("expression item header", in.at);
    const int code = in.i16();
    o.i16(code);
    int kids = 0, vals = 0;
    switch (code) {
        case 1: case 4: {
            const std::uint32_t n = in.count();
            if (n == 0 || n > 8) bad("expression input count", in.at);
            o.i16(std::int32_t(n));
            o.i16(0);
            vals = 3;
            kids = code == 4 ? int(n) : 0;
            break;
        }
        case 2: vals = 4; break;
        case 3: kids = 2; break;
        case 7: vals = 5; break;
        default: throw DecodeError(ErrorCode::Unsupported, "expression item code " + std::to_string(code) + " not decoded", in.at);
    }
    f64s(in, o, std::size_t(vals));
    for (int i = 0; i < kids; ++i) expr_item(in, o, depth + 1);
}

// One map of an edge curve after its code and range: identity (1), affine (2), scale (3), table (5).
// Returns whether the map needs the closing range.
template <class F64>
bool edge_map(In& in, Out& o, int m, F64 f64) {
    if (m == 0) return false;
    o.f64(f64()); o.f64(f64());
    if (m == 3) { o.f64(in.f64()); o.f64(in.f64()); }
    else if (m == 5) {
        const std::uint32_t N = in.count();
        if (N > 100000) bad("edge-curve table size", in.at);
        o.count(N);
        f64s(in, o, 2 * std::size_t(N));
    } else if (m != 2) bad("edge-curve map " + std::to_string(m), in.at);
    return true;
}

// Edge curves (2011): header h (bit 3: one map; bit 2: map 0, code omitted; bit 0: last tolerance 0,
// omitted). A map 1 may omit its identity range: exactly one assignment must consume the body.
void edge_curve(In& in, Out& o) {
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
                const double tol = m == 0 && last && (h & 1) ? 0.0 : in.f64();
                t.i16(m); t.f64(a0); t.f64(a1); t.f64(tol);
                if (m == 1) {
                    if (opt >> j & 1) { t.f64(a0); t.f64(a1); t.f64(tol); continue; }
                    t.f64(in.f64()); t.f64(in.f64());
                    t.f64(last && (h & 1) && !tail ? 0.0 : in.f64());
                    continue;
                }
                tail |= edge_map(in, t, m, [&] { return in.f64(); });
            }
            if (tail) { t.f64(in.f64()); t.f64(in.f64()); t.f64(h & 1 ? 0.0 : in.f64()); }
            if (in.left() == 8) in.raw(8);  // two i32 in raw-integer streams
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

// Edge curves (2012+): header h, byte x = 0x80 | maps << 2 | bit 0: closing t0 = first t0 | bit 1:
// closing t1 = first t1 (one curve with map 0: the first range is that curve's range).
void edge_curve_derived(In& in, Out& o, const std::vector<std::uint32_t>& refs, Ctx& cx, std::uint32_t self) {
    const std::uint8_t h = in.u8(), x = in.u8();
    const std::uint32_t n = x >> 2 & 0x1f;
    if (!(x & 0x80) || n == 0 || ((h & 8) && n != 1)) bad("edge-curve derivation byte " + std::to_string(x), in.at);
    if (n > 16) bad("edge-curve map count", in.at);
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
        const double tol = m == 0 && last && (h & 1) ? 0.0 : in.f64();
        o.i16(m); o.f64(a0); o.f64(a1); o.f64(tol);
        if (m == 1) {  // x bits 0 / 1: the identity's range is the first
            const double b0 = x & 1 ? a0 : in.f64(), b1 = x & 2 ? a1 : in.f64();
            o.f64(b0); o.f64(b1); o.f64(last && (h & 1) && !tail ? 0.0 : in.f64());
            continue;
        }
        tail |= edge_map(in, o, m, [&] { return in.f64(); });
    }
    if (tail) { o.f64(x & 1 ? f0 : in.f64()); o.f64(x & 2 ? f1 : in.f64()); o.f64(h & 1 ? 0.0 : in.f64()); }
    for (int i = 0; i < 4; ++i) o.u8(0);
    // points take start / end from this range, or from the first curve's when this one is the infinite
    // default (|t| >= 1e6) or, for an edge curve of one edge, wider than it (a full period around a half)
    auto it = cx.range.find(refs.at(0));
    const bool infinite = std::abs(f0) >= 1e6 || std::abs(f1) >= 1e6;
    bool wider = false;
    if (it != cx.range.end()) {
        auto u = cx.edge_uses.find(self);
        wider = u != cx.edge_uses.end() && u->second == 1 && std::abs(f1 - f0) > std::abs(it->second[1] - it->second[0]) + 1e-9;
    }
    cx.range[self] = (infinite || wider) && it != cx.range.end() ? it->second : Range{f0, f1};
    cx.first[self] = refs.at(0);
}

// 14 (2012+): refs, code byte; per ref a parameter: 2 start, 1 end of the curve's range, 0 stored f64
// and, unless the curves sharing it include a line or circle, int segment. refs[0] bits 2-3 (refs[1] too when bit 7);
// then refs[1] bits 4-5 or (bit 7) refs[2] bits 0-1. The other refs' fields are not decoded (their
// parameters follow from the edge curve's maps).
void point_derived(In& in, Out& o, Canonical& c, Ctx& cx, std::size_t at) {
    const auto refs = in.refs();
    if (refs.empty()) bad("point without curves", at);
    const std::uint8_t cd = in.u8();
    std::vector<std::pair<std::uint32_t, double>> qs;
    auto at_end = [&](std::uint32_t ref, int code) {
        const auto& rg = cx.curve_range(ref, at);
        return code == 2 ? rg[0] : rg[1];
    };
    const int c0 = cd >> 2 & 3;
    if (c0 == 3) bad("point parameter code 3", at);
    const double q0 = c0 == 0 ? in.f64() : at_end(refs[0], c0);
    auto line_or_circle = [&](std::uint32_t ref) { const std::uint16_t t = cx.curve_type(ref); return t == 150 || t == 151; };
    std::int32_t span = 1;
    if (c0 == 0 && !line_or_circle(refs[0]) && !((cd & 0x80) && refs.size() > 1 && line_or_circle(refs[1]))) span = in.i32();
    qs.push_back({refs[0], q0});
    std::size_t next = 1;
    if (cd & 0x80) { if (refs.size() > 1) qs.push_back({refs[1], c0 == 0 ? q0 : at_end(refs[1], c0)}); next = 2; }
    const int c1 = cd & 0x80 ? (cd & 3) : (cd >> 4 & 3);
    if (next < refs.size() && (c1 == 1 || c1 == 2) && cx.range.count(refs[next])) qs.push_back({refs[next], at_end(refs[next], c1)});
    in.k = in.b.size;
    o.count(std::uint32_t(qs.size()));
    for (std::size_t i = 0; i < qs.size(); ++i) {
        c.refs.push_back(qs[i].first);
        o.f64(qs[i].second);
        o.i32(i < ((cd & 0x80) ? 2u : 1u) ? span : 1);
    }
}

// 150 (2012+): ref, u8 f. Directions: bit 0 (0, 1), bit 1 (1, 0), bit 6 (0, -1), bit 7 (-1, 0), else
// stored. +v / +u: [c, t0, t1], o = (c, 0) / (0, c); bit 4: c = box start of the other parameter;
// bit 2: t0 = box start; bit 3: t1 = box end if t0 derived, else -t0. -v / -u: [t0, t1, o]; bit 2: t =
// minus the box range; bit 3: t1 = -t0; bits 4 / 5: o = box start / end. General: [t0, t1, o, d];
// bit 2: no range (the infinite default); bit 3: t1 = -t0.
void line_derived(In& in, Out& o, std::uint32_t surface, Ctx& cx, std::uint32_t self, std::size_t at) {
    const std::uint8_t f = in.u8();
    auto box = [&]() -> const Box& { return cx.surface_box(surface, at); };
    double o0, o1, d0, d1, t0, t1;
    if (f & 3) {
        const bool vert = f & 1;
        if (f & 0x20) throw DecodeError(ErrorCode::Unsupported, "line form " + std::to_string(f) + " not decoded", at);
        const double cc = f & 16 ? box()[vert ? 0 : 2] : in.f64();
        t0 = f & 4 ? box()[vert ? 2 : 0] : in.f64();
        t1 = f & 8 ? (f & 4 ? box()[vert ? 3 : 1] : -t0) : in.f64();
        o0 = vert ? cc : 0; o1 = vert ? 0 : cc; d0 = vert ? 0 : 1; d1 = vert ? 1 : 0;
    } else if (f & 0xc0) {
        const bool vert = f & 64;
        if (f & 4) { const auto& b = box(); t0 = -b[vert ? 3 : 1]; t1 = -b[vert ? 2 : 0]; }
        else { t0 = in.f64(); t1 = f & 8 ? -t0 : in.f64(); }
        if (f & 0x30) {
            const double cc = box()[(vert ? 0 : 2) + (f & 0x20 ? 1 : 0)];
            o0 = vert ? cc : 0; o1 = vert ? 0 : cc;
        } else { o0 = in.f64(); o1 = in.f64(); }
        d0 = vert ? 0 : -1; d1 = vert ? -1 : 0;
    } else {
        if (f & 0x30) throw DecodeError(ErrorCode::Unsupported, "line form " + std::to_string(f) + " not decoded", at);
        if (f & 4) { t0 = -1e7; t1 = 1e7; } else { t0 = in.f64(); t1 = f & 8 ? -t0 : in.f64(); }
        o0 = in.f64(); o1 = in.f64();
        d0 = in.f64(); d1 = in.f64();
    }
    for (double v : {o0, o1, d0, d1, t0, t1}) o.f64(v);
    cx.range[self] = {t0, t1};
}

// 200 (2012+): int v; bits 0-2: origin component zero; 3: U is the unit vector of its one free
// component; 4-6: U component zero; 7, 8-10: the same for V; 11: scale (1, 1); 12: default box.
void plane_derived(In& in, Out& o, Ctx& cx, std::uint32_t self, std::size_t at) {
    const std::int64_t v = in.sint();
    if (v < 0 || v >> 13) throw DecodeError(ErrorCode::Unsupported, "plane header form not decoded", at);
    double w[9];
    auto vec = [&](double* x, int zero, bool unit) {
        int free = 0;
        for (int j = 0; j < 3; ++j) free += !(zero >> j & 1);
        if (unit && free != 1) throw DecodeError(ErrorCode::Unsupported, "plane header form not decoded", at);
        for (int j = 0; j < 3; ++j) x[j] = zero >> j & 1 ? 0.0 : unit ? 1.0 : in.f64();
    };
    vec(w, int(v & 7), false);
    vec(w + 3, int(v >> 4 & 7), v >> 3 & 1);
    vec(w + 6, int(v >> 8 & 7), v >> 7 & 1);
    for (double x : w) o.f64(x);
    if (v >> 11 & 1) { o.f64(1); o.f64(1); } else { o.f64(in.f64()); o.f64(in.f64()); }
    Box bx{-1e7, 1e7, -1e7, 1e7};
    if (!(v >> 12 & 1)) for (double& x : bx) x = in.f64();
    for (double x : bx) o.f64(x);
    auto nrm = [](const double* x) { return x[0] * x[0] + x[1] * x[1] + x[2] * x[2]; };
    if (std::abs(nrm(w + 3) - 1) > 1e-9 || std::abs(nrm(w + 6) - 1) > 1e-9 || std::abs(w[3] * w[6] + w[4] * w[7] + w[5] * w[8]) > 1e-9)
        throw DecodeError(ErrorCode::Unsupported, "plane header form not decoded (frame not orthonormal)", at);
    cx.box[self] = bx;
}

Canonical canonical_in(const Record& r, const CompactLayout& L, Ctx* cx) {
    In in{r.body, r.id, L.compact_ints, r.offset};
    Out o;
    Canonical c;
    c.type = r.type;
    switch (r.type) {
        case 13: case 513: c.refs = in.refs(); break;
        case 14:
            if (L.derived) { point_derived(in, o, c, *cx, r.offset); break; }
            c.refs = in.refs();
            o.count(std::uint32_t(c.refs.size()));
            for (std::size_t i = 0; i < c.refs.size(); ++i) { o.f64(in.f64()); o.i32(in.i32()); }
            break;
        case 510: c.refs = in.refs(1); o.u8(in.u8()); break;
        case 511: c.refs = in.refs(L.edge_curve_ref ? 5 : 4); o.u8(in.u8()); break;
        case 512: c.refs = in.refs(); o.i16(in.i16()); break;
        case 520: {
            c.refs = in.refs();
            const std::size_t n = c.refs.size();
            if (L.derived && in.left() == (4 + n + 7) / 8 + 2 && in.b[in.k] == 0) in.u8();  // optional 0 byte
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
            int a, b;
            switch (in.u8()) {
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
            o.i16(-1);
            for (const int bit : bt) o.i16(pm(bit));
            o.i16(cc);
            if (cc) {  // 4 f64, 4 ints (absent in 2017+: 25 bytes follow; written as 0), u8, 24 bytes
                o.raw(in.raw(32));
                const bool ints = in.left() != 25;
                for (int i = 0; i < 4; ++i) o.i32(ints ? in.i32() : 0);
                o.u8(in.u8());
                o.raw(in.raw(24));
            }
            break;
        }
        case 524: c.refs = in.refs(1); o.i16(in.i16()); break;
        case 525: c.refs = in.refs(); o.i32(in.i32()); o.i16(in.i16()); break;
        case 500: case 502: c.refs = in.refs(); in.k = in.b.size; break;
        case 191: case 192: case 193:
            c.refs = in.refs();
            if (L.derived) edge_curve_derived(in, o, c.refs, *cx, r.id);
            else edge_curve(in, o);
            break;
        case 150: {
            c.refs = in.refs(1);
            if (L.derived) { line_derived(in, o, c.refs[0], *cx, r.id, r.offset); break; }
            const int form = in.i16();  // 0: o d t0 t1; 1: o.x t0 t1, d = (0, 1); 2: o.y t0 t1, d = (1, 0)
            double v[6];
            if (form == 0) { for (double& x : v) x = in.f64(); }
            else if (form == 1) { v[0] = in.f64(); v[1] = 0; v[2] = 0; v[3] = 1; v[4] = in.f64(); v[5] = in.f64(); }
            else if (form == 2) { v[0] = 0; v[1] = in.f64(); v[2] = 1; v[3] = 0; v[4] = in.f64(); v[5] = in.f64(); }
            else bad("line form " + std::to_string(form), r.offset);
            for (double x : v) o.f64(x);
            break;
        }
        case 151: {
            // c, [two ints 1 1 (some objects)], r t0 t1, u8 (1: s = 1), [s], offset
            c.refs = in.refs(1);
            const double c0 = in.f64(), c1 = in.f64();
            auto fits = [&](std::size_t k) {
                return k + 24 < in.b.size && in.b[k + 24] <= 1 && in.b.size - k == 25 + (in.b[k + 24] == 0 ? 8u : 0u) + 8;
            };
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
            // knot vector, i16 n, i16 x, P D1, int has2, [D2], t0 t1, i16 -2, [2018Fall+: byte 0]
            c.refs = in.refs(1);
            const int n = knot_vector(in, o).n;
            if (in.i16() != n) bad("Hermite value count", r.offset);
            if (load_u32(o.p.data() + 8, false) != std::uint32_t(in.i16())) bad("Hermite x field", r.offset);
            o.i32(n);
            f64s(in, o, 4 * std::size_t(n));
            const int has2 = in.i16();
            o.i16(has2);
            if (has2) f64s(in, o, 2 * std::size_t(n));
            const double t0 = in.f64(), t1 = in.f64();
            o.f64(t0); o.f64(t1);
            if (cx) cx->range[r.id] = {t0, t1};
            if (in.i16() != -2) bad("Hermite closing value", r.offset);
            if (in.left() == 1 && in.b[in.k] == 0) in.u8();
            break;
        }
        case 130: {  // no refs: knot vector, int rational, poles[3], [weights], tail
            const Knots kv = knot_vector(in, o);
            const int rat = in.i16();
            o.i16(rat);
            poles(in, o, kv.poles, 3, rat != 0, r.offset, "B-spline curve");
            bspline_tail(in, o);
            if (cx) cx->range[r.id] = kv.range;
            break;
        }
        case 170: {
            // ref, i16 0, knot vector, poles[2], tail; 2012+: form byte 0x1e and f64 L (the knot range's
            // length) instead of the i16, the tail implied (1, n - 1, 0 L 1 0, 0)
            c.refs = in.refs(1);
            double span = 0;
            if (L.derived) {
                if (in.u8() != 0x1e) throw DecodeError(ErrorCode::Unsupported, "B-spline p-curve form not decoded", r.offset);
                span = in.f64();
            } else if (in.i16() != 0) {
                bad("B-spline p-curve form", r.offset);
            }
            o.i16(0);
            const Knots kv = knot_vector(in, o);
            poles(in, o, kv.poles, 2, false, r.offset, "B-spline p-curve");
            if (L.derived) { o.i16(1); o.i16(kv.n - 1); o.f64(0); o.f64(span); o.f64(1); o.f64(0); o.u8(0); }
            else bspline_tail(in, o);
            if (cx) cx->range[r.id] = kv.range;
            break;
        }
        case 231: {
            // knot vectors u v, int rational, poles (v fastest), [weights], a block (not decoded; 2012+:
            // ends with two ints -2)
            const Knots u = knot_vector(in, o), v = knot_vector(in, o);
            const int rat = in.i16();
            o.i16(rat);
            poles(in, o, u.poles * v.poles, 3, rat != 0, r.offset, "NURBS");
            const ByteView blk = in.raw(in.left());
            if (L.derived && (blk.size < 2 || blk[blk.size - 1] != 0x07 || blk[blk.size - 2] != 0x07)) bad("NURBS closing block", r.offset);
            o.raw(blk);
            if (cx) cx->box[r.id] = {u.range[0], u.range[1], v.range[0], v.range[1]};
            break;
        }
        case 200: {
            if (L.derived) { plane_derived(in, o, *cx, r.id, r.offset); break; }
            // o, u8 m (bits 0-5: U V components zero; 6: scale (1, 1); 7: default box), the rest
            f64s(in, o, 3);
            const std::uint8_t m = in.u8();
            for (int j = 0; j < 6; ++j) o.f64(m >> j & 1 ? 0.0 : in.f64());
            if (m & 64) { o.f64(1); o.f64(1); } else f64s(in, o, 2);
            if (m & 128) { o.f64(-1e7); o.f64(1e7); o.f64(-1e7); o.f64(1e7); } else f64s(in, o, 4);
            break;
        }
        case 201: {
            // o, u8 m (bits 0-5: X Y components zero), r u0 u1 v0 v1, u8 (bits: su = 1, sv = 1, reference = u0), the rest
            f64s(in, o, 3);
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
        case 216: case 243: {
            const Box b = sections(in, o, r.type == 216 ? 10 : 6);
            if (cx) cx->box[r.id] = b;
            break;
        }
        case 10: case 100: case 202: case 203: case 204:
            o.raw(in.raw(in.left()));
            if (cx) {
                // boxes in the p-curves' parameters: cone as stored; sphere and torus store angles, their
                // parameters are arc lengths (u = R a, v = r b)
                auto at = [&](int i) { return load_f64_le(o.p.data() + 8 * i); };
                if (r.type == 202 && o.p.size() == 184) cx->box[r.id] = {at(14), at(15), at(16), at(17)};
                if (r.type == 203 && o.p.size() >= 17 * 8) cx->box[r.id] = {at(13) * at(12), at(14) * at(12), at(15) * at(12), at(16) * at(12)};
                if (r.type == 204 && o.p.size() >= 18 * 8) cx->box[r.id] = {at(14) * at(12), at(15) * at(12), at(16) * at(13), at(17) * at(13)};
            }
            break;
        case 211:  // ref, o X Y Z, box, int, int, k, 3 f64, int, f64
            c.refs = in.refs(1);
            f64s(in, o, 16);
            o.i32(in.i32()); o.i32(in.i32());
            f64s(in, o, 4);
            o.i32(in.i32());
            o.f64(in.f64());
            break;
        case 155:  // ref, f64 t0 t1 a b, items u v
            c.refs = in.refs(1);
            f64s(in, o, 4);
            expr_item(in, o, 0);
            expr_item(in, o, 0);
            break;
        case 12: case 210: case 212: c.refs = in.refs(1); o.raw(in.raw(in.left())); break;
        case 214: c.refs = in.refs(2); o.raw(in.raw(in.left())); break;
        case 241: case 242: case 244: case 251: c.refs = in.refs(); o.raw(in.raw(in.left())); break;
        default:
            throw DecodeError(ErrorCode::Unsupported, "compact class " + std::to_string(r.type) + " (code " + std::to_string(r.code) + ") not decoded", r.offset);
    }
    if (in.left()) bad(std::to_string(in.left()) + " bytes past the class grammar", r.offset);
    c.payload = std::move(o.p);
    return c;
}

// Transcodes record i into out[i], or records why it failed.
template <class F>
void transcode(const Record& r, CanonicalRecord& out, F&& f) {
    try {
        out.c = f();
        out.error.clear();
    } catch (const DecodeError& e) {
        out.error = "class " + std::to_string(r.type) + " record " + std::to_string(r.id) + ": " + e.what();
    }
}

// 2012+ resolution order: values derived from other records need those records first.
int pass(std::uint16_t t) {
    if ((t >= 200 && t < 300) || t == 10 || t == 12 || (t >= 100 && t < 150)) return 0;  // surfaces, points and curves in space
    if (t >= 150 && t < 191) return 1;                                                   // p-curves
    if (is_edge_curve(t)) return 2;
    if (t == 14 || t == 13) return 3;
    return 4;
}

}  // namespace

std::optional<CompactLayout> compact_layout(const Stream& s) {
    if (s.dialect != Dialect::Compact || s.features.size() < 2) return std::nullopt;
    if (s.level >= 83 && s.level <= 85) {
        CompactLayout L;
        L.compact_ints = (s.features[0] & 8) != 0;
        L.fall2009 = true;
        return L;
    }
    const std::uint8_t f0 = s.features[0], f1 = s.features[1];
    if ((f0 == 0xb5 || f0 == 0xbf) && (f1 == 0x75 || f1 == 0x77)) return CompactLayout{(f0 & 8) != 0, (f1 & 2) == 0, false};
    if ((f0 == 0x06 || f0 == 0x08) && s.features.size() >= 3 && f1 == 0xb5 && s.features[2] == 0x75)
        return CompactLayout{true, true, f0 == 0x08};
    return std::nullopt;
}

Canonical canonical(const Record& r, const CompactLayout& L) {
    if (L.fall2009) return canonical_2009(r, L.compact_ints);
    if (L.derived) throw DecodeError(ErrorCode::Unsupported, "record with derived values needs its stream (canonical_records)", r.offset);
    return canonical_in(r, L, nullptr);
}

std::vector<CanonicalRecord> canonical_records(const Stream& s, const CompactLayout& L) {
    const auto& recs = s.records;
    std::vector<CanonicalRecord> out(recs.size());
    if (L.fall2009) {
        for (std::size_t i = 0; i < recs.size(); ++i) transcode(recs[i], out[i], [&] { return canonical_2009(recs[i], L.compact_ints); });
        return out;
    }
    Ctx cx;
    Ctx* dx = L.derived ? &cx : nullptr;
    for (const Record& r : recs) cx.type[r.id] = r.type;
    if (L.derived)
        for (const Record& r : recs) {
            if (r.type != 511 || r.body.size == 0) continue;
            try { In in{r.body, r.id, true, r.offset}; ++cx.edge_uses[in.ref()]; } catch (const DecodeError&) {}
        }
    auto run_pass = [&](int p, bool failed_only) {
        for (std::size_t i = 0; i < recs.size(); ++i)
            if (pass(recs[i].type) == p && (!failed_only || !out[i].error.empty()))
                transcode(recs[i], out[i], [&] { return canonical_in(recs[i], L, dx); });
    };
    for (int p = 0; p < 5; ++p) {
        run_pass(p, false);
        // carriers (212, 241, 242, 244) and 251 take their first reference's box, offsets (214) their base's
        for (bool more = p == 0; more;) {
            more = false;
            for (std::size_t i = 0; i < recs.size(); ++i) {
                const std::uint16_t t = recs[i].type;
                const std::size_t j = t == 214 ? 1 : 0;
                if ((t != 212 && t != 214 && t != 241 && t != 242 && t != 244 && t != 251) || !out[i].error.empty() || out[i].c.refs.size() <= j || cx.box.count(recs[i].id)) continue;
                if (auto it = cx.box.find(out[i].c.refs[j]); it != cx.box.end()) { cx.box[recs[i].id] = it->second; more = true; }
            }
        }
        // extrusions (210): u along the profile (its range), v from the payload (d[3] v0 v1); p-curves on
        // them are tried again once known. Unbounded ones (|v| >= 1e6) get no box: their v start is not established.
        if (p == 1 && L.derived) {
            bool added = false;
            for (std::size_t i = 0; i < recs.size(); ++i) {
                const Record& r = recs[i];
                if (r.type != 210 || !out[i].error.empty() || out[i].c.refs.empty() || out[i].c.payload.size() < 40 || cx.box.count(r.id)) continue;
                std::uint32_t prof = out[i].c.refs[0];
                if (auto t = cx.type.find(prof); t != cx.type.end() && is_edge_curve(t->second))
                    for (std::size_t m = 0; m < recs.size(); ++m)
                        if (recs[m].id == prof && out[m].error.empty() && !out[m].c.refs.empty()) prof = out[m].c.refs[0];
                auto rg = cx.range.find(prof);
                if (rg == cx.range.end()) continue;
                const std::uint8_t* q = out[i].c.payload.data();
                const double v0 = load_f64_le(q + 24), v1 = load_f64_le(q + 32);
                if (std::abs(v0) >= 1e6 || std::abs(v1) >= 1e6) continue;
                cx.box[r.id] = {rg->second[0], rg->second[1], v0, v1};
                added = true;
            }
            if (added) run_pass(1, true);
        }
    }
    return out;
}

}  // namespace cadcore::formats::catia::cgm
