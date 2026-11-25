#include "formats/catia/cgm_compact.hpp"

#include <cstring>
#include <optional>
#include <string>

#include "cadcore/error.hpp"

namespace cadcore::formats::catia::cgm {
namespace {

[[noreturn]] void bad(const std::string& what, std::size_t at) { fail(ErrorCode::Malformed, "CGM compact record: " + what, at); }

// Reads a compact body.
struct In {
    ByteView b;
    std::uint32_t self;
    bool ci;
    std::size_t at, k = 0;

    void need(std::size_t n) const { if (k + n > b.size) fail(ErrorCode::Truncated, "CGM compact record ends early", at + k); }
    std::int64_t sint() {
        need(1);
        const std::uint8_t x = b[k++];
        const unsigned lo = x & 3;
        if (lo == 1) return x >> 2;
        if (lo == 3) return -std::int64_t(x >> 2) - 1;
        const unsigned n = x >> 2;
        if (n > 8) bad("integer width", at + k - 1);
        need(n);
        std::uint64_t v = 0;
        for (unsigned i = 0; i < n; ++i) v |= std::uint64_t(b[k++]) << (8 * i);
        return lo == 0 ? std::int64_t(v) : -std::int64_t(v) - 1;
    }
    std::uint32_t ref() { return std::uint32_t(std::int64_t(self) - sint()); }
    std::uint32_t count() { return varint(b, k); }
    std::vector<std::uint32_t> refs(std::size_t n) { std::vector<std::uint32_t> r(n); for (auto& x : r) x = ref(); return r; }
    std::vector<std::uint32_t> refs() { return refs(count()); }
    std::uint8_t u8() { need(1); return b[k++]; }
    double f64() { need(8); double d; std::memcpy(&d, b.data + k, 8); k += 8; return d; }
    std::int32_t i16() { if (ci) return std::int32_t(sint()); need(2); auto v = std::int16_t(load_u16(b.data + k, false)); k += 2; return v; }
    std::int32_t i32() { if (ci) return std::int32_t(sint()); need(4); auto v = std::int32_t(load_u32(b.data + k, false)); k += 4; return v; }
    std::vector<int> bits(std::size_t n) {
        const std::size_t nb = (n + 7) / 8;
        need(nb);
        std::vector<int> v(n);
        for (std::size_t j = 0; j < n; ++j) v[j] = b[k + j / 8] >> (j % 8) & 1;
        k += nb;
        return v;
    }
    ByteView raw(std::size_t n) { need(n); ByteView v(b.data + k, n); k += n; return v; }
    std::size_t left() const { return b.size - k; }
};

// Writes a payload in the varint dialect's form.
struct Out {
    std::vector<std::uint8_t> p;
    void u8(std::uint8_t v) { p.push_back(v); }
    void i16(std::int32_t v) { for (int i = 0; i < 2; ++i) p.push_back(std::uint8_t(std::uint32_t(v) >> (8 * i))); }
    void i32(std::int32_t v) { for (int i = 0; i < 4; ++i) p.push_back(std::uint8_t(std::uint32_t(v) >> (8 * i))); }
    void f64(double d) { std::uint8_t t[8]; std::memcpy(t, &d, 8); p.insert(p.end(), t, t + 8); }
    void raw(ByteView v) { p.insert(p.end(), v.data, v.data + v.size); }
    void count(std::uint32_t n) {  // varint (cgm.hpp)
        if (n < 0x80) { p.push_back(std::uint8_t(0x80 | n)); return; }
        std::uint8_t mask = 0;
        std::vector<std::uint8_t> bytes;
        for (int j = 0; j < 4; ++j)
            if (std::uint8_t x = std::uint8_t(n >> (8 * j))) { mask |= std::uint8_t(1u << j); bytes.push_back(x); }
        p.push_back(std::uint8_t(mask << 3));
        p.insert(p.end(), bytes.begin(), bytes.end());
    }
};

int pm(int bit) { return bit ? 1 : -1; }

// Knot vector: int degree, int n, u8 z (4: x = 0 omitted; 0: int x follows), f64 knots[n], u8 1
// (standard multiplicities: ends degree + 1, interior degree - 2). Written in the varint dialect's
// form: i32 degree, 0, x, n, 0, knots, i32 multiplicities. Returns the number of poles.
std::size_t knot_vector(In& in, Out& o, int& n_out) {
    const int deg = in.i32(), n = in.i32();
    const std::uint8_t z = in.u8();
    if (z != 4 && z != 0) bad("knot-vector form byte", in.at);
    const int x = z ? 0 : in.i32();
    if (n < 2 || n > 100000 || deg < 1 || deg > 30) bad("knot vector size", in.at);
    if (n > 2 && deg < 3) throw DecodeError(ErrorCode::Unsupported, "interior knot multiplicity of degree " + std::to_string(deg) + " not established", in.at);
    o.i32(deg); o.i32(0); o.i32(x); o.i32(n); o.i32(0);
    for (int i = 0; i < n; ++i) o.f64(in.f64());
    if (in.u8() != 1) throw DecodeError(ErrorCode::Unsupported, "knot multiplicities not standard", in.at);
    std::size_t sum = 0;
    for (int i = 0; i < n; ++i) { const int m = i == 0 || i == n - 1 ? deg + 1 : deg - 2; o.i32(m); sum += std::size_t(m); }
    n_out = n;
    return sum - std::size_t(deg) - 1;
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

}  // namespace

std::optional<CompactLayout> compact_layout(const Stream& s) {
    if (s.dialect != Dialect::Compact || s.features.size() < 2) return std::nullopt;
    const std::uint8_t f0 = s.features[0], f1 = s.features[1];
    // 2011Spring streams: feature bytes b5|bf 75|77 (bit 3 of byte 0: compact integers; bit 1 of
    // byte 1: edges without their edge-curve reference).
    if ((f0 == 0xb5 || f0 == 0xbf) && (f1 == 0x75 || f1 == 0x77)) return CompactLayout{(f0 & 8) != 0, (f1 & 2) == 0};
    return std::nullopt;
}

Canonical canonical(const Record& r, const CompactLayout& L) {
    In in{r.body, r.id, L.compact_ints, r.offset};
    Out o;
    Canonical c;
    c.type = r.type;
    switch (r.type) {
        case 13: case 513: c.refs = in.refs(); break;
        case 14:
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
        case 191: case 192: case 193: c.refs = in.refs(); edge_curve(in, o, c.refs.size()); break;
        case 150: {
            c.refs = in.refs(1);
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
            const double cx = in.f64(), cy = in.f64();
            // r t0 t1, flag byte (1: scale 1), [scale], offset; some objects of raw-integer streams
            // carry two i32 (1, 1) after the centre
            auto fits = [&](std::size_t k) {
                return k + 24 < in.b.size && in.b[k + 24] <= 1 && in.b.size - k == 25 + (in.b[k + 24] == 0 ? 8u : 0u) + 8;
            };
            bool plain = fits(in.k), ints = fits(in.k + 8);
            if (plain && ints && load_u32(in.b.data + in.k, false) == 1 && load_u32(in.b.data + in.k + 4, false) == 1) plain = false;
            if (plain == ints) bad("circle layout", r.offset);
            if (ints) in.raw(8);
            o.f64(cx); o.f64(cy); o.i32(1); o.i32(1);
            const double rad = in.f64(), t0 = in.f64(), t1 = in.f64();
            const std::uint8_t f = in.u8();
            const double s = f == 1 ? 1.0 : in.f64();
            o.f64(rad); o.f64(t0); o.f64(t1); o.f64(s); o.f64(in.f64());
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
            o.f64(in.f64()); o.f64(in.f64());
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
            o.raw(in.raw(in.left()));  // trailing block (152 or 154 bytes; not decoded)
            break;
        }
        case 200: {
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
            break;
        }
        case 10: case 100: case 202: case 203: case 204: case 216: case 243: o.raw(in.raw(in.left())); break;
        // the decoder uses the references and the leading doubles of these (cgm_model.hpp)
        case 12: case 155: case 210: case 211: case 212: c.refs = in.refs(1); o.raw(in.raw(in.left())); break;
        case 214: c.refs = in.refs(2); o.raw(in.raw(in.left())); break;
        case 241: case 242: case 244: case 251: c.refs = in.refs(); o.raw(in.raw(in.left())); break;
        default:
            throw DecodeError(ErrorCode::Unsupported, "compact class " + std::to_string(r.type) + " (code " + std::to_string(r.code) + ") not decoded", r.offset);
    }
    if (in.left()) bad(std::to_string(in.left()) + " bytes past the class grammar", r.offset);
    c.payload = std::move(o.p);
    return c;
}

}  // namespace cadcore::formats::catia::cgm
