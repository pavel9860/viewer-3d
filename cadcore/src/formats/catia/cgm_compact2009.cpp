// CGM compact dialect, 2009Fall streams (directory levels 83..85) -> canonical records.
//
// Grammar established over every record of the corpus's 2009Fall parts (219,050 records, each consumed
// exactly) and against records of the same parts saved by V5R10 / V5R17 (research/cgm/FINDINGS.md):
// the varint dialect's body grammar (cgm_model.hpp), with
//   references      varint count, varint absolute ids (as in the varint dialect)
//   scalar integers compact signed integers when the header's feature byte 0 has bit 3 (0b), else the
//                   raw i16 / i32 of the varint dialect (01); the i16 sequences of shells and loops stay
//                   raw i16 in both
// and these differences per class:
//   150 line        int form: 0 o d t0 t1; 1 o.x t0 t1 (d = (0, 1)); 2 o.y t0 t1 (d = (1, 0))
//   170 B-spline    i16 0, knot vector, poles[2], tail (cgm_compact_io.hpp) closing with an int -2
//   160 Hermite     after the value count (compact streams only) an int equal to the knot vector's fifth
//                   field; after t0 t1 an int -2
//   191-193         after the maps an i32; edge curves listing their curves twice (A B A B) add two ints
//                   (not decoded)
//   521 loop        the closing int c always follows; c != 0: 4 f64, 4 ints, u8, 24 bytes
//   knot vector     ints degree, a, x, n, b (b: the number of knots before parameter 0), f64 knots[n],
//                   ints multiplicities[n]
//   211, 216        as in the other compact layouts (cgm_compact.cpp)
#include "formats/catia/cgm_compact.hpp"

#include <string>

#include "cadcore/error.hpp"
#include "formats/catia/cgm_compact_io.hpp"

namespace cadcore::formats::catia::cgm {
namespace {
using namespace compact_io;

std::vector<std::uint32_t> abs_refs(In& in) {
    std::vector<std::uint32_t> r(in.count());
    for (auto& x : r) x = varint(in.b, in.k);
    return r;
}

// Knot vector -> i32 degree, a, x, n, b, f64 knots[n], i32 multiplicities[n].
struct Knots { int n = 0, b = 0; std::size_t poles = 0; };
Knots knot_vector(In& in, Out& o, bool counted = false) {
    if (counted) {  // 216: the section count precedes the knot vector
        const std::int32_t c = in.i32();
        o.i32(c);
    }
    const int deg = in.i32(), a = in.i32(), x = in.i32(), n = in.i32(), b = in.i32();
    if (n < 2 || n > 100000 || deg < 1 || deg > 30) bad("knot vector size", in.at);
    o.i32(deg); o.i32(a); o.i32(x); o.i32(n); o.i32(b);
    for (int i = 0; i < n; ++i) o.f64(in.f64());
    std::size_t sum = 0;
    for (int i = 0; i < n; ++i) {
        const std::int32_t m = in.i32();
        if (m < 1 || m > deg + 1) bad("knot multiplicity", in.at);
        o.i32(m);
        sum += std::size_t(m);
    }
    if (sum < std::size_t(deg) + 2) bad("knot multiplicities", in.at);
    return {n, b, sum - std::size_t(deg) - 1};
}

void f64s(In& in, Out& o, std::size_t n) { for (std::size_t i = 0; i < n; ++i) o.f64(in.f64()); }

}  // namespace

Canonical canonical_2009(const Record& r, bool compact_ints) {
    In in{r.body, r.id, compact_ints, r.offset};
    Out o;
    Canonical c;
    c.type = r.type;
    if (r.type == 0) throw DecodeError(ErrorCode::Unsupported, "compact type code " + std::to_string(r.code) + " not mapped", r.offset);
    c.refs = abs_refs(in);
    // raw i16 (the shell and loop sequences are raw in both integer forms)
    auto raw16 = [&] { in.need(2); const auto v = std::int16_t(load_u16(in.b.data + in.k, false)); in.k += 2; return std::int32_t(v); };
    switch (r.type) {
        case 13: case 513: break;
        case 14: {
            const std::uint32_t n = in.count();
            o.count(n);
            for (std::uint32_t i = 0; i < n; ++i) { o.f64(in.f64()); o.i32(in.i32()); }
            break;
        }
        case 510: case 511: o.u8(in.u8()); break;
        case 512: o.i16(in.i16()); break;
        case 520: {
            const std::uint32_t n = in.count();
            o.count(n);
            o.i16(in.i16());
            for (std::uint32_t i = 0; i < n; ++i) o.i16(raw16());
            o.i16(in.i16());
            break;
        }
        case 521: {
            const std::uint32_t n = in.count();
            o.count(n);
            for (int i = 0; i < 3; ++i) o.i16(in.i16());  // a, b, the first value of the sequence
            for (std::uint32_t i = 0; i < 3 * n; ++i) o.i16(raw16());
            const std::int32_t cl = in.i16();
            o.i16(cl);
            if (cl) {
                f64s(in, o, 4);
                for (int i = 0; i < 4; ++i) o.i32(in.i32());
                o.u8(in.u8());
                o.raw(in.raw(24));
            }
            break;
        }
        case 191: case 192: case 193: {
            const std::uint32_t n = in.count();
            if (n > 16) bad("edge-curve map count", r.offset);
            o.count(n);
            bool tail = false;
            for (std::uint32_t j = 0; j < n; ++j) {
                const std::int32_t m = in.i16();
                o.i16(m);
                f64s(in, o, 3);
                if (m == 0) continue;
                if (m == 1) { f64s(in, o, 3); continue; }
                tail = true;
                f64s(in, o, 2);
                if (m == 3) f64s(in, o, 2);
                else if (m == 5) {
                    const std::uint32_t N = in.count();
                    if (N > 100000) bad("edge-curve table size", r.offset);
                    o.count(N);
                    f64s(in, o, 2 * std::size_t(N));
                } else if (m != 2) bad("edge-curve map " + std::to_string(m), r.offset);
            }
            if (tail) f64s(in, o, 3);
            o.i32(in.i32());
            const std::size_t h = c.refs.size() / 2;
            if (c.refs.size() == 2 * h && n + 1 == h && std::equal(c.refs.begin(), c.refs.begin() + std::ptrdiff_t(h), c.refs.begin() + std::ptrdiff_t(h))) {
                in.i16();
                in.i16();
            }
            break;
        }
        case 150: {
            const std::int32_t form = in.i16();
            if (form == 0) { f64s(in, o, 6); break; }
            if (form != 1 && form != 2) bad("line form " + std::to_string(form), r.offset);
            const double c0 = in.f64(), t0 = in.f64(), t1 = in.f64();
            if (form == 1) { o.f64(c0); o.f64(0); o.f64(0); o.f64(1); }
            else { o.f64(0); o.f64(c0); o.f64(1); o.f64(0); }
            o.f64(t0); o.f64(t1);
            break;
        }
        case 151:
            f64s(in, o, 2);
            o.i32(in.i32());
            o.i32(in.i32());
            f64s(in, o, 5);
            break;
        case 160: {
            const auto [n, b, np] = knot_vector(in, o);
            (void)np;
            const std::int32_t n2 = in.i32();
            if (n2 != n) bad("Hermite value count differs from the knot count", r.offset);
            o.i32(n2);
            if (compact_ints && in.i32() != b) bad("Hermite knot-offset repeat", r.offset);
            f64s(in, o, 4 * std::size_t(n));
            const std::int32_t has2 = in.i16();
            o.i16(has2);
            if (has2) f64s(in, o, 2 * std::size_t(n));
            f64s(in, o, 2);
            if (in.i16() != -2) bad("Hermite closing value", r.offset);
            break;
        }
        case 170: {
            if (in.i16() != 0) bad("B-spline p-curve form", r.offset);
            o.i16(0);
            const std::size_t np = knot_vector(in, o).poles;
            if (2 * np > in.left() / 8) bad("B-spline p-curve poles past the record", r.offset);
            f64s(in, o, 2 * np);
            bspline_pcurve_tail(in, o);
            break;
        }
        case 231:
            knot_vector(in, o);
            knot_vector(in, o);
            o.i16(in.i16());
            o.raw(in.raw(in.left()));  // poles, weights, a block (cgm_model.hpp reads poles and weights)
            break;
        case 524: o.i16(in.i16()); break;
        case 525: o.i32(in.i32()); o.i16(in.i16()); break;
        // the decoder uses the references and the leading doubles of these (cgm_model.hpp)
        case 216: case 243: {
            const auto [n, b, np] = knot_vector(in, o, true);
            (void)b;
            (void)np;
            f64s(in, o, (r.type == 216 ? 30 : 18) * std::size_t(n) + 4);
            o.raw(in.raw(in.left()));
            break;
        }
        case 211:
            f64s(in, o, 16);
            o.i32(in.i32()); o.i32(in.i32());
            f64s(in, o, 4);
            o.i32(in.i32());
            f64s(in, o, 1);
            break;
        case 10: case 12: case 100: case 155: case 200: case 201: case 202: case 203: case 204: case 210:
        case 212: case 214: case 241: case 242: case 244: case 251: case 500: case 502:
            o.raw(in.raw(in.left()));
            break;
        default:
            throw DecodeError(ErrorCode::Unsupported, "compact class " + std::to_string(r.type) + " (code " + std::to_string(r.code) + ") not decoded", r.offset);
    }
    if (in.left()) bad(std::to_string(in.left()) + " bytes past the class grammar", r.offset);
    c.payload = std::move(o.p);
    return c;
}

}  // namespace cadcore::formats::catia::cgm
