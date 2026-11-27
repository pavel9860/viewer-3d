// CGM compact record grammar: hand-built bodies, one per stored form, against their exact canonical
// (varint-dialect) payloads. No corpus needed.
#include <cstring>
#include <initializer_list>
#include <memory>
#include <vector>

#include "cadcore/error.hpp"
#include "formats/catia/cgm.hpp"
#include "formats/catia/cgm_compact.hpp"
#include "harness.hpp"

using namespace cadcore;
namespace cgm = cadcore::formats::catia::cgm;

namespace {

// Little-endian byte writer for bodies (compact or raw integers) and expected payloads.
struct W {
    std::vector<std::uint8_t> b;
    bool ci = false;  // integers as compact signed integers
    W& u8(std::uint8_t v) { b.push_back(v); return *this; }
    W& le(std::uint64_t v, int n) { for (int i = 0; i < n; ++i) b.push_back(std::uint8_t(v >> (8 * i))); return *this; }
    W& sint(std::int64_t v) {
        const bool neg = v < 0;
        const std::uint64_t m = neg ? std::uint64_t(-v - 1) : std::uint64_t(v);
        if (m < 64) return u8(std::uint8_t(m << 2 | (neg ? 3 : 1)));
        int n = 1;
        while (n < 3 && (m >> (8 * n)) >= 16) ++n;
        u8(std::uint8_t((m >> (8 * n)) << 4 | std::uint64_t(n) << 2 | (neg ? 2 : 0)));
        return le(m, n);
    }
    W& i16(std::int32_t v) { return ci ? sint(v) : le(std::uint32_t(v), 2); }
    W& i32(std::int32_t v) { return ci ? sint(v) : le(std::uint32_t(v), 4); }
    W& f64(double d) { std::uint64_t v; std::memcpy(&v, &d, 8); return le(v, 8); }
    W& f64s(std::initializer_list<double> ds) { for (double d : ds) f64(d); return *this; }
    W& zeros(int n) { for (int i = 0; i < n; ++i) f64(0); return *this; }
    W& count(std::uint32_t n) { return u8(std::uint8_t(0x80 | n)); }  // varint, n < 128
};
W raw() { return W{}; }
W compact() { W w; w.ci = true; return w; }

// Canonical form of one record (layouts without derived values).
cgm::Canonical canon(std::uint16_t type, std::uint32_t id, const W& body, cgm::CompactLayout L) {
    cgm::Record r;
    r.type = type;
    r.id = id;
    r.body = ByteView(body.b.data(), body.b.size());
    return cgm::canonical(r, L);
}
const cgm::CompactLayout kRaw{false, true, false};
const cgm::CompactLayout kCompact{true, true, false};

// A 2012+ stream from (type, id, body) records.
struct Recs {
    std::vector<std::unique_ptr<W>> bodies;
    cgm::Stream s;
    void add(std::uint16_t type, std::uint32_t id, W body) {
        bodies.push_back(std::make_unique<W>(std::move(body)));
        cgm::Record r;
        r.type = type;
        r.id = id;
        r.body = ByteView(bodies.back()->b.data(), bodies.back()->b.size());
        s.records.push_back(r);
    }
};

}  // namespace

TEST_CASE("CGM compact: signed integers and relative references") {
    // 2011Spring edge (compact integers, no edge-curve reference): v0 v1 p0 p1 as self - target, u8
    const W edge = compact().sint(2).sint(-4).sint(-2).sint(-3).u8(0x21);
    const auto c = canon(511, 100, edge, cgm::CompactLayout{true, false, false});
    CHECK((c.refs == std::vector<std::uint32_t>{98, 104, 102, 103}));
    CHECK((c.payload == std::vector<std::uint8_t>{0x21}));
    // multi-byte compact integers (1 to 3 bytes after the first)
    const W big = compact().sint(0x1234).sint(-0x1235).sint(0x12345).sint(64).u8(0);
    const auto b = canon(511, 0x100000, big, cgm::CompactLayout{true, false, false});
    CHECK((b.refs == std::vector<std::uint32_t>{0x100000 - 0x1234, 0x100000 + 0x1235, 0x100000 - 0x12345, 0x100000 - 64}));
}

TEST_CASE("CGM compact: shells and loops (sense bits, loop blocks)") {
    // shell: one face, kind 3 (bits 1100), sense +1, closing -1 (raw i16)
    W shell = raw().count(1).sint(-1);
    shell.u8(0x13).i16(-1);
    const auto s = canon(520, 100, shell, kRaw);
    CHECK((s.payload == raw().count(1).i16(3).i16(1).i16(-1).b));
    // loop: (pcurve, edge) + surface; header 0x41 (a = b = 1); bits for sense, y, next x; closing 1 with
    // the block: 4 f64, 4 ints, u8, 24 bytes (2014/2015) or without the ints (2017+)
    for (bool ints : {true, false}) {
        W l = raw().count(3).sint(1).sint(2).sint(3).count(1).u8(0x41).u8(0b101).i16(1).zeros(4);
        if (ints) l.i32(1).i32(1).i32(1).i32(1);
        l.u8(7);
        for (int i = 0; i < 24; ++i) l.u8(std::uint8_t(i));
        const auto c = canon(521, 100, l, kRaw);
        W want = raw().count(1).i16(1).i16(1).i16(-1).i16(1).i16(-1).i16(1).i16(1).zeros(4);
        for (int i = 0; i < 4; ++i) want.i32(ints ? 1 : 0);
        want.u8(7);
        for (int i = 0; i < 24; ++i) want.u8(std::uint8_t(i));
        CHECK(c.payload == want.b);
    }
}

TEST_CASE("CGM compact: knot vectors (forms of x and multiplicities)") {
    // 130 3D B-spline (no refs): degree 1, knots 0 1 2, x stored, multiplicities listed after a byte 0
    W b = raw().i32(1).i32(3).u8(0).i32(0).f64s({0, 1, 2}).u8(0).i32(2).i32(1).i32(2);
    b.i16(0).zeros(9);                                           // not rational, 3 poles
    b.i16(1).i16(2).f64s({0, 2, 1, 0}).u8(0).i16(-2);           // tail
    W want = raw().i32(1).i32(0).i32(0).i32(3).i32(0).f64s({0, 1, 2}).i32(2).i32(1).i32(2);
    want.i16(0).zeros(9).i16(1).i16(2).f64s({0, 2, 1, 0}).u8(0);
    CHECK(canon(130, 1, b, kRaw).payload == want.b);
    // 170 B-spline p-curve, compact integers: degree 3, x omitted (form 4), standard multiplicities (byte 1)
    W p = compact().sint(1).i16(0).i32(3).i32(2).u8(4).f64s({0, 5}).u8(1).zeros(8);
    p.i16(1).i16(1).f64s({0, 5, 1, 0}).u8(0).i16(-2);
    W pw = raw().i16(0).i32(3).i32(0).i32(0).i32(2).i32(0).f64s({0, 5}).i32(4).i32(4).zeros(8);
    pw.i16(1).i16(1).f64s({0, 5, 1, 0}).u8(0);
    const auto c = canon(170, 10, p, kCompact);
    CHECK((c.refs == std::vector<std::uint32_t>{9}));
    CHECK(c.payload == pw.b);
    // tail with c != 0: the extra fields are kept as written
    W e = raw().sint(1).i16(0).i32(3).i32(2).u8(0x0c).f64s({0, 5}).zeros(8);
    e.i16(1).i16(1).f64s({0, 5, 1, 0}).u8(1).u8(0xaa).u8(0xbb).i16(-2);
    CHECK(canon(170, 10, e, kRaw).payload.back() == 0xbb);
}

TEST_CASE("CGM compact: lines, circles, Hermite curves") {
    // 150 forms: 1 o.x t0 t1 (d = (0, 1)); 2 o.y t0 t1 (d = (1, 0))
    CHECK(canon(150, 5, raw().sint(1).i16(1).f64s({3, 0, 2}), kRaw).payload == raw().f64s({3, 0, 0, 1, 0, 2}).b);
    CHECK(canon(150, 5, raw().sint(1).i16(2).f64s({4, 1, 2}), kRaw).payload == raw().f64s({0, 4, 1, 0, 1, 2}).b);
    // 151 with and without the two ints; flag 1: scale 1
    const W want = raw().f64s({1, 2}).i32(1).i32(1).f64s({3, 0, 6, 1, 0});
    CHECK(canon(151, 5, raw().sint(1).f64s({1, 2, 3, 0, 6}).u8(1).f64(0), kRaw).payload == want.b);
    CHECK(canon(151, 5, raw().sint(1).f64s({1, 2}).i32(1).i32(1).f64s({3, 0, 6}).u8(1).f64(0), kRaw).payload == want.b);
    // 160 Hermite (degree 5, two knots), cubic (no D2); 2018Fall+ adds a byte 0 after the closing -2
    for (bool tail : {false, true}) {
        W h = raw().sint(1).i32(5).i32(2).u8(0x0c).f64s({0, 1}).i16(2).i16(0).zeros(8).i16(0).f64s({0, 1}).i16(-2);
        if (tail) h.u8(0);
        W hw = raw().i32(5).i32(0).i32(0).i32(2).i32(0).f64s({0, 1}).i32(6).i32(6).i32(2).zeros(8).i16(0).f64s({0, 1});
        CHECK(canon(160, 5, h, kRaw).payload == hw.b);
    }
}

TEST_CASE("CGM compact: expression p-curves (155) and section surfaces (216, 243)") {
    // u = item 1 (linear, one input); v = item 3 of items 2 and 7
    W e = raw().sint(1).f64s({0, 1, 0, 0});
    e.u8(0x41).i16(1).count(1).f64s({0, 1, 0});
    e.u8(0x41).i16(3).u8(0x41).i16(2).f64s({1, 0, 1, 0}).u8(0x41).i16(7).f64s({1, 1, 1, 0, 0});
    W ew = raw().f64s({0, 1, 0, 0}).i16(1).i16(1).i16(0).f64s({0, 1, 0});
    ew.i16(3).i16(2).f64s({1, 0, 1, 0}).i16(7).f64s({1, 1, 1, 0, 0});
    CHECK(canon(155, 5, e, kRaw).payload == ew.b);
    bool threw = false;
    try { canon(155, 5, raw().sint(1).f64s({0, 1, 0, 0}).u8(0x41).i16(9), kRaw); } catch (const DecodeError&) { threw = true; }
    CHECK(threw);  // unknown item code
    // 216 (10 values per section) and 243 (6): n, knot vector, P D1 D2, box, rest kept
    for (int m : {10, 6}) {
        W s = raw().i32(2).i32(5).i32(2).u8(0x0c).f64s({0, 4}).zeros(3 * m * 2).f64s({0, 4, -1, 2}).u8(9);
        W sw = raw().i32(2).i32(5).i32(0).i32(0).i32(2).i32(0).f64s({0, 4}).i32(6).i32(6).zeros(3 * m * 2).f64s({0, 4, -1, 2}).u8(9);
        CHECK(canon(m == 10 ? 216 : 243, 5, s, kRaw).payload == sw.b);
    }
}

TEST_CASE("CGM compact 2009Fall: absolute references, line forms, B-spline closing") {
    W line = compact().count(1).u8(0x85).i16(2).f64s({4, 1, 2});  // refs: varint count, varint id 5
    cgm::Record r;
    r.type = 150;
    r.id = 9;
    r.body = ByteView(line.b.data(), line.b.size());
    const auto c = cgm::canonical_2009(r, true);
    CHECK((c.refs == std::vector<std::uint32_t>{5}));
    CHECK(c.payload == raw().f64s({0, 4, 1, 0, 1, 2}).b);
    W bs = compact().count(1).u8(0x85).i16(0).i32(1).i32(0).i32(0).i32(2).i32(0).f64s({0, 1}).i32(2).i32(2).zeros(4);
    bs.i16(1).i16(1).f64s({0, 1, 1, 0}).u8(0).i32(-2);
    r.type = 170;
    r.body = ByteView(bs.b.data(), bs.b.size());
    CHECK(cgm::canonical_2009(r, true).payload == raw().i16(0).i32(1).i32(0).i32(0).i32(2).i32(0).f64s({0, 1}).i32(2).i32(2).zeros(4).i16(1).i16(1).f64s({0, 1, 1, 0}).u8(0).b);
}

TEST_CASE("CGM compact 2012+: values derived from other records") {
    Recs t;
    // plane 100: origin 0, U = x, V = y (unit flags), scale 1 1, box u 0..10, v 0..5
    t.add(200, 100, compact().sint(7 | 8 | 6 << 4 | 0x80 | 5 << 8 | 1 << 11).f64s({0, 10, 0, 5}));
    // line 101 along +u at the box's v start, t0 stored 2, t1 = the box's u end (bit 3)
    t.add(150, 101, compact().sint(1).u8(0x02 | 0x10 | 0x08).f64(2));
    // edge curve 102 of the line: one map 0, range from the line (x bits 0, 1), last tolerance 0
    t.add(192, 102, compact().count(1).sint(1).u8(0x0d).u8(0x80 | 1 << 2 | 3));
    // point 103 shared by 102 and 101 (bit 7) at the start
    t.add(14, 103, compact().count(2).sint(1).sint(2).u8(0x88));
    // Hermite 105 on the plane, its edge curve 106, point 107 stored: offset 0.5 in segment 2 (an int
    // follows: the curve is not a line or circle)
    t.add(160, 105, compact().sint(5).i32(5).i32(3).u8(0x0c).f64s({0, 1, 2}).i16(3).i16(0).zeros(12).i16(0).f64s({0, 2}).i16(-2));
    t.add(192, 106, compact().count(1).sint(1).u8(0x0d).u8(0x80 | 1 << 2 | 3));
    t.add(14, 107, compact().count(1).sint(1).u8(0x00).f64(0.5).sint(2));
    // edge curve 108 of the line with a stored range wider than the line's (0..20); point 109 shared by
    // 108 and 101 at the start takes the line's start (2)
    t.add(192, 108, compact().count(1).sint(7).u8(0x0d).u8(0x80 | 1 << 2).f64s({0, 20}));
    t.add(14, 109, compact().count(2).sint(1).sint(8).u8(0x88));
    const auto out = cgm::canonical_records(t.s, cgm::CompactLayout{true, true, true});
    for (std::size_t i = 0; i < out.size(); ++i) {
        CAPTURE(t.s.records[i].id);
        CAPTURE(out[i].error);
        CHECK(out[i].error.empty());
    }
    auto payload = [&](std::uint32_t id) { for (std::size_t i = 0; i < out.size(); ++i) if (t.s.records[i].id == id) return out[i].c.payload; return std::vector<std::uint8_t>{}; };
    CHECK(payload(100) == raw().f64s({0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 10, 0, 5}).b);
    CHECK(payload(101) == raw().f64s({0, 0, 1, 0, 2, 10}).b);
    CHECK(payload(102) == raw().count(1).i16(0).f64s({2, 10, 0}).le(0, 4).b);
    CHECK(payload(103) == raw().count(2).f64(2).i32(1).f64(2).i32(1).b);
    CHECK(payload(107) == raw().count(1).f64(0.5).i32(2).b);
    CHECK(payload(109) == raw().count(2).f64(2).i32(1).f64(2).i32(1).b);
}
