// CGM compact dialect (2011 and later): record bodies -> the canonical record form (the varint
// dialect's body grammar, cgm_model.hpp), so one decoder serves both dialects.
//
// Compact bodies, as established against records of the same parts saved by older releases
// (research/cgm/FINDINGS.md, "Compact dialect bodies"):
//   signed integer  first byte b: low bits 01 -> +(b >> 2); 11 -> -((b >> 2) + 1); 00 -> n = bits 2-3
//                   bytes follow (little-endian), b's high nibble above them: v = (b >> 4) << 8n | bytes;
//                   10 -> the same, negative (-v - 1)
//   reference       a signed integer d, target = own id - d
//   reference list  varint count (as in the varint dialect) for variable lists, none for fixed arity
//   small integers  (the varint dialect's i16 / i32 fields) signed integers when the stream's header
//                   feature byte 0 has bit 3, else raw i16 / i32 as in the varint dialect
//   +-1 sequences   shells and loops pack their i16 senses as bits (1 = +1, 0 = -1)
// Per class (old class number; refs as in the varint dialect unless noted):
//   13, 513         refs
//   14 point        refs, per curve: f64 parameter, int
//   510 vertex      ref, u8
//   511 edge        [edge curve,] v0 v1 p0 p1, u8 flags; the edge-curve reference is absent when
//                   header feature byte 1 has bit 1
//   512 face        refs, int sense
//   520 shell       refs; bits: kind (4 bits: 1, 3, or 8 = -1), one sense bit per face; int
//   521 loop        refs, varint n, header byte (a, b: 0x41 (1,1), 0x50 (1,-1), 0xc1 (-1,1),
//                   0xd0 (-1,-1)), bits for the i16 sequence after the first x (implied -1): per
//                   coedge sense, y, next x, then the closing i16; int c; c != 0: block of 4 f64,
//                   4 ints (omitted in 2017+ streams: 25 bytes follow the f64), u8, 24 bytes
//   524, 525        ref(s), the i16 / i32 fields as ints
//   191-193 edge curves  refs, header byte h (bit 3: one map, count omitted; bit 2: map 0, map
//                   code omitted; bit 0: the record's last tolerance is zero and omitted); per map:
//                   int code, f64 t0 t1 tol, map data (map 1 may omit its identity range; resolved
//                   by the record length); closing range when some map is 2, 3 or 5; no 4 bytes
//   150 line        ref, int form: 0 o d t0 t1; 1 o.x t0 t1 (d = (0, 1)); 2 o.y t0 t1 (d = (1, 0))
//   151 circle      ref, c, [i32 i32 (raw-integer streams, some objects)], r t0 t1, u8 (1: s = 1),
//                   [s], offset
//   knot vector     int degree, int n, u8 z (4: x = 0, omitted; 0: int x follows), f64 knots[n], u8 1
//                   (standard multiplicities: ends degree + 1, interior degree - 2)
//   160 Hermite     ref, knot vector, i16 n, i16 x, P D1, int has2, [D2], t0 t1, i16 -2 (2018Fall+:
//                   then a byte 0)
//   130 B-spline    (code 22, no refs) knot vector, int rational, poles[np][3], [weights], tail as 170
//   curve
//   170 B-spline    ref, i16 0, knot vector, poles[np][2] (f64), int a, int b, f64 0 L 1 0, u8 c, [c != 0:
//   p-curve         further fields, not decoded], int -2
//   231 NURBS       knot vectors u v, int rational, poles (v fastest), [weights], a 152- or 154-byte
//                   block (not decoded)
//   200 plane       o, u8 m (bits 0-5: U V component is zero; 6: scale (1, 1); 7: default box), rest
//   201 cylinder    o, u8 m (bits 0-5 zero X Y components), r u0 u1 v0 v1, u8 (bits: su = 1, sv = 1,
//                   reference parameter = u0), the others
//   216, 243       int n, knot vector, sections P D1 D2 [n][10 | 6], box, rest unchanged
//   10, 100, 202, 203, 204 unchanged; 12, 155, 210, 211, 212: ref + unchanged; 214: two refs +
//                   unchanged; 241, 242, 244, 251: refs + unchanged (the decoder reads their references
//                   and leading doubles only)
//
// 2012+ streams (feature bytes 08 b5 75): compact integers, edges keep their edge-curve reference, and
// values that equal parameters of other records are omitted (verified by closure: the p-curves of
// every checked edge curve meet in 3D, and every checked vertex parameter lands on all its curves):
//   knot vector     form byte bit 3: standard multiplicities, the multiplicity byte is omitted
//   170 B-spline    ref, form byte 0x1e, f64 L (the knot range's length), knot vector, poles; the tail
//   p-curve         takes its usual values (a 1, b n - 1, 0 L 1 0, c 0) and the closing int is omitted
//   520 shell       some shells carry one 0 byte before the bits (the record length decides)
//   200 plane       int v; bits 0-2: origin component zero; 3: U is the unit vector of its one free
//                   component; 4-6: U component zero; 7, 8-10: the same for V; 11: scale (1, 1);
//                   12: default box; then the stored doubles (o U V scale box)
//   150 line        ref, u8 f. Directions: bit 0 (0, 1), bit 1 (1, 0), bit 6 (0, -1), bit 7 (-1, 0),
//                   none: stored. +v / +u: [c, t0, t1] with o = (c, 0) / (0, c); bit 4: c = box start
//                   of the other parameter; bit 2: t0 = box start; bit 3: t1 = box end if t0 is
//                   derived, else -t0. -v / -u: [t0, t1, o]; bit 2: o = 0; bit 3: t1 = -t0; bit 4 or
//                   5: t = minus the box range. General: [t0, t1, o, d]; bit 2: o = 0; bit 3: t1 = -t0
//                   (box: the surface's u0 u1 v0 v1)
//   191-193         refs, header h, u8 x = 0x80 | map count << 2 | bits (1: closing t0 = first t0,
//                   2: closing t1 = first t1; for one curve with map 0 or 1 the first range is the
//                   p-curve's range), maps as above with compact map codes (no separate count)
//   14 point        refs, u8 code: bit 7 refs 0 and 1 share the parameter; bits 2-3: refs[0]'s
//                   parameter is 2 the start, 1 the end of the edge curve's range, 0 stored f64 [+ int
//                   unless the shared curves include a line or circle]; the other refs' fields are
//                   not decoded (their parameters follow from the edge curve's maps)
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "formats/catia/cgm.hpp"

namespace cadcore::formats::catia::cgm {

// How a compact stream writes its bodies (from the header's feature bytes).
struct CompactLayout {
    bool compact_ints = false;   // small integers as signed compact integers
    bool edge_curve_ref = true;  // edges reference their edge curve
    bool derived = false;        // 2012+: values derived from other records (see canonical_records)
    bool fall2009 = false;       // 2009Fall: the varint dialect's grammar with absolute references
                                 // (cgm_compact2009.cpp)
};

// The layout of a compact stream whose bodies this module decodes (2011Spring: b5|bf 75|77;
// 2011Fall: 06 b5 75, the 2011Spring grammar with compact integers and edge-curve references;
// 2012+: 08 b5 75; 2009Fall: directory levels 83..85, feature byte 0 bit 3 for compact integers);
// nullopt otherwise.
std::optional<CompactLayout> compact_layout(const Stream& s);

struct Canonical {
    std::uint16_t type = 0;
    std::vector<std::uint32_t> refs;
    std::vector<std::uint8_t> payload;  // the varint dialect's payload, little-endian
};

// One compact record in canonical form (layouts without derived values). Throws DecodeError
// (Unsupported for classes not mapped or not decoded) when the body does not fit its class grammar
// exactly.
Canonical canonical(const Record& r, const CompactLayout& layout);

// All records of a compact stream in canonical form, in stream order; `error` is set (and the
// canonical form empty) for a record that did not decode. Layouts with derived values (2012+)
// resolve them in passes: surfaces (parameter boxes), p-curves (ranges), edge curves (ranges),
// points, then the rest.
// One record of a 2009Fall stream in canonical form (grammar in cgm_compact2009.cpp).
Canonical canonical_2009(const Record& r, bool compact_ints);

struct CanonicalRecord { Canonical c; std::string error; };
std::vector<CanonicalRecord> canonical_records(const Stream& s, const CompactLayout& layout);

}  // namespace cadcore::formats::catia::cgm
