// CGM compact dialect (2011 and later): record bodies -> the canonical record form (the varint
// dialect's body grammar, cgm_model.hpp), so one decoder serves both dialects.
//
// Compact bodies, as established against records of the same parts saved by older releases
// (research/cgm/FINDINGS.md, "Compact dialect bodies"):
//   signed integer  first byte b: low bits 01 -> +(b >> 2); 11 -> -((b >> 2) + 1); 00 -> (b >> 2)
//                   bytes follow (little-endian, positive); 10 -> the same, negative (-v - 1)
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
//                   4 ints, u8, 24 bytes
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
//   160 Hermite     ref, knot vector, i16 n, i16 x, P D1, int has2, [D2], t0 t1, i16 -2
//   231 NURBS       knot vectors u v, int rational, poles (v fastest), [weights], a 152- or 154-byte
//                   block (not decoded)
//   200 plane       o, u8 m (bits 0-5: U V component is zero; 6: scale (1, 1); 7: default box), rest
//   201 cylinder    o, u8 m (bits 0-5 zero X Y components), r u0 u1 v0 v1, u8 (bits: su = 1, sv = 1,
//                   reference parameter = u0), the others
//   10, 100, 202, 203, 204, 216, 243 unchanged; 12, 155, 210, 211, 212: ref + unchanged; 214: two refs +
//                   unchanged; 241, 242, 244, 251: refs + unchanged (the decoder reads their references
//                   and leading doubles only)
#pragma once

#include <cstdint>
#include <optional>
#include <vector>

#include "formats/catia/cgm.hpp"

namespace cadcore::formats::catia::cgm {

// How a compact stream writes its bodies (from the header's feature bytes).
struct CompactLayout {
    bool compact_ints = false;   // small integers as signed compact integers
    bool edge_curve_ref = true;  // edges reference their edge curve
};

// The layout of a compact stream whose bodies this module decodes; nullopt otherwise (2009Fall,
// 2011Fall and 2012+ streams use further encodings, not decoded yet).
std::optional<CompactLayout> compact_layout(const Stream& s);

struct Canonical {
    std::uint16_t type = 0;
    std::vector<std::uint32_t> refs;
    std::vector<std::uint8_t> payload;  // the varint dialect's payload, little-endian
};

// One compact record in canonical form. Throws DecodeError (Unsupported for classes not mapped or
// not decoded) when the body does not fit its class grammar exactly.
Canonical canonical(const Record& r, const CompactLayout& layout);

}  // namespace cadcore::formats::catia::cgm
