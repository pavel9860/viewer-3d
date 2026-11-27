// CGM compact dialect: record bodies -> the varint dialect's form (cgm_model.hpp), so one decoder
// serves both. Differences from the varint grammar (research/cgm/FINDINGS.md has the evidence):
//   integers     compact signed integers (cgm_compact_io.hpp In::sint) when feature byte 0 has bit 3,
//                else raw i16 / i32; references are own id minus a signed integer
//   +-1 runs     shells and loops pack their senses as bits (1 = +1)
//   511 edge     the edge-curve reference is absent when feature byte 1 has bit 1
//   521 loop     header byte (a, b), sense bits, int c; c != 0: 4 f64, 4 ints (none in 2017+), u8, 24 bytes
//   191-193      header byte instead of the map count; codes and tolerances may be omitted
//   knot vector  int degree, int n, form byte, knots, standard or listed multiplicities
//   150, 200, 201 masked or implied components; 151 optional ints; 160 adds i16 n x and int -2
//   130 / 170    (codes 22 / 33) tail int a, int b, f64 0 L 1 0, u8 c, int -2
//   155          items: u8 0x41, int code, [1, 4: varint inputs], values, children
// 2012+ streams (08 b5 75) omit values that equal parameters of other records (boxes, curve ranges);
// canonical_records resolves them in passes.
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "formats/catia/cgm.hpp"

namespace cadcore::formats::catia::cgm {

struct CompactLayout {
    bool compact_ints = false;
    bool edge_curve_ref = true;
    bool derived = false;   // 2012+
    bool fall2009 = false;  // the varint grammar with absolute references (cgm_compact2009.cpp)
};

// From the header: 2011Spring b5|bf 75|77; 2011Fall 06 b5 75; 2012+ 08 b5 75; 2009Fall levels 83..85.
std::optional<CompactLayout> compact_layout(const Stream& s);

struct Canonical {
    std::uint16_t type = 0;
    std::vector<std::uint32_t> refs;
    std::vector<std::uint8_t> payload;  // the varint dialect's payload, little-endian
};

// One record (layouts without derived values); throws DecodeError when the body does not fit exactly.
Canonical canonical(const Record& r, const CompactLayout& layout);
Canonical canonical_2009(const Record& r, bool compact_ints);

// All records in stream order; `error` is set for a record that did not decode.
struct CanonicalRecord { Canonical c; std::string error; };
std::vector<CanonicalRecord> canonical_records(const Stream& s, const CompactLayout& layout);

}  // namespace cadcore::formats::catia::cgm
