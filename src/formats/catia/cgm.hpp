// CATIA V5 exact geometry: the CGMGeom container stream ("CGM" = Dassault's geometric modeler).
//
// Location: the Data section whose class is `CGMGeom` carries the container GUID (item 3, four
// u32 in the Data stream's byte order); its geometry is the store stream named
// `Cont_<g0>_<g1>_<g2>_<g3>` (lower-case hex, no leading zeros).
//
// Stream header: u8 order (01 little-endian, 00 big-endian), u32 length of the rest, then one or two
// u32-length-prefixed version strings ("  V5R12-HCN:052003", "2015Fall@Prototype", ...), then
// root object ids and a record directory. Four record dialects occur (corpus counts in
// research/cgm/FINDINGS.md):
//
//   varint  (V5R8..V5R19, CXR6/7): directory  e5 0d 06 01 00 u32 level u32 size, then records
//           e5 0d|0f  u8 flags  u16 type  u32 length  u32 id  [i16 extra if 0f]  body[length]
//           chaining exactly over `size` bytes; record bodies use the varint encoding below.
//   compact (2009Fall and later): directory  a8 06 7a u32 level u32 size, then records
//           u8 marker  u8 flags  u8 type  length  id-field  body[length]
//           marker high nibble b -> 1-byte length, a -> 4-byte length; low nibble = length width +
//           id-field width. Type numbers are per file; bodies are not decoded yet.
//   fixed   (CXR4/5, CGA-V2R4): fixed-width records  u32 id  u16 type  b3 00 ...; recognised
//           but not framed yet (the header's padding before the length varies).
//
// Varint integers (record bodies of the varint dialect): a byte >= 0x80 is the value & 0x7f; any
// other byte t is a tag whose bits 3..6 say which of the value's bytes 0..3 follow (absent bytes are
// zero): `30 84 02` = 0x028400. Streams of the oldest varint releases write some counts as raw u32.
#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "cadcore/bytes.hpp"

namespace cadcore::formats::catia::cgm {

enum class Dialect { Varint, Compact, Fixed };
const char* dialect_name(Dialect d);

struct Record {
    std::uint32_t id = 0;      // varint/fixed: object id; compact: the id field as read
    std::uint16_t type = 0;    // varint/fixed: class number; compact: per-file type byte
    std::uint8_t flags = 0;
    std::size_t offset = 0;    // of the record header in the stream
    ByteView body;
};

struct Stream {
    Dialect dialect = Dialect::Varint;
    bool big = false;            // byte order of fixed-width fields
    std::string version;         // first version string of the header
    std::uint32_t level = 0;     // directory level
    std::vector<Record> records;
    std::size_t end = 0;         // first byte after the records (start of the trailer)
};

// The CGMGeom stream of a V5 store; nullopt when the document has none.
std::optional<std::vector<std::uint8_t>> find_stream(ByteView file);

// Splits a CGM stream into records. Throws DecodeError when the stream fits no dialect or a
// record does not chain to the directory's size, and Unsupported for the fixed dialect (dialect
// set, no records).
Stream frame(ByteView stream);

// One varint integer at `k` (advanced). Throws DecodeError past the end.
std::uint32_t varint(ByteView b, std::size_t& k);

}  // namespace cadcore::formats::catia::cgm
