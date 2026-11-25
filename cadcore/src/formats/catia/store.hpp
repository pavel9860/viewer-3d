// CATIA V5 storage container (the "V5_CFV2" / "V5_CFV3" store of CATPart, CATProduct, CGR and
// the binary .3DRep files of 3DXML archives).
//
//   header   "V5_CFVn\0"  u32be catalog_offset  u32be catalog_length  (ff x 8 | u32be block_offset
//            u32be block_length of a separate "CATIA_V5 CB0001" catalog block)          (32 bytes)
//   catalog  at origin + catalog_offset, catalog_length bytes
//
// The catalog names the store's streams. A named stream is
//   u32be name_bytes, name (UTF-16BE, 00 00 terminated), then a record:
//   +0 u32be record_length (counted from +0)   +4 u32be id   +8 u16be kind (2 = stream)
//   +10 u32be   +14 u32be stream_size   ...   then, closing the record,
//   u32be segment_count, segment_count x { u32be offset (from origin), u32be length,
//   u32be length (again), u32be logical_offset, u32be ordinal }.
// Newer writers list segments in stream order with ordinal 0; older ones (R11-R14) number them
// 1..n and may list them in any order. In stream order the logical offsets run 0, l0, l0+l1, ...
// and the lengths add up to stream_size. Writers interleave the segments of several streams in the
// file (e.g. Skeleton and SurfacicReps, ~16 KB at a time), so a stream must be reassembled
// before it is parsed.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "cadcore/bytes.hpp"

namespace cadcore::formats::catia {

struct Segment {
    std::uint32_t offset, length, logical;
};

struct StreamInfo {
    std::string name;  // UTF-16BE name converted to UTF-8
    std::uint32_t id = 0;
    std::uint32_t size = 0;
    std::vector<Segment> segments;
};

struct Store {
    std::size_t origin = 0;   // offset of the "V5_CFVn" header in the file
    int version = 0;          // n of V5_CFVn
    std::vector<StreamInfo> streams;

    const StreamInfo* find(std::string_view name) const;
};

// The store whose header starts at `origin`. Throws on a malformed header or catalog.
Store read_store(ByteView file, std::size_t origin = 0);

// The stream's bytes, segments joined in logical order.
std::vector<std::uint8_t> read_stream(ByteView file, const Store& store, const StreamInfo& stream);

}  // namespace cadcore::formats::catia
