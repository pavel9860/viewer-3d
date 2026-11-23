// The "Data" stream of a CATIA V5 store: one section per container of the document.
//
//   section: u32 id, u32 length8 (8-byte units, header included), u16 1, u16 kind, u32 index, item*
//   item:    u16 tag, u16 flags, u32 length8, body
// Tags: 1 (8 bytes of container data), 2 (the container's class and superclass names, each
// NUL-terminated), 3 (its GUID), 4 / 5 (payloads: object containers "SPEC", "OSMN", "OSMB",
// formats/catia/spec.hpp; framed OSM streams "OSMX", formats/catia/osm.hpp; user ids).
// Fields are in the byte order of the writer's platform: big-endian in files written on
// big-endian UNIX systems (8 files of the corpus), told by the first section's u16 1.
#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "cadcore/bytes.hpp"

namespace cadcore::formats::catia {

struct DataItem {
    std::uint16_t tag = 0;
    std::size_t offset = 0;  // of the body, in the Data stream
    ByteView body;
};

struct DataSection {
    std::size_t offset = 0;  // in the Data stream
    std::uint16_t kind = 0;
    std::uint32_t index = 0;
    std::string cls, super;  // from item 2; empty when the section has none
    std::vector<DataItem> items;
};

struct DataStream {
    std::shared_ptr<const std::vector<std::uint8_t>> bytes;  // the stream, which the views point into
    bool big = false;
    std::vector<DataSection> sections;
};

// The file's Data stream; nullopt when the file is not a V5 store or has none. Throws DecodeError
// when its section layout is malformed.
std::optional<DataStream> read_data_stream(ByteView file);

}  // namespace cadcore::formats::catia
