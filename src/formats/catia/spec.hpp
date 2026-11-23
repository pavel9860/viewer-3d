// CATIA V5 object containers with inline names ("SPEC", "OSMN", "OSMB": R5..R17), the forms
// that preceded the framed OSM streams.
//
// They are payload items (tag 4 or 5) of the store's "Data" stream (formats/catia/data.hpp).
//
// A SPEC container, with names written inline the first time (u8 n then n-1 bytes, or 00 u32 n
// then n bytes; an empty name is not numbered) and referenced afterwards by number as in OSM:
//   "SPEC" <version> <compact> [fe when version >= 19] <highest object id>
//   "OSMB" <version> <compact> <name>* fe <highest object id>
//   "OSMN" <name>* <version> <compact> [82 <n> (<link blob> <link blob>)*n] fe <highest object id>
//   object*                                   until de
//   de, one relation block per object
//   property tables (one per object, in some files one more), fe, then padding to the section end
// Byte order: the u32/f64 fields inside containers follow the writer's
// platform (big-endian in files from big-endian UNIX systems, 8 files of the corpus); the object id
// after ea is little-endian in both.
// Object (version-dependent fields: struct Layout in spec.cpp, tabulated over containers of
// versions 5..26 in ~2,200 files):
//   {00 | 02} <name> <handle> <stamp> [<stamp> from version 13]
//   (fe f6 <GUID 16> | 80 <5 bytes>) <kind u8> [u8 from version 19] [<compact> when kind is odd]
//   ea <u32 id> [<compact n> <name>*n from version 10]
//   startup: <compact> | e5 <blob> <compact> | fd <compact>
//   chunk* fe   (the chunk encoding of OSM records with inline names; chunks close with fe from
//               version 10; the first object (the catalog manager, in all 79 containers of
//               versions 10+) has no closing fe: its list ends where the next object or de begins)
//   flag 02: a parameter value follows: <compact> <type> [<compact> from version 13]
//            (fe | <unit name>) <value by type: 07, 15 number; 12 name; 11, 1c compact>
// Relation block (the relation header grammar of OSM streams before version 35):
//   <head u8>, form f = head & 7; f != 0 and fe next: inline class fe (fe | e5 <blob> <compact>) <name> fe
//   header token* [fe] (compact | e5 <blob> | fd <compact>; without fe the header ends at the first
//   byte below 80, an attribute name: 5 KWE check blocks)
//   (<attribute> (<count> <index < count> <value>* | <value>))* fe   (an object without chunks: the
//   empty list is fe or left out, ~200 blocks; the next head is never fe)
//   value: compact | e5 <blob> <compact> | 7f <a> <n> <id>*n <index>*(n+1)  (path, not followed)
// Property table: <count> (<name> <type> <value>)*: type 2, 3, 5: compact; 11: name; 0, 14: none;
//   13: byte e7+n (n = 1..8; 485 values), a form byte, n bytes.
#pragma once

#include <string>
#include <vector>

#include "cadcore/bytes.hpp"
#include "cadcore/error.hpp"
#include "formats/catia/osm.hpp"

namespace cadcore::formats::catia {

// A container that could not be decoded: its place in the Data stream and the reason.
struct FailedContainer {
    std::size_t offset = 0;  // in the Data stream
    std::size_t size = 0;
    bool has_products = false;  // its bytes name `_PartNumber` (it may hold product references)
    ErrorCode code = ErrorCode::Malformed;
    std::string message;                    // without the code prefix
    std::size_t error_offset = DecodeError::npos;  // where decoding stopped
};

struct SpecStreams {
    std::vector<Stream> streams;          // decoded containers (relations and headers decoded, records parsed)
    std::vector<FailedContainer> failed;  // each container is decoded on its own; residue or an unknown form fails only it
};

// The SPEC containers of the file's Data stream; empty when the file has none. Throws DecodeError
// when the Data stream's section layout itself is malformed.
SpecStreams parse_spec_streams(ByteView file);

}  // namespace cadcore::formats::catia
