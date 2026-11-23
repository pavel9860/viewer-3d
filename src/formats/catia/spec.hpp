// CATIA V5 object containers with inline names ("SPEC", "OSMN", "OSMB": R5..R17), the forms
// that preceded the framed OSM streams.
//
// They sit in the store's "Data" stream, a sequence of sections
//   u32 id, u32 length8                    length in 8-byte units, this header included
//   body: u16 1, u16 kind, u32 index       then items  u16 tag, u16 flags, u32 length8, data
// Items: tag 1 (8 bytes of container data), 2 (the container's class names), 3 (its GUID),
// 4 / 5 payloads: object containers in this grammar ("SPEC", R11..R16; "OSMN", "OSMB", R16..R17;
// read here), framed OSM streams ("OSMX", formats/catia/osm.hpp), user ids.
//
// A SPEC container, with names written inline the first time (u8 n then n-1 bytes, or 00 u32 n
// then n bytes; an empty name is not numbered) and referenced afterwards by number as in OSM:
//   "SPEC" <version> <compact> [fe when version >= 19] <highest object id>
//   "OSMB" <version> <compact> <name>* fe <highest object id>
//   "OSMN" <name>* <version> <compact> fe <highest object id>
//   object*                                   until de
//   de, one relation block per object
//   property tables (one per object, in some files one more), fe, then padding to the section end
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
//   header token* fe (compact | e5 <blob> | fd <compact>)
//   (<attribute> (<count> <index < count> <value>* | <value>))* fe
//   value: compact | e5 <blob> <compact> | 7f <a> <n> <id>*n <index>*(n+1)  (path, not followed)
// Property table: <count> (<name> <type> <value>)*: type 2, 3: compact; 0, 14: none;
//   13: byte e7+n (n = 1..8; 485 values), a form byte, n bytes.
#pragma once

#include <vector>

#include "cadcore/bytes.hpp"
#include "formats/catia/osm.hpp"

namespace cadcore::formats::catia {

// The SPEC containers of the file's Data stream, as object streams (relations and headers
// decoded, records parsed); empty when the file has none. Throws DecodeError on any residue.
std::vector<Stream> parse_spec_streams(ByteView file);

}  // namespace cadcore::formats::catia
