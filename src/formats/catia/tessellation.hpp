// CATIA V5 stored tessellation: the scene-graph stream of a store ("Skeleton" in a binary
// .3DRep, the "MainDataStream" of "SceneGraph" in a CGR).
//
// Nodes are tagged 04 04 fd ff. A framed node is preceded by  ff u32le length 00 <type>  (the
// length counted from the ff byte) and bounds its content exactly. Its header is
//   04 04 fd ff, 5 bytes, f32 radius, 4 bytes, f32 x3 centre, u32, then at +33
//   ff ff <u8 0 or 2> 00 <u8 lists>
// lists = 1: a geometry node; the point list and face list follow. lists = 0 followed by 00: a
// polyline node, whose face list (without points) holds explicit polylines. lists = 0 followed
// by a P1g child count: a group node, whose children follow.
//
//   points     P1g count, P1g normal_count, u8 has_tex, u8 0, u8 0,
//              [u8 dims, f32[dims x count] texture coordinates when has_tex],
//              encoding (CFV2: u8; CFV3: one bit): 0 = 2-bit prefix table over the points
//              (ceil(count/4) bytes, P1g words, f32[words]; opcodes 0 xyz, 3 yz, 2 z, 1 repeat),
//              1 = raw f32[3 count]
//   normals    when normal_count > 0: P1g m, u16[m], ceil(normal_count/2) bytes
//   face list  f32 sag, u8 nLOD, nLOD x {record, f32 lod_sag, attribute}, u8 groups,
//              groups x { u8 attribute_mode (2 = per item, 1 = per group), u8 type,
//                         P1g items, items, attribute(s) }
//                type 0x0a  facet records: 01 <0x40|tri 1|strip 2|fan 4|normal 8>
//                           [P1g per kind] ff u32le n [f32 x3 normal] strip lengths,
//                           fan lengths, n references
//                types 0..2 polylines: items 02 P1g n, n references (indexed), or
//                           03 P1g n u8 1, f32[3 n] (explicit points)
//   attribute  8 bytes <flags> 04 04 ff A B G R (faces; flags 0x10 solid, 0x80 transparent)
//              or <flags> 24 04 ff A B G R (edges)
//
// The two generations differ only in how lengths and references are stored:
//   CFV2  loop counts P1g; references big-endian, 1 byte when count < 256, else 2 (3, 4 for
//         larger lists); the point encoding is a byte.
//   CFV3  fields of bit_length(count) bits, LSB first, from one bit stream that runs through
//         the whole stream and is refilled with a u32le word taken at the byte position where
//         the reader stands when it runs out of bits; the point encoding is one bit.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "cadcore/bytes.hpp"
#include "cadcore/scene.hpp"

namespace cadcore::formats::catia {

enum class Generation { Cfv2, Cfv3 };

struct TessFace {
    std::vector<std::uint32_t> triangles;  // 3 per triangle, degenerate ones dropped
    std::uint8_t flags = 0;                // 0x10 solid, 0x80 transparent
    Color color;
};

struct TessLines {
    Color color;
    std::vector<std::vector<float>> polylines;  // x, y, z per point
};

struct TessBody {
    std::size_t offset = 0;        // of the node in the stream
    std::vector<float> positions;  // x, y, z per point (empty for a polyline node)
    std::vector<TessFace> faces;
    std::vector<TessLines> lines;  // edge and wire polylines, resolved to coordinates
};

// Every geometry and polyline node of the stream, in stream order. Throws on any byte of a
// node that does not fit the grammar, or a node that does not end at its framed length.
std::vector<TessBody> decode_tessellation(ByteView stream, Generation generation);

// Adds one mesh per body with facets (face groups carry the face colours; the body kind comes
// from the solid flag) and one polyline per edge or wire, all attached to `node`.
void add_to_scene(std::vector<TessBody>&& bodies, const std::string& name, Scene& scene, std::uint32_t node);

}  // namespace cadcore::formats::catia
