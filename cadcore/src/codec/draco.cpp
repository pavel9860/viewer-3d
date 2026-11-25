// Draco geometry decoder (meshes and point clouds).
//
// Written from the Draco bitstream specification and checked against streams produced by
// the Draco encoder; where the two disagree (bit order of bit runs, boundary tests, seam
// handling, metadata entry sizes, legacy layouts) the decoder follows the encoded data.
// Module layout: draco/bitstream (bytes, bit runs, rANS/rABS), draco/corner_table,
// draco/connectivity (sequential + EdgeBreaker), draco/attributes (decoders, prediction).
#include "codec/draco.hpp"

#include <cstring>
#include <string>

#include "codec/draco/attributes.hpp"
#include "codec/draco/bitstream.hpp"
#include "codec/draco/connectivity.hpp"

namespace cadcore::codec::draco {
namespace {

constexpr int POINT_CLOUD = 0, TRIANGULAR_MESH = 1;
constexpr int MESH_SEQUENTIAL = 0, MESH_EDGEBREAKER = 1;
constexpr int POINT_CLOUD_SEQUENTIAL = 0, POINT_CLOUD_KD_TREE = 1;
constexpr int METADATA_FLAG_MASK = 0x8000;

// Metadata: entries (name u8-length string, varint-sized value), then named sub-metadata.
void skip_metadata(Stream& s, int depth) {
    if (depth > 1000) bad("metadata nested too deeply");
    const std::uint32_t entries = s.varint32();
    for (std::uint32_t i = 0; i < entries; ++i) {
        s.skip(s.u8());
        const std::uint32_t size = s.varint32();
        if (size == 0 || size > s.remaining()) bad("invalid metadata entry size");
        s.skip(size);
    }
    const std::uint32_t subs = s.varint32();
    if (subs > s.remaining()) bad("metadata sub-entry count exceeds the stream");
    for (std::uint32_t i = 0; i < subs; ++i) {
        s.skip(s.u8());
        skip_metadata(s, depth + 1);
    }
}

void decode_metadata(Stream& s) {
    const std::uint32_t num_att = s.varint32();
    for (std::uint32_t i = 0; i < num_att; ++i) {
        s.varint32();  // attribute unique id
        skip_metadata(s, 0);
    }
    skip_metadata(s, 0);
}

}  // namespace

const Attribute* Mesh::find(std::uint32_t unique_id) const {
    const Attribute* r = nullptr;
    for (const Attribute& a : attributes)
        if (a.unique_id == unique_id) {
            if (r) bad("attribute unique id " + std::to_string(unique_id) + " is not unique");
            r = &a;
        }
    return r;
}

Mesh decode(ByteView data) {
    Stream s(data.data, data.size);
    if (data.size < 11 || std::memcmp(s.bytes(5), "DRACO", 5) != 0) bad("missing DRACO signature");
    const int major = s.u8(), minor = s.u8();
    const int encoder_type = s.u8(), method = s.u8();
    const int flags = s.u16();
    if (encoder_type != POINT_CLOUD && encoder_type != TRIANGULAR_MESH) bad("unknown geometry type " + std::to_string(encoder_type));
    const int max_major = 2, max_minor = encoder_type == POINT_CLOUD ? 3 : 2;
    if (major < 1 || major > max_major || (major == max_major && minor > max_minor))
        unsupported("bitstream version " + std::to_string(major) + "." + std::to_string(minor));
    s.ver = version(major, minor);
    if (s.ver >= version(1, 3) && (flags & METADATA_FLAG_MASK)) decode_metadata(s);

    Connectivity conn;
    if (encoder_type == TRIANGULAR_MESH) {
        if (method == MESH_SEQUENTIAL) decode_sequential_connectivity(s, conn);
        else if (method == MESH_EDGEBREAKER) decode_edgebreaker_connectivity(s, conn);
        else bad("unknown mesh encoding method " + std::to_string(method));
    } else {
        if (method != POINT_CLOUD_SEQUENTIAL && method != POINT_CLOUD_KD_TREE) bad("unknown point cloud encoding method " + std::to_string(method));
        const std::int32_t n = s.i32();
        if (n < 0) bad("negative point count");
        conn.num_points = std::uint32_t(n);
    }

    Mesh m;
    m.num_points = conn.num_points;
    m.attributes = decode_attributes(s, conn, encoder_type == TRIANGULAR_MESH, encoder_type == POINT_CLOUD && method == POINT_CLOUD_KD_TREE);
    m.faces = std::move(conn.faces);
    for (const std::uint32_t p : m.faces)
        if (p >= m.num_points) bad("face index out of range");
    return m;
}

}  // namespace cadcore::codec::draco
