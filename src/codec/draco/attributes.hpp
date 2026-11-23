// Draco attribute decoding: attribute decoders, traversal sequencing, prediction schemes
// and the quantization / octahedral-normal transforms.
#pragma once

#include <cstdint>
#include <vector>

#include "codec/draco.hpp"
#include "codec/draco/bitstream.hpp"
#include "codec/draco/connectivity.hpp"

namespace cadcore::codec::draco {

// Decodes all attribute decoders of a mesh or point cloud. `conn` carries the point
// count and, for EdgeBreaker meshes, the corner tables. Returns the attributes in
// stream order with one value tuple per point.
// `kd_tree`: the point cloud uses kd-tree attribute coding.
std::vector<Attribute> decode_attributes(Stream& s, Connectivity& conn, bool is_mesh, bool kd_tree);

}  // namespace cadcore::codec::draco
