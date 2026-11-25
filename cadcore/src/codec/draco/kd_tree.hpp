// Draco kd-tree point coding ("DynamicIntegerPointsKdTreeDecoder", bitstream 2.3).
#pragma once

#include <cstdint>
#include <vector>

#include "codec/draco/bitstream.hpp"

namespace cadcore::codec::draco {

// Decodes `expected_points` points of `dimension` unsigned components each, in stream
// order, into `out` (expected_points * dimension values). `level` is 0..6.
void decode_kd_tree_points(Stream& s, int level, std::uint32_t dimension, std::uint32_t expected_points, std::vector<std::uint32_t>& out);

}  // namespace cadcore::codec::draco
