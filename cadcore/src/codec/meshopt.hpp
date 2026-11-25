// EXT_meshopt_compression / KHR_meshopt_compression decoders, written from the extension specification.
#pragma once

#include <cstdint>
#include <vector>

#include "cadcore/bytes.hpp"

namespace cadcore::codec::meshopt {

enum class Mode { Attributes, Triangles, Indices };
enum class Filter { None, Octahedral, Quaternion, Exponential, Color };

// Decodes `count` elements of `stride` bytes. Throws DecodeError(Malformed) on invalid data.
std::vector<std::uint8_t> decode(ByteView src, std::size_t count, std::size_t stride, Mode mode, Filter filter);

}  // namespace cadcore::codec::meshopt
