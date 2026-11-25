// Draco geometry bitstream decoder interface (written from the Draco bitstream specification).
#pragma once

#include <cstdint>
#include <vector>

#include "cadcore/bytes.hpp"

namespace cadcore::codec::draco {

struct Attribute {
    std::uint32_t unique_id = 0;
    int type = 0;            // Draco attribute type: 0 POSITION, 1 NORMAL, 2 COLOR, 3 TEX_COORD, 4 GENERIC
    int components = 0;
    std::vector<float> values;  // num_points * components, dequantized / normalised to float
};

struct Mesh {
    std::uint32_t num_points = 0;
    std::vector<std::uint32_t> faces;          // 3 point indices per triangle (empty for point clouds)
    std::vector<Attribute> attributes;         // in stream order
    // The attribute with this unique id; nullptr if absent. Throws if the id is ambiguous
    // (streams older than 1.3 store 16-bit ids that need not be unique).
    const Attribute* find(std::uint32_t unique_id) const;
};

// Throws DecodeError (Malformed / Unsupported) on invalid or unsupported bitstreams.
Mesh decode(ByteView data);

}  // namespace cadcore::codec::draco
