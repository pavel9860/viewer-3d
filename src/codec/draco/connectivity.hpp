// Draco mesh connectivity: sequential index lists and EdgeBreaker (standard, predictive, valence).
#pragma once

#include <cstdint>
#include <vector>

#include "codec/draco/bitstream.hpp"
#include "codec/draco/corner_table.hpp"

namespace cadcore::codec::draco {

struct AttributeData {               // one per EdgeBreaker "attribute data" (seam set)
    AttributeCornerTable table;
    std::int32_t decoder_id = -1;     // attribute decoder that owns this data
    bool connectivity_used = true;    // false when the owning decoder traverses positions
};

struct Connectivity {
    Connectivity() = default;
    Connectivity(const Connectivity&) = delete;
    Connectivity& operator=(const Connectivity&) = delete;
    bool edgebreaker = false;
    std::uint32_t num_points = 0;
    std::vector<std::uint32_t> faces;  // 3 point ids per face (= corner -> point)
    CornerTable table;                 // EdgeBreaker only
    std::vector<AttributeData> attribute_data;
    std::int32_t pos_data_decoder_id = -1;
};

// Both fill `out` in place: attribute tables point into out.table, so it must not move.
void decode_sequential_connectivity(Stream& s, Connectivity& out);
void decode_edgebreaker_connectivity(Stream& s, Connectivity& out);

}  // namespace cadcore::codec::draco
