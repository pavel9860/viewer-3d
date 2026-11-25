// DEFLATE (RFC 1951) decoder and CRC-32 (ISO 3309 / ITU-T V.42, as used by zip and PNG).
#pragma once

#include <cstdint>
#include <vector>

#include "cadcore/bytes.hpp"
#include "core/buffer.hpp"

namespace cadcore::codec {

// Decodes a raw DEFLATE stream. `expected_size` (when known) pre-sizes the output;
// decoding fails with Malformed if more than `max_output` bytes would be produced.
Buffer inflate(ByteView raw, std::size_t expected_size = 0, std::size_t max_output = std::size_t(1) << 32);

// Same result as inflate (or the same error), decoded on up to `threads` threads for large streams:
// parts start at block boundaries found by search, decode with an unknown preceding window, and
// are accepted only when each part ends exactly where the next begins; otherwise it decodes
// sequentially.
Buffer inflate_parallel(ByteView raw, std::size_t expected_size, std::size_t max_output, unsigned threads);

std::uint32_t crc32(const std::uint8_t* data, std::size_t size, std::uint32_t crc = 0);
// CRC of the concatenation A B from crc(A), crc(B) and B's length (GF(2) matrix powers, as zlib).
std::uint32_t crc32_combine(std::uint32_t crc_a, std::uint32_t crc_b, std::size_t len_b);
// crc32 of large data computed in pieces on up to `threads` threads and combined.
std::uint32_t crc32_parallel(const std::uint8_t* data, std::size_t size, unsigned threads);

}  // namespace cadcore::codec
