// Test-only zip writer with a DEFLATE encoder (fixed-Huffman blocks, greedy LZ77),
// written from RFC 1951 and APPNOTE. Produces inputs with known content for the
// inflate and zip reader tests.
#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace cadcore::testgen {

std::vector<std::uint8_t> deflate_fixed(const std::uint8_t* data, std::size_t size);

struct ZipEntryIn {
    std::string name;
    std::string data;
    bool compress = true;  // method 8 (deflate) or 0 (stored)
};

std::vector<std::uint8_t> write_zip(const std::vector<ZipEntryIn>& entries);

}  // namespace cadcore::testgen
