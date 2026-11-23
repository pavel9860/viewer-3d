// Read-only zip access over an in-memory buffer (miniz), shared by 3MF, FCStd, F3D, SKP.
#pragma once

#include "core/buffer.hpp"
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "cadcore/bytes.hpp"

namespace cadcore::formats {

class ZipArchive {
public:
    // threads: used to verify the CRC of large entries in parallel.
    explicit ZipArchive(ByteView data, unsigned threads = 1);
    ~ZipArchive();
    ZipArchive(const ZipArchive&) = delete;
    ZipArchive& operator=(const ZipArchive&) = delete;

    std::vector<std::string> names() const;
    // Entry contents, or nullopt when absent. Names compare case-insensitively
    // (OPC part names are case-insensitive).
    std::optional<Buffer> read(std::string_view name) const;

private:
    unsigned threads_;
    struct Impl;
    Impl* impl_;
};

}  // namespace cadcore::formats
