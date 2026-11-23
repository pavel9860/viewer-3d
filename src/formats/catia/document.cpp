#include "formats/catia/document.hpp"

#include <cstring>

#include "cadcore/catia.hpp"
#include "cadcore/error.hpp"

namespace cadcore::formats::catia {
namespace {

std::uint32_t get32(ByteView f, std::size_t at, bool big) {
    const std::uint8_t* p = f.data + at;
    return big ? std::uint32_t(p[0]) << 24 | std::uint32_t(p[1]) << 16 | std::uint32_t(p[2]) << 8 | p[3]
               : std::uint32_t(p[3]) << 24 | std::uint32_t(p[2]) << 16 | std::uint32_t(p[1]) << 8 | p[0];
}

std::uint16_t get16(ByteView f, std::size_t at, bool big) {
    return big ? std::uint16_t(f[at] << 8 | f[at + 1]) : std::uint16_t(f[at + 1] << 8 | f[at]);
}

}  // namespace

std::optional<Uuid> document_uuid(ByteView f) {
    // The table is written in the byte order of the writer's platform (big-endian on UNIX
    // workstations); the GUID itself is a byte string in both.
    constexpr std::size_t kScan = 2u << 20;  // the table sits near the start of the file
    const std::size_t limit = f.size < kScan ? f.size : kScan;
    std::optional<Uuid> found;
    int count = 0;
    for (const bool big : {false, true}) {
        const std::string_view marker = big ? std::string_view("\x00\x02\x20\x81", 4) : std::string_view("\x02\x00\x81\x20", 4);
        const std::string_view tail_marker = big ? std::string_view("\x00\x03\x00\xf7\x00\x00\x00\x03", 8) : std::string_view("\x03\x00\xf7\x00\x03\x00\x00\x00", 8);
        for (std::size_t a = f.find(marker, 32, limit); a != ByteView::npos; a = f.find(marker, a + 1, limit)) {
            const std::size_t s = a - 32;
            if (s + 40 > f.size) continue;
            if (get16(f, s + 16, big) != 1 || get16(f, s + 18, big) != 0x6c || get32(f, s + 20, big) != 2 || get32(f, s + 28, big) != 0) continue;
            const std::uint32_t f10 = get32(f, s + 36, big);
            const std::size_t names = s + 40;
            const std::size_t z1 = f.find(std::string_view("\0", 1), names, f.size);
            if (z1 == ByteView::npos) continue;
            const std::size_t z2 = f.find(std::string_view("\0", 1), z1 + 1, f.size);
            if (z2 == ByteView::npos) continue;
            const std::size_t n1 = z1 - names, n2 = z2 - z1 - 1;
            if (f10 != (n1 + n2 + 2 + 7) / 8 + 1) continue;
            const std::size_t tail = names + 8 * (std::size_t(f10) - 1);
            if (tail + 24 > f.size || std::memcmp(f.data + tail, tail_marker.data(), 8) != 0) continue;
            bool padded = true;
            for (std::size_t i = z2 + 1; i < tail; ++i) padded = padded && f[i] == 0;
            if (!padded) continue;
            if (n1 == 11 && std::memcmp(f.data + names, "CATProdCont", 11) == 0) {
                Uuid u;
                std::memcpy(u.data(), f.data + tail + 8, 16);
                found = u;
                ++count;
            }
        }
    }
    if (count != 1) return std::nullopt;
    return found;
}

Uuid link_target(ByteView blob, std::size_t at) {
    if (blob.size < 18 || blob[0] != 0x11 || blob[1] != 0x03) fail(ErrorCode::Malformed, "CATIA document link with an unknown layout", at);
    Uuid u;
    std::memcpy(u.data(), blob.data + 2, 16);
    return u;
}

}  // namespace cadcore::formats::catia

namespace cadcore {

std::optional<std::string> catia_document_id(ByteView file) {
    const auto u = formats::catia::document_uuid(file);
    if (!u) return std::nullopt;
    static const char* k = "0123456789abcdef";
    std::string s;
    for (auto b : *u) { s += k[b >> 4]; s += k[b & 15]; }
    return s;
}

}  // namespace cadcore
