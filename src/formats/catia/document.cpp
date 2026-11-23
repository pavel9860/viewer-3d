#include "formats/catia/document.hpp"

#include <cstring>

#include "cadcore/catia.hpp"
#include "cadcore/error.hpp"

namespace cadcore::formats::catia {
namespace {

std::uint32_t le32(ByteView f, std::size_t at) {
    std::uint32_t x;
    std::memcpy(&x, f.data + at, 4);
    return x;
}

}  // namespace

std::optional<Uuid> document_uuid(ByteView f) {
    constexpr std::size_t kScan = 2u << 20;  // the table sits near the start of the file
    const std::size_t limit = f.size < kScan ? f.size : kScan;
    std::optional<Uuid> found;
    int count = 0;
    for (std::size_t a = f.find(std::string_view("\x02\x00\x81\x20", 4), 32, limit); a != ByteView::npos;
         a = f.find(std::string_view("\x02\x00\x81\x20", 4), a + 1, limit)) {
        const std::size_t s = a - 32;
        if (s + 40 > f.size) continue;
        if (le32(f, s + 16) != 0x006c0001 || le32(f, s + 20) != 2 || le32(f, s + 28) != 0) continue;
        const std::uint32_t f10 = le32(f, s + 36);
        const std::size_t names = s + 40;
        const std::size_t z1 = f.find(std::string_view("\0", 1), names, f.size);
        if (z1 == ByteView::npos) continue;
        const std::size_t z2 = f.find(std::string_view("\0", 1), z1 + 1, f.size);
        if (z2 == ByteView::npos) continue;
        const std::size_t n1 = z1 - names, n2 = z2 - z1 - 1;
        if (f10 != (n1 + n2 + 2 + 7) / 8 + 1) continue;
        const std::size_t tail = names + 8 * (std::size_t(f10) - 1);
        if (tail + 24 > f.size || std::memcmp(f.data + tail, "\x03\x00\xf7\x00\x03\x00\x00\x00", 8) != 0) continue;
        for (std::size_t i = z2 + 1; i < tail; ++i)
            if (f[i] != 0) goto next;  // padding must be zero
        if (n1 == 11 && std::memcmp(f.data + names, "CATProdCont", 11) == 0) {
            Uuid u;
            std::memcpy(u.data(), f.data + tail + 8, 16);
            found = u;
            ++count;
        }
    next:;
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
