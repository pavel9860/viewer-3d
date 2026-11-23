#include "formats/catia/store.hpp"

#include <algorithm>
#include <cstring>

#include "cadcore/error.hpp"

namespace cadcore::formats::catia {
namespace {

constexpr std::size_t kHeader = 32;
constexpr std::uint32_t kMaxName = 1024;
constexpr std::uint32_t kMaxSegments = 1u << 20;

std::uint32_t be32(const std::uint8_t* p) { return (std::uint32_t(p[0]) << 24) | (std::uint32_t(p[1]) << 16) | (std::uint32_t(p[2]) << 8) | p[3]; }
std::uint16_t be16(const std::uint8_t* p) { return std::uint16_t((p[0] << 8) | p[1]); }

[[noreturn]] void bad(const std::string& what, std::size_t at) { fail(ErrorCode::Malformed, "CATIA store: " + what, at); }

// UTF-16BE name of `bytes` bytes ending in 00 00; empty when it is not a printable BMP name.
std::string utf16_name(const std::uint8_t* p, std::uint32_t bytes) {
    if (bytes < 4 || bytes % 2 || be16(p + bytes - 2) != 0) return {};
    std::string out;
    for (std::uint32_t i = 0; i + 2 < bytes; i += 2) {
        const std::uint16_t c = be16(p + i);
        if (c < 0x20 || c == 0x7f || (c >= 0xd800 && c <= 0xdfff)) return {};
        if (c < 0x80) out += char(c);
        else if (c < 0x800) { out += char(0xc0 | (c >> 6)); out += char(0x80 | (c & 0x3f)); }
        else { out += char(0xe0 | (c >> 12)); out += char(0x80 | ((c >> 6) & 0x3f)); out += char(0x80 | (c & 0x3f)); }
    }
    return out;
}

// Segment table closing a stream record [r, r + len): the count n sits 20 n + 4 bytes before the end.
// Entries: offset, length, length, logical offset, ordinal. Newer writers list the segments in
// stream order with ordinal 0; older ones (R11-R14) number them 1..n and may list them in any
// order. Returned in stream order: logical offsets 0, l0, l0+l1, ... adding up to `size`.
bool segment_table(const std::uint8_t* r, std::uint32_t len, std::uint32_t size, std::size_t store_end, std::vector<Segment>& out) {
    for (std::uint32_t n = 0; n <= kMaxSegments && 26 + 4 + std::uint64_t(20) * n <= len; ++n) {  // 0: an empty stream
        const std::uint8_t* t = r + len - 4 - 20 * std::size_t(n);
        if (be32(t) != n) continue;
        std::vector<Segment> segs(n);
        std::vector<bool> seen(n, false);
        bool ok = true, numbered = false;
        for (std::uint32_t k = 0; k < n && ok; ++k) {
            const std::uint8_t* e = t + 4 + 20 * std::size_t(k);
            const Segment s{be32(e), be32(e + 4), be32(e + 12)};
            const std::uint32_t ordinal = be32(e + 16);
            if (k == 0) numbered = ordinal != 0;
            const std::uint32_t slot = numbered ? ordinal - 1 : k;
            ok = be32(e + 8) == s.length && std::uint64_t(s.offset) + s.length <= store_end && (ordinal != 0) == numbered && slot < n && !seen[slot];
            if (ok) { seen[slot] = true; segs[slot] = s; }
        }
        std::uint64_t logical = 0;
        for (std::uint32_t k = 0; k < n && ok; ++k) {
            ok = segs[k].logical == logical;
            logical += segs[k].length;
        }
        if (ok && logical == size) {
            out = std::move(segs);
            return true;
        }
    }
    return false;
}

}  // namespace

const StreamInfo* Store::find(std::string_view name) const {
    for (const auto& s : streams)
        if (s.name == name) return &s;
    return nullptr;
}

Store read_store(ByteView file, std::size_t origin) {
    if (origin + kHeader > file.size || std::memcmp(file.data + origin, "V5_CFV", 6) != 0 || file.data[origin + 7] != 0)
        bad("missing V5_CFVn header", origin);
    const std::uint8_t* h = file.data + origin;
    if (h[6] < '0' || h[6] > '9') bad("bad version digit", origin + 6);
    // Bytes 16..23: ff x 8, or (offset, length) of the previous catalog block ("CATIA_V5 CB0001",
    // as the current one at bytes 8..15: catalogs are written alternately). It is not needed to
    // locate streams and is not read; it may lie at or past the end of a file that shrank since.
    if (!std::all_of(h + 16, h + 24, [](std::uint8_t b) { return b == 0xff; })) {
        const std::uint32_t boff = be32(h + 16), blen = be32(h + 20);
        constexpr char kBlock[] = "CATIA_V5 CB0001";
        const bool past_end = origin + std::uint64_t(boff) >= file.size;
        if (!past_end && (blen < sizeof kBlock || origin + std::uint64_t(boff) + blen > file.size || std::memcmp(h + boff, kBlock, sizeof kBlock) != 0))
            bad("header bytes 16..23 are neither ff padding nor a catalog block", origin + 16);
    }
    Store st;
    st.origin = origin;
    st.version = h[6] - '0';
    const std::uint32_t off = be32(h + 8), len = be32(h + 12);
    if (off < kHeader || origin + std::uint64_t(off) + len > file.size) bad("catalog outside the file", origin + 8);
    const std::uint8_t* c = h + off;
    const std::size_t store_end = std::size_t(off) + len;  // segment offsets are relative to the origin
    std::size_t p = 0;
    while (p + 4 <= len) {
        const std::uint32_t nl = be32(c + p);
        if (nl >= 4 && nl <= kMaxName && p + 4 + nl + 26 <= len) {
            std::string name = utf16_name(c + p + 4, nl);
            const std::size_t r = p + 4 + nl;
            const std::uint32_t rl = name.empty() ? 0 : be32(c + r);
            const std::uint16_t kind = rl ? be16(c + r + 8) : 0;
            if (rl >= 26 && r + rl <= len && (kind == 1 || kind == 2)) {
                if (kind == 2) {
                    StreamInfo s{std::move(name), be32(c + r + 4), be32(c + r + 14), {}};
                    if (!segment_table(c + r, rl, s.size, store_end, s.segments))
                        bad("stream '" + s.name + "' has no valid segment table", origin + off + r);
                    st.streams.push_back(std::move(s));
                }
                p = r + rl;
                continue;
            }
        }
        ++p;
    }
    return st;
}

std::vector<std::uint8_t> read_stream(ByteView file, const Store& store, const StreamInfo& stream) {
    std::vector<std::uint8_t> out;
    out.reserve(stream.size);
    for (const auto& s : stream.segments) {
        const std::size_t at = store.origin + s.offset;
        out.insert(out.end(), file.data + at, file.data + at + s.length);
    }
    return out;
}

}  // namespace cadcore::formats::catia
