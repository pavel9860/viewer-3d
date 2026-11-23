#include "testgen/zipwriter.hpp"

#include <array>

#include "codec/inflate.hpp"

namespace cadcore::testgen {
namespace {

class BitWriter {
public:
    void put(std::uint32_t bits, int n) {  // LSB-first
        buf_ |= std::uint64_t(bits) << cnt_;
        cnt_ += n;
        while (cnt_ >= 8) { out.push_back(std::uint8_t(buf_)); buf_ >>= 8; cnt_ -= 8; }
    }
    void put_code(std::uint32_t code, int len) {  // Huffman codes go MSB-first
        std::uint32_t rev = 0;
        for (int i = 0; i < len; ++i) rev |= ((code >> i) & 1) << (len - 1 - i);
        put(rev, len);
    }
    void flush() { if (cnt_) { out.push_back(std::uint8_t(buf_)); buf_ = 0; cnt_ = 0; } }
    std::vector<std::uint8_t> out;

private:
    std::uint64_t buf_ = 0;
    int cnt_ = 0;
};

void put_literal(BitWriter& w, int sym) {
    if (sym < 144) w.put_code(0x30 + sym, 8);
    else if (sym < 256) w.put_code(0x190 + (sym - 144), 9);
    else if (sym < 280) w.put_code(sym - 256, 7);
    else w.put_code(0xC0 + (sym - 280), 8);
}

constexpr std::uint16_t kLenBase[29] = {3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258};
constexpr std::uint8_t kLenExtra[29] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
constexpr std::uint16_t kDistBase[30] = {1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577};
constexpr std::uint8_t kDistExtra[30] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};

void put_match(BitWriter& w, std::size_t len, std::size_t dist) {
    int li = 28;
    while (kLenBase[li] > len) --li;
    put_literal(w, 257 + li);
    w.put(std::uint32_t(len - kLenBase[li]), kLenExtra[li]);
    int di = 29;
    while (kDistBase[di] > dist) --di;
    w.put_code(std::uint32_t(di), 5);
    w.put(std::uint32_t(dist - kDistBase[di]), kDistExtra[di]);
}

void put_u16(std::vector<std::uint8_t>& b, std::uint32_t v) { b.push_back(std::uint8_t(v)); b.push_back(std::uint8_t(v >> 8)); }
void put_u32(std::vector<std::uint8_t>& b, std::uint32_t v) { for (int i = 0; i < 4; ++i) b.push_back(std::uint8_t(v >> (8 * i))); }

}  // namespace

std::vector<std::uint8_t> deflate_fixed(const std::uint8_t* d, std::size_t n) {
    BitWriter w;
    w.put(1, 1);  // final block
    w.put(1, 2);  // fixed Huffman
    constexpr std::size_t kWindow = 32768, kHash = 1 << 15;
    std::vector<std::int64_t> head(kHash, -1), prev(n, -1);
    auto hash = [&](std::size_t i) { return ((d[i] << 10) ^ (d[i + 1] << 5) ^ d[i + 2]) & (kHash - 1); };
    std::size_t i = 0;
    while (i < n) {
        std::size_t best_len = 0, best_dist = 0;
        if (i + 2 < n) {
            const std::size_t h = hash(i);
            std::int64_t cand = head[h];
            for (int chain = 0; cand >= 0 && chain < 64 && i - std::size_t(cand) <= kWindow; ++chain, cand = prev[std::size_t(cand)]) {
                std::size_t l = 0;
                while (l < 258 && i + l < n && d[std::size_t(cand) + l] == d[i + l]) ++l;
                if (l > best_len) { best_len = l; best_dist = i - std::size_t(cand); }
            }
            prev[i] = head[h];
            head[h] = std::int64_t(i);
        }
        if (best_len >= 3) {
            put_match(w, best_len, best_dist);
            for (std::size_t k = 1; k < best_len && i + k + 2 < n; ++k) {
                const std::size_t h = hash(i + k);
                prev[i + k] = head[h];
                head[h] = std::int64_t(i + k);
            }
            i += best_len;
        } else {
            put_literal(w, d[i]);
            ++i;
        }
    }
    put_literal(w, 256);
    w.flush();
    return std::move(w.out);
}

std::vector<std::uint8_t> write_zip(const std::vector<ZipEntryIn>& entries) {
    std::vector<std::uint8_t> out, cd;
    for (const auto& e : entries) {
        const auto* raw = reinterpret_cast<const std::uint8_t*>(e.data.data());
        const std::uint32_t crc = codec::crc32(raw, e.data.size());
        std::vector<std::uint8_t> packed = e.compress ? deflate_fixed(raw, e.data.size()) : std::vector<std::uint8_t>(raw, raw + e.data.size());
        const std::uint32_t offset = std::uint32_t(out.size());
        const std::uint16_t method = e.compress ? 8 : 0;
        put_u32(out, 0x04034b50);
        put_u16(out, 20); put_u16(out, 0); put_u16(out, method); put_u16(out, 0); put_u16(out, 0x21);
        put_u32(out, crc); put_u32(out, std::uint32_t(packed.size())); put_u32(out, std::uint32_t(e.data.size()));
        put_u16(out, std::uint16_t(e.name.size())); put_u16(out, 0);
        out.insert(out.end(), e.name.begin(), e.name.end());
        out.insert(out.end(), packed.begin(), packed.end());
        put_u32(cd, 0x02014b50);
        put_u16(cd, 20); put_u16(cd, 20); put_u16(cd, 0); put_u16(cd, method); put_u16(cd, 0); put_u16(cd, 0x21);
        put_u32(cd, crc); put_u32(cd, std::uint32_t(packed.size())); put_u32(cd, std::uint32_t(e.data.size()));
        put_u16(cd, std::uint16_t(e.name.size())); put_u16(cd, 0); put_u16(cd, 0); put_u16(cd, 0); put_u16(cd, 0);
        put_u32(cd, 0); put_u32(cd, offset);
        cd.insert(cd.end(), e.name.begin(), e.name.end());
    }
    const std::uint32_t cd_offset = std::uint32_t(out.size());
    out.insert(out.end(), cd.begin(), cd.end());
    put_u32(out, 0x06054b50);
    put_u16(out, 0); put_u16(out, 0);
    put_u16(out, std::uint16_t(entries.size())); put_u16(out, std::uint16_t(entries.size()));
    put_u32(out, std::uint32_t(cd.size())); put_u32(out, cd_offset);
    put_u16(out, 0);
    return out;
}

}  // namespace cadcore::testgen
