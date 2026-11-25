// DEFLATE decoder written from RFC 1951.
//
// Bits are read LSB-first. Huffman codes are canonical (section 3.2.2) and are
// decoded with a 10-bit direct lookup table; longer codes fall back to a
// canonical walk over the per-length counts. Block types: 0 stored, 1 fixed
// Huffman, 2 dynamic Huffman; 3 is an error.
#include "core/parallel.hpp"
#include "codec/inflate.hpp"

#include <array>
#include <atomic>
#include <cstring>

#include "cadcore/error.hpp"

namespace cadcore::codec {
namespace {

constexpr int kMaxBits = 15;
constexpr int kFastBits = 12;

class BitReader {
public:
    explicit BitReader(ByteView v) : p_(v.data), end_(v.data + v.size), base_(v.data) {}
    // Starts at bit `bit` of v (bits are read LSB first, as DEFLATE stores them).
    BitReader(ByteView v, std::size_t bit) : p_(v.data + bit / 8), end_(v.data + v.size), base_(v.data) {
        if (bit / 8 > v.size) fail(ErrorCode::Truncated, "bit offset past the end");
        if (bit % 8) bits(int(bit % 8));
    }
    // Bits consumed from the start of the data.
    std::size_t bit_pos() const { return std::size_t(p_ - base_) * 8 - std::size_t(cnt_); }

    std::uint32_t bits(int n) {  // n <= 24
        fill(n);
        const std::uint32_t v = std::uint32_t(buf_ & ((std::uint64_t(1) << n) - 1));
        buf_ >>= n;
        cnt_ -= n;
        return v;
    }
    // Peek up to n bits without consuming; missing bits past the end read as zero.
    std::uint32_t peek(int n) {
        fill_soft(n);
        return std::uint32_t(buf_ & ((std::uint64_t(1) << n) - 1));
    }
    void consume(int n) {
        if (n > cnt_) fail(ErrorCode::Truncated, "DEFLATE stream ends inside a code", offset());
        buf_ >>= n;
        cnt_ -= n;
    }
    void align_to_byte() { consume(cnt_ % 8); }
    // Bytes of a stored block (after alignment).
    const std::uint8_t* take_bytes(std::size_t n) {
        // Return buffered whole bytes to the stream first.
        while (cnt_ >= 8) { --p_; cnt_ -= 8; }
        buf_ = 0;
        cnt_ = 0;
        if (std::size_t(end_ - p_) < n) fail(ErrorCode::Truncated, "stored block runs past the data", offset());
        const std::uint8_t* r = p_;
        p_ += n;
        return r;
    }
    std::size_t offset() const { return std::size_t(p_ - start_hint_); }

    // Fast mode (at least 8 input bytes left): refill to 56+ bits, then take bits unchecked.
    bool fast_ok() const { return end_ - p_ >= 8; }
    void refill_fast() {
        std::uint64_t w;
        std::memcpy(&w, p_, 8);  // little-endian hosts only (x86, wasm)
        buf_ |= w << cnt_;
        p_ += (63 - cnt_) >> 3;
        cnt_ |= 56;
    }
    std::uint32_t peek_fast(int n) const { return std::uint32_t(buf_ & ((std::uint64_t(1) << n) - 1)); }
    void drop(int n) { buf_ >>= n; cnt_ -= n; }
    std::uint32_t take_fast(int n) {
        const std::uint32_t v = peek_fast(n);
        drop(n);
        return v;
    }
    void set_start(const std::uint8_t* s) { start_hint_ = s; }

private:
    void fill(int n) {
        fill_soft(n);
        if (cnt_ < n) fail(ErrorCode::Truncated, "DEFLATE stream ends early", offset());
    }
    void fill_soft(int n) {
        if (cnt_ >= n) return;
        if (end_ - p_ >= 8) {  // whole 8-byte load: top up to 56..63 buffered bits
            std::uint64_t w;
            std::memcpy(&w, p_, 8);  // little-endian hosts only (x86, wasm)
            buf_ |= w << cnt_;
            p_ += (63 - cnt_) >> 3;
            cnt_ |= 56;
            return;
        }
        while (cnt_ < n && p_ < end_) {
            buf_ |= std::uint64_t(*p_++) << cnt_;
            cnt_ += 8;
        }
    }
    const std::uint8_t* p_;
    const std::uint8_t* end_;
    const std::uint8_t* start_hint_ = nullptr;
    const std::uint8_t* base_;
    std::uint64_t buf_ = 0;
    int cnt_ = 0;
};

// Canonical Huffman decoder for one alphabet.
class Huffman {
public:
    void build(const std::uint8_t* lengths, int n) {
        std::array<std::uint16_t, kMaxBits + 1> count{};
        for (int i = 0; i < n; ++i) ++count[lengths[i]];
        count[0] = 0;
        // Over-subscribed codes are invalid; incomplete codes are allowed (RFC 1951 permits
        // a single distance code, and decoders must reject only codes that are actually used).
        int left = 1;
        for (int len = 1; len <= kMaxBits; ++len) {
            left = (left << 1) - count[len];
            if (left < 0) fail(ErrorCode::Malformed, "over-subscribed Huffman code");
        }
        std::array<std::uint16_t, kMaxBits + 2> offs{};
        for (int len = 1; len <= kMaxBits; ++len) offs[len + 1] = std::uint16_t(offs[len] + count[len]);
        symbols_.assign(std::size_t(n), 0);
        for (int i = 0; i < n; ++i)
            if (lengths[i]) symbols_[offs[lengths[i]]++] = std::uint16_t(i);
        count_ = count;
        // Fast table: for each kFastBits-bit pattern, the symbol and length of a code of length <= kFastBits.
        fast_.fill(0);
        int code = 0, first = 0, index = 0;
        for (int len = 1; len <= kFastBits; ++len) {
            for (int k = 0; k < count[len]; ++k) {
                const std::uint16_t sym = symbols_[std::size_t(index + k)];
                const int c = code + k;
                // Codes are stored MSB-first in the stream; reverse into LSB-first table index.
                int rev = 0;
                for (int b = 0; b < len; ++b) rev |= ((c >> b) & 1) << (len - 1 - b);
                for (int fill = rev; fill < (1 << kFastBits); fill += 1 << len) fast_[std::size_t(fill)] = std::uint16_t(sym << 4 | len);
            }
            index += count[len];
            first += count[len];
            code = (code + count[len]) << 1;
        }
        first_slow_index_ = index;
        first_slow_code_ = code;
    }

    // After refill_fast: a code of at most kFastBits from the table, else the checked path.
    int decode_fast(BitReader& br) const {
        const std::uint16_t e = fast_[br.peek_fast(kFastBits)];
        if (e) {
            br.drop(e & 15);
            return e >> 4;
        }
        return decode(br);
    }

    int decode(BitReader& br) const {
        const std::uint32_t look = br.peek(kFastBits);
        const std::uint16_t e = fast_[look];
        if (e) {
            br.consume(e & 15);
            return e >> 4;
        }
        // Slow path: canonical walk for lengths > kFastBits (bit by bit, MSB-first code order).
        int code = 0, first = 0, index = 0;
        for (int len = 1; len <= kMaxBits; ++len) {
            code |= int(br.bits(1));
            const int count = count_[std::size_t(len)];
            if (code - count < first) return symbols_[std::size_t(index + (code - first))];
            index += count;
            first = (first + count) << 1;
            code <<= 1;
        }
        fail(ErrorCode::Malformed, "invalid Huffman code", br.offset());
    }

private:
    std::array<std::uint16_t, 1 << kFastBits> fast_{};
    std::array<std::uint16_t, kMaxBits + 1> count_{};
    std::vector<std::uint16_t> symbols_;
    int first_slow_index_ = 0, first_slow_code_ = 0;
};

constexpr std::uint16_t kLenBase[29] = {3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258};
constexpr std::uint8_t kLenExtra[29] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
constexpr std::uint16_t kDistBase[30] = {1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577};
constexpr std::uint8_t kDistExtra[30] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};
constexpr std::uint8_t kCodeLengthOrder[19] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};

// Code tables of a dynamic block (after its 3 header bits). strict (speculative search): every
// code must also be complete, as encoders write them (RFC 1951 leaves incomplete codes legal).
bool complete(const std::uint8_t* len, int n) {
    long long left = 1;
    int used = 0;
    std::array<int, kMaxBits + 1> count{};
    for (int i = 0; i < n; ++i) { ++count[len[i]]; used += len[i] != 0; }
    for (int l = 1; l <= kMaxBits; ++l) left = (left << 1) - count[std::size_t(l)];
    return left == 0 || used == 1;
}

bool read_dynamic(BitReader& br, Huffman& lit, Huffman& dist, bool strict) {
    const int hlit = int(br.bits(5)) + 257, hdist = int(br.bits(5)) + 1, hclen = int(br.bits(4)) + 4;
    if (strict && (hlit > 286 || hdist > 30)) fail(ErrorCode::Malformed, "code counts out of range", br.offset());
    std::uint8_t cl[19] = {};
    for (int i = 0; i < hclen; ++i) cl[kCodeLengthOrder[i]] = std::uint8_t(br.bits(3));
    if (strict && !complete(cl, 19)) fail(ErrorCode::Malformed, "incomplete code-length code", br.offset());
    Huffman clh;
    clh.build(cl, 19);
    std::uint8_t lens[320] = {};
    int n = 0;
    while (n < hlit + hdist) {
        const int sym = clh.decode(br);
        if (sym < 16) {
            lens[n++] = std::uint8_t(sym);
        } else {
            int rep = 0;
            std::uint8_t val = 0;
            if (sym == 16) {
                if (n == 0) fail(ErrorCode::Malformed, "repeat code with no previous length", br.offset());
                val = lens[n - 1];
                rep = 3 + int(br.bits(2));
            } else if (sym == 17) {
                rep = 3 + int(br.bits(3));
            } else {
                rep = 11 + int(br.bits(7));
            }
            if (n + rep > hlit + hdist) fail(ErrorCode::Malformed, "code lengths overflow", br.offset());
            while (rep--) lens[n++] = val;
        }
    }
    if (lens[256] == 0) fail(ErrorCode::Malformed, "dynamic block without end-of-block code", br.offset());
    if (strict && (!complete(lens, hlit) || !complete(lens + hlit, hdist))) fail(ErrorCode::Malformed, "incomplete code", br.offset());
    lit.build(lens, hlit);
    dist.build(lens + hlit, hdist);
    return true;
}

// Output buffer: grows geometrically, never beyond max_output (plus kSlack bytes of scratch
// after the end, so the fast loop may copy matches 8 bytes at a time).
constexpr std::size_t kSlack = 16;
struct Output {
    Buffer buf;
    std::size_t n = 0, max;
    std::size_t capacity() const { return buf.size() - kSlack; }
    std::uint8_t* room(std::size_t k) {
        if (n + k > max) fail(ErrorCode::Malformed, "DEFLATE output exceeds the declared size");
        if (n + k > capacity()) buf.resize(std::min(max, std::max(n + k, capacity() * 2 + 4096)) + kSlack);
        std::uint8_t* p = buf.data() + n;
        n += k;
        return p;
    }
};

void copy_match(std::uint8_t* dst, std::size_t d, std::size_t len) {
    const std::uint8_t* src = dst - d;
    if (d >= 8) {  // 8-byte steps never read bytes this copy has not written yet
        for (std::size_t i = 0; i < len; i += 8) std::memcpy(dst + i, src + i, 8);  // may write up to 7 bytes of slack
    } else {
        for (std::size_t i = 0; i < len; ++i) dst[i] = src[i];  // overlapping: repeats the last d bytes
    }
}

void inflate_block(BitReader& br, const Huffman& lit, const Huffman& dist, Output& out) {
    while (true) {
        // Fast loop: one refill covers a literal/length code (15), its extra bits (5), a distance
        // code (15) and its extra bits (13): 48 <= 56 bits. Room for the longest match (258).
        while (br.fast_ok() && out.n + 258 <= out.capacity() && out.n + 258 <= out.max) {
            br.refill_fast();
            const int sym = lit.decode_fast(br);
            if (sym < 256) {
                out.buf[out.n++] = std::uint8_t(sym);
                continue;
            }
            if (sym == 256) return;
            const int li = sym - 257;
            if (li >= 29) fail(ErrorCode::Malformed, "invalid length symbol", br.offset());
            const std::size_t len = kLenBase[li] + br.take_fast(kLenExtra[li]);
            const int di = dist.decode_fast(br);
            if (di >= 30) fail(ErrorCode::Malformed, "invalid distance symbol", br.offset());
            const std::size_t d = kDistBase[di] + br.take_fast(kDistExtra[di]);
            if (d > out.n) fail(ErrorCode::Malformed, "distance reaches before the start of the output", br.offset());
            copy_match(out.buf.data() + out.n, d, len);
            out.n += len;
        }
        const int sym = lit.decode(br);
        if (sym < 256) {
            if (out.n < out.capacity()) out.buf[out.n++] = std::uint8_t(sym);
            else *out.room(1) = std::uint8_t(sym);
            continue;
        }
        if (sym == 256) return;
        const int li = sym - 257;
        if (li >= 29) fail(ErrorCode::Malformed, "invalid length symbol", br.offset());
        const std::size_t len = kLenBase[li] + br.bits(kLenExtra[li]);
        const int di = dist.decode(br);
        if (di >= 30) fail(ErrorCode::Malformed, "invalid distance symbol", br.offset());
        const std::size_t d = kDistBase[di] + br.bits(kDistExtra[di]);
        if (d > out.n) fail(ErrorCode::Malformed, "distance reaches before the start of the output", br.offset());
        if (out.n + len <= out.capacity() && out.n + len <= out.max) {
            copy_match(out.buf.data() + out.n, d, len);
            out.n += len;
        } else {
            std::uint8_t* dst = out.room(len);
            const std::uint8_t* src = dst - d;
            for (std::size_t i = 0; i < len; ++i) dst[i] = src[i];
        }
    }
}

}  // namespace

namespace {
const Huffman& fixed_lit_table() {
    static const Huffman h = [] {
        std::uint8_t l[288];
        for (int i = 0; i < 288; ++i) l[i] = std::uint8_t(i < 144 ? 8 : i < 256 ? 9 : i < 280 ? 7 : 8);
        Huffman t;
        t.build(l, 288);
        return t;
    }();
    return h;
}
const Huffman& fixed_dist_table() {
    static const Huffman h = [] {
        std::uint8_t l[30];
        std::memset(l, 5, sizeof l);
        Huffman t;
        t.build(l, 30);
        return t;
    }();
    return h;
}
}  // namespace

Buffer inflate(ByteView raw, std::size_t expected_size, std::size_t max_output) {
    BitReader br(raw);
    br.set_start(raw.data);
    Output out;
    out.max = max_output;
    out.buf.resize(std::min(max_output, expected_size ? expected_size : raw.size * 3) + kSlack);
    const Huffman& fixed_lit = fixed_lit_table();
    const Huffman& fixed_dist = fixed_dist_table();
    bool last = false;
    while (!last) {
        last = br.bits(1) != 0;
        const std::uint32_t type = br.bits(2);
        if (type == 0) {
            br.align_to_byte();
            const std::uint8_t* hdr = br.take_bytes(4);
            const std::uint16_t len = std::uint16_t(hdr[0] | hdr[1] << 8), nlen = std::uint16_t(hdr[2] | hdr[3] << 8);
            if (len != std::uint16_t(~nlen)) fail(ErrorCode::Malformed, "stored block length check failed", br.offset());
            const std::uint8_t* data = br.take_bytes(len);
            if (len) std::memcpy(out.room(len), data, len);
        } else if (type == 1) {
            inflate_block(br, fixed_lit, fixed_dist, out);
        } else if (type == 2) {
            Huffman lit, dist;
            read_dynamic(br, lit, dist, false);
            inflate_block(br, lit, dist, out);
        } else {
            fail(ErrorCode::Malformed, "invalid DEFLATE block type 3", br.offset());
        }
    }
    out.buf.resize(out.n);
    return std::move(out.buf);
}

// CRC-32 (IEEE 802.3, reflected), slicing-by-8: eight bytes per step through eight tables, where
// table k maps a byte to its CRC contribution k bytes further on.
std::uint32_t crc32(const std::uint8_t* data, std::size_t size, std::uint32_t crc) {
    static const auto T = [] {
        std::array<std::array<std::uint32_t, 256>, 8> t{};
        for (std::uint32_t i = 0; i < 256; ++i) {
            std::uint32_t c = i;
            for (int k = 0; k < 8; ++k) c = c & 1 ? 0xEDB88320u ^ (c >> 1) : c >> 1;
            t[0][i] = c;
        }
        for (std::size_t k = 1; k < 8; ++k)
            for (std::size_t i = 0; i < 256; ++i) t[k][i] = t[0][t[k - 1][i] & 0xFF] ^ (t[k - 1][i] >> 8);
        return t;
    }();
    crc = ~crc;
    std::size_t i = 0;
    for (; i + 8 <= size; i += 8) {
        const std::uint32_t lo = crc ^ (std::uint32_t(data[i]) | std::uint32_t(data[i + 1]) << 8 | std::uint32_t(data[i + 2]) << 16 | std::uint32_t(data[i + 3]) << 24);
        crc = T[7][lo & 0xFF] ^ T[6][lo >> 8 & 0xFF] ^ T[5][lo >> 16 & 0xFF] ^ T[4][lo >> 24] ^ T[3][data[i + 4]] ^ T[2][data[i + 5]] ^ T[1][data[i + 6]] ^
              T[0][data[i + 7]];
    }
    for (; i < size; ++i) crc = T[0][(crc ^ data[i]) & 0xFF] ^ (crc >> 8);
    return ~crc;
}

// ---- speculative parallel decoding ---------------------------------------------------------------
namespace {

// Part output: values 0..255 are bytes; 256 + i is byte i of the 32 KB window before the part
// (not known while the part decodes in parallel).
constexpr std::uint16_t kWindowRef = 256;
constexpr std::size_t kWindow = 32768;
constexpr std::size_t kNoStop = std::size_t(-1);

// Speculative output: values written through a pointer into a buffer grown ahead (no per-value
// capacity check in the fast loop).
struct SpecOut {
    std::vector<std::uint16_t, NoInitAllocator<std::uint16_t>> buf;
    std::size_t n = 0;
    std::size_t size() const { return n; }
    void ensure(std::size_t k) {
        if (n + k > buf.size()) buf.resize(std::max(n + k, buf.size() * 2 + 65536));
    }
    void push_back(std::uint16_t v) { ensure(1); buf[n++] = v; }
    void insert_bytes(const std::uint8_t* d, std::size_t len) { ensure(len); for (std::size_t i = 0; i < len; ++i) buf[n + i] = d[i]; n += len; }
    std::uint16_t operator[](std::size_t i) const { return buf[i]; }
};

// Copies len values from d back; positions before the part become window references.
void spec_copy(SpecOut& out, std::size_t d, std::size_t len) {
    const std::size_t n = out.n;
    out.ensure(len);
    out.n += len;
    std::uint16_t* o = out.buf.data() + n;
    if (d <= n) {  // inside the part
        const std::uint16_t* src = o - d;
        if (d >= len) std::memcpy(o, src, len * 2);
        else for (std::size_t k = 0; k < len; ++k) o[k] = src[k];
        return;
    }
    for (std::size_t k = 0; k < len; ++k) {
        const std::ptrdiff_t src = std::ptrdiff_t(n + k) - std::ptrdiff_t(d);
        o[k] = src >= 0 ? out.buf[std::size_t(src)] : std::uint16_t(kWindowRef + std::size_t(std::ptrdiff_t(kWindow) + src));
    }
}

void spec_block(BitReader& br, const Huffman& lit, const Huffman& dist, SpecOut& out, bool first, std::size_t max) {
    while (true) {
        // fast loop, as inflate_block: one refill per symbol
        while (br.fast_ok()) {
            br.refill_fast();
            const int sym = lit.decode_fast(br);
            if (sym < 256) {
                if (out.n < out.buf.size()) out.buf[out.n++] = std::uint16_t(sym);
                else out.push_back(std::uint16_t(sym));
                continue;
            }
            if (sym == 256) return;
            const int li = sym - 257;
            if (li >= 29) fail(ErrorCode::Malformed, "invalid length symbol");
            const std::size_t len = kLenBase[li] + br.take_fast(kLenExtra[li]);
            const int di = dist.decode_fast(br);
            if (di >= 30) fail(ErrorCode::Malformed, "invalid distance symbol");
            const std::size_t d = kDistBase[di] + br.take_fast(kDistExtra[di]);
            if (d > out.size() && (first || d > out.size() + kWindow)) fail(ErrorCode::Malformed, "distance before the window");
            if (out.size() + len > max) fail(ErrorCode::Malformed, "output too large");
            spec_copy(out, d, len);
        }
        const int sym = lit.decode(br);
        if (sym < 256) {
            out.push_back(std::uint16_t(sym));
            continue;
        }
        if (sym == 256) return;
        const int li = sym - 257;
        if (li >= 29) fail(ErrorCode::Malformed, "invalid length symbol");
        const std::size_t len = kLenBase[li] + br.bits(kLenExtra[li]);
        const int di = dist.decode(br);
        if (di >= 30) fail(ErrorCode::Malformed, "invalid distance symbol");
        const std::size_t d = kDistBase[di] + br.bits(kDistExtra[di]);
        if (d > out.size() && (first || d > out.size() + kWindow)) fail(ErrorCode::Malformed, "distance before the window");
        if (out.size() + len > max) fail(ErrorCode::Malformed, "output too large");
        spec_copy(out, d, len);
    }
}

const Huffman& fixed_lit_table();
const Huffman& fixed_dist_table();

struct Part {
    std::size_t start_bit = 0, stop_bit = kNoStop;
    SpecOut out;
    bool ok = false, final = false;
};

// Decodes blocks from start_bit; succeeds when a block ends exactly at stop_bit (or, for the last
// part, at the final block). Any error or overshoot leaves ok false.
void decode_part(ByteView raw, Part& p, bool first, std::size_t max) {
    try {
        BitReader br(raw, p.start_bit);
        br.set_start(raw.data);
        p.out.buf.resize(std::min<std::size_t>(max, (raw.size - p.start_bit / 8) * 4));
        while (true) {
            if (p.stop_bit != kNoStop) {
                const std::size_t b = br.bit_pos();
                if (b == p.stop_bit) { p.ok = true; return; }
                if (b > p.stop_bit) return;
            }
            const bool last = br.bits(1) != 0;
            const std::uint32_t type = br.bits(2);
            if (type == 0) {
                br.align_to_byte();
                const std::uint8_t* hdr = br.take_bytes(4);
                const std::uint16_t len = std::uint16_t(hdr[0] | hdr[1] << 8), nlen = std::uint16_t(hdr[2] | hdr[3] << 8);
                if (len != std::uint16_t(~nlen)) return;
                const std::uint8_t* d = br.take_bytes(len);
                if (p.out.size() + len > max) return;
                p.out.insert_bytes(d, len);
            } else if (type == 1) {
                spec_block(br, fixed_lit_table(), fixed_dist_table(), p.out, first, max);
            } else if (type == 2) {
                Huffman lit, dist;
                read_dynamic(br, lit, dist, false);
                spec_block(br, lit, dist, p.out, first, max);
            } else {
                return;
            }
            if (last) {
                p.final = true;
                p.ok = p.stop_bit == kNoStop;
                return;
            }
        }
    } catch (const DecodeError&) {
    }
}

// First bit in [from, to) where a non-final dynamic block plausibly starts: complete codes, and the
// whole block decodes (with an unknown window) to its end code followed by a valid block type.
std::size_t find_block_start(ByteView raw, std::size_t from, std::size_t to) {
    auto bit = [&](std::size_t b) { return (raw.data[b >> 3] >> (b & 7)) & 1; };
    to = std::min(to, raw.size * 8 - 64);
    SpecOut scratch;
    for (std::size_t b = from; b < to; ++b) {
        if (bit(b) != 0 || bit(b + 1) != 0 || bit(b + 2) != 1) continue;  // BFINAL 0, BTYPE 10
        // HLIT <= 29, HDIST <= 29
        std::uint32_t h = 0;
        for (int k = 0; k < 10; ++k) h |= std::uint32_t(bit(b + 3 + std::size_t(k))) << k;
        if ((h & 31) > 29 || (h >> 5) > 29) continue;
        {  // the code-length code must be complete: read its lengths straight from the bits
            std::uint32_t hclen = 0;
            for (int k = 0; k < 4; ++k) hclen |= std::uint32_t(bit(b + 13 + std::size_t(k))) << k;
            hclen += 4;
            int count[8] = {};
            for (std::uint32_t i = 0; i < hclen; ++i) {
                int l = 0;
                for (int k = 0; k < 3; ++k) l |= bit(b + 17 + 3 * i + std::size_t(k)) << k;
                ++count[l];
            }
            long long left = 1;
            int used = 0;
            for (int l = 1; l < 8; ++l) { left = (left << 1) - count[l]; used += count[l]; }
            if (!(left == 0 || used == 1)) continue;
        }
        try {
            BitReader br(raw, b + 3);
            Huffman lit, dist;
            read_dynamic(br, lit, dist, true);
            scratch.n = 0;
            spec_block(br, lit, dist, scratch, false, std::size_t(1) << 26);
            if (br.bits(3) >> 1 == 3) continue;  // the next block type must be valid
            return b;
        } catch (const DecodeError&) {
        }
    }
    return kNoStop;
}

}  // namespace

Buffer inflate_parallel(ByteView raw, std::size_t expected_size, std::size_t max_output, unsigned threads) {
    constexpr std::size_t kMinPart = 512 << 10;  // compressed bytes per part
    const std::size_t n = std::min<std::size_t>(thread_count(threads), raw.size / kMinPart);
    if (n < 2) return inflate(raw, expected_size, max_output);
    // 1. a block start in each part after the first
    std::vector<Slot<std::size_t>> found(n);
    parallel_for(n - 1, 1, threads, [&](std::size_t b0, std::size_t e0) {
        for (std::size_t k = b0 + 1; k < e0 + 1; ++k) {
            const std::size_t from = raw.size * 8 * k / n;
            found[k].value = find_block_start(raw, from, from + std::min<std::size_t>(raw.size * 8 / n, std::size_t(8) << 20));
        }
    });
    std::vector<std::size_t> starts{0};
    for (std::size_t k = 1; k < n; ++k)
        if (found[k].value != kNoStop && found[k].value > starts.back()) starts.push_back(found[k].value);
    if (starts.size() < 2) return inflate(raw, expected_size, max_output);
    // 2. every part decodes up to the next start, which must be where one of its blocks ends
    std::vector<Slot<Part>> parts(starts.size());
    std::vector<std::function<void()>> tasks;
    for (std::size_t k = 0; k < starts.size(); ++k)
        tasks.push_back([&, k] {
            Part& p = parts[k].value;
            p.start_bit = starts[k];
            p.stop_bit = k + 1 < starts.size() ? starts[k + 1] : kNoStop;
            decode_part(raw, p, k == 0, max_output);
        });
    run_all(tasks, threads);
    std::size_t total = 0;
    for (auto& p : parts) {
        if (!p.value.ok) return inflate(raw, expected_size, max_output);  // a start was not a real boundary
        total += p.value.out.size();
    }
    if (total > max_output) return inflate(raw, expected_size, max_output);
    // 3. window references: each part's last 32 KB in order (it needs only what precedes it), then
    //    the rest of every part in parallel
    Buffer out(total);  // uninitialised: every byte is written below
    std::vector<std::size_t> at(parts.size() + 1, 0);
    for (std::size_t k = 0; k < parts.size(); ++k) at[k + 1] = at[k] + parts[k].value.out.size();
    std::atomic<bool> bad{false};  // set by any part: a reference before the stream start
    auto resolve = [&](std::size_t k, std::size_t i0, std::size_t i1) {
        const auto& v = parts[k].value.out;
        const std::size_t o = at[k];
        for (std::size_t i = i0; i < i1; ++i) {
            const std::uint16_t x = v[i];
            if (x < kWindowRef) { out[o + i] = std::uint8_t(x); continue; }
            const std::size_t w = x - kWindowRef;  // window byte w = out[o - kWindow + w]
            if (o + w < kWindow) { bad.store(true); return; }  // reaches before the stream
            out[o + i] = out[o - kWindow + w];
        }
    };
    for (std::size_t k = 0; k < parts.size(); ++k) {
        const std::size_t sz = parts[k].value.out.size();
        resolve(k, sz > kWindow ? sz - kWindow : 0, sz);
        if (bad) return inflate(raw, expected_size, max_output);
    }
    std::vector<std::function<void()>> heads;
    for (std::size_t k = 0; k < parts.size(); ++k) {
        const std::size_t sz = parts[k].value.out.size();
        if (sz > kWindow) heads.push_back([&, k, sz] { resolve(k, 0, sz - kWindow); });
    }
    run_all(heads, threads);
    if (bad) return inflate(raw, expected_size, max_output);
    return out;
}

namespace {
std::uint32_t gf2_times(const std::uint32_t* mat, std::uint32_t vec) {
    std::uint32_t sum = 0;
    for (int i = 0; vec; vec >>= 1, ++i)
        if (vec & 1) sum ^= mat[i];
    return sum;
}
void gf2_square(std::uint32_t* sq, const std::uint32_t* mat) {
    for (int n = 0; n < 32; ++n) sq[n] = gf2_times(mat, mat[n]);
}
}  // namespace

std::uint32_t crc32_combine(std::uint32_t crc_a, std::uint32_t crc_b, std::size_t len_b) {
    if (len_b == 0) return crc_a;
    // odd: the operator for one zero bit; squaring it gives the operators for 2, 4, 8 ... bits
    std::uint32_t even[32], odd[32];
    odd[0] = 0xEDB88320u;
    for (int n = 1; n < 32; ++n) odd[n] = 1u << (n - 1);
    gf2_square(even, odd);  // 2 bits
    gf2_square(odd, even);  // 4 bits
    do {  // apply len_b zero bytes to crc_a, one bit of len_b per squaring
        gf2_square(even, odd);
        if (len_b & 1) crc_a = gf2_times(even, crc_a);
        len_b >>= 1;
        if (!len_b) break;
        gf2_square(odd, even);
        if (len_b & 1) crc_a = gf2_times(odd, crc_a);
        len_b >>= 1;
    } while (len_b);
    return crc_a ^ crc_b;
}

std::uint32_t crc32_parallel(const std::uint8_t* data, std::size_t size, unsigned threads) {
    constexpr std::size_t kMinPiece = 4 << 20;
    const std::size_t n = std::max<std::size_t>(1, std::min<std::size_t>(thread_count(threads), size / kMinPiece));
    if (n == 1) return crc32(data, size);
    std::vector<Slot<std::uint32_t>> part(n);
    parallel_for(n, 1, threads, [&](std::size_t b, std::size_t e) {
        for (std::size_t k = b; k < e; ++k) part[k].value = crc32(data + size * k / n, size * (k + 1) / n - size * k / n);
    });
    std::uint32_t crc = part[0].value;
    for (std::size_t k = 1; k < n; ++k) crc = crc32_combine(crc, part[k].value, size * (k + 1) / n - size * k / n);
    return crc;
}

}  // namespace cadcore::codec
