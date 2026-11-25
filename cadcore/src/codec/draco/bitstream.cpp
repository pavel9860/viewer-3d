#include "codec/draco/bitstream.hpp"

#include <algorithm>

namespace cadcore::codec::draco {

void bad(const std::string& what) { fail(ErrorCode::Malformed, "Draco: " + what); }
void unsupported(const std::string& what) { fail(ErrorCode::Unsupported, "Draco: " + what); }

// ---- Stream ------------------------------------------------------------------------------------------

std::uint64_t Stream::varint64() {
    std::uint64_t r = 0;
    for (int shift = 0;; shift += 7) {
        if (shift > 63) bad("varint longer than 64 bits");
        const std::uint8_t in = u8();
        r |= std::uint64_t(in & 0x7f) << shift;
        if (!(in & 0x80)) return r;
    }
}

std::uint32_t Stream::varint32() {
    const std::uint64_t v = varint64();
    if (v > 0xffffffffu) bad("varint does not fit 32 bits");
    return std::uint32_t(v);
}

std::uint32_t Stream::bits(int n) {
    if (!in_bits_) bad("bit read outside a bit run");
    std::uint32_t v = 0;
    for (int i = 0; i < n; ++i) {
        const std::uint64_t byte = bit_ >> 3;
        // Reads past the end yield zero bits, as in the reference decoder; end_bits() then
        // reports the overrun.
        const std::uint32_t b = p_ + byte < n_ ? (d_[p_ + byte] >> (bit_ & 7)) & 1u : 0u;
        v |= b << i;
        ++bit_;
    }
    return v;
}

void Stream::end_bits() {
    in_bits_ = false;
    skip(std::size_t((bit_ + 7) / 8));
}

// ---- BitSpan ------------------------------------------------------------------------------------------

std::uint32_t BitSpan::bits(int k) {
    std::uint32_t v = 0;
    for (int i = 0; i < k; ++i) {
        if (pos_ >= std::uint64_t(n_) * 8) bad("EdgeBreaker symbol buffer exhausted");
        v |= std::uint32_t((d_[pos_ >> 3] >> (pos_ & 7)) & 1) << i;
        ++pos_;
    }
    return v;
}

// ---- rANS core -----------------------------------------------------------------------------------------

namespace {

constexpr std::uint32_t kIoBase = 256;
constexpr std::uint32_t kAnsLBase = 4096;  // rABS
constexpr std::uint32_t kP8Precision = 256;

std::uint32_t le(const std::uint8_t* p, int n) {
    std::uint32_t v = 0;
    for (int i = n - 1; i >= 0; --i) v = v << 8 | p[i];
    return v;
}

// Reads the initial state from the end of `buf[0, size)`; `max_prefix` is 3 for rABS, 4 for rANS.
void read_init(const std::uint8_t* buf, std::size_t size, std::uint32_t l_base, int max_prefix, std::int64_t& off, std::uint32_t& state) {
    if (size < 1) bad("empty entropy-coded buffer");
    const int x = buf[size - 1] >> 6;
    if (x >= max_prefix) bad("invalid entropy-coded buffer prefix");
    const std::size_t n = std::size_t(x) + 1;
    if (size < n) bad("entropy-coded buffer too short");
    static constexpr std::uint32_t kMask[4] = {0x3F, 0x3FFF, 0x3FFFFF, 0x3FFFFFFF};
    off = std::int64_t(size - n);
    state = (le(buf + size - n, int(n)) & kMask[x]) + l_base;
    if (std::uint64_t(state) >= std::uint64_t(l_base) * kIoBase) bad("invalid entropy-coder state");
}

struct SymbolCoder {
    std::uint32_t precision = 0, l_base = 0;
    std::vector<std::uint32_t> prob, cum, lut;
    const std::uint8_t* buf = nullptr;
    std::int64_t off = 0;
    std::uint32_t state = 0;

    explicit SymbolCoder(int unique_symbols_bit_length) {
        const int unclamped = (3 * unique_symbols_bit_length) / 2;
        const int bits = std::clamp(unclamped, 12, 20);
        precision = 1u << bits;
        l_base = precision * 4;
    }

    void create(Stream& s) {
        const std::uint32_t num_symbols = s.ver < version(2, 0) ? s.u32() : s.varint32();
        if (num_symbols / 64 > s.remaining()) bad("rANS symbol count exceeds the data");
        prob.assign(num_symbols, 0);
        for (std::uint32_t i = 0; i < num_symbols; ++i) {
            const std::uint8_t prob_data = s.u8();
            const int token = prob_data & 3;
            if (token == 3) {
                const std::uint32_t offset = prob_data >> 2;
                if (i + offset >= num_symbols) bad("rANS zero run past the symbol count");
                i += offset;  // prob[] already zero
            } else {
                std::uint32_t p = prob_data >> 2;
                for (int b = 0; b < token; ++b) p |= std::uint32_t(s.u8()) << (8 * (b + 1) - 2);
                prob[i] = p;
            }
        }
        if (num_symbols == 0) return;
        cum.resize(num_symbols);
        lut.assign(precision, 0);
        std::uint64_t c = 0;
        for (std::uint32_t i = 0; i < num_symbols; ++i) {
            cum[i] = std::uint32_t(c);
            const std::uint64_t next = c + prob[i];
            if (next > precision) bad("rANS probabilities exceed the precision");
            for (std::uint64_t j = c; j < next; ++j) lut[std::size_t(j)] = i;
            c = next;
        }
        if (c != precision) bad("rANS probabilities do not sum to the precision");
    }

    void start(Stream& s) {
        const std::uint64_t n = s.ver < version(2, 0) ? s.u64() : s.varint64();
        if (n > s.remaining()) bad("rANS data exceeds the stream");
        buf = s.bytes(std::size_t(n));
        read_init(buf, std::size_t(n), l_base, 4, off, state);
    }

    std::uint32_t read() {
        while (state < l_base && off > 0) state = state * kIoBase + buf[--off];
        const std::uint32_t quo = state / precision, rem = state % precision;
        const std::uint32_t sym = lut[rem];
        state = quo * prob[sym] + rem - cum[sym];
        return sym;
    }
};

}  // namespace

// ---- BitDecoder -----------------------------------------------------------------------------------------

void BitDecoder::start(Stream& s) {
    prob_zero_ = s.u8();
    const std::uint32_t n = s.ver < version(2, 2) ? s.u32() : s.varint32();
    if (n > s.remaining()) bad("rABS data exceeds the stream");
    buf_ = s.bytes(n);
    read_init(buf_, n, kAnsLBase, 3, off_, state_);
}

bool BitDecoder::next() {
    if (!buf_) bad("rABS decoder used before it was started");
    const std::uint32_t p = kP8Precision - prob_zero_;
    if (state_ < kAnsLBase && off_ > 0) state_ = state_ * kIoBase + buf_[--off_];
    const std::uint32_t x = state_, quot = x / kP8Precision, rem = x % kP8Precision, xn = quot * p;
    const bool val = rem < p;
    state_ = val ? xn + rem : x - xn - p;
    return val;
}

// ---- DecodeSymbols ---------------------------------------------------------------------------------------

void decode_symbols(Stream& s, std::uint32_t num_values, int num_components, std::uint32_t* out) {
    if (num_values == 0) return;
    const int scheme = s.u8();
    if (scheme == 0) {  // tagged: rANS-coded bit lengths, then the value bits as one bit run
        if (num_components <= 0) bad("tagged symbols need a component count");
        SymbolCoder tags(5);
        tags.create(s);
        if (tags.prob.empty()) bad("tagged symbols without a symbol table");
        tags.start(s);
        s.start_bits();
        std::uint32_t k = 0;
        for (std::uint32_t i = 0; i < num_values; i += std::uint32_t(num_components)) {
            const std::uint32_t bit_length = tags.read();
            if (bit_length > 32) bad("tagged symbol longer than 32 bits");
            for (int j = 0; j < num_components && k < num_values; ++j) out[k++] = s.bits(int(bit_length));
        }
        s.end_bits();
    } else if (scheme == 1) {  // raw
        const int max_bit_length = s.u8();
        if (max_bit_length < 1 || max_bit_length > 18) bad("raw symbols: invalid maximum bit length");
        SymbolCoder coder(max_bit_length);
        coder.create(s);
        if (coder.prob.empty()) bad("raw symbols without a symbol table");
        coder.start(s);
        for (std::uint32_t i = 0; i < num_values; ++i) out[i] = coder.read();
    } else {
        bad("unknown symbol coding scheme " + std::to_string(scheme));
    }
}

}  // namespace cadcore::codec::draco
