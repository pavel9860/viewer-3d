// Bounds-checked binary reading over an immutable byte view.
#pragma once

#include <cstdint>
#include <cstring>
#include <string_view>

#include "cadcore/error.hpp"

namespace cadcore {

// First occurrence of needle[0, m) in hay[0, n), or n when absent. Scans with memchr (vectorised by
// the C library) for the needle's least common byte (binary data is dominated by 00, ff and small
// values), then compares the whole needle.
std::size_t find_bytes(const std::uint8_t* hay, std::size_t n, const std::uint8_t* needle, std::size_t m);

struct ByteView {
    const std::uint8_t* data = nullptr;
    std::size_t size = 0;

    ByteView() = default;
    ByteView(const void* d, std::size_t n) : data(static_cast<const std::uint8_t*>(d)), size(n) {}

    std::uint8_t operator[](std::size_t i) const { return data[i]; }
    bool empty() const { return size == 0; }
    std::string_view text() const { return {reinterpret_cast<const char*>(data), size}; }
    ByteView sub(std::size_t pos, std::size_t len) const;
    bool starts_with(std::string_view s, std::size_t pos = 0) const {
        return pos + s.size() <= size && std::memcmp(data + pos, s.data(), s.size()) == 0;
    }
    // First occurrence of `needle` in [from, to); returns npos when absent.
    std::size_t find(std::string_view needle, std::size_t from = 0, std::size_t to = npos) const;

    static constexpr std::size_t npos = static_cast<std::size_t>(-1);
};

// Sequential little/big-endian reader. Every read checks bounds and throws
// DecodeError(Truncated) with the failing offset.
class Reader {
public:
    explicit Reader(ByteView v, std::size_t pos = 0) : v_(v), pos_(pos) {}

    std::size_t pos() const { return pos_; }
    std::size_t remaining() const { return pos_ <= v_.size ? v_.size - pos_ : 0; }
    void seek(std::size_t p) {
        if (p > v_.size) fail(ErrorCode::Truncated, "seek past end", p);
        pos_ = p;
    }
    void skip(std::size_t n) { need(n); pos_ += n; }
    void need(std::size_t n) const {
        if (n > remaining()) fail(ErrorCode::Truncated, "unexpected end of data", pos_);
    }

    std::uint8_t u8() { need(1); return v_.data[pos_++]; }
    std::uint16_t u16le() { return load_le<std::uint16_t>(); }
    std::uint32_t u32le() { return load_le<std::uint32_t>(); }
    std::uint64_t u64le() { return load_le<std::uint64_t>(); }
    std::uint16_t u16be() { return load_be<std::uint16_t>(); }
    std::uint32_t u32be() { return load_be<std::uint32_t>(); }
    float f32le() { std::uint32_t u = u32le(); float f; std::memcpy(&f, &u, 4); return f; }
    double f64le() { std::uint64_t u = u64le(); double d; std::memcpy(&d, &u, 8); return d; }
    float f32be() { std::uint32_t u = u32be(); float f; std::memcpy(&f, &u, 4); return f; }
    double f64be() { std::uint64_t u = load_be<std::uint64_t>(); double d; std::memcpy(&d, &u, 8); return d; }
    const std::uint8_t* bytes(std::size_t n) { need(n); const std::uint8_t* p = v_.data + pos_; pos_ += n; return p; }

private:
    template <class T> T load_le() {
        need(sizeof(T));
        T v = 0;
        for (std::size_t i = 0; i < sizeof(T); ++i) v |= T(v_.data[pos_ + i]) << (8 * i);
        pos_ += sizeof(T);
        return v;
    }
    template <class T> T load_be() {
        need(sizeof(T));
        T v = 0;
        for (std::size_t i = 0; i < sizeof(T); ++i) v = T(v << 8) | T(v_.data[pos_ + i]);
        pos_ += sizeof(T);
        return v;
    }

    ByteView v_;
    std::size_t pos_;
};

}  // namespace cadcore
