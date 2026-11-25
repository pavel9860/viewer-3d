// Draco byte/bit input and entropy decoders (rANS symbols, rABS bits).
//
// Wire facts (checked against streams written by the Draco encoder, which the published
// specification gets wrong in places):
//  * multi-byte values are little-endian; varints are LEB128;
//  * a bit run stores bits least-significant first in each byte and assembles values
//    least-significant bit first; a run always ends on the next whole byte;
//  * rANS / rABS buffers are read backwards from their end.
#pragma once

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "cadcore/error.hpp"

namespace cadcore::codec::draco {

[[noreturn]] void bad(const std::string& what);
[[noreturn]] void unsupported(const std::string& what);

constexpr std::uint16_t version(int major, int minor) { return std::uint16_t(major << 8 | minor); }

class Stream {
public:
    Stream(const std::uint8_t* d, std::size_t n) : d_(d), n_(n) {}
    std::uint16_t ver = 0;  // bitstream version (major << 8 | minor), set after the header

    std::size_t pos() const { return p_; }
    std::size_t remaining() const { return n_ - p_; }
    const std::uint8_t* head() const { return d_ + p_; }
    std::uint8_t u8() { need(1); return d_[p_++]; }
    std::int8_t i8() { return std::int8_t(u8()); }
    std::uint16_t u16() { need(2); const std::uint16_t v = std::uint16_t(d_[p_] | d_[p_ + 1] << 8); p_ += 2; return v; }
    std::uint32_t u32() { need(4); std::uint32_t v; std::memcpy(&v, d_ + p_, 4); p_ += 4; return v; }
    std::int32_t i32() { return std::int32_t(u32()); }
    std::uint64_t u64() { need(8); std::uint64_t v; std::memcpy(&v, d_ + p_, 8); p_ += 8; return v; }
    float f32() { const std::uint32_t u = u32(); float f; std::memcpy(&f, &u, 4); return f; }
    std::uint64_t varint64();
    std::uint32_t varint32();
    const std::uint8_t* bytes(std::size_t n) { need(n); const std::uint8_t* r = d_ + p_; p_ += n; return r; }
    void skip(std::size_t n) { bytes(n); }

    // Bit run: start_bits(), any number of bits(), end_bits() (advances to the next byte).
    void start_bits() { bit_ = 0; in_bits_ = true; }
    std::uint32_t bits(int n);
    void end_bits();

private:
    void need(std::size_t k) const {
        if (k > n_ - p_) fail(ErrorCode::Truncated, "Draco stream ends early", p_);
    }
    const std::uint8_t* d_;
    std::size_t n_, p_ = 0;
    std::uint64_t bit_ = 0;  // bits consumed in the current run
    bool in_bits_ = false;
};

// Size-prefixed bit buffer (EdgeBreaker symbols): reads bits LSB first from a fixed span.
class BitSpan {
public:
    BitSpan() = default;
    BitSpan(const std::uint8_t* d, std::size_t n) : d_(d), n_(n) {}
    std::uint32_t bits(int k);

private:
    const std::uint8_t* d_ = nullptr;
    std::size_t n_ = 0;
    std::uint64_t pos_ = 0;
};

// Binary rABS decoder over a buffer that is read from its end ("RAnsBitDecoder").
class BitDecoder {
public:
    void start(Stream& s);  // prob_zero u8, size, data
    bool next();
    bool started() const { return buf_ != nullptr; }

private:
    const std::uint8_t* buf_ = nullptr;
    std::int64_t off_ = 0;
    std::uint32_t state_ = 0;
    std::uint8_t prob_zero_ = 0;
};

// Symbol decoding ("DecodeSymbols"): tagged or raw rANS-coded unsigned values.
void decode_symbols(Stream& s, std::uint32_t num_values, int num_components, std::uint32_t* out);

}  // namespace cadcore::codec::draco
