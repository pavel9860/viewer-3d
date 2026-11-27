// Compact record bodies in, canonical (varint-dialect) payloads out; shared by the compact transcoders.
#pragma once

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "cadcore/bytes.hpp"
#include "cadcore/error.hpp"
#include "formats/catia/cgm.hpp"

namespace cadcore::formats::catia::cgm::compact_io {

[[noreturn]] inline void bad(const std::string& what, std::size_t at) { fail(ErrorCode::Malformed, "CGM compact record: " + what, at); }

inline double load_f64_le(const std::uint8_t* p) { double d; std::memcpy(&d, p, 8); return d; }

struct In {
    ByteView b;
    std::uint32_t self;
    bool ci;  // small integers as compact signed integers (else raw i16 / i32)
    std::size_t at, k = 0;

    void need(std::size_t n) const { if (k + n > b.size) fail(ErrorCode::Truncated, "CGM compact record ends early", at + k); }
    // low bits 01: +(x >> 2); 11: -(x >> 2) - 1; 00 / 10: (x >> 2 & 3) bytes follow, x's high nibble on top, + / -
    std::int64_t sint() {
        need(1);
        const std::uint8_t x = b[k++];
        const unsigned lo = x & 3;
        if (lo == 1) return x >> 2;
        if (lo == 3) return -std::int64_t(x >> 2) - 1;
        const unsigned n = x >> 2 & 3;
        need(n);
        std::uint64_t v = 0;
        for (unsigned i = 0; i < n; ++i) v |= std::uint64_t(b[k++]) << (8 * i);
        v |= std::uint64_t(x >> 4) << (8 * n);
        return lo == 0 ? std::int64_t(v) : -std::int64_t(v) - 1;
    }
    std::uint32_t ref() { return std::uint32_t(std::int64_t(self) - sint()); }
    std::uint32_t count() { return varint(b, k); }
    std::vector<std::uint32_t> refs(std::size_t n) { std::vector<std::uint32_t> r(n); for (auto& x : r) x = ref(); return r; }
    std::vector<std::uint32_t> refs() { return refs(count()); }
    std::uint8_t u8() { need(1); return b[k++]; }
    double f64() { need(8); const double d = load_f64_le(b.data + k); k += 8; return d; }
    std::int32_t i16() { if (ci) return std::int32_t(sint()); need(2); auto v = std::int16_t(load_u16(b.data + k, false)); k += 2; return v; }
    std::int32_t i32() { if (ci) return std::int32_t(sint()); need(4); auto v = std::int32_t(load_u32(b.data + k, false)); k += 4; return v; }
    std::vector<int> bits(std::size_t n) {
        const std::size_t nb = (n + 7) / 8;
        need(nb);
        std::vector<int> v(n);
        for (std::size_t j = 0; j < n; ++j) v[j] = b[k + j / 8] >> (j % 8) & 1;
        k += nb;
        return v;
    }
    ByteView raw(std::size_t n) { need(n); ByteView v(b.data + k, n); k += n; return v; }
    std::size_t left() const { return b.size - k; }
};

struct Out {
    std::vector<std::uint8_t> p;
    void u8(std::uint8_t v) { p.push_back(v); }
    void i16(std::int32_t v) { for (int i = 0; i < 2; ++i) p.push_back(std::uint8_t(std::uint32_t(v) >> (8 * i))); }
    void i32(std::int32_t v) { for (int i = 0; i < 4; ++i) p.push_back(std::uint8_t(std::uint32_t(v) >> (8 * i))); }
    void f64(double d) { std::uint8_t t[8]; std::memcpy(t, &d, 8); p.insert(p.end(), t, t + 8); }
    void raw(ByteView v) { p.insert(p.end(), v.data, v.data + v.size); }
    void count(std::uint32_t n) {  // varint (cgm.hpp)
        if (n < 0x80) { p.push_back(std::uint8_t(0x80 | n)); return; }
        std::uint8_t mask = 0;
        std::vector<std::uint8_t> bytes;
        for (int j = 0; j < 4; ++j)
            if (std::uint8_t x = std::uint8_t(n >> (8 * j))) { mask |= std::uint8_t(1u << j); bytes.push_back(x); }
        p.push_back(std::uint8_t(mask << 3));
        p.insert(p.end(), bytes.begin(), bytes.end());
    }
};

inline void f64s(In& in, Out& o, std::size_t n) { for (std::size_t i = 0; i < n; ++i) o.f64(in.f64()); }

// B-spline (130, 170) tail as stored: int a, int b, f64 0 L 1 0, u8 c, [c != 0: kept as written], int -2.
// The canonical form drops the closing -2.
inline void bspline_tail(In& in, Out& o) {
    o.i16(in.i16());
    o.i16(in.i16());
    f64s(in, o, 4);
    const std::uint8_t c = in.u8();
    o.u8(c);
    const std::size_t close = in.ci ? 1 : 2;
    if (in.left() < close) bad("B-spline p-curve closing value", in.at);
    if (c) o.raw(in.raw(in.left() - close));
    if (in.i16() != -2) bad("B-spline p-curve closing value", in.at);
}

}  // namespace cadcore::formats::catia::cgm::compact_io
