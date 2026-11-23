// Allocation-free cursor for line-oriented ASCII formats (STL, OBJ, PLY, OFF, DXF, STEP).
#pragma once

#include <cstdint>
#include <string_view>
#include <system_error>

#include "cadcore/error.hpp"
#include <charconv>
#include <cstring>
#include <type_traits>
#include <functional>
#include <vector>

namespace cadcore::text {

// Decimal number at [b, e): sets out, returns the end of the number or nullptr.
// Fast path (Clinger): up to 19 significant digits and a power of ten whose magnitude is at most
// 22, with the digits below 2^53, convert exactly as one correctly rounded multiply or divide.
// A float goes through that double; double rounding can only differ from direct rounding when
// the double is exactly halfway between two floats, and then (and for everything outside the
// fast path) std::from_chars decides. Results are identical to std::from_chars.
template <class T> const char* parse_real(const char* b, const char* e, T& out) {
    static constexpr double pow10[23] = {1e0,  1e1,  1e2,  1e3,  1e4,  1e5,  1e6,  1e7,  1e8,  1e9,  1e10, 1e11,
                                         1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
    const char* p = b;
    const bool neg = p < e && *p == '-';
    if (neg) ++p;
    std::uint64_t m = 0;
    int digits = 0, exp10 = 0;
    const char* d0 = p;
    while (p < e && unsigned(*p - '0') < 10) { if (digits < 19) { m = m * 10 + unsigned(*p - '0'); if (m) ++digits; } else ++exp10; ++p; }
    const bool int_part = p != d0;
    bool frac_part = false;
    if (p < e && *p == '.') {
        ++p;
        const char* f0 = p;
        while (p < e && unsigned(*p - '0') < 10) { if (digits < 19) { m = m * 10 + unsigned(*p - '0'); if (m) ++digits; --exp10; } ++p; }
        frac_part = p != f0;
    }
    if (!int_part && !frac_part) {  // inf, nan, or not a number
        const auto r = std::from_chars(b, e, out, std::chars_format::general);
        return r.ec == std::errc() ? r.ptr : nullptr;
    }
    if (p < e && (*p == 'e' || *p == 'E')) {
        const char* q = p + 1;
        const bool eneg = q < e && *q == '-';
        if (q < e && (*q == '-' || *q == '+')) ++q;
        if (q < e && unsigned(*q - '0') < 10) {
            int x = 0;
            while (q < e && unsigned(*q - '0') < 10) { if (x < 100000) x = x * 10 + (*q - '0'); ++q; }
            exp10 += eneg ? -x : x;
            p = q;
        }
    }
    if (m <= (std::uint64_t(1) << 53) && exp10 >= -22 && exp10 <= 22) {
        double v = double(m);
        v = exp10 < 0 ? v / pow10[-exp10] : v * pow10[exp10];
        if (neg) v = -v;
        if constexpr (sizeof(T) == sizeof(float)) {
            std::uint64_t bits;
            std::memcpy(&bits, &v, 8);
            const int be = int(bits >> 52 & 0x7ff);
            const bool midpoint = (bits & 0x1fffffff) == 0x10000000;
            if (!midpoint && be > 1023 - 126 && be < 1023 + 127) { out = T(v); return p; }
        } else {
            out = T(v);
            return p;
        }
    }
    const auto r = std::from_chars(b, e, out, std::chars_format::general);
    return r.ec == std::errc() ? r.ptr : nullptr;
}

inline bool is_space(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\f' || c == '\v'; }

class Cursor {
public:
    explicit Cursor(std::string_view s, std::size_t base_offset = 0) : s_(s), p_(0), base_(base_offset) {}

    // When enabled, '#' starts a comment that runs to the end of the line (OBJ, OFF, commented STL).
    void set_hash_comments(bool on) { hash_ = on; }

    bool eof() const { return p_ >= s_.size(); }
    std::size_t offset() const { return base_ + p_; }
    char peek() const { return eof() ? '\0' : s_[p_]; }

    void skip_spaces() { while (p_ < s_.size() && is_space(s_[p_])) ++p_; }
    void skip_ws_and_newlines() {
        while (p_ < s_.size()) {
            if (is_space(s_[p_]) || s_[p_] == '\n') ++p_;
            else if (hash_ && s_[p_] == '#') next_line();
            else break;
        }
    }
    bool at_eol() { skip_spaces(); return eof() || s_[p_] == '\n' || (hash_ && s_[p_] == '#'); }
    void next_line() {
        while (p_ < s_.size() && s_[p_] != '\n') ++p_;
        if (p_ < s_.size()) ++p_;
    }
    // Rest of the current line (without the newline, trailing spaces trimmed); advances past it.
    std::string_view rest_of_line() {
        skip_spaces();
        std::size_t start = p_;
        while (p_ < s_.size() && s_[p_] != '\n' && !(hash_ && s_[p_] == '#')) ++p_;
        std::size_t end = p_;
        while (end > start && is_space(s_[end - 1])) --end;
        next_line();
        return s_.substr(start, end - start);
    }
    // Next whitespace-delimited token on the current line; empty at end of line.
    std::string_view token() {
        skip_spaces();
        std::size_t start = p_;
        if (hash_ && p_ < s_.size() && s_[p_] == '#') return s_.substr(start, 0);
        while (p_ < s_.size() && !is_space(s_[p_]) && s_[p_] != '\n' && !(hash_ && s_[p_] == '#')) ++p_;
        return s_.substr(start, p_ - start);
    }
    // The next token on the current line, without consuming it.
    std::string_view token_peek() {
        std::size_t save = p_;
        std::string_view t = token();
        p_ = save;
        return t;
    }
    // Next token across line breaks.
    std::string_view token_any() {
        skip_ws_and_newlines();
        return token();
    }

    template <class T> bool try_number(T& out) {
        skip_spaces();
        const char* b = s_.data() + p_;
        const char* e = s_.data() + s_.size();
        if (b < e && *b == '+') ++b;
        const char* end;
        if constexpr (std::is_floating_point_v<T>) end = parse_real(b, e, out);
        else { auto r = std::from_chars(b, e, out); end = r.ec == std::errc() ? r.ptr : nullptr; }
        if (!end) return false;
        p_ = static_cast<std::size_t>(end - s_.data());
        return true;
    }
    template <class T> T number(const char* what) {
        T v{};
        if (!try_number(v)) fail(ErrorCode::Malformed, std::string("expected number for ") + what, offset());
        return v;
    }
    bool try_int(long long& out) {
        skip_spaces();
        const char* b = s_.data() + p_;
        const char* e = s_.data() + s_.size();
        if (b < e && *b == '+') ++b;
        auto r = std::from_chars(b, e, out);
        if (r.ec != std::errc()) return false;
        p_ = static_cast<std::size_t>(r.ptr - s_.data());
        return true;
    }
    long long integer(const char* what) {
        long long v = 0;
        if (!try_int(v)) fail(ErrorCode::Malformed, std::string("expected integer for ") + what, offset());
        return v;
    }
    long long integer_any(const char* what) {
        skip_ws_and_newlines();
        return integer(what);
    }
    bool consume(char c) {
        if (peek() == c) { ++p_; return true; }
        return false;
    }
    std::size_t position() const { return p_; }
    std::string_view text() const { return s_; }
    void set_position(std::size_t p) { p_ = p; }

private:
    std::string_view s_;
    std::size_t p_;
    std::size_t base_;
    bool hash_ = false;
};

// Text of a buffer as UTF-8: strips a UTF-8 BOM, transcodes UTF-16 LE/BE (with BOM).
// `storage` holds the transcoded text when a copy is needed.
std::string_view as_utf8(std::string_view raw, std::string& storage);

// Start offsets of the lines of t[from, end) (after each '\n'), found on up to `threads` threads.
std::vector<std::size_t> line_starts(std::string_view t, std::size_t from, unsigned threads);

// Parallel parse of `count` records written one per line, the first at line start `pos` of a text
// whose line starts are `lines` (line_starts). Part k gets records [first, first + n) starting at
// its line, and parse(k, line_offset, first, n) returns the offset after its last record. Returns
// the offset after record count - 1, or npos when the layout is not one record per line (a part
// ends elsewhere than the next part's line, or throws): the caller then parses sequentially.
// Requires count >= 2 * threads; parts = 1 returns npos at once.
std::size_t parallel_line_records(const std::vector<std::size_t>& lines, std::size_t pos, std::size_t count, std::size_t text_size, unsigned threads,
                                  std::size_t& parts, const std::function<std::size_t(std::size_t, std::size_t, std::size_t, std::size_t)>& parse);

inline std::string_view trim(std::string_view s) {
    auto ws = [](char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; };
    while (!s.empty() && ws(s.front())) s.remove_prefix(1);
    while (!s.empty() && ws(s.back())) s.remove_suffix(1);
    return s;
}

inline bool iequals_prefix(std::string_view s, std::string_view prefix) {
    if (s.size() < prefix.size()) return false;
    for (std::size_t i = 0; i < prefix.size(); ++i) {
        char a = s[i], b = prefix[i];
        if (a >= 'A' && a <= 'Z') a = char(a - 'A' + 'a');
        if (b >= 'A' && b <= 'Z') b = char(b - 'A' + 'a');
        if (a != b) return false;
    }
    return true;
}

}  // namespace cadcore::text
