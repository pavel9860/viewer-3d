#include "formats/dxf/groups.hpp"

#include <charconv>
#include <exception>
#include <functional>
#include <cmath>
#include <string>

#include "cadcore/error.hpp"
#include "core/parallel.hpp"
#include "core/text.hpp"

namespace cadcore::formats::autocad {
namespace {

// Lines of a text, each without its terminator (\n, \r\n or a lone \r).
struct Lines {
    std::string_view t;
    std::size_t p = 0;
    bool next(std::string_view& out) {
        if (p >= t.size()) return false;
        const std::size_t start = p;
        while (p < t.size() && t[p] != '\n' && t[p] != '\r') ++p;
        out = t.substr(start, p - start);
        if (p < t.size() && t[p] == '\r') ++p;
        if (p < t.size() && t[p] == '\n') ++p;
        return true;
    }
};

std::string_view without_bom(std::string_view t) {
    if (t.size() >= 3 && std::uint8_t(t[0]) == 0xEF && std::uint8_t(t[1]) == 0xBB && std::uint8_t(t[2]) == 0xBF) t.remove_prefix(3);
    return t;
}

constexpr std::string_view kSentinel{"AutoCAD Binary DXF\r\n\x1a\0", 22};

bool in(int c, int lo, int hi) { return c >= lo && c <= hi; }

using text::trim;

bool parse_double(std::string_view s, double& out) {
    s = trim(s);
    if (!s.empty() && s.front() == '+') s.remove_prefix(1);
    return text::parse_real(s.data(), s.data() + s.size(), out) == s.data() + s.size();
}

bool parse_int(std::string_view s, long long& out) {
    s = trim(s);
    if (!s.empty() && s.front() == '+') s.remove_prefix(1);
    auto r = std::from_chars(s.data(), s.data() + s.size(), out);
    if (r.ec == std::errc() && r.ptr == s.data() + s.size()) return true;
    // Some writers put a decimal number into an integer group ("1.0"); accept integral values.
    double d;
    if (parse_double(s, d) && std::isfinite(d) && d == std::floor(d) && std::fabs(d) < 9e18) {
        out = static_cast<long long>(d);
        return true;
    }
    return false;
}

// Groups of t[begin, end) (begin at a code line) into out; stops after "0 EOF" (returns true).
template <class Out> bool read_ascii_range(std::string_view t, std::size_t base, std::size_t begin, std::size_t end, bool first, Out& out) {
    Lines lines{t.substr(0, end), begin};
    std::size_t& p = lines.p;
    std::string_view code_line, value;
    while (true) {
        const std::size_t at = p;
        if (!lines.next(code_line)) return false;
        if (trim(code_line).empty() && p >= t.size()) return false;  // trailing blank line
        long long code;
        if (!parse_int(code_line, code)) {
            // AutoCAD 1.x wrote "NAME,1" records with comma-separated values: not group-code DXF
            if (first && out.empty() && code_line.find(',') != std::string_view::npos)
                fail(ErrorCode::Unsupported, "AutoCAD 1.x interchange file (comma-separated records, before group-code DXF)", base + at);
            fail(ErrorCode::Malformed, "DXF group code is not an integer", base + at);
        }
        if (!lines.next(value)) fail(ErrorCode::Truncated, "DXF group code without a value", base + at);
        Group g;
        g.code = int(code);
        g.offset = base + at;
        g.str = value;
        switch (value_type(g.code)) {
            case ValueType::Double:
                if (!parse_double(value, g.num)) fail(ErrorCode::Malformed, "DXF group " + std::to_string(code) + ": not a number", base + at);
                break;
            case ValueType::Int:
            case ValueType::Bool:
                if (!parse_int(value, g.ival)) fail(ErrorCode::Malformed, "DXF group " + std::to_string(code) + ": not an integer", base + at);
                g.num = double(g.ival);
                break;
            case ValueType::Binary:
                g.str = trim(value);
                g.hex = true;
                break;
            default:
                g.str = trim(value);
                break;
        }
        out.push_back(g);
        if (g.code == 0 && g.str == "EOF") return true;
    }
}

bool is_terminator(std::string_view t, std::size_t i) { return t[i] == '\n' || (t[i] == '\r' && (i + 1 == t.size() || t[i + 1] != '\n')); }

// Large files: pieces cut at line starts with an even line number (code lines), read in
// parallel and joined in order up to the first "0 EOF". A piece's error counts only when no
// EOF comes before it, as in a sequential read.
GroupStream read_ascii(ByteView data, unsigned threads) {
    GroupStream gs;
    const std::string_view t = without_bom(data.text());
    const std::size_t base = data.size - t.size();
    constexpr std::size_t kMinPiece = 2 << 20;
    const std::size_t parts = std::max<std::size_t>(1, std::min<std::size_t>(thread_count(threads), t.size() / kMinPiece));
    if (parts == 1) {
        gs.saw_eof = read_ascii_range(t, base, 0, t.size(), true, gs.groups);
        return gs;
    }
    // line terminators before each piece's nominal start: global line numbers
    std::vector<std::size_t> lo(parts + 1);
    for (std::size_t k = 0; k <= parts; ++k) lo[k] = t.size() * k / parts;
    std::vector<Slot<std::size_t>> count(parts);
    parallel_for(parts, 1, threads, [&](std::size_t b0, std::size_t e0) {
        for (std::size_t k = b0; k < e0; ++k) {
            std::size_t c = 0;
            for (std::size_t i = lo[k]; i < lo[k + 1]; ++i) c += is_terminator(t, i);
            count[k].value = c;
        }
    });
    std::vector<std::size_t> cut{0}, line_no{0};  // piece starts and their line numbers
    std::size_t lines_before = 0;                   // terminators in [0, lo[k])
    for (std::size_t k = 1; k < parts; ++k) {
        lines_before += count[k - 1].value;
        // the line starting after the j-th terminator from lo[k] has number lines_before + j
        std::size_t number = lines_before, at = t.size();
        for (std::size_t i = lo[k]; i < t.size(); ++i)
            if (is_terminator(t, i)) {
                ++number;
                if (number % 2 == 0) { at = i + 1; break; }
            }
        if (at <= cut.back()) { at = cut.back(); number = line_no.back(); }
        cut.push_back(at);
        line_no.push_back(number);
    }
    lines_before += count[parts - 1].value;
    cut.push_back(t.size());
    line_no.push_back(lines_before + 1);  // a last line without a terminator
    // groups of piece k: its lines / 2; all pieces write in place into one array
    std::vector<std::size_t> first_group(parts + 1, 0);
    for (std::size_t k = 0; k < parts; ++k) first_group[k + 1] = first_group[k] + (line_no[k + 1] - line_no[k] + 1) / 2;
    gs.groups.resize(first_group[parts]);
    struct Slots {  // writes into gs.groups[from, to)
        Group* at;
        Group* end;
        std::size_t n = 0;
        bool empty() const { return n == 0; }
        void push_back(const Group& g) {
            if (at + n == end) fail(ErrorCode::Malformed, "DXF: more groups than lines");
            at[n++] = g;
        }
    };
    std::vector<Slot<std::size_t>> written(parts);
    std::vector<Slot<bool>> eof(parts);
    std::vector<std::exception_ptr> error(parts);
    std::vector<std::function<void()>> tasks;
    for (std::size_t k = 0; k < parts; ++k)
        tasks.push_back([&, k] {
            Slots out{gs.groups.data() + first_group[k], gs.groups.data() + first_group[k + 1]};
            try {
                eof[k].value = read_ascii_range(t, base, cut[k], cut[k + 1], k == 0, out);
            } catch (...) {
                error[k] = std::current_exception();
            }
            written[k].value = out.n;
        });
    run_all(tasks, threads);
    // join in order up to the first EOF: pieces that wrote fewer groups than lines/2 (blank
    // trailing line, EOF) are closed up
    std::size_t n = 0;
    for (std::size_t k = 0; k < parts; ++k) {
        if (error[k]) std::rethrow_exception(error[k]);
        if (n != first_group[k]) std::move(gs.groups.begin() + std::ptrdiff_t(first_group[k]), gs.groups.begin() + std::ptrdiff_t(first_group[k] + written[k].value), gs.groups.begin() + std::ptrdiff_t(n));
        n += written[k].value;
        if (eof[k].value) { gs.saw_eof = true; break; }
    }
    gs.groups.resize(n);
    return gs;
}

GroupStream read_binary(ByteView data) {
    GroupStream gs;
    gs.binary = true;
    // R12 and older: one-byte codes. The first group is always "0 SECTION": with one-byte codes
    // its string starts right after the code byte.
    const bool short_codes = data.size > kSentinel.size() + 1 && data[kSentinel.size() + 1] == 'S';
    Reader r(data, kSentinel.size());
    while (r.remaining() > 0) {
        Group g;
        g.offset = r.pos();
        if (short_codes) {
            const std::uint8_t c = r.u8();
            g.code = c == 255 ? int(std::int16_t(r.u16le())) : int(c);
        } else {
            g.code = int(std::int16_t(r.u16le()));
        }
        switch (value_type(g.code)) {
            case ValueType::String: {
                const std::size_t start = r.pos();
                while (r.u8() != 0) {}
                g.str = data.sub(start, r.pos() - 1 - start).text();
                break;
            }
            case ValueType::Double: g.num = r.f64le(); break;
            case ValueType::Bool: g.ival = r.u8(); g.num = double(g.ival); break;
            case ValueType::Int:
                if (in(g.code, 90, 99) || in(g.code, 420, 429) || in(g.code, 440, 459) || g.code == 1071) g.ival = std::int32_t(r.u32le());
                else if (in(g.code, 160, 169)) g.ival = static_cast<long long>(r.u64le());
                else g.ival = std::int16_t(r.u16le());
                g.num = double(g.ival);
                break;
            case ValueType::Binary: {
                const std::size_t n = r.u8();
                g.str = ByteView(r.bytes(n), n).text();
                break;
            }
        }
        gs.groups.push_back(g);
        if (g.code == 0 && g.str == "EOF") { gs.saw_eof = true; break; }
    }
    return gs;
}

}  // namespace

ValueType value_type(int c) {
    if (in(c, 0, 9) || c == 100 || c == 102 || c == 105 || in(c, 300, 309) || in(c, 320, 369) || in(c, 390, 399) ||
        in(c, 410, 419) || in(c, 430, 439) || in(c, 470, 481) || c == 999 || in(c, 1000, 1003) || in(c, 1005, 1009))
        return ValueType::String;
    if (in(c, 10, 59) || in(c, 110, 149) || in(c, 210, 239) || in(c, 460, 469) || in(c, 1010, 1059)) return ValueType::Double;
    if (in(c, 290, 299)) return ValueType::Bool;
    if (in(c, 310, 319) || c == 1004) return ValueType::Binary;
    if (in(c, 60, 99) || in(c, 160, 179) || in(c, 270, 289) || in(c, 370, 389) || in(c, 400, 409) || in(c, 420, 459) || in(c, 1060, 1071))
        return ValueType::Int;
    return ValueType::String;  // unassigned codes: keep the text
}

bool is_binary_dxf(ByteView data) { return data.starts_with(kSentinel); }

bool looks_like_ascii_dxf(ByteView data) {
    Lines lines{without_bom(data.text().substr(0, std::min<std::size_t>(data.size, 4096)))};
    std::string_view code, value;
    for (int i = 0; i < 16 && lines.next(code) && lines.next(value); ++i) {
        if (trim(code) == "999") continue;  // comment
        return trim(code) == "0" && trim(value) == "SECTION";
    }
    return false;
}

GroupStream read_groups(ByteView data, unsigned threads) { return is_binary_dxf(data) ? read_binary(data) : read_ascii(data, threads); }

}  // namespace cadcore::formats::autocad
