#include "codec/json.hpp"

#include <charconv>
#include <cmath>

#include "cadcore/error.hpp"

namespace cadcore::json {

const Value* Value::get(std::string_view key) const {
    if (type != Type::Object) return nullptr;
    for (const auto& [k, v] : object)
        if (k == key) return &v;
    return nullptr;
}

namespace {

class Parser {
public:
    explicit Parser(std::string_view s) : s_(s) {}

    Value document() {
        ws();
        Value v = value(0);
        ws();
        if (p_ != s_.size()) error("unexpected content after the JSON value");
        return v;
    }

private:
    [[noreturn]] void error(const char* what) { fail(ErrorCode::Malformed, std::string("JSON: ") + what, p_); }
    void ws() { while (p_ < s_.size() && (s_[p_] == ' ' || s_[p_] == '\t' || s_[p_] == '\n' || s_[p_] == '\r')) ++p_; }
    bool literal(std::string_view w) {
        if (s_.compare(p_, w.size(), w) == 0) { p_ += w.size(); return true; }
        return false;
    }

    Value value(int depth) {
        if (depth > 512) error("nesting deeper than 512");
        if (p_ >= s_.size()) error("unexpected end of input");
        Value v;
        const char c = s_[p_];
        if (c == '{') {
            v.type = Type::Object;
            ++p_;
            ws();
            if (p_ < s_.size() && s_[p_] == '}') { ++p_; return v; }
            while (true) {
                ws();
                if (p_ >= s_.size() || s_[p_] != '"') error("expected a member name");
                std::string key = string();
                ws();
                if (p_ >= s_.size() || s_[p_] != ':') error("expected ':'");
                ++p_;
                ws();
                v.object.emplace_back(std::move(key), value(depth + 1));
                ws();
                if (p_ < s_.size() && s_[p_] == ',') { ++p_; continue; }
                if (p_ < s_.size() && s_[p_] == '}') { ++p_; return v; }
                error("expected ',' or '}'");
            }
        }
        if (c == '[') {
            v.type = Type::Array;
            ++p_;
            ws();
            if (p_ < s_.size() && s_[p_] == ']') { ++p_; return v; }
            while (true) {
                ws();
                v.array.push_back(value(depth + 1));
                ws();
                if (p_ < s_.size() && s_[p_] == ',') { ++p_; continue; }
                if (p_ < s_.size() && s_[p_] == ']') { ++p_; return v; }
                error("expected ',' or ']'");
            }
        }
        if (c == '"') { v.type = Type::String; v.string = string(); return v; }
        if (literal("true")) { v.type = Type::Bool; v.boolean = true; return v; }
        if (literal("false")) { v.type = Type::Bool; return v; }
        if (literal("null")) return v;
        if (c == '-' || (c >= '0' && c <= '9')) return number();
        error("unexpected character");
    }

    Value number() {
        const std::size_t start = p_;
        if (s_[p_] == '-') ++p_;
        if (p_ >= s_.size()) error("bad number");
        if (s_[p_] == '0') ++p_;
        else if (s_[p_] >= '1' && s_[p_] <= '9') while (p_ < s_.size() && s_[p_] >= '0' && s_[p_] <= '9') ++p_;
        else error("bad number");
        bool integral = true;
        if (p_ < s_.size() && s_[p_] == '.') {
            integral = false;
            ++p_;
            if (p_ >= s_.size() || s_[p_] < '0' || s_[p_] > '9') error("bad number fraction");
            while (p_ < s_.size() && s_[p_] >= '0' && s_[p_] <= '9') ++p_;
        }
        if (p_ < s_.size() && (s_[p_] == 'e' || s_[p_] == 'E')) {
            integral = false;
            ++p_;
            if (p_ < s_.size() && (s_[p_] == '+' || s_[p_] == '-')) ++p_;
            if (p_ >= s_.size() || s_[p_] < '0' || s_[p_] > '9') error("bad number exponent");
            while (p_ < s_.size() && s_[p_] >= '0' && s_[p_] <= '9') ++p_;
        }
        Value v;
        v.type = Type::Number;
        const char* b = s_.data() + start;
        const char* e = s_.data() + p_;
        auto r = std::from_chars(b, e, v.number);
        if (r.ec == std::errc::result_out_of_range) v.number = (*b == '-') ? -HUGE_VAL : HUGE_VAL;
        else if (r.ec != std::errc() || r.ptr != e) error("bad number");
        if (integral) {
            std::int64_t i = 0;
            auto ri = std::from_chars(b, e, i);
            if (ri.ec == std::errc() && ri.ptr == e) { v.is_integer = true; v.integer = i; }
        }
        return v;
    }

    unsigned hex4() {
        if (p_ + 4 > s_.size()) error("truncated \\u escape");
        unsigned v = 0;
        for (int i = 0; i < 4; ++i) {
            const char c = s_[p_++];
            const int d = c >= '0' && c <= '9' ? c - '0' : c >= 'a' && c <= 'f' ? c - 'a' + 10 : c >= 'A' && c <= 'F' ? c - 'A' + 10 : -1;
            if (d < 0) error("bad \\u escape");
            v = v * 16 + unsigned(d);
        }
        return v;
    }

    std::string string() {
        ++p_;  // opening quote
        std::string out;
        while (true) {
            if (p_ >= s_.size()) error("unterminated string");
            const char c = s_[p_++];
            if (c == '"') return out;
            if (static_cast<unsigned char>(c) < 0x20) error("control character in string");
            if (c != '\\') { out += c; continue; }
            if (p_ >= s_.size()) error("unterminated escape");
            const char e = s_[p_++];
            switch (e) {
                case '"': out += '"'; break;
                case '\\': out += '\\'; break;
                case '/': out += '/'; break;
                case 'b': out += '\b'; break;
                case 'f': out += '\f'; break;
                case 'n': out += '\n'; break;
                case 'r': out += '\r'; break;
                case 't': out += '\t'; break;
                case 'u': {
                    unsigned cp = hex4();
                    if (cp >= 0xD800 && cp < 0xDC00) {
                        if (p_ + 6 > s_.size() || s_[p_] != '\\' || s_[p_ + 1] != 'u') error("unpaired surrogate");
                        p_ += 2;
                        const unsigned lo = hex4();
                        if (lo < 0xDC00 || lo >= 0xE000) error("bad surrogate pair");
                        cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
                    } else if (cp >= 0xDC00 && cp < 0xE000) {
                        error("unpaired surrogate");
                    }
                    if (cp < 0x80) out += char(cp);
                    else if (cp < 0x800) { out += char(0xC0 | cp >> 6); out += char(0x80 | (cp & 63)); }
                    else if (cp < 0x10000) { out += char(0xE0 | cp >> 12); out += char(0x80 | (cp >> 6 & 63)); out += char(0x80 | (cp & 63)); }
                    else { out += char(0xF0 | cp >> 18); out += char(0x80 | (cp >> 12 & 63)); out += char(0x80 | (cp >> 6 & 63)); out += char(0x80 | (cp & 63)); }
                    break;
                }
                default: error("bad escape");
            }
        }
    }

    std::string_view s_;
    std::size_t p_ = 0;
};

}  // namespace

Value parse(std::string_view text) { return Parser(text).document(); }

}  // namespace cadcore::json
