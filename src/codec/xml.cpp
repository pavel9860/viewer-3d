#include "codec/xml.hpp"

#include <algorithm>
#include <cstring>

#include "cadcore/error.hpp"
#include "core/parallel.hpp"
#include "core/text.hpp"

namespace cadcore::xml {

std::string_view local_name(std::string_view q) {
    const auto c = q.find(':');
    return c == std::string_view::npos ? q : q.substr(c + 1);
}

std::string_view Element::local() const { return local_name(name); }

std::string_view Element::attr(std::string_view n) const {
    for (const auto& a : attributes)
        if (a.name == n || local_name(a.name) == n) return a.value;
    return {};
}

bool Element::has_attr(std::string_view n) const {
    for (const auto& a : attributes)
        if (a.name == n || local_name(a.name) == n) return true;
    return false;
}

const Element* Element::child(std::string_view ln) const {
    for (const Element* c : children)
        if (c->local() == ln) return c;
    return nullptr;
}

namespace {

bool is_ws(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }
bool is_name_char(char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == ':' || c == '-' || c == '.' ||
           (static_cast<unsigned char>(c) >= 0x80);
}

void append_utf8(char*& out, std::uint32_t cp) {
    if (cp < 0x80) { *out++ = char(cp); }
    else if (cp < 0x800) { *out++ = char(0xC0 | cp >> 6); *out++ = char(0x80 | (cp & 63)); }
    else if (cp < 0x10000) { *out++ = char(0xE0 | cp >> 12); *out++ = char(0x80 | (cp >> 6 & 63)); *out++ = char(0x80 | (cp & 63)); }
    else { *out++ = char(0xF0 | cp >> 18); *out++ = char(0x80 | (cp >> 12 & 63)); *out++ = char(0x80 | (cp >> 6 & 63)); *out++ = char(0x80 | (cp & 63)); }
}

class Parser {
public:
    // Parses b in place: entity references are decoded into the same bytes (w = b's writable start).
    Parser(std::string_view b, char* w, std::deque<Element>& nodes, std::size_t base = 0, const std::vector<std::string_view>* raw = nullptr)
        : b_(b), w_(w), nodes_(nodes), p_(0), base_(base), raw_(raw) {}

    // Fragment mode: the buffer is element content; calls fn for each top-level child element.
    // Nodes are reused from child to child (no allocation per element once warm).
    template <class F> void fragment(F&& fn) {
        while (true) {
            const std::size_t lt = b_.find('<', p_);
            if (lt == std::string_view::npos) return;  // trailing text is ignored, as in element content
            p_ = lt;
            if (starts("<!--")) { skip_until("-->"); continue; }
            if (starts("<![CDATA[")) { skip_until("]]>"); continue; }
            if (starts("<?")) { skip_until("?>"); continue; }
            if (starts("</")) error("end tag without a start tag");
            used_ = 0;
            fn(*element(0));
        }
    }

    Element* parse() {
        skip_misc();
        if (p_ >= b_.size() || b_[p_] != '<') error("expected the root element");
        Element* root = element(0);
        skip_misc();
        if (p_ != b_.size()) error("content after the root element");
        return root;
    }

private:
    [[noreturn]] void error(const char* what) { fail(ErrorCode::Malformed, std::string("XML: ") + what, base_ + p_); }

    Element* new_node() {
        if (used_ < nodes_.size()) {
            Element& e = nodes_[used_++];
            e.attributes.clear();
            e.children.clear();
            e.text = {};
            e.raw = {};
            e.is_raw = false;
            e.raw_flat = false;
            return &e;
        }
        nodes_.emplace_back();
        ++used_;
        return &nodes_.back();
    }
    bool keep_raw(std::string_view n) const {
        if (!raw_) return false;
        const std::string_view l = local_name(n);
        for (auto r : *raw_)
            if (r == l) return true;
        return false;
    }
    // After the start tag of `el` (p_ past '>'): finds its end tag without building nodes and
    // stores the content as raw. Its well-formedness is checked when for_each_element parses it.
    void skip_raw(Element* el) {
        const std::size_t begin = p_;
        // Fast path, one memchr pass over the '<' (which cannot occur in attribute values): the first
        // "</" closes the element when no comment, CDATA section or PI comes before it; children
        // before it can only be self-closing (a nested same-name one then does not change depth).
        for (std::size_t lt = begin;;) {
            const void* hit = std::memchr(b_.data() + lt, '<', b_.size() - lt);
            if (!hit) break;
            lt = std::size_t(static_cast<const char*>(hit) - b_.data());
            const char next = lt + 1 < b_.size() ? b_[lt + 1] : '\0';
            if (next == '!' || next == '?') break;  // slow path
            if (next != '/') { ++lt; continue; }
            p_ = lt + 2;
            if (name() != el->name) error("mismatched end tag");
            skip_ws();
            if (p_ >= b_.size() || b_[p_] != '>') error("expected '>'");
            ++p_;
            el->raw = b_.substr(begin, lt - begin);
            el->raw_offset = base_ + begin;
            el->is_raw = true;
            el->raw_flat = true;
            return;
        }
        p_ = begin;
        int depth = 0;
        while (true) {
            const std::size_t lt = b_.find('<', p_);
            if (lt == std::string_view::npos) error("unterminated element");
            p_ = lt;
            if (starts("<!--")) { skip_until("-->"); continue; }
            if (starts("<![CDATA[")) { skip_until("]]>"); continue; }
            if (starts("<?")) { skip_until("?>"); continue; }
            if (starts("</")) {
                const std::size_t end = p_;
                p_ += 2;
                const std::string_view n = name();
                skip_ws();
                if (p_ >= b_.size() || b_[p_] != '>') error("expected '>'");
                ++p_;
                if (depth == 0) {
                    if (n != el->name) error("mismatched end tag");
                    el->raw = b_.substr(begin, end - begin);
                    el->raw_offset = base_ + begin;
                    el->is_raw = true;
                    return;
                }
                --depth;
                continue;
            }
            // a start tag: find its '>' outside quoted attribute values
            ++p_;
            char q = 0;
            for (; p_ < b_.size(); ++p_) {
                const char ch = b_[p_];
                if (q) { if (ch == q) q = 0; }
                else if (ch == '"' || ch == '\'') q = ch;
                else if (ch == '>') break;
            }
            if (p_ >= b_.size()) error("unterminated start tag");
            if (b_[p_ - 1] != '/') ++depth;
            ++p_;
        }
    }

    bool starts(const char* s) const { return b_.compare(p_, std::strlen(s), s) == 0; }
    void skip_ws() { while (p_ < b_.size() && is_ws(b_[p_])) ++p_; }
    void skip_until(const char* end) {
        const auto e = b_.find(end, p_);
        if (e == std::string_view::npos) error("unterminated construct");
        p_ = e + std::strlen(end);
    }
    // Prolog, comments, PIs, DOCTYPE and whitespace outside the root.
    void skip_misc() {
        while (true) {
            skip_ws();
            if (starts("<?")) skip_until("?>");
            else if (starts("<!--")) skip_until("-->");
            else if (starts("<!DOCTYPE")) skip_doctype();
            else return;
        }
    }
    void skip_doctype() {
        int depth = 0;
        for (; p_ < b_.size(); ++p_) {
            if (b_[p_] == '[') ++depth;
            else if (b_[p_] == ']') --depth;
            else if (b_[p_] == '>' && depth <= 0) { ++p_; return; }
        }
        error("unterminated DOCTYPE");
    }
    std::string_view name() {
        const std::size_t s = p_;
        while (p_ < b_.size() && is_name_char(b_[p_])) ++p_;
        if (p_ == s) error("expected a name");
        return b_.substr(s, p_ - s);
    }
    // Decodes entities in [s, e) in place; returns the decoded view.
    std::string_view decode(std::size_t s, std::size_t e) {
        char* out = w_ + s;
        char* const begin = out;
        for (std::size_t i = s; i < e;) {
            if (b_[i] != '&') { *out++ = w_[i++]; continue; }
            const std::size_t semi = b_.find(';', i);
            if (semi == std::string_view::npos || semi > e) { p_ = i; error("unterminated entity reference"); }
            const std::string_view ent(b_.data() + i + 1, semi - i - 1);
            if (ent == "lt") *out++ = '<';
            else if (ent == "gt") *out++ = '>';
            else if (ent == "amp") *out++ = '&';
            else if (ent == "quot") *out++ = '"';
            else if (ent == "apos") *out++ = '\'';
            else if (ent.size() > 1 && ent[0] == '#') {
                std::uint32_t cp = 0;
                const bool hex = ent[1] == 'x' || ent[1] == 'X';
                for (char c : ent.substr(hex ? 2 : 1)) {
                    const int v = c >= '0' && c <= '9' ? c - '0' : hex && c >= 'a' && c <= 'f' ? c - 'a' + 10 : hex && c >= 'A' && c <= 'F' ? c - 'A' + 10 : -1;
                    if (v < 0) { p_ = i; error("bad character reference"); }
                    cp = cp * (hex ? 16 : 10) + std::uint32_t(v);
                    if (cp > 0x10FFFF) { p_ = i; error("character reference out of range"); }
                }
                append_utf8(out, cp);  // at most as long as the reference itself
            } else {
                p_ = i;
                error("unknown entity");
            }
            i = semi + 1;
        }
        return std::string_view(begin, std::size_t(out - begin));
    }
    Element* element(int depth) {
        if (depth > 4096) error("element nesting deeper than 4096");
        ++p_;  // '<'
        Element* el = new_node();
        el->name = name();
        while (true) {
            skip_ws();
            if (p_ >= b_.size()) error("unterminated start tag");
            if (b_[p_] == '/') {
                if (p_ + 1 >= b_.size() || b_[p_ + 1] != '>') error("expected '/>'");
                p_ += 2;
                return el;
            }
            if (b_[p_] == '>') {
                ++p_;
                if (keep_raw(el->name)) { skip_raw(el); return el; }
                break;
            }
            Attribute a;
            a.name = name();
            skip_ws();
            if (p_ >= b_.size() || b_[p_] != '=') error("expected '=' after attribute name");
            ++p_;
            skip_ws();
            if (p_ >= b_.size() || (b_[p_] != '"' && b_[p_] != '\'')) error("expected a quoted attribute value");
            const char q = b_[p_++];
            const std::size_t e = b_.find(q, p_);
            if (e == std::string_view::npos) error("unterminated attribute value");
            a.value = decode(p_, e);
            p_ = e + 1;
            el->attributes.push_back(a);
        }
        bool have_text = false;
        while (true) {
            const std::size_t lt = b_.find('<', p_);
            if (lt == std::string_view::npos) error("unterminated element");
            if (lt > p_ && !have_text) {
                std::string_view t = decode(p_, lt);
                std::size_t a = 0, z = t.size();
                while (a < z && is_ws(t[a])) ++a;
                while (z > a && is_ws(t[z - 1])) --z;
                if (z > a) { el->text = t.substr(a, z - a); have_text = true; }
            }
            p_ = lt;
            if (starts("</")) {
                p_ += 2;
                if (name() != el->name) error("mismatched end tag");
                skip_ws();
                if (p_ >= b_.size() || b_[p_] != '>') error("expected '>'");
                ++p_;
                return el;
            }
            if (starts("<!--")) { skip_until("-->"); continue; }
            if (starts("<![CDATA[")) {
                const std::size_t s = p_ + 9;
                skip_until("]]>");
                if (!have_text) { el->text = b_.substr(s, p_ - 3 - s); have_text = true; }
                continue;
            }
            if (starts("<?")) { skip_until("?>"); continue; }
            el->children.push_back(element(depth + 1));
        }
    }

    std::string_view b_;
    char* w_;
    std::deque<Element>& nodes_;
    std::size_t p_;
    std::size_t base_;
    const std::vector<std::string_view>* raw_;
    std::size_t used_ = 0;  // nodes in use (fragment mode reuses the rest)
};

}  // namespace

Document::Document(ByteView bytes, std::vector<std::string_view> raw_elements)
    : Document(Buffer(bytes.data, bytes.data + bytes.size), std::move(raw_elements)) {}

Document::Document(Buffer&& bytes, std::vector<std::string_view> raw_elements) : buf_(std::move(bytes)) {
    std::string utf16;
    const std::string_view view(reinterpret_cast<const char*>(buf_.data()), buf_.size());
    const std::string_view utf8 = text::as_utf8(view, utf16);
    if (utf8.data() != view.data() + (view.size() - utf8.size()) || utf8.size() != view.size()) {
        if (utf8.data() == utf16.data()) buf_.assign(utf16.begin(), utf16.end());           // transcoded
        else buf_.erase(buf_.begin(), buf_.begin() + std::ptrdiff_t(view.size() - utf8.size()));  // BOM
    }
    char* w = reinterpret_cast<char*>(buf_.data());
    root_ = Parser(std::string_view(w, buf_.size()), w, nodes_, 0, raw_elements.empty() ? nullptr : &raw_elements).parse();
}

void Document::for_each_element(const Element& parent, unsigned threads, std::size_t& pieces, const std::function<void(std::size_t, const Element&)>& fn) {
    const std::string_view c = parent.raw;
    char* const w = reinterpret_cast<char*>(buf_.data()) + (c.data() - reinterpret_cast<const char*>(buf_.data()));
    constexpr std::size_t kMinPiece = 1 << 20;
    const std::size_t n = parent.raw_flat ? std::max<std::size_t>(1, std::min<std::size_t>(thread_count(threads), c.size() / kMinPiece)) : 1;
    std::vector<std::size_t> cut{0};
    for (std::size_t k = 1; k < n; ++k) {
        const std::size_t at = c.find('<', std::max(cut.back(), c.size() * k / n));
        cut.push_back(at == std::string_view::npos ? c.size() : at);
    }
    cut.push_back(c.size());
    pieces = n;
    std::vector<std::function<void()>> tasks;
    for (std::size_t k = 0; k < n; ++k)
        tasks.push_back([&, k] {
            std::deque<Element> nodes;
            Parser(c.substr(cut[k], cut[k + 1] - cut[k]), w + cut[k], nodes, parent.raw_offset + cut[k]).fragment([&](const Element& e) { fn(k, e); });
        });
    run_all(tasks, threads);
}

}  // namespace cadcore::xml
