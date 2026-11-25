// Non-validating XML 1.0 parser (W3C XML 1.0, 5th ed.) building a compact read-only DOM.
//
// Handles: prolog / processing instructions, comments, CDATA, DOCTYPE (skipped,
// including an internal subset), elements, attributes (single or double quotes),
// predefined entities and character references (&#n; &#xh;). Text and attribute
// values are decoded in place into the document's own copy, so names and values
// are string_views valid for the lifetime of the Document. UTF-8 input, or UTF-16
// with a BOM (transcoded). Namespace prefixes are kept; local() strips them.
#pragma once

#include <deque>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

#include "cadcore/bytes.hpp"
#include "core/buffer.hpp"

namespace cadcore::xml {

struct Attribute {
    std::string_view name, value;
};

struct Element {
    std::string_view name;
    std::vector<Attribute> attributes;
    std::vector<Element*> children;
    std::string_view text;  // concatenation of direct text/CDATA content when it is a single run, else the first run
    // For elements the Document was told to keep raw: their content, undecoded and unparsed
    // (children stays empty); read it with for_each_element. raw_offset: its byte offset.
    std::string_view raw;
    std::size_t raw_offset = 0;
    bool is_raw = false;
    bool raw_flat = false;  // raw content is only self-closing elements, text and whitespace

    // Attribute value by exact or local (prefix-stripped) name; empty view when absent.
    std::string_view attr(std::string_view n) const;
    bool has_attr(std::string_view n) const;
    // First child element with the given local name, or nullptr.
    const Element* child(std::string_view local_name) const;
    std::string_view local() const;
};

std::string_view local_name(std::string_view qualified);

class Document {
public:
    // Throws DecodeError(Malformed) with the byte offset of the first syntax error. Elements whose
    // local name is in `raw_elements` keep their content raw (Element::raw), for large regular
    // content such as 3MF vertex and triangle lists. The document parses in place, in its own
    // buffer (the moved-in bytes, or a copy).
    explicit Document(ByteView bytes, std::vector<std::string_view> raw_elements = {});
    explicit Document(Buffer&& bytes, std::vector<std::string_view> raw_elements = {});
    const Element& root() const { return *root_; }

    // Parses the raw content of `raw_parent` (an element of this document) with the same grammar
    // and calls fn for each child element in document order (valid only during the call). The
    // content is decoded in place, so each raw element can be read once. Content made only of
    // self-closing children is split into up to `threads` pieces parsed in parallel; fn is then
    // called concurrently, with `piece` = 0.. in content order; `pieces` returns their number.
    void for_each_element(const Element& raw_parent, unsigned threads, std::size_t& pieces,
                          const std::function<void(std::size_t piece, const Element&)>& fn);

    Document(const Document&) = delete;
    Document& operator=(const Document&) = delete;

private:
    Buffer buf_;
    std::deque<Element> nodes_;
    Element* root_ = nullptr;
};

}  // namespace cadcore::xml
