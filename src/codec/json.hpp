// Strict JSON (RFC 8259) parser producing a small read-only DOM.
//
// Numbers keep both a double and, when the literal is an integer that fits,
// an exact int64. Strings are unescaped to UTF-8 (\uXXXX with surrogate pairs).
// Object members keep their order; duplicate keys keep the first occurrence for
// lookup. Rejects comments, trailing commas, NaN/Infinity, control characters in
// strings, and nesting deeper than 512.
#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cadcore/bytes.hpp"

namespace cadcore::json {

enum class Type : std::uint8_t { Null, Bool, Number, String, Array, Object };

class Value {
public:
    Type type = Type::Null;
    bool boolean = false;
    double number = 0;
    bool is_integer = false;
    std::int64_t integer = 0;
    std::string string;
    std::vector<Value> array;
    std::vector<std::pair<std::string, Value>> object;

    bool is(Type t) const { return type == t; }
    // Member lookup; nullptr when absent or when this is not an object.
    const Value* get(std::string_view key) const;
    std::size_t size() const { return type == Type::Array ? array.size() : type == Type::Object ? object.size() : 0; }
};

Value parse(std::string_view text);  // throws DecodeError(Malformed) with the byte offset

}  // namespace cadcore::json
