// Error model shared by every decoder.
#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>

namespace cadcore {

enum class ErrorCode {
    UnknownFormat,   // no decoder recognises the input
    Truncated,       // the input ends inside a structure
    Malformed,       // a structure violates its grammar
    Unsupported,     // valid input using a feature this build does not decode
    MissingResource  // a referenced external file (MTL, .bin, CATPart) is unavailable
};

const char* to_string(ErrorCode code) noexcept;

class DecodeError : public std::runtime_error {
public:
    DecodeError(ErrorCode code, const std::string& message, std::size_t offset = npos);

    ErrorCode code() const noexcept { return code_; }
    // Byte offset where decoding stopped, or npos when not applicable.
    std::size_t offset() const noexcept { return offset_; }
    // The message without the code prefix and offset suffix of what(), for wrapping with context.
    const std::string& message() const noexcept { return message_; }

    static constexpr std::size_t npos = static_cast<std::size_t>(-1);

private:
    ErrorCode code_;
    std::size_t offset_;
    std::string message_;
};

[[noreturn]] void fail(ErrorCode code, const std::string& message, std::size_t offset = DecodeError::npos);

}  // namespace cadcore
