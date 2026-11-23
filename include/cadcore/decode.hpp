// Public entry point: detect the format of a byte buffer and decode it into a Scene.
#pragma once

#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "cadcore/bytes.hpp"
#include "cadcore/scene.hpp"

namespace cadcore {

// Loads a file referenced by the input (OBJ -> MTL, glTF -> .bin/.png,
// CATProduct -> CATPart). The path is as written in the referencing file.
using ResourceLoader = std::function<std::optional<std::vector<std::uint8_t>>(std::string_view path)>;

// Lists the files a decoder may look through to find a referenced document by content
// (CATIA links name documents by UUID, not by path). Paths are as load_resource accepts them.
using ResourceLister = std::function<std::vector<std::string>()>;

// Finds a document by content id (CATIA: catia_document_id()); returns a path for load_resource.
using DocumentLocator = std::function<std::optional<std::string>(std::string_view document_id)>;

struct DecodeOptions {
    // Threads a decode may use (0 = all cores; never more than 6). With more than one,
    // load_resource, list_resources and locate_document may be called from several threads at
    // once and must be thread-safe.
    unsigned threads = 0;
    ResourceLoader load_resource;      // may be empty: external references then raise MissingResource warnings
    ResourceLister list_resources;     // used to find linked CATIA documents when locate_document is empty
    DocumentLocator locate_document;   // optional, preferred over scanning list_resources
};

struct FormatInfo {
    std::string_view id;           // "stl", "obj", ...
    std::string_view description;
    std::vector<std::string_view> extensions;  // lower case, without dot
};

class Decoder {
public:
    virtual ~Decoder() = default;
    virtual const FormatInfo& info() const = 0;
    // Confidence that `data` is this format: 0 = no, 1 = extension only, 2 = signature match.
    virtual int sniff(ByteView data, std::string_view extension) const = 0;
    virtual Scene decode(ByteView data, const DecodeOptions& options) const = 0;
};

// All decoders compiled into this build.
const std::vector<const Decoder*>& decoders();

// Picks the decoder with the highest sniff score (signature beats extension).
const Decoder* detect(ByteView data, std::string_view filename);

// detect + decode. Throws DecodeError.
Scene decode(ByteView data, std::string_view filename, const DecodeOptions& options = {});

}  // namespace cadcore
