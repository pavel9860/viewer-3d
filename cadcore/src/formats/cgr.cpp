// CATIA CGR (.cgr): a V5_CFV2 store holding only the display data of a part or product.
//
// The catalog lists, in order, the containers of the document (MAIN, Describe, Header,
// SceneGraph), each followed by its own MainDataStream. The tessellation is the scene graph:
// the MainDataStream listed right after the "SceneGraph" container, decoded with the CFV2
// grammar of catia/tessellation.hpp (byte references).
#include <cstring>

#include "cadcore/error.hpp"
#include "formats/catia/store.hpp"
#include "formats/catia/tessellation.hpp"
#include "formats/registry.hpp"

namespace cadcore::formats {
namespace {

int sniff(ByteView d, std::string_view ext) {
    const bool magic = d.size >= 8 && std::memcmp(d.data, "V5_CFV2", 8) == 0;
    if (magic && ext == "cgr") return 2;
    return ext == "cgr" ? 1 : 0;
}

Scene decode(ByteView d, const DecodeOptions&) {
    const catia::Store store = catia::read_store(d, 0);
    if (store.version != 2) fail(ErrorCode::Unsupported, "CGR: store version V5_CFV" + std::to_string(store.version));
    const catia::StreamInfo* graph = nullptr;
    for (std::size_t i = 1; i < store.streams.size(); ++i)
        if (store.streams[i].name == "MainDataStream" && store.streams[i - 1].name == "SceneGraph") {
            if (graph) fail(ErrorCode::Malformed, "CGR: more than one SceneGraph data stream");
            graph = &store.streams[i];
        }
    if (!graph) fail(ErrorCode::Malformed, "CGR: no SceneGraph data stream in the catalog");
    const auto stream = catia::read_stream(d, store, *graph);
    auto bodies = catia::decode_tessellation(ByteView(stream.data(), stream.size()), catia::Generation::Cfv2);
    Scene scene;
    scene.format = "cgr";
    scene.nodes.push_back(Node{});
    if (bodies.empty()) scene.warnings.push_back("CGR holds no tessellation");
    catia::add_to_scene(std::move(bodies), "cgr", scene, 0);
    return scene;
}

}  // namespace

const Decoder& cgr() {
    static const SimpleDecoder d({"cgr", "CATIA CGR (display tessellation)", {"cgr"}}, sniff, decode);
    return d;
}

}  // namespace cadcore::formats
