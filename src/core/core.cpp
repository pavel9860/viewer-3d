#include "cadcore/traverse.hpp"
#include <algorithm>
#include <cctype>
#include <cstring>

#include "cadcore/bytes.hpp"
#include "cadcore/decode.hpp"
#include "cadcore/error.hpp"
#include "cadcore/scene.hpp"
#include "formats/registry.hpp"

namespace cadcore {

const char* to_string(ErrorCode code) noexcept {
    switch (code) {
        case ErrorCode::UnknownFormat: return "unknown-format";
        case ErrorCode::Truncated: return "truncated";
        case ErrorCode::Malformed: return "malformed";
        case ErrorCode::Unsupported: return "unsupported";
        case ErrorCode::MissingResource: return "missing-resource";
    }
    return "?";
}

DecodeError::DecodeError(ErrorCode code, const std::string& message, std::size_t offset)
    : std::runtime_error(std::string(to_string(code)) + ": " + message +
                         (offset == npos ? std::string() : " (at byte " + std::to_string(offset) + ")")),
      code_(code),
      offset_(offset) {}

void fail(ErrorCode code, const std::string& message, std::size_t offset) { throw DecodeError(code, message, offset); }

ByteView ByteView::sub(std::size_t pos, std::size_t len) const {
    if (pos > size || len > size - pos) fail(ErrorCode::Truncated, "sub-range outside the buffer", pos);
    return {data + pos, len};
}

std::size_t ByteView::find(std::string_view needle, std::size_t from, std::size_t to) const {
    if (to > size) to = size;
    if (needle.empty() || from >= to || needle.size() > to - from) return npos;
    const std::size_t at = find_bytes(data + from, to - from, reinterpret_cast<const std::uint8_t*>(needle.data()), needle.size());
    return at == to - from ? npos : from + at;
}

std::size_t find_bytes(const std::uint8_t* hay, std::size_t n, const std::uint8_t* needle, std::size_t m) {
    if (m == 0) return 0;
    if (m > n) return n;
    auto rarity = [](std::uint8_t b) { return b == 0 || b == 0xff ? 0 : b < 8 ? 1 : 2; };
    std::size_t k = 0;
    for (std::size_t i = 1; i < m; ++i)
        if (rarity(needle[i]) > rarity(needle[k])) k = i;
    const std::uint8_t* p = hay + k;
    const std::uint8_t* const last = hay + (n - m) + k;  // last position the anchor byte may take
    while (p <= last) {
        p = static_cast<const std::uint8_t*>(std::memchr(p, needle[k], std::size_t(last - p) + 1));
        if (!p) return n;
        if (std::memcmp(p - k, needle, m) == 0) return std::size_t(p - k - hay);
        ++p;
    }
    return n;
}

std::size_t Scene::triangle_count() const {
    std::size_t n = 0;
    for (const auto& m : meshes) n += m.triangle_count();
    return n;
}

std::size_t Scene::vertex_count() const {
    std::size_t n = 0;
    for (const auto& m : meshes) n += m.vertex_count();
    return n;
}

void Scene::make_flat_root(const std::string& name) {
    Node root;
    root.name = name;
    for (std::uint32_t i = 0; i < meshes.size(); ++i) root.meshes.push_back(i);
    for (std::uint32_t i = 0; i < polylines.size(); ++i) root.polylines.push_back(i);
    nodes.clear();
    nodes.push_back(std::move(root));
}

const std::vector<const Decoder*>& decoders() {
    static const std::vector<const Decoder*> all = formats::all_decoders();
    return all;
}

static std::string extension_of(std::string_view filename) {
    auto dot = filename.rfind('.');
    auto slash = filename.find_last_of("/\\");
    if (dot == std::string_view::npos || (slash != std::string_view::npos && dot < slash)) return {};
    std::string ext(filename.substr(dot + 1));
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return std::tolower(c); });
    return ext;
}

const Decoder* detect(ByteView data, std::string_view filename) {
    const std::string ext = extension_of(filename);
    const Decoder* best = nullptr;
    int best_score = 0;
    for (const Decoder* d : decoders()) {
        int s = d->sniff(data, ext);
        if (s > best_score) {
            best = d;
            best_score = s;
        }
    }
    return best;
}

Scene decode(ByteView data, std::string_view filename, const DecodeOptions& options) {
    if (data.size == 0) fail(ErrorCode::Malformed, "empty file (0 bytes)");
    const Decoder* d = detect(data, filename);
    if (!d) fail(ErrorCode::UnknownFormat, "no decoder recognises '" + std::string(filename) + "'");
    Scene s = d->decode(data, options);
    if (s.nodes.empty()) s.make_flat_root();
    return s;
}

Matrix4 multiply(const Matrix4& a, const Matrix4& b) {
    Matrix4 r{};
    for (int c = 0; c < 4; ++c)
        for (int row = 0; row < 4; ++row) {
            double s = 0;
            for (int k = 0; k < 4; ++k) s += a[k * 4 + row] * b[c * 4 + k];
            r[c * 4 + row] = s;
        }
    return r;
}

std::array<double, 3> apply(const Matrix4& t, const float* p) {
    std::array<double, 3> r{};
    for (std::size_t k = 0; k < 3; ++k) r[k] = t[k] * p[0] + t[4 + k] * p[1] + t[8 + k] * p[2] + t[12 + k];
    return r;
}

void for_each_node(const Scene& scene, const std::function<void(const Node&, const Matrix4&)>& visit) {
    if (scene.nodes.empty()) {
        Node root;
        for (std::uint32_t i = 0; i < scene.meshes.size(); ++i) root.meshes.push_back(i);
        for (std::uint32_t i = 0; i < scene.polylines.size(); ++i) root.polylines.push_back(i);
        visit(root, kIdentity);
        return;
    }
    struct Item { std::uint32_t node; Matrix4 world; int depth; };
    std::vector<Item> stack{{0, scene.nodes[0].transform, 0}};
    while (!stack.empty()) {
        const Item it = stack.back();
        stack.pop_back();
        if (it.depth > 512) fail(ErrorCode::Malformed, "node tree deeper than 512 (cycle?)");
        const Node& n = scene.nodes.at(it.node);
        visit(n, it.world);
        for (auto c : n.children) stack.push_back({c, multiply(it.world, scene.nodes.at(c).transform), it.depth + 1});
    }
}

}  // namespace cadcore
