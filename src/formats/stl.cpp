// STL — binary and ASCII.
//
// Binary: 80-byte header, u32le triangle count n, then n * 50-byte records
//   (normal f32x3, 3 vertices f32x3, u16 attribute). A file is binary when its
//   size is exactly 84 + 50 n; ASCII files can also start with "solid", so the
//   size law takes precedence.
// Colour: Materialise writes "COLOR=r g b a" in the header, and a facet attribute
//   with bit 15 clear carries its own RGB555 colour. VisCAM/SolidView set bit 15
//   to mark a valid per-facet RGB555 colour (blue in the low bits).
// ASCII: solid / facet normal / outer loop / vertex x3+ / endloop / endfacet / endsolid.
//   Several solids per file become several meshes; loops with more than three
//   vertices are fan-triangulated.
#include <cstring>

#include <functional>
#include <optional>

#include "core/parallel.hpp"
#include "core/text.hpp"
#include "formats/registry.hpp"

namespace cadcore::formats {
namespace {

bool binary_size_law(ByteView d) {
    if (d.size < 84) return false;
    std::uint32_t n = Reader(d, 80).u32le();
    return std::uint64_t(84) + std::uint64_t(50) * n == d.size;
}

bool looks_ascii(ByteView d) {
    std::string tmp;
    const std::string_view head = text::as_utf8(d.text().substr(0, std::min<std::size_t>(d.size, 8192)), tmp);
    text::Cursor c(head);
    c.set_hash_comments(true);
    c.skip_ws_and_newlines();
    if (!text::iequals_prefix(head.substr(c.position()), "solid")) return false;
    return head.find("facet") != std::string_view::npos || head.find("endsolid") != std::string_view::npos;
}

int sniff(ByteView d, std::string_view ext) {
    if (binary_size_law(d) || looks_ascii(d)) return 2;
    return ext == "stl" ? 1 : 0;
}

Color rgb555(std::uint16_t v, bool blue_low) {
    std::uint8_t lo = std::uint8_t((v & 31) * 255 / 31), mid = std::uint8_t(((v >> 5) & 31) * 255 / 31),
                 hi = std::uint8_t(((v >> 10) & 31) * 255 / 31);
    return blue_low ? Color{hi, mid, lo, 255} : Color{lo, mid, hi, 255};
}

Scene decode_binary(ByteView d, unsigned threads) {
    const std::uint32_t n = Reader(d, 80).u32le();
    const std::uint64_t need = std::uint64_t(84) + std::uint64_t(50) * n;
    Scene scene;
    scene.format = "stl-binary";
    if (need > d.size) fail(ErrorCode::Truncated, "binary STL declares " + std::to_string(n) + " triangles", d.size);
    if (need < d.size) scene.warnings.push_back("binary STL has " + std::to_string(d.size - need) + " trailing bytes");

    const std::string_view header = d.text().substr(0, 80);
    const auto color_at = header.find("COLOR=");
    const bool materialise = color_at != std::string_view::npos && color_at + 10 <= 80;
    Mesh mesh;
    mesh.name = std::string(header.substr(0, header.find('\0')));
    if (materialise) mesh.color = Color{d[color_at + 6], d[color_at + 7], d[color_at + 8], d[color_at + 9]};
    mesh.positions.resize(std::size_t(n) * 9);
    mesh.indices.resize(std::size_t(n) * 3);
    std::vector<std::uint16_t> attrs(n);
    // records: 12 bytes normal (recomputed from the vertices by consumers), 36 bytes vertices,
    // u16 attribute; all n lie inside the data (checked above)
    parallel_for(n, 1 << 16, threads, [&](std::size_t b0, std::size_t e0) {
        for (std::size_t i = b0; i < e0; ++i) {
            const std::uint8_t* rec = d.data + 84 + i * 50;
            std::memcpy(&mesh.positions[i * 9], rec + 12, 36);  // little-endian hosts only (x86, wasm)
            attrs[i] = std::uint16_t(rec[48] | rec[49] << 8);
            const auto v = std::uint32_t(i * 3);
            mesh.indices[i * 3] = v;
            mesh.indices[i * 3 + 1] = v + 1;
            mesh.indices[i * 3 + 2] = v + 2;
        }
    });
    bool any_viscam = false, any_materialise = false;
    for (auto at : attrs) {
        any_viscam |= (at & 0x8000) != 0;
        any_materialise |= materialise && (at & 0x8000) == 0;
    }
    if (any_materialise || (!materialise && any_viscam)) {
        mesh.colors.resize(std::size_t(n) * 12);
        for (std::uint32_t i = 0; i < n; ++i) {
            Color c = mesh.color;
            if (materialise && !(attrs[i] & 0x8000)) c = rgb555(attrs[i], false);
            if (!materialise && (attrs[i] & 0x8000)) c = rgb555(attrs[i], true);
            for (int k = 0; k < 3; ++k) std::memcpy(&mesh.colors[(std::size_t(i) * 3 + k) * 4], &c, 4);
        }
    }
    scene.meshes.push_back(std::move(mesh));
    return scene;
}

// What one range of an ASCII STL says, in order: solid starts and ends, and runs of facets.
struct Op {
    enum Kind { Solid, Facets, End } kind;
    std::string name;
    std::vector<float> positions;
    std::vector<std::uint32_t> indices;  // into this op's positions
    std::size_t skipped = 0;             // facets with fewer than 3 vertices
};

// Parses [begin, end) of the text, which starts and ends at a keyword line.
std::vector<Op> parse_range(std::string_view text, std::size_t begin, std::size_t end) {
    std::vector<Op> ops;
    text::Cursor c(text.substr(begin, end - begin), begin);
    c.set_hash_comments(true);
    float loop[3 * 64];
    auto facets = [&]() -> Op& {
        if (ops.empty() || ops.back().kind != Op::Facets) {
            ops.push_back({Op::Facets, {}, {}, {}, 0});
            const std::size_t facets_guess = (end - begin) / 200;  // ~250 bytes per written facet
            ops.back().positions.reserve(facets_guess * 9);
            ops.back().indices.reserve(facets_guess * 3);
        }
        return ops.back();
    };
    while (true) {
        std::string_view kw = c.token_any();
        if (kw.empty()) break;
        if (kw == "solid") {
            ops.push_back({Op::Solid, std::string(c.rest_of_line()), {}, {}, 0});
        } else if (kw == "facet") {
            c.next_line();
            if (c.token_any() != "outer") fail(ErrorCode::Malformed, "expected 'outer loop'", c.offset());
            c.next_line();
            int nv = 0;
            while (true) {
                std::string_view t = c.token_any();
                if (t == "vertex") {
                    if (nv == 64) fail(ErrorCode::Malformed, "facet loop with more than 64 vertices", c.offset());
                    for (int k = 0; k < 3; ++k) loop[nv * 3 + k] = c.number<float>("vertex coordinate");
                    ++nv;
                } else if (t == "endloop") {
                    break;
                } else {
                    fail(ErrorCode::Malformed, "unexpected token in facet loop", c.offset());
                }
            }
            if (c.token_any() != "endfacet") fail(ErrorCode::Malformed, "expected 'endfacet'", c.offset());
            Op& f = facets();
            if (nv < 3) { ++f.skipped; continue; }
            const auto base = static_cast<std::uint32_t>(f.positions.size() / 3);
            f.positions.insert(f.positions.end(), loop, loop + nv * 3);
            for (int k = 1; k + 1 < nv; ++k) f.indices.insert(f.indices.end(), {base, base + std::uint32_t(k), base + std::uint32_t(k) + 1});
        } else if (kw == "endsolid") {
            c.next_line();
            ops.push_back({Op::End, {}, {}, {}, 0});
        } else {
            fail(ErrorCode::Malformed, "unexpected token '" + std::string(kw.substr(0, 32)) + "'", c.offset());
        }
    }
    return ops;
}

// Start of the first line at or after `pos` whose first token is a top-level keyword.
std::size_t next_keyword_line(std::string_view t, std::size_t pos) {
    while (true) {
        pos = t.find('\n', pos);
        if (pos == std::string_view::npos) return t.size();
        ++pos;
        std::size_t q = pos;
        while (q < t.size() && text::is_space(t[q])) ++q;
        const std::string_view rest = t.substr(q, 8);
        if (rest.substr(0, 5) == "facet" || rest.substr(0, 5) == "solid" || rest.substr(0, 8) == "endsolid") return pos;
    }
}

// Large files are split at keyword lines and the pieces parsed in parallel; the pieces' results
// are then applied in file order, so the scene is the same as a sequential parse.
Scene decode_ascii(ByteView d, unsigned threads) {
    Scene scene;
    scene.format = "stl-ascii";
    std::string utf8;
    const std::string_view text = text::as_utf8(d.text(), utf8);
    constexpr std::size_t kMinChunk = 1 << 20;
    const std::size_t parts = std::max<std::size_t>(1, std::min<std::size_t>(thread_count(threads), text.size() / kMinChunk));
    std::vector<std::size_t> cut{0};
    for (std::size_t k = 1; k < parts; ++k) cut.push_back(std::max(cut.back(), next_keyword_line(text, text.size() * k / parts)));
    cut.push_back(text.size());
    std::vector<std::vector<Op>> results(parts);
    std::vector<std::function<void()>> tasks;
    for (std::size_t k = 0; k < parts; ++k) tasks.push_back([&, k] { results[k] = parse_range(text, cut[k], cut[k + 1]); });
    run_all(tasks, threads);
    // Apply the ops in order: first decide, for every facet run, its mesh and offsets; then
    // size each mesh once and copy all runs into place in parallel.
    struct Place { Op* op; std::size_t mesh, pos_at, idx_at; std::uint32_t base; };
    std::vector<Place> places;
    std::vector<std::pair<std::size_t, std::size_t>> sizes;  // per mesh: floats, indices
    std::optional<std::size_t> mesh;
    for (auto& ops : results)
        for (auto& op : ops) {
            if (op.kind == Op::Solid) {
                scene.meshes.emplace_back();
                scene.meshes.back().name = std::move(op.name);
                sizes.push_back({0, 0});
                mesh = scene.meshes.size() - 1;
            } else if (op.kind == Op::End) {
                mesh.reset();
            } else {
                if (!mesh) {  // tolerate a missing "solid" line
                    scene.meshes.emplace_back();
                    sizes.push_back({0, 0});
                    mesh = scene.meshes.size() - 1;
                }
                auto& [nf, ni] = sizes[*mesh];
                places.push_back({&op, *mesh, nf, ni, std::uint32_t(nf / 3)});
                nf += op.positions.size();
                ni += op.indices.size();
                for (std::size_t k = 0; k < op.skipped; ++k) scene.warnings.push_back("facet with fewer than 3 vertices skipped");
            }
        }
    for (std::size_t m = 0; m < scene.meshes.size(); ++m) {
        scene.meshes[m].positions.resize(sizes[m].first);
        scene.meshes[m].indices.resize(sizes[m].second);
    }
    std::vector<std::function<void()>> copies;
    for (const Place& p : places)
        copies.push_back([&scene, p] {
            Mesh& m = scene.meshes[p.mesh];
            std::copy(p.op->positions.begin(), p.op->positions.end(), m.positions.begin() + std::ptrdiff_t(p.pos_at));
            for (std::size_t k = 0; k < p.op->indices.size(); ++k) m.indices[p.idx_at + k] = p.base + p.op->indices[k];
        });
    run_all(copies, threads);
    if (scene.meshes.empty()) fail(ErrorCode::Malformed, "ASCII STL without solids");
    return scene;
}

Scene decode(ByteView d, const DecodeOptions& o) {
    if (binary_size_law(d)) return decode_binary(d, o.threads);
    if (looks_ascii(d)) return decode_ascii(d, o.threads);
    if (d.size >= 84) return decode_binary(d, o.threads);  // size law broken: binary with trailing data, or truncated
    fail(ErrorCode::Truncated, "too short for a binary STL and not ASCII STL");
}

}  // namespace

const Decoder& stl() {
    static const SimpleDecoder d({"stl", "Stereolithography (binary and ASCII)", {"stl"}}, sniff, decode);
    return d;
}

}  // namespace cadcore::formats
