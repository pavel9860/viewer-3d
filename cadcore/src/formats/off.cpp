// OFF (Object File Format, Geomview) — ASCII and "OFF BINARY".
//
// Header keyword [ST][C][N][4]OFF (the prefixes announce per-vertex texture,
// colour, normal and homogeneous coordinates), then "nv nf ne" (may share the
// keyword's line). Vertices: x y z [w] [nx ny nz] [colour] [s t].
// Faces: n i0 .. in-1 [colour], colour = r g b [a] as 0..1 floats or 0..255 ints.
// Binary: big-endian i32 counts, f32 vertex data, faces as i32 n, n * i32, i32 k, k * f32.
#include <string>
#include <vector>

#include "core/parallel.hpp"
#include "core/text.hpp"
#include "formats/registry.hpp"

namespace cadcore::formats {
namespace {

struct Header {
    bool st = false, c = false, n = false, four = false, binary = false;
    std::size_t keyword_len = 0;
};

bool parse_keyword(std::string_view s, Header& h) {
    std::size_t i = 0;
    if (s.substr(i, 2) == "ST") { h.st = true; i += 2; }
    if (s.substr(i, 1) == "C") { h.c = true; ++i; }
    if (s.substr(i, 1) == "N") { h.n = true; ++i; }
    if (s.substr(i, 1) == "4") { h.four = true; ++i; }
    if (s.substr(i, 3) != "OFF") return false;
    i += 3;
    if (i != s.size()) return false;
    h.keyword_len = i;
    return true;
}

int sniff(ByteView d, std::string_view ext) {
    text::Cursor c(d.text().substr(0, std::min<std::size_t>(d.size, 256)));
    c.skip_ws_and_newlines();
    Header h;
    if (parse_keyword(c.token(), h)) return 2;
    return ext == "off" ? 1 : 0;
}

Color color_from(const double* v, int n) {
    bool unit = true;
    for (int i = 0; i < n; ++i) unit &= v[i] <= 1.0;
    auto cv = [&](double x) { return std::uint8_t(unit ? (x <= 0 ? 0 : x >= 1 ? 255 : x * 255 + 0.5) : (x < 0 ? 0 : x > 255 ? 255 : x)); };
    return {cv(v[0]), cv(v[1]), cv(v[2]), n >= 4 ? cv(v[3]) : std::uint8_t(255)};
}

void add_face_color(Mesh& m, std::uint32_t tri_before, Color c) {
    std::uint32_t tri_after = std::uint32_t(m.triangle_count());
    if (!m.groups.empty() && m.groups.back().color == c && m.groups.back().first + m.groups.back().count == tri_before)
        m.groups.back().count += tri_after - tri_before;
    else
        m.groups.push_back({tri_before, tri_after - tri_before, c});
}

Scene decode_binary(ByteView d, std::size_t pos, const Header& h) {
    Reader r(d, pos);
    const std::uint32_t nv = r.u32be(), nf = r.u32be();
    r.u32be();  // edges
    if (std::uint64_t(nv) * 12 > d.size) fail(ErrorCode::Malformed, "vertex count exceeds file size");
    Mesh m;
    const int dim = h.four ? 4 : 3;
    for (std::uint32_t i = 0; i < nv; ++i) {
        float v[4] = {0, 0, 0, 1};
        for (int k = 0; k < dim; ++k) v[k] = r.f32be();
        if (h.four && v[3] != 0) { v[0] /= v[3]; v[1] /= v[3]; v[2] /= v[3]; }
        m.positions.insert(m.positions.end(), {v[0], v[1], v[2]});
        if (h.n) r.skip(12);
        if (h.c) r.skip(16);
        if (h.st) r.skip(8);
    }
    for (std::uint32_t f = 0; f < nf; ++f) {
        std::uint32_t n = r.u32be();
        if (n > nv) fail(ErrorCode::Malformed, "face larger than vertex count", r.pos());
        std::vector<std::uint32_t> idx(n);
        for (auto& i : idx) {
            i = r.u32be();
            if (i >= nv) fail(ErrorCode::Malformed, "face index out of range", r.pos());
        }
        std::uint32_t before = std::uint32_t(m.triangle_count());
        for (std::uint32_t k = 1; k + 1 < n; ++k) m.indices.insert(m.indices.end(), {idx[0], idx[k], idx[k + 1]});
        std::uint32_t nc = r.u32be();
        if (nc > 4) fail(ErrorCode::Malformed, "face colour with more than 4 components", r.pos());
        double col[4];
        for (std::uint32_t k = 0; k < nc; ++k) col[k] = r.f32be();
        if (nc >= 3) add_face_color(m, before, color_from(col, int(nc)));
    }
    Scene s;
    s.format = "off-binary";
    s.meshes.push_back(std::move(m));
    return s;
}

Scene decode(ByteView d, const DecodeOptions& options) {
    std::string utf8;
    text::Cursor c(text::as_utf8(d.text(), utf8));
    c.set_hash_comments(true);
    c.skip_ws_and_newlines();
    Header h;
    if (!parse_keyword(c.token(), h)) fail(ErrorCode::Malformed, "missing OFF keyword", c.offset());
    std::string_view after = c.token();
    if (after == "BINARY") {
        c.next_line();
        return decode_binary(d, c.position(), h);
    }
    if (!after.empty()) c.set_position(c.position() - after.size());
    auto skip_comments = [](text::Cursor& k) {
        k.skip_ws_and_newlines();
        while (k.peek() == '#') { k.next_line(); k.skip_ws_and_newlines(); }
    };
    skip_comments(c);
    const long long nv = c.integer_any("vertex count");
    const long long nf = c.integer_any("face count");
    long long ne = 0;
    c.try_int(ne);
    if (nv < 0 || nf < 0 || std::uint64_t(nv) > d.size) fail(ErrorCode::Malformed, "bad OFF counts");
    Mesh m;
    const std::string_view text = c.text();
    std::vector<std::size_t> lines;
    const unsigned threads = options.threads;
    auto cursor_at = [&](std::size_t at) {
        text::Cursor k(text);
        k.set_hash_comments(true);
        k.set_position(at);
        return k;
    };
    // vertex records -> positions (+ colours)
    struct Verts { std::vector<float> pos; std::vector<std::uint8_t> col; };
    auto vertices = [&](text::Cursor& k, std::size_t n, Verts& o) {
        o.pos.reserve(n * 3);
        for (std::size_t i = 0; i < n; ++i) {
            skip_comments(k);
            double v[4] = {0, 0, 0, 1};
            for (int j = 0; j < (h.four ? 4 : 3); ++j) v[j] = k.number<double>("vertex coordinate");
            if (h.four && v[3] != 0) { v[0] /= v[3]; v[1] /= v[3]; v[2] /= v[3]; }
            o.pos.insert(o.pos.end(), {float(v[0]), float(v[1]), float(v[2])});
            double rest[16];
            int nr = 0;
            while (nr < 16 && k.try_number(rest[nr])) ++nr;
            if (h.c) {
                const int off = h.n ? 3 : 0;
                const int ncol = nr - off - (h.st ? 2 : 0);
                const Color col = ncol >= 3 ? color_from(rest + off, ncol >= 4 ? 4 : 3) : Color{200, 200, 200, 255};
                o.col.insert(o.col.end(), {col.r, col.g, col.b, col.a});
            }
            k.next_line();
        }
    };
    // face records -> triangles, and per face (triangles, colour, has colour)
    struct FaceCol { std::uint32_t triangles; Color color; bool colored; };
    struct Faces { std::vector<std::uint32_t> idx; std::vector<FaceCol> faces; };
    auto face_records = [&](text::Cursor& k, std::size_t n, Faces& o) {
        std::vector<std::uint32_t> idx;
        o.idx.reserve(n * 3);
        o.faces.reserve(n);
        for (std::size_t f = 0; f < n; ++f) {
            skip_comments(k);
            const long long sz = k.integer("face size");
            if (sz < 0 || sz > nv) fail(ErrorCode::Malformed, "bad face size", k.offset());
            idx.resize(std::size_t(sz));
            for (auto& i : idx) {
                const long long v = k.integer("face index");
                if (v < 0 || v >= nv) fail(ErrorCode::Malformed, "face index out of range", k.offset());
                i = std::uint32_t(v);
            }
            const std::size_t before = o.idx.size();
            for (std::size_t j = 1; j + 1 < idx.size(); ++j) o.idx.insert(o.idx.end(), {idx[0], idx[j], idx[j + 1]});
            double col[4];
            int ncol = 0;
            while (ncol < 4 && k.try_number(col[ncol])) ++ncol;
            o.faces.push_back({std::uint32_t((o.idx.size() - before) / 3), ncol >= 3 ? color_from(col, ncol) : Color{200, 200, 200, 255}, ncol >= 3});
            k.next_line();
        }
    };
    // Records one per line are parsed in parallel pieces; anything else sequentially.
    auto run = [&](std::size_t count, auto& outs, auto&& parse_part) {
        std::size_t parts = 1;
        std::size_t end = std::string_view::npos;
        if (count >= 8192) {
            if (lines.empty()) lines = text::line_starts(text, 0, threads);
            outs.assign(kMaxThreads, {});
            end = text::parallel_line_records(lines, c.position(), count, text.size(), threads, parts,
                                              [&](std::size_t k, std::size_t at, std::size_t, std::size_t n) {
                                                  auto kc = cursor_at(at);
                                                  outs[k].value = {};
                                                  parse_part(kc, n, outs[k].value);
                                                  return kc.position();
                                              });
        }
        if (end == std::string_view::npos) {
            parts = 1;
            outs.assign(1, {});
            parse_part(c, count, outs[0].value);
        } else {
            c.set_position(end);
        }
        return parts;
    };
    if (c.at_eol()) c.next_line();  // records start on the next line (the sequential parse skips this newline too)
    std::vector<Slot<Verts>> vo;
    const std::size_t vparts = run(std::size_t(nv), vo, vertices);
    append_all(m.positions, vparts, [&](std::size_t k) -> auto& { return vo[k].value.pos; }, threads);
    append_all(m.colors, vparts, [&](std::size_t k) -> auto& { return vo[k].value.col; }, threads);
    std::vector<Slot<Faces>> fo;
    const std::size_t fparts = run(std::size_t(nf), fo, face_records);
    append_all(m.indices, fparts, [&](std::size_t k) -> auto& { return fo[k].value.idx; }, threads);
    std::vector<std::pair<std::uint32_t, Color>> face_cols;  // (first triangle, colour) per face
    bool any_color = false;
    std::uint32_t tri = 0;
    for (std::size_t k = 0; k < fparts; ++k)
        for (const FaceCol& f : fo[k].value.faces) {
            face_cols.push_back({tri, f.color});
            any_color |= f.colored;
            tri += f.triangles;
        }
    if (any_color) {
        for (std::size_t f = 0; f < face_cols.size(); ++f) {
            std::uint32_t end = f + 1 < face_cols.size() ? face_cols[f + 1].first : std::uint32_t(m.triangle_count());
            Color col = face_cols[f].second;
            if (!m.groups.empty() && m.groups.back().color == col)
                m.groups.back().count += end - face_cols[f].first;
            else
                m.groups.push_back({face_cols[f].first, end - face_cols[f].first, col});
        }
    }
    if (m.groups.size() == 1 && m.groups[0].count == m.triangle_count()) {
        m.color = m.groups[0].color;
        m.groups.clear();
    }
    Scene s;
    s.format = "off-ascii";
    s.meshes.push_back(std::move(m));
    return s;
}

}  // namespace

const Decoder& off() {
    static const SimpleDecoder d({"off", "Object File Format (Geomview)", {"off"}}, sniff, decode);
    return d;
}

}  // namespace cadcore::formats
