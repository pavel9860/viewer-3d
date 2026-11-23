// CATIA V5 .CATPart — the display tessellation stored in the file.
//
// Port of the Python reference decoder (catia/catv5); the byte grammar is
// documented in catia/README.md. Summary:
//   nested store   "V5_CFV2\0" + u32be catalog offset + u32be catalog length + ff*8;
//                  body = [origin + 32, origin + catalog offset)
//   point list     node "04 04 fd ff", token "ff ff 02 00 01" within 64 bytes, then
//                  P1g count, P1g edge bound, u8 has_tex, u8, u8, [u8 dims, f32*count*dims],
//                  u8 encoding (0 = 2-bit prefix table, 1 = raw xyz)
//   face list      f32 sag, u8 nLOD, nLOD x {record, f32, attribute}, u8 kind, "02 0a",
//                  P1g nFaces, nFaces facet records back to back, nFaces attributes
//                  "<flags> 04 04 ff A B G R" (flags 0x10 solid, 0x80 transparent)
//   facet record   01 <0x40|tri 1|strip 2|fan 4|normal 8> [P1g per kind] ff u32le n
//                  [f32x3 normal] [P1g per strip] [P1g per fan] n refs (u8 | u16be)
#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "core/parallel.hpp"
#include "formats/registry.hpp"

namespace cadcore::formats {
namespace {

constexpr char kMagic[] = "V5_CFV2";  // followed by \0
constexpr std::uint8_t kNode[] = {0x04, 0x04, 0xfd, 0xff};
constexpr std::uint8_t kToken[] = {0xff, 0xff, 0x02, 0x00, 0x01};
constexpr std::uint32_t kMaxElements = 50'000'000;
constexpr std::uint32_t kMaxRefs = 50'000'000;

int sniff(ByteView d, std::string_view ext) {
    if (d.size >= 8 && std::memcmp(d.data, kMagic, 8) == 0 && (ext == "catpart" || ext.empty())) return 2;
    if (d.size >= 8 && std::memcmp(d.data, kMagic, 8) == 0) return ext == "catpart" ? 2 : 1;
    return 0;
}

// --- low level --------------------------------------------------------------------------------

struct Span {
    const std::uint8_t* p;
    std::size_t n;
    std::uint8_t at(std::size_t i) const { return p[i]; }
};

bool read_count(Span s, std::size_t& pos, std::size_t limit, std::uint32_t& out) {
    if (pos >= limit) return false;
    std::uint8_t v = s.p[pos];
    if (v != 0xff) { out = v; ++pos; return true; }
    if (pos + 5 > limit) return false;
    std::memcpy(&out, s.p + pos + 1, 4);  // little-endian hosts only (x86, wasm)
    pos += 5;
    return true;
}
float f32(Span s, std::size_t pos) { float f; std::memcpy(&f, s.p + pos, 4); return f; }
std::uint32_t u32le(Span s, std::size_t pos) { std::uint32_t v; std::memcpy(&v, s.p + pos, 4); return v; }
std::uint32_t u32be(Span s, std::size_t pos) {
    return (std::uint32_t(s.p[pos]) << 24) | (std::uint32_t(s.p[pos + 1]) << 16) | (std::uint32_t(s.p[pos + 2]) << 8) | s.p[pos + 3];
}
std::size_t find(Span s, const std::uint8_t* needle, std::size_t nlen, std::size_t from, std::size_t to) {
    if (to > s.n) to = s.n;
    if (from >= to || nlen > to - from) return std::string::npos;
    const std::size_t at = find_bytes(s.p + from, to - from, needle, nlen);
    return at == to - from ? std::string::npos : from + at;
}

// --- containers -------------------------------------------------------------------------------

struct Store { std::size_t origin, lo, hi; };

std::vector<Store> nested_stores(Span s) {
    std::vector<Store> out;
    std::size_t pos = 8;
    static const std::uint8_t magic[8] = {'V', '5', '_', 'C', 'F', 'V', '2', 0};
    while ((pos = find(s, magic, 8, pos, s.n)) != std::string::npos) {
        if (pos + 32 <= s.n && std::all_of(s.p + pos + 16, s.p + pos + 24, [](std::uint8_t b) { return b == 0xff; })) {
            std::size_t off = u32be(s, pos + 8), len = u32be(s, pos + 12);
            if (off > 32 && pos + off + len <= s.n) out.push_back({pos, pos + 32, pos + off});
        }
        pos += 8;
    }
    return out;
}

// --- point lists ------------------------------------------------------------------------------

struct PointList {
    std::size_t node, end;
    std::uint32_t count;
    std::vector<float> xyz;
};

std::optional<PointList> decode_points(Span s, std::size_t node, std::size_t limit) {
    std::size_t tok = find(s, kToken, 5, node + 4, std::min(node + 4 + 64, limit));
    if (tok == std::string::npos) return std::nullopt;
    std::size_t p = tok + 5;
    std::uint32_t ec, eb;
    if (!read_count(s, p, limit, ec) || !read_count(s, p, limit, eb)) return std::nullopt;
    if (ec == 0 || ec > kMaxElements || p + 4 > limit) return std::nullopt;
    const std::uint8_t has_tex = s.at(p);
    if (has_tex > 1) return std::nullopt;
    p += 3;
    if (has_tex) {
        std::uint32_t dims = s.at(p++);
        if (dims < 1 || dims > 4 || p + std::size_t(4) * dims * ec > limit) return std::nullopt;
        p += std::size_t(4) * dims * ec;
    }
    if (p >= limit) return std::nullopt;
    const std::uint8_t enc = s.at(p++);
    PointList pl{node, 0, ec, {}};
    pl.xyz.resize(std::size_t(ec) * 3);
    if (enc == 0) {
        const std::size_t tl = (ec + 3) / 4;
        if (p + tl >= limit) return std::nullopt;
        const std::size_t table = p;
        p += tl;
        std::uint32_t words;
        if (!read_count(s, p, limit, words)) return std::nullopt;
        static constexpr int kWords[4] = {3, 0, 1, 2};
        std::uint64_t need = 0;
        for (std::uint32_t i = 0; i < ec; ++i) need += kWords[(s.at(table + (i >> 2)) >> ((i & 3) * 2)) & 3];
        if (need != words || p + std::size_t(4) * words > limit) return std::nullopt;
        float x = 0, y = 0, z = 0;
        std::size_t q = p;
        for (std::uint32_t i = 0; i < ec; ++i) {
            switch ((s.at(table + (i >> 2)) >> ((i & 3) * 2)) & 3) {
                case 0: x = f32(s, q); y = f32(s, q + 4); z = f32(s, q + 8); q += 12; break;
                case 3: y = f32(s, q); z = f32(s, q + 4); q += 8; break;
                case 2: z = f32(s, q); q += 4; break;
                default: break;  // opcode 1: same position as the previous element
            }
            pl.xyz[i * 3] = x; pl.xyz[i * 3 + 1] = y; pl.xyz[i * 3 + 2] = z;
        }
        pl.end = q;
    } else if (enc == 1) {
        if (p + std::size_t(12) * ec > limit) return std::nullopt;
        std::memcpy(pl.xyz.data(), s.p + p, std::size_t(12) * ec);
        pl.end = p + std::size_t(12) * ec;
    } else {
        return std::nullopt;
    }
    for (float v : pl.xyz)
        if (!std::isfinite(v) || std::fabs(v) >= 1e7f) return std::nullopt;
    return pl;
}

// --- facet records ----------------------------------------------------------------------------

struct Header {
    std::size_t offset, refs_at;
    std::uint8_t marker;
    std::uint32_t triangles, n;
    std::vector<std::uint32_t> strips, fans;
};

// Loop counts inside a record: single P1g bytes (V5R25+ files) or u16be (older releases, e.g. V5R17).
// The writer's choice is not recorded; find_face_list resolves it by exact parsing of the whole list.
enum class CountEnc { P1g, U16BE };

bool read_loop_count(Span s, std::size_t& p, std::size_t limit, CountEnc enc, std::uint32_t& out) {
    if (enc == CountEnc::P1g) return read_count(s, p, limit, out);
    if (p + 2 > limit) return false;
    out = (std::uint32_t(s.at(p)) << 8) | s.at(p + 1);
    p += 2;
    return true;
}

std::optional<Header> parse_header(Span s, std::size_t pos, std::size_t limit, CountEnc enc) {
    if (pos + 7 > limit || s.at(pos) != 0x01 || (s.at(pos + 1) & 0xf0) != 0x40) return std::nullopt;
    Header h{pos, 0, s.at(pos + 1), 0, 0, {}, {}};
    std::size_t p = pos + 2;
    std::uint32_t c[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k)
        if (h.marker & (1 << k))
            if (!read_count(s, p, limit, c[k])) return std::nullopt;
    if (p + 5 > limit || s.at(p) != 0xff) return std::nullopt;
    h.n = u32le(s, p + 1);
    p += 5;
    if (h.n > kMaxRefs || c[1] > h.n || c[2] > h.n) return std::nullopt;
    if (h.marker & 0x08) {
        if (p + 12 > limit) return std::nullopt;
        p += 12;
    }
    h.triangles = c[0];
    std::uint64_t sum = std::uint64_t(3) * c[0];
    for (int k = 1; k <= 2; ++k) {
        auto& v = k == 1 ? h.strips : h.fans;
        v.resize(c[k]);
        for (auto& x : v) {
            if (!read_loop_count(s, p, limit, enc, x)) return std::nullopt;
            sum += x;
        }
    }
    if (sum != h.n) return std::nullopt;
    h.refs_at = p;
    return h;
}

struct Record {
    Header h;
    std::vector<std::uint32_t> refs;
    std::size_t end;
};

std::optional<Record> read_record(Span s, const Header& h, int width, std::size_t limit, std::uint32_t ec) {
    const std::size_t end = h.refs_at + std::size_t(width) * h.n;
    if (end > limit) return std::nullopt;
    Record r{h, std::vector<std::uint32_t>(h.n), end};
    // Big-endian references of 1..4 bytes: u8, u16be, u24be (element spaces above 65535), u32be;
    // all inside [refs_at, end), checked above.
    const std::uint8_t* p = s.p + h.refs_at;
    std::uint32_t* out = r.refs.data();
    switch (width) {
        case 1: for (std::uint32_t i = 0; i < h.n; ++i) out[i] = p[i]; break;
        case 2: for (std::uint32_t i = 0; i < h.n; ++i) out[i] = std::uint32_t(p[2 * i]) << 8 | p[2 * i + 1]; break;
        default:
            for (std::uint32_t i = 0; i < h.n; ++i) {
                std::uint32_t v = 0;
                for (int k = 0; k < width; ++k) v = (v << 8) | p[std::size_t(width) * i + std::size_t(k)];
                out[i] = v;
            }
    }
    for (std::uint32_t i = 0; i < h.n; ++i)
        if (out[i] >= ec) return std::nullopt;
    return r;
}

bool is_attr(Span s, std::size_t pos, std::size_t limit) {
    return pos + 8 <= limit && s.at(pos + 1) == 0x04 && s.at(pos + 2) == 0x04 && s.at(pos + 3) == 0xff;
}

void triangulate(const Record& r, std::vector<std::uint32_t>& out, std::uint32_t base) {
    const auto& f = r.refs;
    std::size_t k = 0;
    auto emit = [&](std::uint32_t a, std::uint32_t b, std::uint32_t c) {
        if (a != b && b != c && a != c) out.insert(out.end(), {a + base, b + base, c + base});
    };
    for (std::uint32_t t = 0; t < r.h.triangles; ++t, k += 3) emit(f[k], f[k + 1], f[k + 2]);
    for (auto c : r.h.strips) {
        for (std::uint32_t t = 0; t + 2 < c; ++t)
            t % 2 == 0 ? emit(f[k + t], f[k + t + 1], f[k + t + 2]) : emit(f[k + t + 1], f[k + t], f[k + t + 2]);
        k += c;
    }
    for (auto c : r.h.fans) {
        for (std::uint32_t t = 1; t + 1 < c; ++t) emit(f[k], f[k + t], f[k + t + 1]);
        k += c;
    }
}

// --- face lists -------------------------------------------------------------------------------

struct Face {
    Record rec;
    std::uint8_t flags;
    Color color;
};

std::optional<std::vector<Record>> parse_records(Span s, std::size_t pos, std::size_t limit, std::uint32_t ec, std::uint32_t count, CountEnc enc,
                                                  std::size_t& end) {
    std::vector<Record> out;
    out.reserve(count);
    std::optional<Header> h = parse_header(s, pos, limit, enc);
    for (std::uint32_t k = 0; k < count; ++k) {
        if (!h) return std::nullopt;
        const bool last = k + 1 == count;
        std::optional<Record> chosen;
        std::optional<Header> chosen_next;
        int options = 0;
        // Narrowest width able to address the element space, up to 4 bytes; a record may also use a
        // wider width. A width is an option when what follows it is valid (checked first: cheap) and
        // all its references are inside the element space.
        const int min_width = ec <= 0x100 ? 1 : ec <= 0x10000 ? 2 : ec <= 0x1000000 ? 3 : 4;
        for (int w = h->n == 0 ? 1 : min_width; w <= (h->n == 0 ? 1 : 4); ++w) {
            const std::size_t rec_end = h->refs_at + std::size_t(w) * h->n;
            if (rec_end > limit) continue;
            std::optional<Header> next;
            if (last ? !is_attr(s, rec_end, limit) : !(next = parse_header(s, rec_end, limit, enc))) continue;
            auto r = read_record(s, *h, w, limit, ec);
            if (!r) continue;
            chosen = std::move(r);
            chosen_next = std::move(next);
            ++options;
        }
        if (options != 1) return std::nullopt;
        pos = chosen->end;
        out.push_back(std::move(*chosen));
        h = std::move(chosen_next);
    }
    end = pos;
    return out;
}

bool same_records(const std::vector<Record>& a, const std::vector<Record>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i].end != b[i].end || a[i].refs != b[i].refs || a[i].h.strips != b[i].h.strips || a[i].h.fans != b[i].h.fans ||
            a[i].h.triangles != b[i].h.triangles)
            return false;
    return true;
}

struct FaceList {
    std::vector<Face> faces;
};

// `file_enc`: the loop-count encoding once a list of this file has shown it (only lists with
// strips or fans can: the writer uses one encoding throughout); later lists are parsed with it
// alone, and with both only when they do not parse with it.
std::optional<FaceList> find_face_list(Span s, std::size_t start, std::size_t limit, std::uint32_t ec, std::optional<CountEnc>& file_enc) {
    static const std::uint8_t marker[2] = {0x02, 0x0a};
    std::size_t pos = start;
    while ((pos = find(s, marker, 2, pos, limit)) != std::string::npos) {
        std::size_t rp = pos + 2;
        std::uint32_t count = 0;
        if (read_count(s, rp, limit, count) && count >= 1 && count <= 10'000'000) {
            std::size_t attr = 0;
            std::optional<std::vector<Record>> recs;
            if (file_enc) recs = parse_records(s, rp, limit, ec, count, *file_enc, attr);
            if (!recs) {
                std::size_t a1 = 0, a2 = 0;
                auto by_p1g = parse_records(s, rp, limit, ec, count, CountEnc::P1g, a1);
                // without strips and fans no loop count is read, so the u16 reading is the same list
                const bool loops = by_p1g && std::any_of(by_p1g->begin(), by_p1g->end(), [](const Record& r) { return !r.h.strips.empty() || !r.h.fans.empty(); });
                auto by_u16 = by_p1g && !loops ? std::optional<std::vector<Record>>{} : parse_records(s, rp, limit, ec, count, CountEnc::U16BE, a2);
                if (by_p1g && by_u16) {
                    // Lists without strips or fans read identically under both encodings; anything else is a real conflict.
                    if (a1 != a2 || !same_records(*by_p1g, *by_u16)) fail(ErrorCode::Malformed, "face list parses differently under both loop-count encodings", rp);
                } else if ((by_p1g && loops) || (!by_p1g && by_u16)) {  // only lists with loops show the encoding
                    file_enc = by_p1g ? CountEnc::P1g : CountEnc::U16BE;
                }
                if (by_p1g) { recs = std::move(by_p1g); attr = a1; }
                else if (by_u16) { recs = std::move(by_u16); attr = a2; }
            }
            if (recs) {
                FaceList fl;
                bool complete = true;
                for (auto& r : *recs) {
                    if (!is_attr(s, attr, limit)) { complete = false; break; }
                    fl.faces.push_back({std::move(r), s.at(attr), Color{s.at(attr + 7), s.at(attr + 6), s.at(attr + 5), s.at(attr + 4)}});
                    attr += 8;
                }
                if (complete) return fl;
            }
        }
        ++pos;
    }
    return std::nullopt;
}

Scene decode(ByteView d, const DecodeOptions& options) {
    if (d.size < 8 || std::memcmp(d.data, kMagic, 8) != 0) fail(ErrorCode::Malformed, "missing V5_CFV2 magic");
    const Span s{d.data, d.size};
    const auto stores = nested_stores(s);
    auto innermost = [&](std::size_t off) -> const Store* {
        const Store* best = nullptr;
        for (const auto& st : stores)
            if (st.lo <= off && off < st.hi && (!best || st.hi - st.lo < best->hi - best->lo)) best = &st;
        return best;
    };
    const unsigned threads = options.threads;
    constexpr std::size_t kChunk = 4 << 20;  // bytes per scan task

    // 1. Point-list nodes: every kNode inside a store, scanned in chunks (a match may start in a
    //    chunk and end in the next), then decoded independently.
    const std::size_t nchunks = (s.n + kChunk - 1) / kChunk;
    std::vector<std::vector<std::size_t>> hits(nchunks);
    parallel_for(nchunks, 1, threads, [&](std::size_t b0, std::size_t e0) {
        for (std::size_t c = b0; c < e0; ++c) {
            const std::size_t lo = c * kChunk, hi = std::min(s.n, lo + kChunk + 3);
            for (std::size_t pos = lo; (pos = find(s, kNode, 4, pos, hi)) != std::string::npos && pos < lo + kChunk; pos += 4)
                if (const Store* st = innermost(pos); st && pos + 4 <= st->hi) hits[c].push_back(pos);
        }
    });
    std::vector<std::size_t> nodes;
    for (auto& h : hits) nodes.insert(nodes.end(), h.begin(), h.end());
    std::vector<std::optional<PointList>> decoded(nodes.size());
    parallel_for(nodes.size(), 64, threads, [&](std::size_t b0, std::size_t e0) {
        for (std::size_t i = b0; i < e0; ++i) decoded[i] = decode_points(s, nodes[i], innermost(nodes[i])->hi);
    });
    std::vector<PointList> lists;
    for (auto& pl : decoded)
        if (pl) lists.push_back(std::move(*pl));

    // 2. Face list of each point list (up to the next list), converted to a mesh.
    std::vector<std::optional<Mesh>> meshes(lists.size());
    parallel_for(lists.size(), 16, threads, [&](std::size_t b0, std::size_t e0) {
        std::optional<CountEnc> file_enc;  // learned per range: one writer, one encoding
        for (std::size_t i = b0; i < e0; ++i) {
            auto& pl = lists[i];
            std::size_t limit = innermost(pl.node)->hi;
            if (i + 1 < lists.size()) limit = std::min(limit, lists[i + 1].node);
            auto fl = find_face_list(s, pl.end, limit, pl.count, file_enc);
            if (!fl) continue;
            Mesh m;
            m.name = "body @" + std::to_string(pl.node);
            m.positions = std::move(pl.xyz);  // each list belongs to one range: not read again
            bool any_solid = false, any_surface = false;
            for (const auto& f : fl->faces) {
                const auto first = std::uint32_t(m.triangle_count());
                triangulate(f.rec, m.indices, 0);
                m.groups.push_back({first, std::uint32_t(m.triangle_count()) - first, f.color});
                if (f.rec.h.n) ((f.flags & 0x10) ? any_solid : any_surface) = true;
            }
            m.kind = any_solid && !any_surface ? BodyKind::Solid : any_surface && !any_solid ? BodyKind::Surface : BodyKind::Unknown;
            if (!m.groups.empty()) m.color = m.groups[0].color;
            meshes[i] = std::move(m);
        }
    });

    Scene scene;
    scene.format = "catpart";
    std::size_t without_faces = 0;
    for (auto& m : meshes) {
        if (m) scene.meshes.push_back(std::move(*m));
        else ++without_faces;
    }
    if (without_faces) scene.warnings.push_back(std::to_string(without_faces) + " point lists without a face list");
    if (scene.meshes.empty()) scene.warnings.push_back("no stored tessellation found (geometry only as exact B-rep)");
    return scene;
}

}  // namespace

const Decoder& catpart() {
    static const SimpleDecoder d({"catpart", "CATIA V5 Part (stored tessellation)", {"catpart"}}, sniff, decode);
    return d;
}

}  // namespace cadcore::formats
