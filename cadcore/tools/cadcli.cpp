// cadcli — command line front end.
//
//   cadcli info FILE...            decode, print one JSON line per file
//   cadcli gen OUT_DIR             write every canonical test part in every format
//   cadcli corpus DIR... [--json]  decode every file with a known extension, report per format;
//                                  files under */invalid/* or *MUSTFAIL* are expected to fail
//   cadcli bench FILE [N]          decode N times, print MB/s
//   cadcli compare A B TOL         match world triangles (or points) of A and B within TOL
//   cadcli triangles FILE          print world triangles, one per line
//   cadcli ply FILE OUT.ply        write the assembled world-space triangles with face colours
//   cadcli dxfcheck FILE...        decoded bounding box against AutoCAD's stored extents
//   cadcli svg FILE OUT.svg [PX]   top view (x right, y up) of triangles, polylines and points, for review
//   cadcli benchsuite DIR... [--write BASE.tsv] [--check BASE.tsv]
//                                  min-of-5 decode time of every file; --write records a baseline,
//                                  --check fails when a file got > 1.5x (+2 ms) slower than its baseline
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include "cadcore/catia.hpp"
#include "cadcore/compare.hpp"
#include "cadcore/decode.hpp"
#include "cadcore/stats.hpp"
#include "formats/catia/cgm_model.hpp"
#include "formats/catia/osm.hpp"
#include "formats/catia/spec.hpp"
#include "formats/catia/store.hpp"
#include "formats/dxf/document.hpp"
#include "formats/dxf/groups.hpp"
#include "formats/dxf/geometry.hpp"
#include "testgen/testgen.hpp"

namespace fs = std::filesystem;
using namespace cadcore;

static std::vector<std::uint8_t> read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary | std::ios::ate);
    if (!f) return {};
    std::vector<std::uint8_t> b(std::size_t(f.tellg()));
    f.seekg(0);
    f.read(reinterpret_cast<char*>(b.data()), std::streamsize(b.size()));
    return b;
}

static std::string json_escape(const std::string& s) {
    std::string o;
    for (char c : s) {
        if (c == '"' || c == '\\') { o += '\\'; o += c; }
        else if (static_cast<unsigned char>(c) < 0x20) o += ' ';
        else o += c;
    }
    return o;
}

// Folders whose CATIA documents may be linked from the file being decoded (default: its own).
static std::vector<fs::path> g_document_roots;
static bool g_corpus_roots = false;  // corpus mode: roots are the corpus arguments

// CATIA document id -> path, built once over g_document_roots (first 2 MB of each file).
static const std::map<std::string, fs::path>& document_index() {
    static std::map<std::string, fs::path> index;
    static std::vector<fs::path> built_for;
    if (built_for == g_document_roots) return index;
    index.clear();
    built_for = g_document_roots;
    std::vector<char> head(2u << 20);
    for (const fs::path& root : g_document_roots) {
        std::error_code ec;
        for (fs::recursive_directory_iterator it(root, fs::directory_options::skip_permission_denied, ec), end; it != end; it.increment(ec)) {
            if (ec) break;
            if (!it->is_regular_file(ec)) continue;
            std::string e = it->path().extension().string();
            std::transform(e.begin(), e.end(), e.begin(), ::tolower);
            if (e != ".catpart" && e != ".catproduct") continue;
            std::ifstream f(it->path(), std::ios::binary);
            f.read(head.data(), std::streamsize(head.size()));
            const auto id = catia_document_id(ByteView(head.data(), std::size_t(f.gcount())));
            if (id) index.emplace(*id, fs::absolute(it->path(), ec).lexically_normal());  // first found wins
        }
    }
    return index;
}

static DecodeOptions options_for(const fs::path& file) {
    DecodeOptions o;
    if (const char* t = std::getenv("CADCORE_THREADS")) o.threads = unsigned(std::atoi(t));  // 0 or unset: all cores (max 6)
    const fs::path dir = file.parent_path();
    o.load_resource = [dir](std::string_view rel) -> std::optional<std::vector<std::uint8_t>> {
        std::string r(rel);
        std::replace(r.begin(), r.end(), '\\', '/');
        // percent-decoding for URIs such as "my%20file.bin"
        std::string dec;
        for (std::size_t i = 0; i < r.size(); ++i) {
            if (r[i] == '%' && i + 2 < r.size()) { dec += char(std::stoi(r.substr(i + 1, 2), nullptr, 16)); i += 2; }
            else dec += r[i];
        }
        fs::path p = fs::path(dec).is_absolute() ? fs::path(dec) : dir / dec;
        std::error_code ec;
        if (!fs::is_regular_file(p, ec)) return std::nullopt;
        return read_file(p);
    };
    if (!g_corpus_roots) g_document_roots = {dir.empty() ? fs::path(".") : dir};
    o.locate_document = [](std::string_view id) -> std::optional<std::string> {
        const auto& idx = document_index();
        const auto it = idx.find(std::string(id));
        if (it == idx.end()) return std::nullopt;
        return it->second.string();  // absolute or root-relative; load_resource accepts both
    };
    return o;
}

struct Outcome {
    bool ok = false;
    std::string error, format;
    std::size_t triangles = 0, meshes = 0, warnings = 0;
    std::string first_warning;
    double ms = 0;
};

static Outcome run(const fs::path& p) {
    Outcome o;
    auto bytes = read_file(p);
    auto t0 = std::chrono::steady_clock::now();
    try {
        Scene s = decode(ByteView(bytes.data(), bytes.size()), p.filename().string(), options_for(p));
        o.ok = true;
        o.format = s.format;
        o.triangles = measure(s).triangles;
        o.meshes = s.meshes.size();
        o.warnings = s.warnings.size();
        if (!s.warnings.empty()) o.first_warning = s.warnings.front();
    } catch (const DecodeError& e) {
        o.error = e.what();
    } catch (const std::exception& e) {
        o.error = std::string("INTERNAL: ") + e.what();
    }
    o.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return o;
}

static int cmd_info(int argc, char** argv) {
    for (int i = 2; i < argc; ++i) {
        fs::path p = argv[i];
        auto bytes = read_file(p);
        try {
            Scene s = decode(ByteView(bytes.data(), bytes.size()), p.filename().string(), options_for(p));
            SceneStats st = measure(s);
            std::printf("{\"file\":\"%s\",\"format\":\"%s\",\"meshes\":%zu,\"nodes\":%zu,\"instances\":%zu,\"triangles\":%zu,\"area\":%.6g,"
                        "\"volume\":%.6g,\"bbox\":[%.6g,%.6g,%.6g,%.6g,%.6g,%.6g],\"unit_to_mm\":%g,\"warnings\":[",
                        json_escape(p.string()).c_str(), s.format.c_str(), s.meshes.size(), s.nodes.size(), st.instances, st.triangles, st.area,
                        st.signed_volume, st.min[0], st.min[1], st.min[2], st.max[0], st.max[1], st.max[2], s.unit_to_mm);
            for (std::size_t k = 0; k < s.warnings.size(); ++k) std::printf("%s\"%s\"", k ? "," : "", json_escape(s.warnings[k]).c_str());
            std::printf("]}\n");
        } catch (const DecodeError& e) {
            std::printf("{\"file\":\"%s\",\"error\":\"%s\"}\n", json_escape(p.string()).c_str(), json_escape(e.what()).c_str());
        }
    }
    return 0;
}

static int cmd_gen(int argc, char** argv) {
    if (argc < 3) return 2;
    fs::path out = argv[2];
    int n = 0;
    for (const auto& m : testgen::canonical_models())
        for (const auto& v : testgen::variants()) {
            if (v.id.rfind("catpart", 0) == 0 && !m.instances.empty()) continue;
            if (v.id == "catpart-r17-counts" && m.name == "dense_cylinder35k") continue;  // not representable with u16 counts
            fs::path dir = out / v.id;
            fs::create_directories(dir);
            for (const auto& [name, bytes] : v.write(m)) {
                std::ofstream(dir / name, std::ios::binary).write(reinterpret_cast<const char*>(bytes.data()), std::streamsize(bytes.size()));
                ++n;
            }
        }
    std::printf("wrote %d files to %s\n", n, out.string().c_str());
    return 0;
}

// corpus: every file with a supported extension must decode with no warnings, unless
// expectations.tsv (found next to a scanned directory) says otherwise. Conformance-negative
// suites (paths containing /invalid/ or MUSTFAIL) are reported separately: a viewer may
// accept a non-conformant file, so those counts are informational.
static std::map<std::string, std::pair<std::string, std::string>> load_expectations(const fs::path& dir) {
    std::map<std::string, std::pair<std::string, std::string>> out;
    std::ifstream f(dir / "expectations.tsv");
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto t1 = line.find('\t'), t2 = line.find('\t', t1 + 1);
        if (t1 == std::string::npos || t2 == std::string::npos) continue;
        out[(dir / line.substr(0, t1)).lexically_normal().string()] = {line.substr(t1 + 1, t2 - t1 - 1), line.substr(t2 + 1)};
    }
    return out;
}

static int cmd_corpus(int argc, char** argv) {
    std::set<std::string> exts;
    for (const Decoder* d : decoders())
        for (auto e : d->info().extensions) exts.insert(std::string(e));
    std::vector<fs::path> files;
    std::map<std::string, std::pair<std::string, std::string>> expect;
    g_document_roots.clear();
    g_corpus_roots = true;
    for (int i = 2; i < argc; ++i) g_document_roots.push_back(fs::path(argv[i]).lexically_normal());
    for (int i = 2; i < argc; ++i) {
        for (auto& [k, v] : load_expectations(argv[i])) expect[k] = v;
        for (auto it = fs::recursive_directory_iterator(argv[i], fs::directory_options::skip_permission_denied); it != fs::recursive_directory_iterator(); ++it) {
            if (!it->is_regular_file()) continue;
            std::string e = it->path().extension().string();
            if (e.empty()) continue;
            e = e.substr(1);
            std::transform(e.begin(), e.end(), e.begin(), ::tolower);
            if (exts.count(e)) files.push_back(it->path().lexically_normal());
        }
    }
    std::sort(files.begin(), files.end());
    struct Tally { int clean = 0, warned = 0, failed = 0, rejected_as_expected = 0, negative_rejected = 0, negative_accepted = 0, internal = 0; double ms = 0; std::uintmax_t bytes = 0; };
    std::map<std::string, Tally> by_ext;
    int unexpected = 0;
    std::set<std::string> seen_expect;
    for (const auto& p : files) {
        std::string e = p.extension().string().substr(1);
        std::transform(e.begin(), e.end(), e.begin(), ::tolower);
        const std::string ps = p.string();
        const bool negative = ps.find("/invalid/") != std::string::npos || ps.find("MUSTFAIL") != std::string::npos;
        auto ex = expect.find(ps);
        const std::string want = ex == expect.end() ? "" : ex->second.first;
        if (ex != expect.end()) seen_expect.insert(ps);
        Outcome o = run(p);
        Tally& t = by_ext[e];
        t.ms += o.ms;
        t.bytes += fs::file_size(p);
        if (o.error.rfind("INTERNAL", 0) == 0) { ++t.internal; ++unexpected; std::printf("INTERNAL %s: %s\n", ps.c_str(), o.error.c_str()); continue; }
        if (negative && want.empty()) { o.ok ? ++t.negative_accepted : ++t.negative_rejected; continue; }
        if (want == "reject") {
            if (o.ok) { ++unexpected; ++t.failed; std::printf("ACCEPTED-BUT-INVALID %s (%s)\n", ps.c_str(), ex->second.second.c_str()); }
            else ++t.rejected_as_expected;
            continue;
        }
        if (want == "warn") {
            if (o.ok && o.warnings) ++t.rejected_as_expected;
            else { ++unexpected; ++t.failed; std::printf("EXPECTED-WARNING-MISSING %s\n", ps.c_str()); }
            continue;
        }
        if (!o.ok) { ++t.failed; ++unexpected; std::printf("FAIL %s: %s\n", ps.c_str(), o.error.c_str()); continue; }
        if (o.warnings) { ++t.warned; ++unexpected; std::printf("WARN %s: %s\n", ps.c_str(), o.first_warning.c_str()); continue; }
        ++t.clean;
    }
    for (auto& [k, v] : expect)
        if (!seen_expect.count(k)) { ++unexpected; std::printf("EXPECTATION-WITHOUT-FILE %s\n", k.c_str()); }
    std::printf("\n%-8s %6s %6s %6s %9s %10s %10s %8s %9s\n", "ext", "clean", "warn", "fail", "expected", "neg:reject", "neg:accept", "internal", "MB/s");
    Tally total;
    for (auto& [e, t] : by_ext) {
        std::printf("%-8s %6d %6d %6d %9d %10d %10d %8d %9.1f\n", e.c_str(), t.clean, t.warned, t.failed, t.rejected_as_expected, t.negative_rejected,
                    t.negative_accepted, t.internal, t.ms > 0 ? t.bytes / 1e6 / (t.ms / 1e3) : 0.0);
        total.clean += t.clean; total.warned += t.warned; total.failed += t.failed; total.rejected_as_expected += t.rejected_as_expected;
        total.negative_rejected += t.negative_rejected; total.negative_accepted += t.negative_accepted; total.internal += t.internal;
    }
    std::printf("%-8s %6d %6d %6d %9d %10d %10d %8d\n", "TOTAL", total.clean, total.warned, total.failed, total.rejected_as_expected, total.negative_rejected,
                total.negative_accepted, total.internal);
    std::printf("unexpected outcomes: %d\n", unexpected);
    return unexpected ? 1 : 0;
}

static int cmd_bench(int argc, char** argv) {
    if (argc < 3) return 2;
    fs::path p = argv[2];
    int n = argc > 3 ? std::atoi(argv[3]) : 10;
    auto bytes = read_file(p);
    auto opts = options_for(p);
    auto t0 = std::chrono::steady_clock::now();
    std::size_t tris = 0;
    for (int i = 0; i < n; ++i) tris = decode(ByteView(bytes.data(), bytes.size()), p.filename().string(), opts).triangle_count();
    double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / n;
    std::printf("%s: %.2f ms/decode, %.1f MB/s, %zu triangles\n", p.filename().string().c_str(), s * 1e3, bytes.size() / 1e6 / s, tris);
    return 0;
}

static int cmd_compare(int argc, char** argv) {
    if (argc < 5) return 2;
    const fs::path pa = argv[2], pb = argv[3];
    const double tol = std::atof(argv[4]);
    auto ba = read_file(pa), bb = read_file(pb);
    const Scene a = decode(ByteView(ba.data(), ba.size()), pa.filename().string(), options_for(pa));
    const Scene b = decode(ByteView(bb.data(), bb.size()), pb.filename().string(), options_for(pb));
    const MatchReport t = match_triangles(world_triangles(a), world_triangles(b), tol);
    const MatchReport p = match_points(world_points(a), world_points(b), tol);
    std::printf("{\"triangles\":[%zu,%zu],\"unmatched\":[%zu,%zu],\"max_dev\":%.6g,\"points\":[%zu,%zu],\"unmatched_points\":[%zu,%zu],"
                "\"max_dev_points\":%.6g,\"identical\":%s}\n",
                t.a, t.b, t.unmatched_a, t.unmatched_b, t.max_deviation, p.a, p.b, p.unmatched_a, p.unmatched_b, p.max_deviation,
                t.identical() && p.identical() ? "true" : "false");
    return t.identical() && p.identical() ? 0 : 1;
}

static int cmd_triangles(int argc, char** argv) {
    if (argc < 3) return 2;
    const fs::path p = argv[2];
    auto bytes = read_file(p);
    const Scene s = decode(ByteView(bytes.data(), bytes.size()), p.filename().string(), options_for(p));
    for (const auto& t : world_triangles(s))
        std::printf("%.6g %.6g %.6g | %.6g %.6g %.6g | %.6g %.6g %.6g\n", t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8]);
    return 0;
}

static int cmd_ply(int argc, char** argv) {
    if (argc < 4) return 2;
    const fs::path p = argv[2];
    auto bytes = read_file(p);
    const Scene s = decode(ByteView(bytes.data(), bytes.size()), p.filename().string(), options_for(p));
    std::vector<float> xyz;
    std::vector<std::uint8_t> rgb;
    for_each_node(s, [&](const Node& n, const Matrix4& world) {
        for (auto mi : n.meshes) {
            const Mesh& m = s.meshes.at(mi);
            if (m.topology != Topology::Triangles) continue;
            for (std::size_t t = 0; t < m.triangle_count(); ++t) {
                Color c = m.color;
                for (const auto& g : m.groups)
                    if (t >= g.first && t < g.first + g.count) c = g.color;
                for (int k = 0; k < 3; ++k) {
                    const auto w = apply(world, &m.positions[std::size_t(m.indices[t * 3 + std::size_t(k)]) * 3]);
                    xyz.insert(xyz.end(), {float(w[0]), float(w[1]), float(w[2])});
                }
                rgb.insert(rgb.end(), {c.r, c.g, c.b, c.a});
            }
        }
    });
    const std::size_t nt = rgb.size() / 4;
    std::ofstream o(argv[3], std::ios::binary);
    o << "ply\nformat binary_little_endian 1.0\nelement vertex " << nt * 3 << "\nproperty float x\nproperty float y\nproperty float z\n"
      << "element face " << nt << "\nproperty list uchar int vertex_indices\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha\nend_header\n";
    o.write(reinterpret_cast<const char*>(xyz.data()), std::streamsize(xyz.size() * 4));
    for (std::size_t t = 0; t < nt; ++t) {
        const std::uint8_t three = 3;
        o.write(reinterpret_cast<const char*>(&three), 1);
        const std::int32_t idx[3] = {std::int32_t(t * 3), std::int32_t(t * 3 + 1), std::int32_t(t * 3 + 2)};
        o.write(reinterpret_cast<const char*>(idx), 12);
        o.write(reinterpret_cast<const char*>(&rgb[t * 4]), 4);
    }
    std::printf("%zu triangles -> %s\n", nt, argv[3]);
    return 0;
}

// svg: orthographic top view of the world-space scene (triangles filled, polylines stroked,
// points as dots) on a dark background, `px` wide.
static int cmd_svg(int argc, char** argv) {
    if (argc < 4) return 2;
    const fs::path p = argv[2];
    const double px = argc > 4 ? std::atof(argv[4]) : 900;
    auto bytes = read_file(p);
    const Scene s = decode(ByteView(bytes.data(), bytes.size()), p.filename().string(), options_for(p));
    struct Shape { int kind; std::vector<double> xy; Color c; };  // 0 triangle, 1 polyline, 2 point
    std::vector<Shape> shapes;
    double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
    auto xy = [](const Matrix4& w, const float* v) { const auto p = apply(w, v); return std::array<double, 2>{p[0], p[1]}; };
    for_each_node(s, [&](const Node& n, const Matrix4& world) {
        for (auto mi : n.meshes) {
            const Mesh& m = s.meshes.at(mi);
            if (m.topology == Topology::Points) {
                for (std::size_t i = 0; i < m.vertex_count(); ++i) {
                    const auto p = xy(world, &m.positions[3 * i]);
                    Shape sh{2, {p[0], p[1]}, m.color};
                    if (m.colors.size() >= 4 * (i + 1)) sh.c = {m.colors[4 * i], m.colors[4 * i + 1], m.colors[4 * i + 2], m.colors[4 * i + 3]};
                    shapes.push_back(sh);
                }
                continue;
            }
            for (std::size_t t = 0; t < m.triangle_count(); ++t) {
                Shape sh{0, {}, m.color};
                for (const auto& g : m.groups)
                    if (t >= g.first && t < g.first + g.count) sh.c = g.color;
                for (int k = 0; k < 3; ++k) {
                    const auto p = xy(world, &m.positions[3 * std::size_t(m.indices[3 * t + std::size_t(k)])]);
                    sh.xy.insert(sh.xy.end(), {p[0], p[1]});
                }
                shapes.push_back(std::move(sh));
            }
        }
        for (auto pi : n.polylines) {
            const Polyline& pl = s.polylines.at(pi);
            Shape sh{1, {}, pl.color};
            for (std::size_t i = 0; i + 2 < pl.points.size(); i += 3) {
                const auto p = xy(world, &pl.points[i]);
                sh.xy.insert(sh.xy.end(), {p[0], p[1]});
            }
            shapes.push_back(std::move(sh));
        }
    });
    for (const auto& sh : shapes)
        for (std::size_t i = 0; i + 1 < sh.xy.size(); i += 2) {
            lo[0] = std::min(lo[0], sh.xy[i]); hi[0] = std::max(hi[0], sh.xy[i]);
            lo[1] = std::min(lo[1], sh.xy[i + 1]); hi[1] = std::max(hi[1], sh.xy[i + 1]);
        }
    if (shapes.empty()) { lo[0] = lo[1] = 0; hi[0] = hi[1] = 1; }
    const double w = std::max(hi[0] - lo[0], 1e-9), h = std::max(hi[1] - lo[1], 1e-9);
    const double scale = (px - 20) / std::max(w, h);
    const double W = w * scale + 20, H = h * scale + 20;
    auto X = [&](double x) { return 10 + (x - lo[0]) * scale; };
    auto Y = [&](double y) { return H - 10 - (y - lo[1]) * scale; };
    auto rgb = [](Color c) { char b[32]; std::snprintf(b, sizeof b, "rgb(%d,%d,%d)", c.r, c.g, c.b); return std::string(b); };
    std::ofstream o(argv[3]);
    o << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << H << "\" viewBox=\"0 0 " << W << " " << H << "\">\n";
    o << "<rect width=\"100%\" height=\"100%\" fill=\"rgb(33,40,48)\"/>\n";
    char buf[64];
    for (const auto& sh : shapes) {
        if (sh.kind == 2) {
            o << "<circle cx=\"" << X(sh.xy[0]) << "\" cy=\"" << Y(sh.xy[1]) << "\" r=\"1.5\" fill=\"" << rgb(sh.c) << "\"/>\n";
            continue;
        }
        o << (sh.kind == 0 ? "<polygon" : "<polyline") << " points=\"";
        for (std::size_t i = 0; i + 1 < sh.xy.size(); i += 2) {
            std::snprintf(buf, sizeof buf, "%.2f,%.2f ", X(sh.xy[i]), Y(sh.xy[i + 1]));
            o << buf;
        }
        if (sh.kind == 0) o << "\" fill=\"" << rgb(sh.c) << "\" fill-opacity=\"" << sh.c.a / 255.0 << "\" stroke=\"none\"/>\n";
        else o << "\" fill=\"none\" stroke=\"" << rgb(sh.c) << "\" stroke-width=\"0.8\"/>\n";
    }
    o << "</svg>\n";
    std::printf("%zu shapes, %zu warnings -> %s\n", shapes.size(), s.warnings.size(), argv[3]);
    return 0;
}

// dxfcheck: the decoded drawing's world bounding box against AutoCAD's own extents ($EXTMIN /
// $EXTMAX). MATCH: equal; INSIDE: within them (expected when content such as text is not drawn);
// OUTSIDE: geometry beyond what AutoCAD computed (misplaced or extra geometry).
static int cmd_dxfcheck(int argc, char** argv) {
    bool verbose = false;
    int outside = 0, match = 0, inside = 0, skipped = 0, unreliable = 0;
    for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == "-v") { verbose = true; continue; }  // also print both boxes
        const fs::path p = argv[i];
        const auto bytes = read_file(p);
        const ByteView v(bytes.data(), bytes.size());
        formats::autocad::Document doc;
        formats::autocad::GroupStream gs;
        Scene s;
        try {
            gs = formats::autocad::read_groups(v);
            doc = formats::autocad::parse_document(gs);
            s = decode(v, p.filename().string(), options_for(p));
        } catch (const std::exception&) { ++skipped; continue; }
        if (!doc.has_extents) { ++skipped; continue; }
        // Are the stored extents those of the file's own geometry? Exact extents of the simplest
        // model-space entities (LINE, POINT, CIRCLE, ARC without an extrusion direction), from their
        // groups alone: if these reach outside, the writer did not compute the extents.
        double alo[3] = {1e300, 1e300, 1e300}, ahi[3] = {-1e300, -1e300, -1e300};
        auto agrow = [&](double x, double y, double z) {
            const double w[3] = {x, y, z};
            for (int k = 0; k < 3; ++k) { alo[k] = std::min(alo[k], w[k]); ahi[k] = std::max(ahi[k], w[k]); }
        };
        bool infinite = false;
        for (const auto& e : doc.entities) {
            if (e.integer(67) == 1) continue;
            if (e.type == "RAY" || e.type == "XLINE") { infinite = true; continue; }
            if (e.num(39) != 0) continue;  // thickness: an extruded body
            const bool tilted = e.has(210) && !(e.num(210) == 0 && e.num(220) == 0 && e.num(230, 1) == 1);
            if (tilted && e.type == "CIRCLE") {
                // centre in the object coordinate system (DXF arbitrary axis algorithm), then the exact
                // box of a circle with unit normal n: centre +- r sqrt(1 - n_k^2)
                const formats::autocad::Ocs o(e.point(210, {0, 0, 1}));
                const formats::autocad::Vec3 c = o.to_wcs(e.point(10)), n = o.az;
                const double r = e.num(40), nk[3] = {n.x, n.y, n.z}, ck[3] = {c.x, c.y, c.z};
                double lo3[3], hi3[3];
                for (int k = 0; k < 3; ++k) {
                    const double h = r * std::sqrt(std::max(0.0, 1 - nk[k] * nk[k]));
                    lo3[k] = ck[k] - h;
                    hi3[k] = ck[k] + h;
                }
                agrow(lo3[0], lo3[1], lo3[2]);
                agrow(hi3[0], hi3[1], hi3[2]);
                continue;
            }
            if (tilted && e.type != "LINE" && e.type != "POINT") continue;  // LINE / POINT are in world coordinates anyway
            if (e.type == "LINE") { agrow(e.num(10), e.num(20), e.num(30)); agrow(e.num(11), e.num(21), e.num(31)); }
            else if (e.type == "POINT") agrow(e.num(10), e.num(20), e.num(30));
            else if (e.type == "CIRCLE" || e.type == "ARC") {
                const double cx = e.num(10), cy = e.num(20), cz = e.num(30), r = e.num(40);
                const double pi = 3.14159265358979323846;
                double a0 = e.type == "ARC" ? e.num(50) : 0, a1 = e.type == "ARC" ? e.num(51) : 360;
                if (e.type == "ARC") {
                    a0 = std::fmod(std::fmod(a0, 360) + 360, 360);
                    a1 = std::fmod(std::fmod(a1, 360) + 360, 360);
                    if (a1 <= a0) a1 += 360;
                }
                agrow(cx + r * std::cos(a0 * pi / 180), cy + r * std::sin(a0 * pi / 180), cz);
                agrow(cx + r * std::cos(a1 * pi / 180), cy + r * std::sin(a1 * pi / 180), cz);
                for (int q = 0; q < 8; ++q) {  // quadrant points inside the sweep
                    const double a = q * 90.0;
                    if (a >= a0 && a <= a1) agrow(cx + r * std::cos(a * pi / 180), cy + r * std::sin(a * pi / 180), cz);
                }
            }
        }
        double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
        auto grow = [&](const std::array<double, 3>& w) { for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], w[k]); hi[k] = std::max(hi[k], w[k]); } };
        for_each_node(s, [&](const Node& n, const Matrix4& world) {
            for (auto mi : n.meshes) {
                const Mesh& m = s.meshes[mi];
                for (std::size_t k = 0; k + 2 < m.positions.size(); k += 3) grow(apply(world, &m.positions[k]));
            }
            for (auto pi : n.polylines) {
                const Polyline& pl = s.polylines[pi];
                for (std::size_t k = 0; k + 2 < pl.points.size(); k += 3) grow(apply(world, &pl.points[k]));
            }
        });
        const double e_lo[3] = {doc.extmin.x, doc.extmin.y, doc.extmin.z}, e_hi[3] = {doc.extmax.x, doc.extmax.y, doc.extmax.z};
        if (lo[0] > hi[0]) { std::printf("EMPTY    %s (extents set, nothing drawn; %zu warnings)\n", p.string().c_str(), s.warnings.size()); ++inside; continue; }
        double diag = 0;
        for (int k = 0; k < 3; ++k) diag = std::max(diag, e_hi[k] - e_lo[k]);
        const double tol = 5e-4 * std::max(diag, 1.0);  // AutoCAD's own curve extents are approximate (~1e-4)
        double out = 0, gap = 0;  // largest overshoot beyond, and largest shortfall within, the extents
        const int axes = doc.extents_z ? 3 : 2;  // R10 and older store extents without z
        for (int k = 0; k < axes; ++k) {
            out = std::max({out, e_lo[k] - lo[k], hi[k] - e_hi[k]});
            gap = std::max({gap, lo[k] - e_lo[k], e_hi[k] - hi[k]});
        }
        double aout = 0;
        if (alo[0] <= ahi[0])
            for (int k = 0; k < (doc.extents_z ? 3 : 2); ++k) aout = std::max({aout, e_lo[k] - alo[k], ahi[k] - e_hi[k]});
        if (aout > tol) {
            std::printf("UNRELIABLE %s  (stored extents exclude the file's own lines/arcs by %.4g: not computed by the writer)\n", p.string().c_str(), aout);
            ++unreliable;
            continue;
        }
        if (infinite && out > tol) {
            std::printf("INFINITE %s  (RAY/XLINE drawn as segments; not in AutoCAD's extents)\n", p.string().c_str());
            ++skipped;
            continue;
        }
        const char* verdict = out > tol ? "OUTSIDE" : gap > tol ? "INSIDE " : "MATCH  ";
        (out > tol ? outside : gap > tol ? inside : match)++;
        std::printf("%s %s  over %.4g  under %.4g  (extents size %.4g; %zu warnings)\n", verdict, p.string().c_str(), out, gap, diag, s.warnings.size());
        if (verbose)
            std::printf("    decoded  x %.6g..%.6g  y %.6g..%.6g  z %.6g..%.6g\n    extents  x %.6g..%.6g  y %.6g..%.6g  z %.6g..%.6g\n", lo[0], hi[0], lo[1], hi[1],
                        lo[2], hi[2], e_lo[0], e_hi[0], e_lo[1], e_hi[1], e_lo[2], e_hi[2]);
    }
    std::printf("\nMATCH %d  INSIDE %d  OUTSIDE %d  UNRELIABLE extents %d  skipped %d (no extents / not decodable / infinite lines)\n", match, inside, outside, unreliable, skipped);
    return outside ? 1 : 0;
}

static int cmd_benchsuite(int argc, char** argv) {
    std::vector<std::string> dirs;
    std::string write_to, check_against;
    double speedup = 0;  // --speedup K: every format must decode its baseline files K times faster in total
    for (int i = 2; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--write" && i + 1 < argc) write_to = argv[++i];
        else if (a == "--check" && i + 1 < argc) check_against = argv[++i];
        else if (a == "--speedup" && i + 1 < argc) speedup = std::atof(argv[++i]);
        else dirs.push_back(a);
    }
    std::set<std::string> exts;
    for (const Decoder* d : decoders())
        for (auto e : d->info().extensions) exts.insert(std::string(e));
    std::vector<fs::path> files;
    g_document_roots.clear();
    g_corpus_roots = true;
    for (const auto& d : dirs) {
        g_document_roots.push_back(fs::path(d).lexically_normal());
        for (auto it = fs::recursive_directory_iterator(d, fs::directory_options::skip_permission_denied); it != fs::recursive_directory_iterator(); ++it) {
            if (!it->is_regular_file()) continue;
            std::string e = it->path().extension().string();
            if (e.size() < 2) continue;
            e = e.substr(1);
            std::transform(e.begin(), e.end(), e.begin(), ::tolower);
            if (exts.count(e)) files.push_back(it->path().lexically_normal());
        }
    }
    std::sort(files.begin(), files.end());
    std::map<std::string, double> base;
    if (!check_against.empty()) {
        std::ifstream f(check_against);
        if (!f) { std::fprintf(stderr, "cannot read %s\n", check_against.c_str()); return 2; }
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            const auto t1 = line.find('\t'), t2 = line.find('\t', t1 + 1);
            base[line.substr(0, t1)] = std::stod(line.substr(t2 + 1));
        }
    }
    struct Agg { double bytes = 0, seconds = 0, base_ms = 0, now_ms = 0; int n = 0; };
    std::map<std::string, Agg> by_ext;
    std::ofstream out;
    if (!write_to.empty()) {
        out.open(write_to);
        out << "# path\tbytes\tmin_ms (cadcli benchsuite, min of 5 decodes)\n";
    }
    int regressions = 0, compared = 0;
    for (const auto& p : files) {
        const auto bytes = read_file(p);
        const auto opts = options_for(p);
        double best = 1e30;
        bool ok = true;
        for (int rep = 0; rep < 5 && ok; ++rep) {
            const auto t0 = std::chrono::steady_clock::now();
            try {
                (void)decode(ByteView(bytes.data(), bytes.size()), p.filename().string(), opts);
            } catch (const std::exception&) {
                ok = false;  // failures are the corpus command's business
            }
            best = std::min(best, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        }
        if (!ok) continue;
        std::string e = p.extension().string().substr(1);
        std::transform(e.begin(), e.end(), e.begin(), ::tolower);
        Agg& a = by_ext[e];
        a.bytes += double(bytes.size());
        a.seconds += best;
        ++a.n;
        const double ms = best * 1e3;
        if (out) out << p.string() << '\t' << bytes.size() << '\t' << ms << '\n';
        if (auto it = base.find(p.string()); it != base.end()) {
            ++compared;
            a.base_ms += it->second;
            a.now_ms += ms;
            if (ms > 1.5 * it->second + 2.0) {
                ++regressions;
                std::printf("SLOWER %s: %.2f ms (baseline %.2f ms)\n", p.string().c_str(), ms, it->second);
            }
        }
    }
    int slow_formats = 0;
    std::printf("\n%-10s %6s %10s %12s %10s\n", "ext", "files", "MB", "MB/s", "speedup");
    for (const auto& [e, a] : by_ext) {
        const double f = a.now_ms > 0 && a.base_ms > 0 ? a.base_ms / a.now_ms : 0;
        const bool short_ = speedup > 0 && f < speedup;
        slow_formats += short_;
        std::printf("%-10s %6d %10.1f %12.1f %9.2fx%s\n", e.c_str(), a.n, a.bytes / 1e6, a.bytes / 1e6 / a.seconds, f, short_ ? "  BELOW TARGET" : "");
    }
    if (speedup > 0) std::printf("speedup target %.1fx: %d formats below\n", speedup, slow_formats);
    if (!check_against.empty()) std::printf("compared %d files with the baseline: %d slower\n", compared, regressions);
    return regressions || slow_formats ? 1 : 0;
}

// Scene node tree: names, children, meshes, polylines, non-identity transforms.
static int cmd_tree(int argc, char** argv) {
    if (argc < 3) return 2;
    const fs::path p = argv[2];
    auto bytes = read_file(p);
    try {
        Scene s = decode(ByteView(bytes.data(), bytes.size()), p.filename().string(), options_for(p));
        std::vector<bool> on_path(s.nodes.size(), false);
        auto walk = [&](auto&& self, std::uint32_t n, int depth) -> void {
            const Node& node = s.nodes[n];
            std::printf("%*s[%u] %s", depth * 2, "", n, node.name.c_str());
            if (!node.meshes.empty()) {
                std::size_t tris = 0;
                for (auto m : node.meshes) tris += s.meshes[m].triangle_count();
                std::printf("  meshes %zu (%zu triangles)", node.meshes.size(), tris);
            }
            if (!node.polylines.empty()) std::printf("  polylines %zu", node.polylines.size());
            if (node.transform != kIdentity) std::printf("  t=(%g %g %g)", node.transform[12], node.transform[13], node.transform[14]);
            std::printf("\n");
            for (const auto& [k, v] : node.properties) std::printf("%*s  %s = %s\n", depth * 2, "", k.c_str(), v.c_str());
            if (on_path[n]) { std::printf("%*s(cycle)\n", depth * 2 + 2, ""); return; }
            on_path[n] = true;
            for (auto c : node.children) self(self, c, depth + 1);
            on_path[n] = false;
        };
        walk(walk, 0, 0);
        for (const auto& w : s.warnings) std::printf("warning: %s\n", w.c_str());
        for (const auto& w : s.notes) std::printf("note: %s\n", w.c_str());
    } catch (const DecodeError& e) {
        std::printf("error: %s\n", e.what());
        return 1;
    }
    return 0;
}

// CATIA object streams: every object with its class, name, attributes, relation header and relations.
// Streams of a CATIA store: name, size, segment count and the first bytes.
static int cmd_streams(int argc, char** argv) {
    if (argc < 3) return 2;
    auto bytes = read_file(argv[2]);
    const ByteView f(bytes.data(), bytes.size());
    try {
        const auto store = formats::catia::read_store(f);
        std::printf("V5_CFV%d, %zu streams\n", store.version, store.streams.size());
        for (const auto& st : store.streams) {
            const auto data = formats::catia::read_stream(f, store, st);
            std::string head;
            char hx[4];
            for (std::size_t i = 0; i < std::min<std::size_t>(24, data.size()); ++i) { std::snprintf(hx, sizeof hx, "%02x ", data[i]); head += hx; }
            std::printf("  %-32s id %u size %u segs %zu | %s\n", st.name.c_str(), st.id, st.size, st.segments.size(), head.c_str());
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s\n", e.what());
        return 1;
    }
    return 0;
}

// cgm FILE [OUT.cgmir]: the CGMGeom stream's dialect, records, decoded B-rep and left-out faces;
// writes the flat IR when OUT is given.
static int cmd_cgm(int argc, char** argv) {
    if (argc < 3) return 2;
    auto bytes = read_file(argv[2]);
    namespace cgm = formats::catia::cgm;
    try {
        const auto st = cgm::find_stream(ByteView(bytes.data(), bytes.size()));
        if (!st) { std::printf("{\"file\":\"%s\",\"cgm\":false}\n", argv[2]); return 0; }
        const cgm::Stream s = cgm::frame(ByteView(st->data(), st->size()));
        std::vector<std::vector<std::uint32_t>> shown;
        std::string shown_error;
        try { shown = cgm::see_body_lists(ByteView(bytes.data(), bytes.size())); } catch (const DecodeError& e) { shown_error = e.what(); }
        const cgm::Ir ir = cgm::decode(s, shown);
        std::size_t want = 0;
        for (const auto& r : s.records) want += r.type == 512;
        std::printf("{\"file\":\"%s\",\"dialect\":\"%s\",\"version\":\"%s\",\"records\":%zu,\"faces\":%zu,\"faces_decoded\":%zu,\"shells\":%zu,\"surfaces\":%zu,\"pcurves\":%zu,\"edges\":%zu,\"see_body\":\"%s\",\"skipped\":[",
                    argv[2], cgm::dialect_name(s.dialect), s.version.c_str(), s.records.size(), want, ir.faces.size(), ir.shells.size(), ir.surfaces.size(), ir.pcurves.size(), ir.edges.size(), shown_error.empty() ? (shown.empty() ? "none" : "ok") : "malformed");
        for (std::size_t i = 0; i < ir.skipped.size(); ++i) {
            std::string why = ir.skipped[i].why;
            for (char& c : why) if (c == '"' || c == '\\') c = '\'';
            std::printf("%s{\"id\":%u,\"type\":%u,\"why\":\"%s\"}", i ? "," : "", ir.skipped[i].id, ir.skipped[i].type, why.c_str());
        }
        std::printf("]}\n");
        if (argc > 3) {
            const auto b = cgm::serialize(ir);
            FILE* f = std::fopen(argv[3], "wb");
            if (!f || std::fwrite(b.data(), 1, b.size(), f) != b.size()) { std::fprintf(stderr, "cannot write %s\n", argv[3]); if (f) std::fclose(f); return 1; }
            std::fclose(f);
        }
    } catch (const std::exception& e) {
        std::string why = e.what();
        for (char& c : why) if (c == '"' || c == '\\') c = '\'';
        std::printf("{\"file\":\"%s\",\"error\":\"%s\"}\n", argv[2], why.c_str());
        return 1;
    }
    return 0;
}

static int cmd_objects(int argc, char** argv) {
    if (argc < 3) return 2;
    auto bytes = read_file(argv[2]);
    const std::string filter = argc > 3 ? argv[3] : "";
    try {
        auto streams = formats::catia::parse_streams(ByteView(bytes.data(), bytes.size()));
        auto spec = formats::catia::parse_spec_streams(ByteView(bytes.data(), bytes.size()));
        for (const auto& f : spec.failed) std::printf("container @%zu (%zu bytes) not decoded: %s (at byte %zu)\n", f.offset, f.size, f.message.c_str(), f.error_offset);
        for (auto& sp : spec.streams) streams.push_back(std::move(sp));
        for (auto& st : streams) {
            formats::catia::parse_records(st);
            std::printf("stream @%zu: %zu objects%s\n", st.offset, st.objects.size(), st.spec_version ? (" (SPEC v" + std::to_string(st.spec_version) + ")").c_str() : "");
            for (const auto& o : st.objects) {
                const std::string cls = st.name(o.cls);
                if (!filter.empty() && cls.find(filter) == std::string::npos) continue;
                std::printf("  #%u %s '%s' flag %u @%zu\n", o.id, cls.c_str(), st.name(o.name).c_str(), o.flag, o.offset);
                for (const auto& c : o.chunks) {
                    std::string v = formats::catia::value_text(st, c.value);
                    if (v.empty()) v = "-";
                    if (v.size() > 200) v = v.substr(0, 200) + "... (" + std::to_string(v.size()) + " chars)";
                    std::printf("      %s = %s\n", st.name(c.attr).c_str(), v.c_str());
                }
                std::string h;
                for (const auto& t : formats::catia::relation_header(o)) {
                    if (t.is_blob || t.is_guid) {
                        h += t.is_guid ? " guid[" : " blob[";
                        for (std::size_t k = 0; k < t.blob.size && k < 20; ++k) { char b[4]; std::snprintf(b, sizeof b, "%02x", t.blob[k]); h += b; }
                        h += "]";
                    } else {
                        h += " " + std::to_string(t.number);
                    }
                }
                if (!h.empty()) std::printf("      header:%s\n", h.c_str());
                for (const auto& r : formats::catia::relations(st, o)) {
                    std::string ids;
                    for (auto id : r.objects) ids += " " + std::to_string(id);
                    std::printf("      -> %s:%s\n", st.name(r.attr).c_str(), ids.c_str());
                }
            }
        }
    } catch (const DecodeError& e) {
        std::printf("error: %s\n", e.what());
        return 1;
    }
    return 0;
}

int main(int argc, char** argv) {
    const std::string cmd = argc > 1 ? argv[1] : "";
    if (cmd == "info") return cmd_info(argc, argv);
    if (cmd == "gen") return cmd_gen(argc, argv);
    if (cmd == "corpus") return cmd_corpus(argc, argv);
    if (cmd == "bench") return cmd_bench(argc, argv);
    if (cmd == "compare") return cmd_compare(argc, argv);
    if (cmd == "triangles") return cmd_triangles(argc, argv);
    if (cmd == "ply") return cmd_ply(argc, argv);
    if (cmd == "svg") return cmd_svg(argc, argv);
    if (cmd == "dxfcheck") return cmd_dxfcheck(argc, argv);
    if (cmd == "benchsuite") return cmd_benchsuite(argc, argv);
    if (cmd == "tree") return cmd_tree(argc, argv);
    if (cmd == "objects") return cmd_objects(argc, argv);
    if (cmd == "streams") return cmd_streams(argc, argv);
    if (cmd == "cgm") return cmd_cgm(argc, argv);
    std::fprintf(stderr, "usage: cadcli info FILE... | gen OUT_DIR | corpus DIR... | bench FILE [N] | compare A B TOL | triangles FILE | ply FILE OUT | svg FILE OUT [PX] | benchsuite DIR... [--write F] [--check F] | tree FILE | objects CATIA_FILE [CLASS] | streams CATIA_FILE | cgm CATPART [OUT.cgmir]\n");
    return 2;
}
