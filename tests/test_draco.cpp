// Draco decoder against ground truth: the uncompressed OBJ sources that the Draco project
// encoded into its test streams (samples/draco-testdata). Every decoded corner must match a
// source corner in position, normal and texture coordinate within the encoder's quantization.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "cadcore/error.hpp"
#include "codec/draco.hpp"
#include "harness.hpp"

using namespace cadcore;

namespace {

const std::string kDir = std::string(CADCORE_SAMPLES_DIR) + "/draco-testdata/";

std::vector<std::uint8_t> read(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("missing ground-truth file " + path);
    return std::vector<std::uint8_t>((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

struct Corner {
    double p[3] = {0, 0, 0}, n[3] = {0, 0, 0}, uv[2] = {0, 0};
    bool has_n = false, has_uv = false;
};
using Tri = std::array<Corner, 3>;

// Test-local OBJ reader for the triangle-only sources (v, vt, vn, f with a/b/c indices).
std::vector<Tri> read_obj(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("missing ground-truth file " + path);
    std::vector<std::array<double, 3>> v, vn;
    std::vector<std::array<double, 2>> vt;
    std::vector<Tri> tris;
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream ls(line);
        std::string kw;
        ls >> kw;
        if (kw == "v") { std::array<double, 3> a{}; ls >> a[0] >> a[1] >> a[2]; v.push_back(a); }
        else if (kw == "vn") { std::array<double, 3> a{}; ls >> a[0] >> a[1] >> a[2]; vn.push_back(a); }
        else if (kw == "vt") { std::array<double, 2> a{}; ls >> a[0] >> a[1]; vt.push_back(a); }
        else if (kw == "f") {
            std::vector<Corner> poly;
            std::string tok;
            while (ls >> tok) {
                Corner c;
                int idx[3] = {0, 0, 0};
                std::size_t start = 0;
                for (int k = 0; k < 3; ++k) {
                    const std::size_t slash = tok.find('/', start);
                    const std::string part = tok.substr(start, slash == std::string::npos ? std::string::npos : slash - start);
                    if (!part.empty()) idx[k] = std::stoi(part);
                    if (slash == std::string::npos) break;
                    start = slash + 1;
                }
                for (int i = 0; i < 3; ++i) c.p[i] = v.at(std::size_t(idx[0] - 1))[std::size_t(i)];
                if (idx[1]) { c.has_uv = true; for (int i = 0; i < 2; ++i) c.uv[i] = vt.at(std::size_t(idx[1] - 1))[std::size_t(i)]; }
                if (idx[2]) {
                    c.has_n = true;
                    const auto& nn = vn.at(std::size_t(idx[2] - 1));
                    const double len = std::sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
                    for (int i = 0; i < 3; ++i) c.n[i] = nn[std::size_t(i)] / len;
                }
                poly.push_back(c);
            }
            if (poly.size() != 3) throw std::runtime_error("ground-truth OBJ must be triangulated: " + path);
            tris.push_back({poly[0], poly[1], poly[2]});
        }
    }
    return tris;
}

std::vector<Tri> read_drc(const std::string& path) {
    const auto bytes = read(path);
    const codec::draco::Mesh m = codec::draco::decode(ByteView(bytes.data(), bytes.size()));
    const codec::draco::Attribute *pos = nullptr, *nrm = nullptr, *uv = nullptr;
    for (const auto& a : m.attributes) {
        if (a.type == 0 && !pos) pos = &a;
        if (a.type == 1 && !nrm) nrm = &a;
        if (a.type == 3 && !uv) uv = &a;
    }
    if (!pos) throw std::runtime_error("no position attribute in " + path);
    std::vector<Tri> tris;
    for (std::size_t f = 0; f < m.faces.size() / 3; ++f) {
        Tri t;
        for (int k = 0; k < 3; ++k) {
            const std::size_t p = m.faces[f * 3 + std::size_t(k)];
            Corner& c = t[std::size_t(k)];
            for (int i = 0; i < 3; ++i) c.p[i] = pos->values[p * 3 + std::size_t(i)];
            if (nrm) { c.has_n = true; for (int i = 0; i < 3; ++i) c.n[i] = nrm->values[p * 3 + std::size_t(i)]; }
            if (uv) { c.has_uv = true; for (int i = 0; i < 2; ++i) c.uv[i] = uv->values[p * 2 + std::size_t(i)]; }
        }
        tris.push_back(t);
    }
    return tris;
}

double dist3(const double* a, const double* b) {
    return std::sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]));
}

struct Deviation {
    std::size_t unmatched = 0;
    double pos = 0, normal_deg = 0, uv = 0;
};

// Brute-force one-to-one matching (the sources are small or spatially hashed below).
Deviation compare(const std::vector<Tri>& truth, const std::vector<Tri>& decoded, double pos_tol) {
    Deviation d;
    std::map<std::array<long long, 3>, std::vector<std::size_t>> grid;
    const double cell = pos_tol * 4;
    auto key = [&](const Tri& t) {
        std::array<long long, 3> k{};
        for (int i = 0; i < 3; ++i) k[std::size_t(i)] = (long long)std::floor((t[0].p[i] + t[1].p[i] + t[2].p[i]) / 3 / cell);
        return k;
    };
    for (std::size_t i = 0; i < decoded.size(); ++i) grid[key(decoded[i])].push_back(i);
    std::vector<bool> used(decoded.size(), false);
    for (const Tri& t : truth) {
        const auto k = key(t);
        double best = -1;
        std::size_t bj = 0;
        int brot = 0;
        for (long long dx = -1; dx <= 1; ++dx)
            for (long long dy = -1; dy <= 1; ++dy)
                for (long long dz = -1; dz <= 1; ++dz) {
                    auto it = grid.find({k[0] + dx, k[1] + dy, k[2] + dz});
                    if (it == grid.end()) continue;
                    for (const std::size_t j : it->second) {
                        if (used[j]) continue;
                        for (int rot = 0; rot < 3; ++rot) {
                            double worst = 0;
                            for (int c = 0; c < 3; ++c) worst = std::max(worst, dist3(t[std::size_t(c)].p, decoded[j][std::size_t((c + rot) % 3)].p));
                            if (worst <= pos_tol && (best < 0 || worst < best)) { best = worst; bj = j; brot = rot; }
                        }
                    }
                }
        if (best < 0) { ++d.unmatched; continue; }
        used[bj] = true;
        d.pos = std::max(d.pos, best);
        for (int c = 0; c < 3; ++c) {
            const Corner& a = t[std::size_t(c)];
            const Corner& b = decoded[bj][std::size_t((c + brot) % 3)];
            if (a.has_n != b.has_n || a.has_uv != b.has_uv) { ++d.unmatched; break; }
            if (a.has_n) {
                const double dot = std::clamp(a.n[0] * b.n[0] + a.n[1] * b.n[1] + a.n[2] * b.n[2], -1.0, 1.0);
                d.normal_deg = std::max(d.normal_deg, std::acos(dot) * 180.0 / 3.14159265358979);
            }
            if (a.has_uv) d.uv = std::max(d.uv, std::max(std::abs(a.uv[0] - b.uv[0]), std::abs(a.uv[1] - b.uv[1])));
        }
    }
    for (bool u : used) d.unmatched += u ? 0 : 1;
    return d;
}

void check_pair(const std::string& obj, const std::string& drc, double pos_tol, double max_normal_deg, double max_uv) {
    CAPTURE(drc);
    const auto truth = read_obj(kDir + obj);
    std::vector<Tri> decoded;
    REQUIRE_NOTHROW(decoded = read_drc(kDir + drc));
    REQUIRE(decoded.size() == truth.size());
    const Deviation d = compare(truth, decoded, pos_tol);
    std::printf("  %-40s pos %.3g  normal %.3g deg  uv %.3g  unmatched %zu\n", drc.c_str(), d.pos, d.normal_deg, d.uv, d.unmatched);
    CAPTURE(d.pos);
    CAPTURE(d.normal_deg);
    CAPTURE(d.uv);
    CHECK(d.unmatched == 0);
    CHECK(d.pos <= pos_tol);
    CHECK(d.normal_deg <= max_normal_deg);
    CHECK(d.uv <= max_uv);
}

}  // namespace

TEST_CASE("Draco test_nm streams of every bitstream version match test_nm.obj") {
    for (const char* f : {"test_nm.obj.edgebreaker.0.9.1.drc", "test_nm.obj.edgebreaker.0.10.0.drc", "test_nm.obj.edgebreaker.1.0.0.drc",
                          "test_nm.obj.edgebreaker.1.1.0.drc", "test_nm.obj.edgebreaker.cl4.2.2.drc", "test_nm.obj.edgebreaker.cl10.2.2.drc",
                          "test_nm.obj.sequential.0.9.1.drc", "test_nm.obj.sequential.0.10.0.drc", "test_nm.obj.sequential.1.0.0.drc",
                          "test_nm.obj.sequential.1.1.0.drc", "test_nm.obj.sequential.cl3.2.2.drc", "test_nm_quant.0.9.0.drc"})
        check_pair("test_nm.obj", f, 2e-3, 2.0, 0);
}

TEST_CASE("Draco cube_att streams match cube_att.obj (positions, texture coordinates, normals)") {
    for (const char* f : {"cube_att.drc", "cube_att.obj.edgebreaker.cl4.2.2.drc", "cube_att.obj.edgebreaker.cl10.2.2.drc",
                          "cube_att.obj.sequential.cl3.2.2.drc", "cube_att_sub_o_no_metadata.drc"})
        check_pair("cube_att.obj", f, 1e-3, 1.0, 1e-3);
    check_pair("cube_att_sub_o_2.obj", "cube_att_sub_o_2.drc", 1e-3, 1.0, 1e-3);
}

TEST_CASE("Draco bunny stream matches bunny_norm.obj") { check_pair("bunny_norm.obj", "bunny_gltf.drc", 1e-4, 2.0, 0); }

TEST_CASE("Draco kd-tree point cloud equals the sequential encoding of the same cloud") {
    auto rows = [](const std::string& f) {
        const auto bytes = read(kDir + f);
        const codec::draco::Mesh m = codec::draco::decode(ByteView(bytes.data(), bytes.size()));
        std::vector<std::vector<float>> r(m.num_points);
        for (const auto& a : m.attributes)
            for (std::uint32_t i = 0; i < m.num_points; ++i)
                for (int c = 0; c < a.components; ++c) r[i].push_back(a.values[std::size_t(i) * std::size_t(a.components) + std::size_t(c)]);
        std::sort(r.begin(), r.end());
        return r;
    };
    std::vector<std::vector<float>> seq, kd;
    REQUIRE_NOTHROW(seq = rows("pc_color.drc"));
    REQUIRE_NOTHROW(kd = rows("pc_kd_color.drc"));
    CHECK(seq.size() == 7733);
    CHECK(seq == kd);
}
