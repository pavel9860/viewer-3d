#include "assembly_gt.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>

#include "harness.hpp"

using namespace cadcore;
namespace fs = std::filesystem;

namespace gt {
namespace {

struct Row {
    std::string name;
    int parent = 0, child = 0;
    double m[12] = {};
};

std::vector<Row> read_rows(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("missing fixture " + path);
    std::vector<Row> rows;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        Row r;
        std::getline(ls, r.name, '\t');
        ls >> r.parent >> r.child;
        for (double& v : r.m) ls >> v;
        rows.push_back(r);
    }
    return rows;
}

struct Occurrence {
    std::uint32_t node, parent_target, child_target;
};

}  // namespace

std::vector<std::uint8_t> read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error("missing ground-truth file " + p.string());
    return std::vector<std::uint8_t>((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

DecodeOptions folder_options(const fs::path& dir) {
    DecodeOptions o;
    o.load_resource = [dir](std::string_view rel) -> std::optional<std::vector<std::uint8_t>> {
        const fs::path p = dir / std::string(rel);
        std::error_code ec;
        if (!fs::is_regular_file(p, ec)) return std::nullopt;
        return read_file(p);
    };
    o.list_resources = [dir]() {
        std::vector<std::string> out;
        for (const auto& e : fs::recursive_directory_iterator(dir))
            if (e.is_regular_file()) out.push_back(fs::relative(e.path(), dir).generic_string());
        return out;
    };
    return o;
}

void check_assembly(const fs::path& p, const std::string& fixture, const std::string& format,
                    const std::function<bool(const std::string&)>& allowed) {
    const std::string product = p.string();
    CAPTURE(product);
    const auto bytes = read_file(p);
    Scene s;
    REQUIRE_NOTHROW(s = decode(ByteView(bytes.data(), bytes.size()), p.filename().string(), folder_options(p.parent_path())));
    CHECK(s.format == format);
    for (const auto& w : s.warnings) {
        CAPTURE(w);
        CHECK(allowed(w));
    }

    // Walk: product/reference node -> instance nodes -> target node -> ...
    std::vector<Occurrence> occ;
    std::vector<std::uint32_t> stack{s.nodes.at(0).children.at(0)};
    std::vector<bool> seen(s.nodes.size(), false);
    while (!stack.empty()) {
        const std::uint32_t target = stack.back();
        stack.pop_back();
        if (seen[target]) continue;
        seen[target] = true;
        for (const std::uint32_t inst : s.nodes[target].children) {
            if (s.nodes[inst].children.empty()) continue;  // the reference's own representation
            REQUIRE(s.nodes[inst].children.size() == 1);
            const std::uint32_t child = s.nodes[inst].children[0];
            occ.push_back({inst, target, child});
            stack.push_back(child);
        }
    }
    const auto rows = read_rows(std::string(CADCORE_FIXTURES_DIR) + "/" + fixture);
    CHECK(occ.size() == rows.size());

    std::multimap<std::string, std::size_t> by_name;
    for (std::size_t i = 0; i < occ.size(); ++i) by_name.emplace(s.nodes[occ[i].node].name, i);
    std::map<int, std::uint32_t> pd_to_node;
    std::map<std::uint32_t, int> node_to_pd;
    auto consistent = [&](int pd, std::uint32_t node) {
        const auto a = pd_to_node.find(pd);
        const auto b = node_to_pd.find(node);
        return (a == pd_to_node.end() || a->second == node) && (b == node_to_pd.end() || b->second == pd);
    };
    std::vector<bool> used(occ.size(), false);
    std::size_t unmatched = 0;
    double worst = 0;
    for (const Row& r : rows) {
        bool found = false;
        auto [lo, hi] = by_name.equal_range(r.name);
        for (auto it = lo; it != hi && !found; ++it) {
            const Occurrence& o = occ[it->second];
            if (used[it->second] || !consistent(r.parent, o.parent_target) || !consistent(r.child, o.child_target)) continue;
            const auto& t = s.nodes[o.node].transform;
            const double ours[12] = {t[0], t[1], t[2], t[4], t[5], t[6], t[8], t[9], t[10], t[12], t[13], t[14]};
            double d = 0;
            for (int k = 0; k < 12; ++k) d = std::max(d, std::abs(ours[k] - r.m[k]) / (1.0 + std::abs(r.m[k])));
            if (d > 1e-9) continue;
            worst = std::max(worst, d);
            used[it->second] = true;
            pd_to_node[r.parent] = o.parent_target;
            node_to_pd[o.parent_target] = r.parent;
            pd_to_node[r.child] = o.child_target;
            node_to_pd[o.child_target] = r.child;
            found = true;
        }
        if (!found) ++unmatched;
    }
    CAPTURE(worst);
    CHECK(unmatched == 0);
    std::printf("  %-60s %zu occurrences, worst relative placement difference %.3g\n", product.c_str(), occ.size(), worst);
}

}  // namespace gt
