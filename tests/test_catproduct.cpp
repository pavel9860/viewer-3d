// CATProduct decoder against ground truth: the occurrence tables of STEP exports of the same
// assemblies (see assembly_gt.hpp).
#include <filesystem>
#include <string>

#include "assembly_gt.hpp"
#include "harness.hpp"

namespace fs = std::filesystem;

namespace {

// The only accepted losses are linked parts that store no tessellation (exact B-rep only),
// a CATPart limitation tracked by the corpus run, not an assembly one.
bool brep_only(const std::string& w) { return w.find("no stored tessellation found") != std::string::npos; }

}  // namespace

TEST_CASE("CATProduct assemblies match the STEP exports of the same products") {
    const fs::path catia(CADCORE_CATIA_DIR);
    gt::check_assembly(catia / "project K/Board.CATProduct", "Board.occurrences.tsv", "catproduct", brep_only);
    gt::check_assembly(catia / "local/9f754ead", "board_b.occurrences.tsv", "catproduct", brep_only);
    gt::check_assembly(catia / "local/dee2bd00", "part-0ac6.occurrences.tsv", "catproduct", brep_only);
    gt::check_assembly(fs::path(CADCORE_SAMPLES_DIR) / "local/0e91cca9", "AsmA.occurrences.tsv", "catproduct", brep_only);
}

namespace {

// Indented outline of the decoded tree: node names (file names for linked documents) and
// non-zero translations.
std::string outline(const cadcore::Scene& s) {
    std::string out;
    auto walk = [&](auto&& self, std::uint32_t n, int depth) -> void {
        std::string name = s.nodes[n].name;
        if (const auto slash = name.rfind('/'); slash != std::string::npos) name = name.substr(slash + 1);
        out += std::string(std::size_t(depth) * 2, ' ') + name;
        const auto& t = s.nodes[n].transform;
        if (t[12] != 0 || t[13] != 0 || t[14] != 0) out += " @" + std::to_string(t[12]) + "," + std::to_string(t[13]) + "," + std::to_string(t[14]);
        out += "\n";
        for (auto c : s.nodes[n].children) self(self, c, depth + 1);
    };
    walk(walk, 0, 0);
    return out;
}

}  // namespace

// CATIA R11-era products (SPEC containers in the Data stream; formats/catia/spec.hpp). No STEP
// export of these CAA samples exists, so this is a regression check of the decoded structure.
// Its consistency evidence: every linked instance resolves (by document UUID) to the CATPart
// or CATProduct file whose name the instance carries (Jack_Frame -> Jack_Frame.CATPart, ...),
// and the placements are the _Position values stored in the files.
TEST_CASE("CATProduct R11 SPEC containers decode to their assembly trees") {
    const fs::path caa = fs::path(CADCORE_SAMPLES_DIR) / "github/catia-v5/zhoujason457-stack__CAA_PrjSrc";
    auto tree = [&](const std::string& rel) {
        const fs::path p = caa / rel;
        CAPTURE(rel);
        const auto bytes = gt::read_file(p);
        cadcore::Scene s;
        REQUIRE_NOTHROW(s = cadcore::decode(cadcore::ByteView(bytes.data(), bytes.size()), p.filename().string(), gt::folder_options(p.parent_path())));
        CHECK(s.format == "catproduct");
        return outline(s);
    };
    CHECK(tree("CAAAssemblyUI.edu/CNext/resources/graphic/Fix.CATProduct") ==  // SPEC version 9
          "Fix\n  Fix\n    Jack_Frame\n      Jack_Frame.CATPart\n    Jack_Branch_3\n      Jack_Branch_3.CATPart\n");
    CHECK(tree("CAACATIAApplicationFrm.edu/InputData/Product2.CATProduct") ==  // version 19, linked sub-product
          "Product2\n  Product2\n    Product3.1\n      Product3\n        Part1.1.1\n          CAACafSearch.CATPart\n    Part1.1\n      Part1.CATPart\n");
    CHECK(tree("CAAArrangmentInterfaces.edu/CNext/resources/graphic/CAAArrEduIn.CATProduct") ==  // local references, placements
          "Product1\n  Product1\n    Item Reservation1.1 @-100.000000,-4400.000000,0.000000\n      Item Reservation1\n"
          "    Item Reservation2.1 @-4500.000000,-1700.000000,0.000000\n      Item Reservation2\n");
}
