// CATIA CGMGeom: varint integers, framing, compact transcoding and B-rep decode of corpus parts
// (skipped when the samples folder is absent).
#include <cmath>
#include <functional>
#include <optional>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <vector>

#include "formats/catia/cgm.hpp"
#include "formats/catia/cgm_compact.hpp"
#include "formats/catia/cgm_model.hpp"
#include "harness.hpp"

using namespace cadcore;
namespace cgm = cadcore::formats::catia::cgm;
namespace fs = std::filesystem;

namespace {
std::vector<std::uint8_t> load(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), {}};
}
// A corpus part's framed stream, or nullopt when the samples are absent.
std::optional<cgm::Stream> stream(const char* rel, std::vector<std::uint8_t>& keep) {
    const fs::path p = fs::path(CADCORE_SAMPLES_DIR) / rel;
    if (!fs::exists(p)) return std::nullopt;
    const auto file = load(p);
    const auto st = cgm::find_stream(ByteView(file.data(), file.size()));
    REQUIRE(st.has_value());
    keep.assign(st->begin(), st->end());
    return cgm::frame(ByteView(keep.data(), keep.size()));
}
cgm::Ir decode(const cgm::Stream& s, bool pending) {
    cgm::DecodeOptions opt;
    opt.pending_kinds = pending;
    return cgm::decode(s, {}, opt);
}
std::uint32_t var(std::vector<std::uint8_t> b) {
    std::size_t k = 0;
    return cgm::varint(ByteView(b.data(), b.size()), k);
}
}  // namespace

TEST_CASE("CGM varint: small values and byte-presence tags") {
    CHECK(var({0x80}) == 0);
    CHECK(var({0x85}) == 5);
    CHECK(var({0x08, 0xda}) == 0xda);
    CHECK(var({0x10, 0x01}) == 0x100);           // byte 1 only
    CHECK(var({0x18, 0x01, 0x01}) == 0x101);
    CHECK(var({0x30, 0x84, 0x02}) == 0x028400);  // bytes 1 and 2
    CHECK(var({0x38, 0x93, 0x9f, 0x01}) == 0x019f93);
}

TEST_CASE("CGM corpus parts: framing and B-rep decode") {
    std::vector<std::uint8_t> keep;
    const auto st = stream("caa/caadoc/CAAElecDeviceItf.edu/CNext/resources/graphic/CAAElecConnect/Part4.CATPart", keep);
    if (!st) return;
    const cgm::Stream& fr = *st;
    CHECK(fr.dialect == cgm::Dialect::Varint);
    CHECK(fr.records.size() == 221);
    const cgm::Ir ir = cgm::decode(fr);
    CHECK(ir.faces.size() == 8);
    CHECK(ir.skipped.empty());
    std::size_t solids = 0;
    for (const auto& sh : ir.shells) solids += sh.solid;
    CHECK(solids == 1);
    const auto bin = cgm::serialize(ir);
    CHECK(bin.size() > 8);
    CHECK(std::string(reinterpret_cast<const char*>(bin.data()), 5) == "CGMIR");
}

TEST_CASE("CGM compact records: canonical form") {
    // 2011Spring edge (compact integers, no edge-curve reference): v0 v1 p0 p1 as self - target
    const std::vector<std::uint8_t> edge = {0x09, 0x0f, 0x07, 0x0b, 0x21};
    cgm::Record r;
    r.id = 100;
    r.type = 511;
    r.body = ByteView(edge.data(), edge.size());
    const cgm::Canonical c = cgm::canonical(r, cgm::CompactLayout{true, false});
    REQUIRE(c.refs.size() == 4);
    CHECK(c.refs[0] == 98);   // 09: +2
    CHECK(c.refs[1] == 104);  // 0f: -4
    CHECK(c.refs[2] == 102);  // 07: -2
    CHECK(c.refs[3] == 103);  // 0b: -3
    CHECK((c.payload == std::vector<std::uint8_t>{0x21}));
    // shell: one face, bits kind 3 (11), sense +1; closing -1 (raw i16 when integers are not compact)
    const std::vector<std::uint8_t> shell = {0x81, 0x03, 0x13, 0xff, 0xff};
    r.type = 520;
    r.body = ByteView(shell.data(), shell.size());
    const cgm::Canonical s = cgm::canonical(r, cgm::CompactLayout{false, true});
    REQUIRE(s.refs.size() == 1);
    CHECK(s.refs[0] == 101);
    CHECK((s.payload == std::vector<std::uint8_t>{0x81, 3, 0, 1, 0, 0xff, 0xff}));
}

TEST_CASE("CGM compact layouts") {
    const struct { const char* path; std::uint8_t f0; bool ints, ref, derived, fall2009; } parts[] = {
        {"caa/caadoc/CAACommonLayoutItf.edu/CNext/resources/graphic/CAACloPspEduIn/GASKET-RF-150-4in.CATPart", 0, false, true, false, false},
        {"caa/caadoc/CAAAssemblyUI.edu/CNext/resources/graphic/Nut.CATPart", 0x06, true, true, false, false},
        {"caa/caadoc/CAAElecHarnessItf.edu/CNext/resources/graphic/CAA_EHIMultiBranchableSample/PN-06122000-0836.CATPart", 0x0b, true, true, false, true},
        {"caa/caadoc/CAAStructureInterfaces.edu/CNext/resources/graphic/WT15x54.CATPart", 0x01, false, true, false, true},
    };
    for (const auto& pt : parts) {
        std::vector<std::uint8_t> keep;
        const auto st = stream(pt.path, keep);
        if (!st) continue;
        CHECK(st->dialect == cgm::Dialect::Compact);
        const auto layout = cgm::compact_layout(*st);
        REQUIRE(layout.has_value());
        if (pt.f0) CHECK(st->features.at(0) == pt.f0);
        CHECK(layout->fall2009 == pt.fall2009);
        if (pt.fall2009) {
            CHECK(layout->compact_ints == pt.ints);
            for (const auto& c : cgm::canonical_records(*st, *layout)) CHECK(c.error.empty());
        } else if (pt.f0) {
            CHECK(layout->compact_ints == pt.ints);
            CHECK(layout->edge_curve_ref == pt.ref);
            CHECK(layout->derived == pt.derived);
        }
    }
}

// Corpus parts: faces decoded (all of them when `faces` equals `decoded`), with and without pending kinds.
TEST_CASE("CGM corpus parts: faces per stored form") {
    const struct { const char* path; bool pending; std::size_t faces, decoded; } parts[] = {
        // 2011Spring, 2011Fall, 2009Fall compact integers
        {"caa/caadoc/CAACommonLayoutItf.edu/CNext/resources/graphic/CAACloPspEduIn/GASKET-RF-150-4in.CATPart", false, 8, 8},
        {"caa/caadoc/CAAAssemblyUI.edu/CNext/resources/graphic/Nut.CATPart", false, 40, 40},
        {"caa/caadoc/CAAElecHarnessItf.edu/CNext/resources/graphic/CAA_EHIMultiBranchableSample/PN-06122000-0836.CATPart", false, 22, 22},
        // 2017Fall+ loop blocks without ints; 2018Fall+ Hermite closing byte; carriers over fillets
        {"github/catia-v5/HOLYNOTBAD__holydrone/Hardware/modules/misc/telemetry-cuav-air-link-4g-1.snapshot.3/Telemetry CUAV Air Link 4G.CATPart", false, 396, 396},
        {"nist/nist-ftc-ctc-pmi-cad-models/NIST-FTC-CTC-PMI-CAD-models/CATIA V5-6R2019/nist_ftc_07_asme1_rd_ct5290_ra.CATPart", false, 331, 331},
        {"local/d9dcb20f", true, 2286, 2286},
    };
    for (const auto& pt : parts) {
        std::vector<std::uint8_t> keep;
        const auto st = stream(pt.path, keep);
        if (!st) continue;
        const cgm::Ir ir = decode(*st, pt.pending);
        CHECK(ir.faces.size() == pt.decoded);
        CHECK(ir.faces.size() + ir.skipped.size() >= pt.faces);
        if (pt.faces == pt.decoded) CHECK(ir.skipped.empty());
    }
}

// Newer kinds: `decoded` faces; `check` tests their parameters.
struct KindPart { const char* path; std::size_t faces, decoded, count = 0; };  // count: kind instances (0: some)
void kind_part(const KindPart& pt, bool surface, std::uint32_t kind, const std::function<void(const cgm::Ir&, std::size_t)>& check) {
    std::vector<std::uint8_t> keep;
    const auto st = stream(pt.path, keep);
    if (!st) return;
    const cgm::Ir ir = decode(*st, false);
    CHECK(ir.faces.size() == pt.decoded);
    CHECK(ir.faces.size() + ir.skipped.size() >= pt.faces);
    std::size_t n = 0;
    const std::size_t count = surface ? ir.surfaces.size() : ir.pcurves.size();
    for (std::size_t i = 0; i < count; ++i)
        if ((surface ? ir.surfaces[i].kind : ir.pcurves[i].kind) == kind) { ++n; check(ir, i); }
    CHECK((pt.count ? n == pt.count : n > 0));
}

TEST_CASE("CGM surfaces: fillets (216), revolutions (211), ruled (243, 251), extrusions of 3D curves (210)") {
    const char* g = "caa/caadoc/CAAElecHarnessItf.edu/CNext/resources/graphic/";
    const std::string plug = std::string(g) + "CAA_EhiFLEXComputeBundleSegment_Sample/plug.CATPart";
    const std::string dryer = std::string(g) + "CAA_EHIMultiSegmentSample/hair-dryer.CATPart";
    kind_part({plug.c_str(), 571, 571, 28}, true, cgm::SurfFillet, [](const cgm::Ir& ir, std::size_t i) {
        const auto& p = ir.surfaces[i].p;
        const std::size_t k = std::size_t(p.at(0));
        REQUIRE(p.size() == 1 + 31 * k + 4);
        // first section: both contact points at the same distance from the centre
        auto d = [&](int a) { double t = 0; for (int j = 0; j < 3; ++j) { const double x = p[1 + k + a + j] - p[1 + k + 6 + j]; t += x * x; } return t; };
        CHECK(std::abs(d(0) - d(3)) < 1e-9 * (1 + d(0)));
    });
    kind_part({dryer.c_str(), 65, 65, 6}, true, cgm::SurfRevolution, [](const cgm::Ir& ir, std::size_t i) {
        const auto& p = ir.surfaces[i].p;
        REQUIRE(p.size() >= 15);
        CHECK((p[13] >= 1 && p[13] <= 3));
    });
    auto ruled = [](const cgm::Ir& ir, std::size_t i) {
        const auto& p = ir.surfaces[i].p;
        const std::size_t k = std::size_t(p.at(0));
        REQUIRE(p.size() == 1 + 19 * k + 4);
        CHECK(p[1] == p[1 + 19 * k]);          // u0 is the first knot
        CHECK(p[2 + 19 * k] > p[1 + 19 * k]);
    };
    for (const KindPart& pt : {KindPart{"local/5f11d4d4", 39, 39},
                                  KindPart{"local/04190b04", 717, 716},
                                  KindPart{"caa/caadoc/Doc/online/CAAScdStrUseCases/samples/EditCutout/Plate_0053.CATPart", 20, 20}})
        kind_part(pt, true, cgm::SurfRuled, ruled);
    // extrusion of a 3D curve: d, then a profile of form 2 (line) or 3 (B-spline: degree nk knots np rational poles)
    kind_part({"local/5fdf6e0d", 18, 15, 4}, true, cgm::SurfExtrusionCurve, [](const cgm::Ir& ir, std::size_t i) {
        const auto& p = ir.surfaces[i].p;
        REQUIRE(p.size() > 4);
        CHECK(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] > 0);
        CHECK((p[3] == 2 || p[3] == 3));
        if (p[3] == 3) {
            const std::size_t deg = std::size_t(p[4]), nk = std::size_t(p[5]), np = std::size_t(p[6 + nk]);
            CHECK(nk == np + deg + 1);
            CHECK(p.size() == 8 + nk + 3 * np + (p[7 + nk] != 0 ? np : 0));
        }
    });
}

TEST_CASE("CGM p-curves: B-splines (170), expressions (155)") {
    // B-spline coedge ranges lie in their knot ranges (vertex parameters are offsets within a knot
    // segment, the first segment starting at the first knot)
    auto bspline = [](const cgm::Ir& ir, std::size_t i) {
        const auto& p = ir.pcurves[i].p;
        const std::size_t deg = std::size_t(p.at(0)), nk = std::size_t(p.at(1)), np = std::size_t(p.at(2 + nk));
        REQUIRE(nk == np + deg + 1);
        REQUIRE(p.size() == 4 + nk + 2 * np);
        const double lo = p[2 + deg], hi = p[2 + nk - deg - 1], tol = 1e-4 * (hi - lo);
        for (const auto& ce : ir.coedges)
            if (ce.pcurve == i) {
                CHECK(std::min(ce.t0, ce.t1) >= lo - tol);
                CHECK(std::max(ce.t0, ce.t1) <= hi + tol);
            }
    };
    for (const KindPart& pt : {KindPart{"local/d04dc2b9", 299, 299},
                                  KindPart{"local/5838914f", 507, 507},
                                  KindPart{"github/catia-v5/Guilherme115__Motor-V6/05_Admissao_Escape/Duto.CATPart", 379, 379}})
        kind_part(pt, false, cgm::PcBSpline, bspline);
    // expression trees: two well-formed trees after the 4 head values
    auto expr = [](const cgm::Ir& ir, std::size_t i) {
        const auto& p = ir.pcurves[i].p;
        std::size_t k = 4;
        std::function<void()> node = [&] {
            REQUIRE(k + 2 <= p.size());
            const int code = int(p[k]);
            CHECK((code == 1 || code == 2 || code == 3 || code == 4 || code == 7));
            k += 2 + std::size_t(p[k + 1]);
            REQUIRE(k < p.size());
            const int kids = int(p[k++]);
            for (int j = 0; j < kids; ++j) node();
        };
        node();
        node();
        CHECK(k == p.size());
    };
    for (const KindPart& pt : {KindPart{"caa/caadoc/CAAAssemblyUI.edu/CNext/resources/graphic/Jack_Frame.CATPart", 216, 216},
                                  KindPart{"local/703cd616", 3872, 3872},
                                  KindPart{"local/b33cd8bf", 328, 328}})
        kind_part(pt, false, cgm::PcExpr, expr);
}

TEST_CASE("CGM points: parameters are offsets within a Hermite segment") {
    // segment numbers are shifted by the knot vector's offset field (-7 on trimmed curves of this part)
    std::vector<std::uint8_t> keep;
    const auto st = stream("nist/nist-ctc-pmi-cad-models/NIST-CTC-PMI/CATIA V5R21/nist_ctc_02_asme1_ct5210_rc.CATPart", keep);
    if (!st) return;
    const cgm::Ir ir = cgm::decode(*st);
    std::size_t hermite = 0, outside = 0;
    for (const auto& c : ir.coedges) {
        const auto& pc = ir.pcurves[c.pcurve];
        if (pc.kind != cgm::PcHermite) continue;
        const auto n = std::size_t(pc.p[0]);
        const double k0 = pc.p[1], k1 = pc.p[n], tol = 1e-9 * (1 + std::abs(k1));
        ++hermite;
        if (std::min(c.t0, c.t1) < k0 - tol || std::max(c.t0, c.t1) > k1 + tol) ++outside;
    }
    CHECK(hermite > 100);
    CHECK(outside == 0);
    for (const auto& k : ir.skipped) CHECK(k.why.find("point segment") == std::string::npos);
}
