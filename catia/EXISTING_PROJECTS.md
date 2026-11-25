# Existing projects that open CATIA V5 files (survey 2026-11-21)

## Open source, reads the binary directly

| Project | License | Language | What it gets out of a .CATPart / .CATProduct | Browser? |
|---|---|---|---|---|
| **[catpart-parser](https://github.com/MadManJohnSmith/catpart-parser)** (v1.0, 2026-11-18, active daily) | Apache-2.0 | Python reference + Rust→WASM port (665 KB `.wasm`) | **Triangle mesh** from the embedded tessellation (~80% of parts in their 1.8k-file corpus), container inventory/metadata (~99.7%), **assemblies**: instance tree + transforms `[tx,ty,tz,qw,qx,qy,qz]`, exploded view. Exact B-rep: no (facet subset only) | **Yes**: a three.js viewer, drag-drop CATPart or CATProduct + folder/.zip. [Live demo](https://madmanjohnsmith.github.io/catpart-parser/viewer/viewer.html) |
| **[cadmpeg](https://github.com/cadmpeg/cadmpeg)** (`cadmpeg-codec-catia`, ~39k lines of Rust) | Apache-2.0 code, CC-BY-4.0 docs | Rust CLI/library | **Exact geometry** (L3): planes, cylinders, cones, spheres, tori, NURBS, partial topology → STEP. Formulas/parameters partly. JPEG preview. **No tessellation, no product structure** | Not yet (Rust; could be compiled to WASM) |
| [catia-sketch-cpp](https://github.com/gcourchinoux/catia-sketch-cpp) | — | C++20 | Prototype: sketch points/constraints only, based on cadmpeg docs | No |

Both main projects are clean-room (byte analysis only) and well documented:
- catpart-parser `docs/format/*.md`: header, tail catalog, container table, nested store, vertex pool, facet records.
- cadmpeg `docs/formats/catia.md` (2,054 lines) + `docs/layouts/catia.toml` (byte offsets): B-rep records `b2 03 19` circle, `b2 03 28` cylinder, `a5 03 34` NURBS surface, and others.

## Automation only (needs a licensed CATIA running on Windows)
[pycatia](https://github.com/evereux/pycatia), DEH-CATIA (C#), MATLAB2CAD. Not useful for a standalone viewer.

## Commercial
HOOPS Exchange / Communicator (Tech Soft 3D; reads CATPart/CATProduct/CGR up to V5-6R2026), CAD Exchanger SDK (has JS), Spatial 3D InterOp, Datakit, 3DViewStation WebViewer (renders on the server), eDrawings for CATIA. These read everything but are paid and closed.

## What our samples confirm
`project K/Extruded.CATPart` (V5R28 SP6): header `V5_CFV2\0` + BE u32 catalog offset (0x78e9) + BE u32 catalog length (0x2528), and offset + length = file size. The tail catalog `CATIA_V5 CB0001` holds UTF-16BE stream names (`CATIA1`, `Format`, `Data`, `CATStorageProperty`, `CATPreview`, …). Containers: `CATProdCont`, `CATPrtCont`, `CGMGeom` (exact B-rep), `CATStdCont` (nested `V5_CFV2` holding the display mesh), `CATCGRCont`, `CATSeeBodyCont`. The feature tree names are plain strings (`MechanicalPart`, `HybridBody`, `xy-plane`, `Open_body`, `MANIFOLD_SOLID_BREP`). This matches both projects' docs.

## Recommendation
Don't start from zero:
1. **Mesh and assemblies**: use catpart-parser's WASM as the loader. It already does what a fast browser viewer needs.
2. **Exact geometry, later**: take cadmpeg's spec or crate (compiled to WASM) for picking real edges and faces, measurements, and STEP export.
3. **Our work**: the CATIA-style UI (spec tree, compass, toolbars, Part / Assembly / GSD workbench layouts, view modes), plus filling in the parts the parsers miss when tested against our 304 CATParts and 47 CATProducts.

## catpart-parser on our samples (2026-11-21)
Python reference at commit of 2026-11-21. Covers all 351 loose files under `samples/catia`; the .rar/.zip archives were not unpacked. Data, renders and scripts are in [eval/](eval/).

| | Files | Result |
|---|---|---|
| CATPart: metadata (`catpart-inspect`) | 304 | 304 OK. Releases V5R25–R28 |
| CATPart: mesh export | 304 | 294 exported (7 of them with 0 triangles). 10 `undecoded-variant`: all about 35 KB, so probably empty parts |
| Mesh quality (of 287 with triangles) | | **97 closed solids**, 17 under 5% open edges, 82 at 5–20%, **98 at 20% or more** |
| CATProduct | 47 | 47 metadata OK. No geometry, as expected because products only reference parts. The Python CLI does not decode the assembly tree. That lives in the Rust/WASM `product` module (not run here) |

Speed: median well under 1 s per part. The largest is 7.7 MB / 31k triangles in 4.6 s (pure Python; the WASM build should be faster).

Visual check ([renders](eval/catpart-parser-sample-renders.png), 12 parts sorted by open-edge ratio): closed parts and those under about 10% open edges look correct (plates with holes, flange, tubes). Parts with 20% or more open edges are **missing faces**: the resistor shows only rings with no body, one bracket comes out in fragments, and one shaft has gaps. These are the "curve-entry" loops that the parser deliberately leaves out.

**Verdict:** good enough to start a viewer on today for about 2/3 of our parts. The remaining third needs either a fix to the parser's curve-entry decode, or a fallback to cadmpeg's exact B-rep with our own tessellation.
