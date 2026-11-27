# CGMGeom decoding — hand-off (2026-11-27)

> Temporary. Delete (and commit the deletion) once a session has continued from it. Part names are
> omitted; find the cases with the kernel audit (failed faces per part) or by the symptoms below.

Scope: CATIA V5 CGMGeom exact B-rep -> flat IR (cadcore, C++) -> kernel `cgm` mesher (Rust, owned by the
kernel session: never edit `kernel/`, commit only cadcore files). No encoding exists. Last cadcore commit:
4ac49cc. Findings with evidence: `FINDINGS.md`; IR kinds for the kernel: `KERNEL_SPEC.md`.

## Code map
- `src/formats/catia/cgm.{hpp,cpp}`: stream location, framing, `compact_class` (code -> class).
- `cgm_compact.{hpp,cpp}`: compact record -> canonical (varint-dialect) record. 2012+ derived values in
  `canonical_records` passes; derived forms in `point_derived`, `line_derived`, `plane_derived`,
  `edge_curve_derived`; `Ctx` holds boxes, ranges, edge-curve maps (`Map::inverse`).
- `cgm_compact2009.cpp`, `cgm_compact_io.hpp` (In/Out, `bspline_tail`).
- `cgm_model.{hpp,cpp}`: decoder -> IR, kinds 1-12 surfaces, 1-5 p-curves (all emitted by default).
- `tools/cadcli.cpp`: `cgm PART [OUT] [--pending]`, `cgmrec PART [--raw]`.
- Tests: `tests/test_cgm.cpp` (corpus), `tests/test_cgm_compact.cpp` (byte-built grammar cases),
  `tests/perf.cpp` (CGM decode exponent 1.10, 863k records/s). `cadcore_tests` ~2.7 s, CGM 0.24 s.

## Status
Corpus (samples/, 1,261 CATParts): 489,893 of 490,991 faces decode; 64 parts have no B-rep stream (29
fixed dialect, not framed). Skipped: code 61 (338, one part), class 152/code 26 ellipse-like p-curve
(~280), unbounded extrusions (228: v start unknown), point-count mismatch 43, code 60 (36), Hermite
value count 30, edge-curve maps 4/15 (~25), code 15 (24).

## Open (vertex/box parameters; the kernel session's audit, classes by its numbering)
1. 2012+ NURBS derived lines use a box from the 231 trailing block, not the knots (a 2015 head part:
   true u = 0 where the block starts 0 12.83 0 12.5 and the knots start at -15.1). Block layout depends on its leading compact ints
   (05050505: box = knots; 051d0505, 05090509, ...: differs). Decode the block, set `cx->box` from it.
3. Circle arcs on the wrong half / overshooting (three user parts). Slot loops in a 2015 corpus
   part: both points of a half-circle edge resolve to one end (edge curve of a line and a circle, map 3
   s = -1; points coded 0x84 shared end and 0x18 start with refs[1] end).
4. Collapsed full-turn circle ranges (zero UV length between distinct vertices).
5. NURBS weights <= 0 or huge (81 faces in the kernel audit): likely 2012+ rational weights.
- Stored point parameters that are offsets from the edge curve's range start (a 2014 corpus part:
  stored 3.0 / 0.0 on a range -19.09..-16.09, edge curve of a line and a B-spline with map 3). No field separates them yet; the
  record's derived second parameter gives the right value through the map.

## How to measure
- Kernel with a fixed binary copy (the kernel session rebuilds it): copy `kernel/target/release/cgm`,
  compare HEAD cadcli vs new per part (failed faces, open edges); BREP_DEBUG=gaps,faces,open; examples
  `cgm_edge IR EDGE`, `cgm_loops IR FACE`, `cgm_face IR FACE`.
- Ground truth: CATIA's stored tessellation (`cadcli ply`), vertex/edge consistency, never other parsers.
- Regression: decoded faces per part must not drop; output byte-identical for refactors.
- Limits: `ulimit -v 8000000`, at most 4 workers.
