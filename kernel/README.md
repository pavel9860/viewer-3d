# kernel

A browser-first B-rep kernel and its research spikes (Rust, WASM). Separate from the format decoders (`../cadcore`, `../catia`).

```
brep/         kernel core: geom (Surf trait, analytic surfaces), nurbs (B-spline, extrusion, revolution), p21, topo (chunked COW arena, typed ids, provenance, tolerances),
              step (STEP -> Body import), check (topology + geometry validator), ssi, exact (predicates), trim (faces in UV),
              make (primitives), mesh (UV triangulation, volume), classify (point in body), boolean (union/intersection/difference on Body);
              `brep check|bool FILE|LIST`
kcore/        no_std shared math (PRNG, hash, own sin/sqrt, synthetic + real mesh workloads); native and wasm32;
              raster (exact CPU rasterizer: 1/16-px coverage, sort-middle over threads, cluster culling; tests/raster.rs,
              examples/raster_bench.rs FILE.mesh: native per-phase timing)
spike/src/
  ssi_cases.rs  hand-made hard cases + randomized analytic stress test
  pred.rs       exact predicate benchmark
  boolean.rs    boolean spike: shared edges, vertex registry, UV arrangement, classification, closure + volume
  p21.rs        minimal STEP Part 21 parser
  step.rs       edge check against CATIA STEP exports (local or full SSI)
  bench.rs      determinism, transfer, incremental recompute, hashing, naming, assembly instancing
  real.rs       real CATIA mesh workloads (per-part processing, hashing, assembly picking)
wasm/         kcore exports for the browser (nightly, real atomics: wasm/build.sh)
web/          mesh.html (viewer: WebGL, else the CPU renderer: raster.js driving kcore::raster in k.wasm on 4 threads; it
              draws the coarsest detail level within 1 px of the full mesh), render_test.html (renderer vs exact ray
              casting over 14 directions x 7 distances x 2 pans per model and per detail level: cracks, spills, occlusion,
              outline, determinism incl. 1 vs 4 threads; 800x600 budget 3.3 ms for the picked level), bench.html (frame times),
              index.html (WASM tests), serve.py (COOP/COEP, POST /result)
tools/        export_parts.py (CATIA -> PLY via cadcli)
```

```bash
cargo build --release -j4
./target/release/brep check spike/catia_steps.tsv   # import + validate 15 CATIA STEP files, ~40 s
./target/release/brep bool FILE.stp 1               # A-B, A^B, AuB of the largest shell with a cylinder (validator + volume identities)
cargo test --release -j4
./target/release/spike ssi                         # 7 hard SSI cases, < 1 s
./target/release/spike stress 60                   # 600 random SSI trials, ~9 s
./target/release/spike pred                        # predicate benchmark
./target/release/spike bool                        # 7 boolean cases, ~5 s
./target/release/spike bench                       # architecture benchmarks, ~10 s
./target/release/spike step spike/catia_steps.tsv  # 27,782 CATIA STEP edges, ~70 s (--pair/--file/--full/--debug)
python3 -I tools/export_parts.py samples/catia DATA ../cadcore/build/cadcli "ASSEMBLY.CATProduct"
./target/release/spike real DATA                   # real data; writes DATA/parts.bin/.idx for the web page
wasm/build.sh && cp DATA/parts.* web/              # then serve web/ (python3 -I web/serve.py) and open it
./target/release/brep mesh spike/catia_steps.tsv --out web/meshes   # meshes for mesh.html / render_test.html / bench.html (needs k.wasm)
BREP_DEBUG=sew,tess ./target/release/brep mesh FILE.stp            # diagnostics: BREP_DEBUG=1 (all) or topics check,bool,time,mesh,tess,sew,orient,open,tri
```

Status (2026-11-25): booleans (A-B, A^B, AuB of each file's largest shell with a cylinder) pass on all 15 CATIA STEP files: 0 checker violations, no dangling cuts, volume identities <= 1.5e-5, 25 ms - 2.6 s per operation, 39 s for all; brep import + check of 15 CATIA STEP files (83k edges, 13.9k faces): 0 violations; SSI 7/7, stress 600/600, booleans 7/7 (closed shells, volume error <= 1.6e-5), STEP edges
27,782 with 0 misses, exact predicates 0 wrong signs; browser results bit-identical to native, 3.8-4.1x on 4 workers.
Resources: at most 4 workers, 8 GB per process.
