# kernel

Research spikes for a browser-first B-rep kernel (Rust, WASM). Separate from the format decoders (`../cadcore`, `../catia`).

```
spike/src/
  geom.rs       vectors, Surf trait, analytic surfaces (plane, cylinder, sphere, torus, cone, Bezier), placement
  nurbs.rs      rational B-spline curves/surfaces, extrusion, revolution (derivatives to order 2)
  ssi.rs        surface-surface intersection: transversal, branch points, tangential, tolerant contacts
  ssi_cases.rs  hand-made hard cases + randomized analytic stress test
  p21.rs        minimal STEP Part 21 parser
  step.rs       edge check against CATIA STEP exports (local or full SSI)
  k.rs          deterministic PRNG
```

```bash
cargo build --release
./target/release/spike ssi                                   # 7 hard cases, < 1 s
./target/release/spike stress 60                             # 600 random trials, ~9 s
./target/release/spike step spike/catia_steps.tsv            # 27,782 CATIA edges, ~70 s
./target/release/spike step spike/catia_steps.tsv --pair bspline-cylinder --debug
./target/release/spike step spike/catia_steps.tsv --file Screw --full   # complete SSI per edge (slow)
```

Status (2026-11-23): SSI cases 7/7, stress 600/600, CATIA STEP edges 27,782 classified with 0 misses.
