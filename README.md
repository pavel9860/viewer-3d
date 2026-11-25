# viewer-3d

A next-generation CAD viewer and kernel for the browser: native CAD file decoding, an exact B-rep geometry kernel, and WebAssembly builds that run in parallel in the browser.

```
cadcore/   C++ format decoders (CATIA V5 CATPart/CATProduct/CGR, 3DXML, STEP, glTF/Draco, DXF, ...) and the cadcli tool
kernel/    Rust B-rep kernel: exact geometry, topology, STEP import, validation, surface intersection, booleans; WASM build
catia/     Python reference decoder for CATIA V5 tessellation (grammar research)
FORMAT_COVERAGE.md   which formats and features are covered
```

Each folder has its own README with build and test commands. Sample CAD files are not part of the repository.
