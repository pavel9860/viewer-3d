# catia — CATIA V5 reader and viewer

```
catia/
  catv5/            Python reference decoder (package)
    binary.py       P1g counts, f32 / u16 readers
    container.py    outer V5_CFV2 container, nested stores
    points.py       point lists: table / raw / texcoord layouts
    scenegraph.py   node scan: point lists vs symbol nodes
    facets.py       facet-record grammar, triangulation
    facelist.py     face lists, per-face colours
    mesh.py         bodies, part mesh, topology statistics
    export.py       PLY writer (per-face RGBA)
    __main__.py     CLI
  tests/            unit tests on synthetic records (stdlib unittest)
  tools/            corpus_eval.py (whole sample set), render_sheet.py (PNG sheets)
  review/           renders waiting for the user's visual review
  eval/             first evaluation of catpart-parser (superseded)
```

```bash
python3 -I -m catv5 PART.CATPart -o part.ply --json        # run from catia/
python3 -I -m unittest discover -s tests -t .
python3 -I tools/corpus_eval.py samples/catia OUT_DIR --ply
```

## Tessellation format (as decoded)

1. **Nested store.** A `V5_CFV2\0` header inside the file: `u32be catalog_offset`, `u32be catalog_len`, `ff×8`. The body is `[origin+32, origin+catalog_offset)`.
2. **Point list node.** It starts with `04 04 fd ff`. Within its first 64 bytes comes the list token `ff ff 02 00 01`, then:
   - `P1g element_count`, `P1g edge_bound`, `u8 has_tex`, `u8`, `u8`;
   - only when `has_tex`: `u8 dims` and `f32[dims·count]`;
   - `u8 encoding`. Encoding 0 is a 2-bit prefix-compression table over the lexicographically sorted points (opcodes 0 = xyz, 3 = yz, 2 = z, 1 = repeat), then `P1g words`, then `f32[words]`. Encoding 1 is raw `f32[3·count]`.
3. **Face list.** It follows the point list's per-element normal block:
   - `f32 sag`: the chordal tolerance in mm;
   - `u8 nLOD`, then nLOD × {`record`, `f32 lod_sag`, `attribute`}: coarse whole-body meshes, not faces;
   - `u8 kind`, `02 0a`, `P1g faces`, then that many facet records back to back;
   - one 8-byte attribute per face: `<flags> 04 04 ff A B G R`.
     - Flag `0x10` = the face belongs to a **solid**; surface faces have it clear.
     - Flag `0x80` = transparent: set on 22 of 22 faces with A < 255, and on no other face.
4. **Facet record.**
   - Layout: `01 <0x40|kinds|normal> [P1g per kind] ff u32 n [f32×3 normal] [P1g per strip] [P1g per fan] refs`.
   - Kinds: `0x01` triangles, `0x02` strips, `0x04` fans, `0x08` normal.
   - Refs are u8 or u16be, fixed by the requirement that each record ends exactly at the next structure.
5. **Not geometry.**
   - The environment-map texture: `P1g 256, P1g 256, P1g 262144`, then texels, then the name `CATIA_EnvironmentMap_Default_256x256…`.
   - The default reference-plane symbol nodes, which reuse the list token.

## Status (2026-11-21, samples in samples/catia, 304 CATParts)

Byte-exact results:
- 610 point lists decode exactly: 578 table, 28 raw, 4 with texture coordinates.
- Every point list has a face list.
- Every face list's sag and LOD preamble resolves exactly; there are 15 LOD meshes in total.
- Every face has its attribute.

Body kind, from the file's solid flag: 536 solids, 73 surfaces, 1 empty.

Indicators only (not ground truth):
- 527 solids are closed.
- 2 open solids (`part-ac39`, `part-e5d6`) have every open edge covered by neighbouring edges within 2·sag (T-junctions and differently sampled curves).
- 7 solids have edges used more than twice:
  - 2 are folds in CATIA's fans;
  - 60 edges are consistently oriented touching lumps;
  - 0 are unexplained.
- All 72 non-empty surfaces are open, as expected.

User-confirmed (visual review, 2026-11-21): the 8 solid edge-case parts in `review/solid_edge_cases.png` are correct. Yellow open bodies are surfaces, not broken solids.

Open:
- **No stored tessellation**: `part-40e3`, `part-af48`. Their geometry exists only as the exact B-rep in `CGMGeom`, which needs a B-rep decoder plus our own tessellator.
- **Nine ~35 KB parts produce no bodies.** Their strings show only an empty geometrical set, but that's an indicator, not byte-level proof.
- **Byte accounting is incomplete.** Not yet attributed: edge polylines (`01 01 …`), the normal block, the environment-map texture (identified, not walked), and other node types.
  - Since 2026-11-22 the grammar of the normal block (`P1g m, u16[m], ceil(normal_count/2) bytes`), the group list after the face list (facets, indexed and explicit polylines, `kind` = group count) and the node framing (`ff u32 length` before each geometry node) are known and implemented exactly for CGR and 3DXML binary reps (`cadcore/src/formats/catia/tessellation.hpp`). The CATPart decoder still locates face lists by search; moving it onto that exact parser is open.
- Not started: CATProduct assemblies, edges and wireframe display, the browser viewer (JS port), the CATIA-style UI.
