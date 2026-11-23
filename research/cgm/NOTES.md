# CATIA V5 exact geometry (CGMGeom) — research notes

Goal: tessellate parts stored only as exact B-rep (333 CATParts in the corpus have no stored
tessellation; 55 CAA products lose geometry through them).

## Locating the stream
The Data section of class `CGMGeom` carries the container GUID (item 3); the geometry is the store
stream `Cont_<guid parts>` of that container (e.g. CAAPstBlockMovable_link.CATPart:
`Cont_80a90eda_3abc_3921380b_67611`, 12,149 bytes).

## Framing seen in an R12 file (CAAPstBlockMovable_link)
- Header: `01 70 2f 00 00 12 00 00 00 "  V5R12-HCN:052003" ...` then records.
- First record `e5 0d 06 01 00 <u32 79> <u32 12019>`: a directory (record count, size), 13 bytes.
- Records: `e5 0d <3-byte type> <u32 length> <u32 id> <body of length-4>`; ids increase from 1.
  - `13 c8 00` (121 bytes): model space — axes (f64 1.0 entries) and a box of ±1e6.
  - `03 96 00` (66): f64 data, four per body (planes?).
  - per edge loop a repeating group `03 fe 01`, `03 0d 00`, `03 ff 01`, `03 c0 00`, `03 0e 00` ×2.
- This framing chains to the stream end only in early-release files; R25–R28 files (most of the
  corpus) break after the first records (other markers, e.g. `e5 01 ...`): to be established.

`cgm_census.py PARTLIST`: locates each part's CGM stream and tabulates record framing and types.
