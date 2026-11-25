"""Decode every .CATPart under a directory and report per-body topology.

Usage::

    python3 -I tools/corpus_eval.py SAMPLES_DIR OUT_DIR [--ply]

Writes ``OUT_DIR/results.jsonl`` (one row per file) and prints a summary. With
``--ply`` each part's mesh is written to ``OUT_DIR/ply/<hash>.ply``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import traceback
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from catv5 import decode_part  # noqa: E402
from catv5.__main__ import summarize  # noqa: E402
from catv5.export import write_ply  # noqa: E402


def evaluate(job: tuple[str, str, str, bool]) -> dict:
    root, path, out, want_ply = job
    f = Path(path)
    row = {"file": str(f.relative_to(root)), "size": f.stat().st_size}
    t0 = time.perf_counter()
    try:
        mesh = decode_part(f.read_bytes())
        row.update(summarize(mesh))
        if want_ply:
            p = Path(out) / "ply" / (hashlib.md5(row["file"].encode()).hexdigest()[:10] + ".ply")
            with open(p, "wb") as fh:
                write_ply(mesh, fh)
            row["ply"] = str(p)
    except Exception:
        row["error"] = traceback.format_exc(limit=-3)
    row["seconds"] = round(time.perf_counter() - t0, 3)
    return row


def classify(row: dict) -> str:
    if "error" in row:
        return "error"
    if not row["bodies"]:
        return "no-bodies"
    if all(b["closed"] for b in row["bodies"]):
        return "all-bodies-closed"
    return "has-open-or-nonmanifold-body"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--ply", action="store_true")
    args = ap.parse_args()
    (args.out / "ply").mkdir(parents=True, exist_ok=True)
    files = sorted(p for p in args.root.rglob("*") if p.suffix.lower() == ".catpart")
    jobs = [(str(args.root), str(p), str(args.out), args.ply) for p in files]
    rows = []
    with ProcessPoolExecutor() as ex, open(args.out / "results.jsonl", "w") as fh:
        for row in ex.map(evaluate, jobs, chunksize=4):
            rows.append(row)
            fh.write(json.dumps(row) + "\n")
    files_by_class = Counter(classify(r) for r in rows)
    bodies = [b for r in rows if "bodies" in r for b in r["bodies"]]
    print(f"files: {len(rows)}  {dict(files_by_class)}")
    by_kind = Counter(b["kind"] for b in bodies)
    print(f"bodies: {len(bodies)}  by kind (face solid flag): {dict(by_kind)}")
    solids = [b for b in bodies if b["kind"] == "solid"]
    print(f"solid bodies: closed {sum(b['closed'] for b in solids)}  "
          f"open but boundary covered within 2*sag {sum(1 for b in solids if b['boundary_edges'] and b['solid_boundary_gaps'] == 0)}  "
          f"open with gaps {sum(1 for b in solids if b['solid_boundary_gaps'])}  "
          f"non-manifold {sum(1 for b in solids if b['nonmanifold_edges'])}")
    nmp = Counter()
    for b in bodies:
        for k, v in (b["nonmanifold_patterns"] or {}).items():
            nmp[k] += v
    print(f"non-manifold edges by pattern (all bodies): {dict(nmp)}")
    surfaces = [b for b in bodies if b["kind"] == "surface"]
    print(f"surface bodies: open {sum(1 for b in surfaces if b['boundary_edges'])}  closed {sum(b['closed'] for b in surfaces)}")
    print(f"point lists without faces: {sum(len(r.get('point_lists_without_faces', [])) for r in rows)}  "
          f"face lists without resolved preamble: {sum(len(r.get('face_lists_without_preamble', [])) for r in rows)}  "
          f"symbol nodes: {sum(r.get('symbol_nodes', 0) for r in rows)}  "
          f"LOD meshes: {sum(b['lods'] for b in bodies)}")
    print(f"total triangles: {sum(r.get('triangles', 0) for r in rows)}  "
          f"slowest: {max(rows, key=lambda r: r['seconds'])['seconds']} s")
    for r in rows:
        if "error" in r:
            print("ERROR", r["file"], r["error"].strip().splitlines()[-1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
