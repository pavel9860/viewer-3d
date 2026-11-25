"""Command line: ``python3 -m catv5 FILE.CATPart [-o OUT.ply] [--json]``."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import CatV5Error, body_topology, decode_part
from .validate import boundary_gaps, nonmanifold_patterns
from .export import write_ply


def summarize(mesh) -> dict:
    bodies = []
    for b in mesh.bodies:
        t = body_topology(b)
        sag = b.faces.sag
        gaps = None
        if b.kind == "solid" and t.boundary_edges and sag is not None:
            gaps = len(boundary_gaps(b, 2 * sag))
        bodies.append(
            {
                "node_offset": b.points.node_offset,
                "kind": b.kind,
                "sag": sag,
                "lods": len(b.faces.lods),
                "transparent_faces": sum(f.transparent for f in b.faces.faces),
                "elements": b.points.element_count,
                "encoding": b.points.encoding,
                "faces": len(b.faces.faces),
                "triangles": t.triangles,
                "boundary_edges": t.boundary_edges,
                "nonmanifold_edges": t.nonmanifold_edges,
                "closed": t.closed,
                "solid_boundary_gaps": gaps,
                "nonmanifold_patterns": nonmanifold_patterns(b) if t.nonmanifold_edges else None,
            }
        )
    d = mesh.diagnostics
    return {
        "bodies": bodies,
        "triangles": mesh.triangle_count,
        "stores": d.stores,
        "symbol_nodes": len(d.symbol_nodes),
        "point_lists_without_faces": d.point_lists_without_faces,
        "face_lists_without_preamble": d.face_lists_without_preamble,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="catv5", description="Decode the tessellation of a CATIA V5 .CATPart")
    ap.add_argument("file", type=Path)
    ap.add_argument("-o", "--output", type=Path, help="write a binary PLY with per-face colours")
    ap.add_argument("--json", action="store_true", help="print a JSON summary")
    args = ap.parse_args(argv)
    try:
        mesh = decode_part(args.file.read_bytes())
    except CatV5Error as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    if args.output:
        with open(args.output, "wb") as fh:
            write_ply(mesh, fh)
    s = summarize(mesh)
    if args.json:
        print(json.dumps(s, indent=2))
    else:
        closed = sum(b["closed"] for b in s["bodies"])
        print(f"{args.file.name}: {len(s['bodies'])} bodies ({closed} closed), {s['triangles']} triangles")
    return 0


if __name__ == "__main__":
    sys.exit(main())
