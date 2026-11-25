"""Assemble decoded point lists and face lists into bodies and a part mesh."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from .facelist import FaceList, find_face_list
from .facets import Triangle, triangulate
from .points import PointList
from .scenegraph import SceneScan, SymbolNode, scan


@dataclass(frozen=True)
class Body:
    """One tessellated representation: a point list plus its faces."""

    points: PointList
    faces: FaceList
    triangles: tuple[Triangle, ...]
    #: half-open triangle index range of each face, parallel to ``faces.faces``
    face_ranges: tuple[tuple[int, int], ...]

    @property
    def kind(self) -> str:
        """``solid`` / ``surface`` from the faces' solid flag; ``mixed`` if they disagree; ``empty``."""
        flags = {f.solid for f in self.faces.faces if f.record.header.n}
        if not flags:
            return "empty"
        if flags == {True}:
            return "solid"
        if flags == {False}:
            return "surface"
        return "mixed"


@dataclass(frozen=True)
class TopologyStats:
    triangles: int
    edges: int
    boundary_edges: int
    nonmanifold_edges: int

    @property
    def closed(self) -> bool:
        return self.boundary_edges == 0 and self.nonmanifold_edges == 0


@dataclass
class Diagnostics:
    stores: int = 0
    symbol_nodes: list[SymbolNode] = field(default_factory=list)
    point_lists_without_faces: list[int] = field(default_factory=list)
    face_lists_without_preamble: list[int] = field(default_factory=list)


@dataclass
class PartMesh:
    bodies: list[Body]
    diagnostics: Diagnostics

    @property
    def triangle_count(self) -> int:
        return sum(len(b.triangles) for b in self.bodies)


def build_body(data, pl: PointList, limit: int) -> Body | None:
    fl = find_face_list(data, pl.end, limit, pl.element_count)
    if fl is None:
        return None
    tris: list[Triangle] = []
    ranges = []
    for face in fl.faces:
        start = len(tris)
        tris.extend(triangulate(face.record))
        ranges.append((start, len(tris)))
    return Body(pl, fl, tuple(tris), tuple(ranges))


def decode_part(data) -> PartMesh:
    """Decode every tessellated body of a ``.CATPart`` image (bytes or mmap)."""
    sg: SceneScan = scan(data)
    diag = Diagnostics(stores=len(sg.stores), symbol_nodes=list(sg.symbols))
    bodies: list[Body] = []
    for pl, limit in zip(sg.point_lists, sg.limits):
        body = build_body(data, pl, limit)
        if body is None:
            diag.point_lists_without_faces.append(pl.node_offset)
            continue
        if body.faces.sag is None:
            diag.face_lists_without_preamble.append(pl.node_offset)
        bodies.append(body)
    return PartMesh(bodies, diag)


def body_topology(body: Body) -> TopologyStats:
    """Edge-use statistics with elements welded by exact position.

    Welding is needed because one position appears as several elements, one per
    face normal. A correctly decoded solid body has every edge used exactly twice.
    """
    weld: dict[tuple[float, float, float], int] = {}
    idx = [weld.setdefault(p, len(weld)) for p in body.points.points]
    uses: Counter = Counter()
    n = 0
    for a, b, c in body.triangles:
        a, b, c = idx[a], idx[b], idx[c]
        if a == b or b == c or a == c:
            continue
        n += 1
        for u, v in ((a, b), (b, c), (c, a)):
            uses[(u, v) if u < v else (v, u)] += 1
    return TopologyStats(
        triangles=n,
        edges=len(uses),
        boundary_edges=sum(1 for k in uses.values() if k == 1),
        nonmanifold_edges=sum(1 for k in uses.values() if k > 2),
    )
