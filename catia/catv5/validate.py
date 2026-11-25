"""Geometric checks on decoded bodies.

These are indicators, not ground truth. They test properties that a correctly
decoded solid must satisfy; passing them does not prove that nothing is missing.

``boundary_gaps``: CATIA tessellates each face independently, within the stored
chordal tolerance ``sag``. Neighbouring faces therefore need not share vertices:
- T-junctions occur where one face splits an edge that another does not;
- curved edges are sampled differently on each side.

For a watertight solid, every edge used by one triangle only must be covered,
along its length, by other such edges that lie within ``2 * sag`` of its line
(each side is within ``sag`` of the true curve).
"""

from __future__ import annotations

import math
from collections import Counter

from .mesh import Body

Point = tuple[float, float, float]


def boundary_segments(body: Body) -> list[tuple[Point, Point]]:
    weld: dict[Point, int] = {}
    idx = [weld.setdefault(p, len(weld)) for p in body.points.points]
    pos = {v: k for k, v in weld.items()}
    uses: Counter = Counter()
    for a, b, c in body.triangles:
        a, b, c = idx[a], idx[b], idx[c]
        if a == b or b == c or a == c:
            continue
        for u, v in ((a, b), (b, c), (c, a)):
            uses[(u, v) if u < v else (v, u)] += 1
    return [(pos[u], pos[v]) for (u, v), k in uses.items() if k == 1]


def _uncovered_length(seg, others, tol: float) -> float:
    a, b = seg
    d = [b[i] - a[i] for i in range(3)]
    length = math.sqrt(sum(x * x for x in d))
    if length <= tol:
        return 0.0
    u = [x / length for x in d]
    spans = []
    for c, e in others:
        ts = []
        for p in (c, e):
            w = [p[i] - a[i] for i in range(3)]
            t = sum(w[i] * u[i] for i in range(3))
            if sum(x * x for x in w) - t * t > tol * tol:
                break
            ts.append(t)
        else:
            lo, hi = min(ts), max(ts)
            if hi > 0 and lo < length:
                spans.append((max(lo, 0.0), min(hi, length)))
    spans.sort()
    covered_to, gap = 0.0, 0.0
    for lo, hi in spans:
        if lo > covered_to:
            gap += lo - covered_to
        covered_to = max(covered_to, hi)
    return gap + max(0.0, length - covered_to)


def boundary_gaps(body: Body, tol: float) -> list[tuple[Point, Point, float]]:
    """Boundary segments with more than ``tol`` of their length not covered by others."""
    segs = boundary_segments(body)
    boxes = [tuple(min(s[0][i], s[1][i]) - tol for i in range(3)) + tuple(max(s[0][i], s[1][i]) + tol for i in range(3)) for s in segs]
    out = []
    for i, seg in enumerate(segs):
        bi = boxes[i]
        others = [
            segs[j]
            for j in range(len(segs))
            if j != i and all(boxes[j][k] <= bi[k + 3] and bi[k] <= boxes[j][k + 3] for k in range(3))
        ]
        gap = _uncovered_length(seg, others, tol)
        if gap > tol:
            out.append((seg[0], seg[1], gap))
    return out


def nonmanifold_patterns(body: Body) -> dict[str, int]:
    """Classify edges used by more than two triangles (positions welded).

    - ``fold``: the extra uses come from a triangle and its exact reverse within
      one face. CATIA emits these where a fan's last rim vertex coincides with its
      first.
    - ``touching``: an even number of uses, half in each direction. These are
      consistently oriented sheets meeting along the edge, i.e. lumps of one solid
      that touch.
    - ``other``: anything else.
    """
    weld: dict[Point, int] = {}
    idx = [weld.setdefault(p, len(weld)) for p in body.points.points]
    tris = []
    for fi, (lo, hi) in enumerate(body.face_ranges):
        for t in body.triangles[lo:hi]:
            w = (idx[t[0]], idx[t[1]], idx[t[2]])
            if len(set(w)) == 3:
                tris.append((w, fi))
    by_set: dict[frozenset, list] = {}
    for w, fi in tris:
        by_set.setdefault(frozenset(w), []).append((w, fi))
    folds = set()
    for group in by_set.values():
        if len(group) == 2 and group[0][1] == group[1][1]:
            a, c = group[0][0], group[1][0]
            if (a[0], a[2], a[1]) in {(c[i], c[(i + 1) % 3], c[(i + 2) % 3]) for i in range(3)}:
                folds.add(frozenset(a))
    def directed_uses(selected):
        uses: dict[tuple[int, int], list[int]] = {}
        for w in selected:
            for u, v in ((w[0], w[1]), (w[1], w[2]), (w[2], w[0])):
                uses.setdefault((min(u, v), max(u, v)), [0, 0])[0 if u < v else 1] += 1
        return uses

    all_uses = directed_uses([w for w, _ in tris])
    kept_uses = directed_uses([w for w, _ in tris if frozenset(w) not in folds])
    out = {"fold": 0, "touching": 0, "other": 0}
    for k, (f, r) in all_uses.items():
        if f + r <= 2:
            continue
        kf, kr = kept_uses.get(k, [0, 0])
        if kf + kr <= 2:
            out["fold"] += 1
        elif kf == kr:
            out["touching"] += 1
        else:
            out["other"] += 1
    return out
