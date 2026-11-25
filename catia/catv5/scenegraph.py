"""Scene-graph scan: classify every node in the nested stores.

Nodes start with ``04 04 fd ff``. A node is a *point list* when it carries the
list token and its declared layout decodes exactly (:mod:`catv5.points`). Nodes
with the token whose declared data does not fit are *symbol nodes* (for
example the default reference-plane glyphs). They are reported, never silently
dropped.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .container import Store, find_nested_stores
from .points import NODE_MAGIC, NotAPointList, PointList, decode_point_list


@dataclass(frozen=True)
class SymbolNode:
    offset: int
    size: int
    reason: str


@dataclass
class SceneScan:
    stores: list[Store]
    point_lists: list[PointList] = field(default_factory=list)
    symbols: list[SymbolNode] = field(default_factory=list)
    #: per point list: exclusive limit for structures that belong to it
    limits: list[int] = field(default_factory=list)


def _innermost_store(stores: list[Store], offset: int) -> Store | None:
    best = None
    for s in stores:
        if s.body_lo <= offset < s.body_hi and (best is None or s.body_hi - s.body_lo < best.body_hi - best.body_lo):
            best = s
    return best


def scan(data) -> SceneScan:
    stores = find_nested_stores(data)
    result = SceneScan(stores)
    seen: set[int] = set()
    for store in stores:
        pos = store.body_lo
        while True:
            pos = data.find(NODE_MAGIC, pos, store.body_hi)
            if pos < 0:
                break
            if pos not in seen:
                seen.add(pos)
                inner = _innermost_store(stores, pos) or store
                nxt = data.find(NODE_MAGIC, pos + 4, inner.body_hi)
                node_end = inner.body_hi if nxt < 0 else nxt
                try:
                    pl = decode_point_list(data, pos, inner.body_hi)
                except NotAPointList as exc:
                    result.symbols.append(SymbolNode(pos, node_end - pos, str(exc)))
                else:
                    if pl is not None:
                        result.point_lists.append(pl)
            pos += len(NODE_MAGIC)
    result.point_lists.sort(key=lambda p: p.node_offset)
    for i, pl in enumerate(result.point_lists):
        inner = _innermost_store(stores, pl.node_offset)
        limit = inner.body_hi if inner else len(data)
        if i + 1 < len(result.point_lists):
            limit = min(limit, result.point_lists[i + 1].node_offset)
        result.limits.append(limit)
    return result
