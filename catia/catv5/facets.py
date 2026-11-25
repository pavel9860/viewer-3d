"""Facet records: the triangulated content of one B-rep face.

Grammar (all counts are P1g counts, see :mod:`catv5.binary`)::

    01  <marker>  [count per primitive kind present]  ff  u32le n
        [f32 * 3 normal]            when marker & 0x08
        [P1g * strip_count]         vertex count of each strip
        [P1g * fan_count]           vertex count of each fan
        n references                u8 or u16be element indices

    marker = 0x40 | kinds | normal
        0x01  independent triangles   (header count = triangles; 3 refs each, listed first)
        0x02  triangle strips         (header count = strips)
        0x04  triangle fans           (header count = fans)
        0x08  a 12-byte unit normal follows n (planar face)

    n == 3 * triangles + sum(strip counts) + sum(fan counts).
    A record with no kinds (marker 0x40/0x48) is an empty face with n == 0.

The reference width is not stored in the record. It is resolved by the
enclosing face list (:mod:`catv5.facelist`), which requires every record to end
exactly where the next structure begins.

Winding: strips alternate (even triangles a,b,c; odd triangles b,a,c), and fans
pivot on their first vertex. Planar faces wind counter-clockwise around their
stored normal.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass

from .binary import read_count, u16be_array
from .errors import FormatError

RECORD_TAG = 0x01
MARKER_BASE = 0x40
KIND_TRIANGLES = 0x01
KIND_STRIPS = 0x02
KIND_FANS = 0x04
HAS_NORMAL = 0x08
HEADER_END = 0xFF
MAX_REFS = 50_000_000

_NORMAL = struct.Struct("<3f")

Triangle = tuple[int, int, int]


@dataclass(frozen=True)
class FacetHeader:
    offset: int
    marker: int
    triangles: int
    strips: tuple[int, ...]
    fans: tuple[int, ...]
    normal: tuple[float, float, float] | None
    n: int
    refs_offset: int


@dataclass(frozen=True)
class FacetRecord:
    header: FacetHeader
    refs: tuple[int, ...]
    width: int  #: bytes per reference: 1 = u8, 2 = u16be, 3 = u24be, 4 = u32be
    end: int


COUNTS_P1G = "p1g"      # V5R25+ writers: loop counts are P1g counts
COUNTS_U16BE = "u16be"  # older writers (e.g. V5R17): loop counts are u16be


def _loop_count(data, p: int, limit: int, enc: str) -> tuple[int, int]:
    if enc == COUNTS_P1G:
        return read_count(data, p, limit)
    if p + 2 > limit:
        raise FormatError("loop count runs past the record", p)
    return (data[p] << 8) | data[p + 1], p + 2


def parse_header(data, pos: int, limit: int, counts: str = COUNTS_P1G) -> FacetHeader | None:
    """Parse a facet header at ``pos``. Returns ``None`` when the bytes are not one.

    ``counts`` selects the loop-count encoding; the face list resolves it (see facelist.py)."""
    if pos + 7 > limit or data[pos] != RECORD_TAG:
        return None
    marker = data[pos + 1]
    if marker & 0xF0 != MARKER_BASE:
        return None
    try:
        p = pos + 2
        kinds = {}
        for kind in (KIND_TRIANGLES, KIND_STRIPS, KIND_FANS):
            if marker & kind:
                kinds[kind], p = read_count(data, p, limit)
        if p + 5 > limit or data[p] != HEADER_END:
            return None
        n = int.from_bytes(data[p + 1 : p + 5], "little")
        p += 5
        if n > MAX_REFS:
            return None
        normal = None
        if marker & HAS_NORMAL:
            if p + 12 > limit:
                return None
            normal = _NORMAL.unpack_from(data, p)
            p += 12
        strips = []
        for _ in range(kinds.get(KIND_STRIPS, 0)):
            c, p = _loop_count(data, p, limit, counts)
            strips.append(c)
        fans = []
        for _ in range(kinds.get(KIND_FANS, 0)):
            c, p = _loop_count(data, p, limit, counts)
            fans.append(c)
    except FormatError:
        return None
    tri = kinds.get(KIND_TRIANGLES, 0)
    if 3 * tri + sum(strips) + sum(fans) != n:
        return None
    return FacetHeader(pos, marker, tri, tuple(strips), tuple(fans), normal, n, p)


def read_record(data, header: FacetHeader, width: int, limit: int, element_count: int) -> FacetRecord | None:
    """Read the references at the given width. ``None`` if they overrun or leave the element space."""
    end = header.refs_offset + width * header.n
    if end > limit:
        return None
    if width == 1:
        refs = tuple(data[header.refs_offset : end])
    elif width == 2:
        refs = u16be_array(data, header.refs_offset, header.n)
    else:  # u24be / u32be: element spaces above 65535
        b = bytes(data[header.refs_offset : end])
        refs = tuple(int.from_bytes(b[i : i + width], "big") for i in range(0, len(b), width))
    if refs and max(refs) >= element_count:
        return None
    return FacetRecord(header, refs, width, end)


def triangulate(record: FacetRecord) -> list[Triangle]:
    """Expand the record's primitives into triangles. Index-degenerate triangles
    (strip joins) are dropped."""
    h, refs = record.header, record.refs
    out: list[Triangle] = []
    k = 0
    for _ in range(h.triangles):
        out.append((refs[k], refs[k + 1], refs[k + 2]))
        k += 3
    for c in h.strips:
        for t in range(c - 2):
            a, b, d = refs[k + t], refs[k + t + 1], refs[k + t + 2]
            out.append((a, b, d) if t % 2 == 0 else (b, a, d))
        k += c
    for c in h.fans:
        pivot = refs[k]
        for t in range(1, c - 1):
            out.append((pivot, refs[k + t], refs[k + t + 1]))
        k += c
    if k != h.n:
        raise FormatError("primitive counts do not consume the references", h.offset)
    return [t for t in out if t[0] != t[1] and t[1] != t[2] and t[0] != t[2]]
