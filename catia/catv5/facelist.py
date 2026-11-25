"""The face list that follows a point list, with its level-of-detail preamble.

Layout after the point list's per-element normal block (whose layouts vary and
are not walked here)::

    f32  sag                       chordal tolerance used to tessellate the faces (mm)
    u8   lod_count
    lod_count * { facet record     coarse whole-body mesh
                  f32  lod_sag
                  attribute }      8 bytes, see below
    u8   list_kind
    02 0a
    P1g  face_count
    face_count facet records, back to back       (see :mod:`catv5.facets`)
    face_count attribute entries:
        <flags> 04 04 ff  A B G R

    flags bit 0x10: the face belongs to a solid. Surface faces leave it clear.
    flags bit 0x80: the face is transparent (exactly when A < 255).

The face list is identified exactly. A candidate ``02 0a <count>`` is accepted
only when:
- exactly ``face_count`` records parse back to back;
- every reference lies inside the element space;
- the byte after the last record starts an attribute entry.

Each record's reference width (u8 or u16be) is the unique width under which the
record ends exactly at the next structure. The preamble is accepted only when
its LOD chain ends exactly at ``list_kind``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .binary import f32, read_count
from .errors import FormatError
from .facets import COUNTS_P1G, COUNTS_U16BE, FacetRecord, parse_header, read_record

LIST_MARKER = b"\x02\x0a"
ATTR_SIG = b"\x04\x04\xff"
ATTR_LEN = 8
MAX_FACES = 10_000_000
MAX_LODS = 32
FLAG_SOLID = 0x10
FLAG_TRANSPARENT = 0x80

RGBA = tuple[int, int, int, int]


@dataclass(frozen=True)
class Face:
    record: FacetRecord
    color: RGBA  #: (r, g, b, a) 0..255
    attr_flags: int

    @property
    def solid(self) -> bool:
        return bool(self.attr_flags & FLAG_SOLID)

    @property
    def transparent(self) -> bool:
        return bool(self.attr_flags & FLAG_TRANSPARENT)


@dataclass(frozen=True)
class Lod:
    record: FacetRecord
    sag: float
    color: RGBA
    attr_flags: int


@dataclass(frozen=True)
class FaceList:
    marker_offset: int
    list_kind: int
    sag: float | None  #: None when the preamble did not resolve
    lods: tuple[Lod, ...]
    faces: tuple[Face, ...]
    end: int  #: absolute offset after the attribute entries


def _is_attr(data, pos: int, limit: int) -> bool:
    return pos + ATTR_LEN <= limit and bytes(data[pos + 1 : pos + 4]) == ATTR_SIG


def _attr(data, pos: int) -> tuple[int, RGBA]:
    a, b, g, r = data[pos + 4 : pos + 8]
    return data[pos], (r, g, b, a)


def _plausible_sag(v: float) -> bool:
    return math.isfinite(v) and 0.0 < v < 1e4


def min_width(element_count: int) -> int:
    """Narrowest big-endian reference width able to address the element space."""
    return 1 if element_count <= 0x100 else 2 if element_count <= 0x10000 else 3 if element_count <= 0x1000000 else 4


def parse_records(data, pos: int, limit: int, element_count: int, count: int, counts: str = COUNTS_P1G) -> tuple[list[FacetRecord], int] | None:
    """Parse ``count`` back-to-back records starting at ``pos`` (see module docstring)."""
    records: list[FacetRecord] = []
    for k in range(count):
        header = parse_header(data, pos, limit, counts)
        if header is None:
            return None
        last = k + 1 == count
        options = []
        for width in ((1,) if header.n == 0 else range(min_width(element_count), 5)):
            rec = read_record(data, header, width, limit, element_count)
            if rec is None:
                continue
            nxt_ok = _is_attr(data, rec.end, limit) if last else parse_header(data, rec.end, limit, counts) is not None
            if nxt_ok:
                options.append(rec)
        if len(options) != 1:
            return None
        records.append(options[0])
        pos = options[0].end
    return records, pos


def _parse_lod_chain(data, pos: int, count: int, stop: int, element_count: int, counts: str = COUNTS_P1G) -> list[Lod] | None:
    lods = []
    for _ in range(count):
        header = parse_header(data, pos, stop, counts)
        if header is None:
            return None
        chosen = None
        for width in ((1,) if header.n == 0 else range(min_width(element_count), 5)):
            rec = read_record(data, header, width, stop, element_count)
            if rec is not None and _is_attr(data, rec.end + 4, stop + 1):
                if chosen is not None:
                    return None
                chosen = rec
        if chosen is None:
            return None
        flags, rgba = _attr(data, chosen.end + 4)
        lods.append(Lod(chosen, f32(data, chosen.end), rgba, flags))
        pos = chosen.end + 4 + ATTR_LEN
    return lods if pos == stop else None


def parse_preamble(data, lo: int, kind_pos: int, element_count: int) -> tuple[float | None, tuple[Lod, ...]]:
    """Resolve ``sag`` and the LOD meshes that end exactly at ``kind_pos`` (the ``list_kind`` byte)."""
    for q in range(lo, kind_pos - 5):
        n = data[q + 4]
        if 1 <= n <= MAX_LODS and data[q + 5] == 0x01 and _plausible_sag(f32(data, q)):
            lods = _parse_lod_chain(data, q + 5, n, kind_pos, element_count) or _parse_lod_chain(data, q + 5, n, kind_pos, element_count, COUNTS_U16BE)
            if lods is not None:
                return f32(data, q), tuple(lods)
    if kind_pos - 5 >= lo and data[kind_pos - 1] == 0 and _plausible_sag(f32(data, kind_pos - 5)):
        return f32(data, kind_pos - 5), ()
    return None, ()


def _resolve_counts(data, pos, limit, element_count, count):
    """Parse under both loop-count encodings. Lists without strips/fans read identically under
    both; a list that parses differently under both is a conflict and raises."""
    a = parse_records(data, pos, limit, element_count, count, COUNTS_P1G)
    b = parse_records(data, pos, limit, element_count, count, COUNTS_U16BE)
    if a is not None and b is not None and (a[1] != b[1] or [(r.refs, r.header.strips, r.header.fans) for r in a[0]] != [(r.refs, r.header.strips, r.header.fans) for r in b[0]]):
        raise FormatError("face list parses differently under both loop-count encodings", pos)
    return a if a is not None else b


def find_face_list(data, start: int, limit: int, element_count: int) -> FaceList | None:
    """Locate and parse the face list in ``[start, limit)``. ``None`` when there is none."""
    pos = start
    while True:
        pos = data.find(LIST_MARKER, pos, limit)
        if pos < 0:
            return None
        try:
            count, rec_pos = read_count(data, pos + 2, limit)
        except Exception:
            count = 0
        if 1 <= count <= MAX_FACES:
            parsed = _resolve_counts(data, rec_pos, limit, element_count, count)
            if parsed is not None:
                records, attr_pos = parsed
                faces = []
                for rec in records:
                    if not _is_attr(data, attr_pos, limit):
                        break
                    flags, rgba = _attr(data, attr_pos)
                    faces.append(Face(rec, rgba, flags))
                    attr_pos += ATTR_LEN
                if len(faces) == count:
                    sag, lods = parse_preamble(data, start, pos - 1, element_count)
                    return FaceList(pos, data[pos - 1], sag, lods, tuple(faces), attr_pos)
        pos += 1
