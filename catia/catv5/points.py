"""Point lists ("element spaces") of the tessellation scene graph.

A point-list node starts with ``04 04 fd ff``; within its first 64 bytes sits the
list token ``ff ff 02 00 01``, followed by::

    P1g  element_count          number of elements (vertices) faces may reference
    P1g  edge_bound             (bound used by the edge tables; not needed for faces)
    u8   has_texcoords
    u8   reserved_a
    u8   reserved_b
    [u8 tex_dims, f32 * element_count * tex_dims]     only when has_texcoords
    u8   encoding               0 = prefix-compressed table, 1 = raw xyz

Encoding 0 (prefix compression over the lexicographically sorted list)::

    ceil(element_count / 4) bytes   2-bit opcode per element, LSB first
    P1g                             number of f32 words that follow
    f32 * words

    opcode 0: x, y, z written      opcode 3: y, z written (x repeats)
    opcode 2: z written            opcode 1: nothing written (x, y, z repeat)

    The state starts at the origin. Opcode 1 therefore yields an element at the
    same position as its predecessor: the same vertex position used with a
    different normal on another face.

Encoding 1: ``f32 * 3 * element_count`` (x, y, z per element).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .binary import f32_array, read_count
from .errors import FormatError

NODE_MAGIC = b"\x04\x04\xfd\xff"
LIST_TOKEN = b"\xff\xff\x02\x00\x01"
TOKEN_WINDOW = 64
ENC_TABLE = 0
ENC_RAW = 1
_OPCODE_WORDS = (3, 0, 1, 2)  # f32 words consumed by opcodes 0..3
MAX_ELEMENTS = 50_000_000
MAX_ABS_COORD = 1e7

Point = tuple[float, float, float]


@dataclass(frozen=True)
class PointList:
    node_offset: int  #: absolute offset of the node magic
    element_count: int
    edge_bound: int
    encoding: int
    points: tuple[Point, ...]
    texcoords: tuple[float, ...] | None
    tex_dims: int
    end: int  #: absolute offset just past the point data


class NotAPointList(Exception):
    """The node carries the list token but its declared data does not fit the node.

    Raised for the small symbol nodes (reference-plane glyphs and similar) that
    reuse the token; they are classified, not decoded.
    """


def decode_point_list(data, node_offset: int, node_end: int) -> PointList | None:
    """Decode the point list of the node at ``node_offset``.

    Returns ``None`` when the node has no list token. Raises :class:`NotAPointList`
    when the token is present but the declared layout does not fit inside the node.
    """
    base = node_offset + len(NODE_MAGIC)
    tok = data.find(LIST_TOKEN, base, min(base + TOKEN_WINDOW, node_end))
    if tok < 0:
        return None
    try:
        return _decode(data, node_offset, tok + len(LIST_TOKEN), node_end)
    except FormatError as exc:
        raise NotAPointList(str(exc)) from None


def _decode(data, node_offset: int, pos: int, end: int) -> PointList:
    ec, pos = read_count(data, pos, end)
    edge_bound, pos = read_count(data, pos, end)
    if not 0 < ec <= MAX_ELEMENTS:
        raise FormatError(f"implausible element count {ec}", pos)
    if pos + 4 > end:
        raise FormatError("header runs past the node", pos)
    has_tex = data[pos]
    if has_tex not in (0, 1):
        raise FormatError(f"texcoord flag {has_tex:#x} is not 0/1", pos)
    pos += 3
    tex: tuple[float, ...] | None = None
    dims = 0
    if has_tex:
        dims = data[pos]
        pos += 1
        if not 1 <= dims <= 4 or pos + 4 * dims * ec > end:
            raise FormatError(f"texcoord block (dims={dims}) does not fit the node", pos)
        tex = f32_array(data, pos, dims * ec)
        pos += 4 * dims * ec
    if pos >= end:
        raise FormatError("encoding byte runs past the node", pos)
    encoding = data[pos]
    pos += 1
    if encoding == ENC_TABLE:
        points, pos = _decode_table(data, pos, end, ec)
    elif encoding == ENC_RAW:
        if pos + 12 * ec > end:
            raise FormatError("raw point block does not fit the node", pos)
        flat = f32_array(data, pos, 3 * ec)
        points = tuple(flat[i : i + 3] for i in range(0, 3 * ec, 3))
        pos += 12 * ec
    else:
        raise FormatError(f"unknown point encoding {encoding:#x}", pos - 1)
    for p in points:
        if not all(math.isfinite(c) and abs(c) < MAX_ABS_COORD for c in p):
            raise FormatError(f"implausible coordinate {p}", pos)
    return PointList(node_offset, ec, edge_bound, encoding, points, tex, dims, pos)


def _decode_table(data, pos: int, end: int, ec: int):
    table_len = (ec + 3) // 4
    if pos + table_len >= end:
        raise FormatError("opcode table runs past the node", pos)
    table = data[pos : pos + table_len]
    pos += table_len
    words, pos = read_count(data, pos, end)
    opcodes = [(table[i >> 2] >> ((i & 3) * 2)) & 3 for i in range(ec)]
    need = sum(_OPCODE_WORDS[o] for o in opcodes)
    if need != words:
        raise FormatError(f"opcode table needs {need} words, header declares {words}", pos)
    if pos + 4 * words > end:
        raise FormatError("point words run past the node", pos)
    vals = f32_array(data, pos, words)
    x = y = z = 0.0
    k = 0
    points = []
    for op in opcodes:
        if op == 0:
            x, y, z = vals[k], vals[k + 1], vals[k + 2]
            k += 3
        elif op == 3:
            y, z = vals[k], vals[k + 1]
            k += 2
        elif op == 2:
            z = vals[k]
            k += 1
        points.append((x, y, z))
    return tuple(points), pos + 4 * words
