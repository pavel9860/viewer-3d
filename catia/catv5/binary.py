"""Low-level binary readers shared by every layer.

CATIA's scene-graph streams use one variable-length integer encoding everywhere
("P1g count"): a single byte ``v``; when ``v == 0xFF`` the value follows as a
little-endian ``u32``.
"""

from __future__ import annotations

import struct

from .errors import FormatError

_F32 = struct.Struct("<f")
_F32x3 = struct.Struct("<3f")
_U32LE = struct.Struct("<I")
_U32BE = struct.Struct(">I")

COUNT_ESCAPE = 0xFF


def read_count(buf, pos: int, limit: int | None = None) -> tuple[int, int]:
    """Read a P1g count at ``pos``. Returns ``(value, next_pos)``."""
    end = len(buf) if limit is None else limit
    if pos >= end:
        raise FormatError("count runs past the buffer", pos)
    v = buf[pos]
    if v != COUNT_ESCAPE:
        return v, pos + 1
    if pos + 5 > end:
        raise FormatError("escaped count runs past the buffer", pos)
    return _U32LE.unpack_from(buf, pos + 1)[0], pos + 5


def u32le(buf, pos: int) -> int:
    return _U32LE.unpack_from(buf, pos)[0]


def u32be(buf, pos: int) -> int:
    return _U32BE.unpack_from(buf, pos)[0]


def f32(buf, pos: int) -> float:
    return _F32.unpack_from(buf, pos)[0]


def f32x3(buf, pos: int) -> tuple[float, float, float]:
    return _F32x3.unpack_from(buf, pos)


def f32_array(buf, pos: int, count: int) -> tuple[float, ...]:
    return struct.unpack_from(f"<{count}f", buf, pos)


def u16be_array(buf, pos: int, count: int) -> tuple[int, ...]:
    return struct.unpack_from(f">{count}H", buf, pos)
