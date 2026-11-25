"""Outer ``V5_CFV2`` container and the nested stores embedded in it.

File layout (big-endian header fields)::

    [0:8)   b"V5_CFV2\\0"
    [8:12)  u32be  catalog offset (relative to the container origin)
    [12:16) u32be  catalog length
    [16:24) ff * 8

A ``.CATPart`` embeds further ``V5_CFV2`` containers ("nested stores") inside
its storage blocks. The display tessellation lives in their bodies, i.e. the
bytes between the 32-byte header and the catalog: ``[origin + 32, origin + catalog_offset)``.
"""

from __future__ import annotations

from dataclasses import dataclass

from .binary import u32be
from .errors import NotCatiaFile

MAGIC = b"V5_CFV2\x00"
HEADER_LEN = 32
_FF8 = b"\xff" * 8


@dataclass(frozen=True)
class Store:
    """One nested ``V5_CFV2`` store."""

    origin: int  #: absolute offset of the store's magic
    body_lo: int  #: first body byte (origin + 32)
    body_hi: int  #: catalog start (origin + catalog offset), exclusive end of the body
    catalog_len: int


def check_magic(data) -> None:
    if bytes(data[:8]) != MAGIC:
        raise NotCatiaFile("not a CATIA V5 file: missing V5_CFV2 magic")


def find_nested_stores(data) -> list[Store]:
    """Every nested store whose header is self-consistent, in file order."""
    check_magic(data)
    size = len(data)
    stores: list[Store] = []
    pos = 8
    while True:
        pos = data.find(MAGIC, pos)
        if pos < 0:
            return stores
        if pos + HEADER_LEN <= size and bytes(data[pos + 16 : pos + 24]) == _FF8:
            cat_off = u32be(data, pos + 8)
            cat_len = u32be(data, pos + 12)
            if HEADER_LEN < cat_off and pos + cat_off + cat_len <= size:
                stores.append(Store(pos, pos + HEADER_LEN, pos + cat_off, cat_len))
        pos += len(MAGIC)
