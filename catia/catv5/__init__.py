"""catv5 — reader for the tessellation stored in CATIA V5 ``.CATPart`` files.

Layers, bottom to top:

- :mod:`.binary`      variable-length counts, f32/u16 readers
- :mod:`.container`   outer ``V5_CFV2`` container and nested stores
- :mod:`.points`      point lists (element spaces) in three encodings
- :mod:`.scenegraph`  node scan: point lists vs symbol nodes
- :mod:`.facets`      facet-record grammar and triangulation
- :mod:`.facelist`    face lists and per-face colours
- :mod:`.mesh`        bodies, part mesh, topology statistics
- :mod:`.export`      PLY writer
"""

from .errors import CatV5Error, FormatError, NotCatiaFile
from .mesh import Body, PartMesh, TopologyStats, body_topology, decode_part

__all__ = [
    "Body",
    "CatV5Error",
    "FormatError",
    "NotCatiaFile",
    "PartMesh",
    "TopologyStats",
    "body_topology",
    "decode_part",
]
