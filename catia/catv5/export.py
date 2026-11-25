"""Writers for decoded meshes."""

from __future__ import annotations

import struct
from typing import BinaryIO

from .mesh import PartMesh


def write_ply(mesh: PartMesh, fh: BinaryIO) -> None:
    """Binary little-endian PLY: all bodies concatenated, per-face RGBA."""
    verts = [p for b in mesh.bodies for p in b.points.points]
    faces = []
    base = 0
    for b in mesh.bodies:
        for face, (lo, hi) in zip(b.faces.faces, b.face_ranges):
            rgba = face.color
            for t in b.triangles[lo:hi]:
                faces.append((t[0] + base, t[1] + base, t[2] + base, rgba))
        base += len(b.points.points)
    header = (
        "ply\nformat binary_little_endian 1.0\ncomment catv5 decoder\n"
        f"element vertex {len(verts)}\nproperty float x\nproperty float y\nproperty float z\n"
        f"element face {len(faces)}\nproperty list uchar int vertex_indices\n"
        "property uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha\nend_header\n"
    )
    fh.write(header.encode("ascii"))
    fh.write(b"".join(struct.pack("<3f", *v) for v in verts))
    fh.write(b"".join(struct.pack("<B3i4B", 3, a, b, c, *rgba) for a, b, c, rgba in faces))
