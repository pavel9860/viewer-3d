"""Render PLY files (from ``catv5 -o`` or ``corpus_eval --ply``) into one PNG contact sheet.

Pure Python z-buffer renderer, isometric view, per-face colours when present.

Usage::

    python3 -I tools/render_sheet.py OUT.png CELL_PX COLUMNS file1.ply file2.ply ...
"""

from __future__ import annotations

import math
import struct
import sys
import zlib

LIGHT = (0.4, 0.5, -0.77)
BG = 245


def read_ply(path):
    data = open(path, "rb").read()
    end = data.index(b"end_header\n") + len(b"end_header\n")
    header = data[:end].decode("ascii")
    nv = int(header.split("element vertex ")[1].split()[0])
    nf = int(header.split("element face ")[1].split()[0])
    vprops = header.split("element face")[0].count("property float")
    has_color = "property uchar red" in header
    pos = end
    verts = []
    for _ in range(nv):
        verts.append(struct.unpack_from("<3f", data, pos))
        pos += 4 * vprops
    faces = []
    for _ in range(nf):
        k = data[pos]
        idx = struct.unpack_from(f"<{k}i", data, pos + 1)
        pos += 1 + 4 * k
        rgb = (180, 190, 230)
        if has_color:
            rgb = tuple(data[pos : pos + 3])
            pos += 4
        for j in range(1, k - 1):
            faces.append((idx[0], idx[j], idx[j + 1], rgb))
    return verts, faces


def render(verts, faces, size):
    img = bytearray([BG] * (size * size * 3))
    zbuf = [math.inf] * (size * size)
    if not verts or not faces:
        return img
    ay, ax = math.radians(-35), math.radians(30)
    cy, sy, cx, sx = math.cos(ay), math.sin(ay), math.cos(ax), math.sin(ax)

    def project(v):
        x, y, z = v
        x, y = x * cy - y * sy, x * sy + y * cy
        y, z = y * cx - z * sx, y * sx + z * cx
        return x, z, y

    used = sorted({i for f in faces for i in f[:3]})
    proj = {i: project(verts[i]) for i in used}
    xs = [p[0] for p in proj.values()]
    ys = [p[1] for p in proj.values()]
    scale = 0.9 * size / max(max(xs) - min(xs), max(ys) - min(ys), 1e-9)
    mx, my = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2
    scr = {i: ((p[0] - mx) * scale + size / 2, size / 2 - (p[1] - my) * scale, p[2]) for i, p in proj.items()}
    for a, b, c, rgb in faces:
        A, B, C = scr[a], scr[b], scr[c]
        ux, uy, uz = B[0] - A[0], B[1] - A[1], B[2] - A[2]
        vx, vy, vz = C[0] - A[0], C[1] - A[1], C[2] - A[2]
        nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        nl = math.sqrt(nx * nx + ny * ny + nz * nz) or 1.0
        shade = 0.35 + 0.65 * abs(nx * LIGHT[0] + ny * LIGHT[1] + nz * LIGHT[2]) / nl
        col = bytes(min(255, int(ch * shade)) for ch in rgb)
        den = (B[1] - C[1]) * (A[0] - C[0]) + (C[0] - B[0]) * (A[1] - C[1])
        if abs(den) < 1e-12:
            continue
        x0, x1 = max(0, int(min(A[0], B[0], C[0]))), min(size - 1, int(max(A[0], B[0], C[0])) + 1)
        y0, y1 = max(0, int(min(A[1], B[1], C[1]))), min(size - 1, int(max(A[1], B[1], C[1])) + 1)
        for y in range(y0, y1 + 1):
            for x in range(x0, x1 + 1):
                w0 = ((B[1] - C[1]) * (x - C[0]) + (C[0] - B[0]) * (y - C[1])) / den
                w1 = ((C[1] - A[1]) * (x - C[0]) + (A[0] - C[0]) * (y - C[1])) / den
                w2 = 1 - w0 - w1
                if w0 < -1e-6 or w1 < -1e-6 or w2 < -1e-6:
                    continue
                z = w0 * A[2] + w1 * B[2] + w2 * C[2]
                i = y * size + x
                if z < zbuf[i]:
                    zbuf[i] = z
                    img[3 * i : 3 * i + 3] = col
    return img


def write_png(path, width, height, rgb):
    raw = b"".join(b"\x00" + bytes(rgb[y * width * 3 : (y + 1) * width * 3]) for y in range(height))

    def chunk(tag, body):
        return struct.pack(">I", len(body)) + tag + body + struct.pack(">I", zlib.crc32(tag + body) & 0xFFFFFFFF)

    with open(path, "wb") as fh:
        fh.write(b"\x89PNG\r\n\x1a\n")
        fh.write(chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)))
        fh.write(chunk(b"IDAT", zlib.compress(raw, 6)))
        fh.write(chunk(b"IEND", b""))


def main(argv):
    out, size, cols, paths = argv[0], int(argv[1]), int(argv[2]), argv[3:]
    rows = (len(paths) + cols - 1) // cols
    width, height = cols * size, rows * size
    sheet = bytearray([255] * (width * height * 3))
    for k, path in enumerate(paths):
        img = render(*read_ply(path), size)
        r, c = divmod(k, cols)
        for y in range(size):
            o = ((r * size + y) * width + c * size) * 3
            sheet[o : o + size * 3] = img[y * size * 3 : (y + 1) * size * 3]
        # 1-px separator
        for y in range(size):
            o = ((r * size + y) * width + c * size) * 3
            sheet[o : o + 3] = b"\xc8\xc8\xc8"
    write_png(out, width, height, sheet)


if __name__ == "__main__":
    main(sys.argv[1:])
