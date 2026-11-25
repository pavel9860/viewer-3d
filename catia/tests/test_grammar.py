"""Unit tests for the byte grammars, on synthetic records (no sample files needed).

Run: python3 -I -m unittest discover -s tests -t .   (from catia/)
"""

import struct
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from catv5.binary import read_count  # noqa: E402
from catv5.facelist import find_face_list  # noqa: E402
from catv5.facets import parse_header, read_record, triangulate  # noqa: E402
from catv5.points import NODE_MAGIC, LIST_TOKEN, NotAPointList, decode_point_list  # noqa: E402


def facet(marker, counts, n, normal=None, loop_counts=(), refs=(), width=1):
    b = bytes([0x01, marker]) + bytes(counts) + b"\xff" + struct.pack("<I", n)
    if normal:
        b += struct.pack("<3f", *normal)
    b += bytes(loop_counts)
    b += bytes(refs) if width == 1 else struct.pack(f">{len(refs)}H", *refs)
    return b


ATTR_GREY = b"\x30\x04\x04\xff\xff\x5f\x5f\x5f"  # solid, opaque grey
ATTR_SURF = b"\xa0\x04\x04\xff\xc6\xd6\xaa\x83"  # surface, transparent (131,170,214,198)


class CountTests(unittest.TestCase):
    def test_plain_and_escaped(self):
        self.assertEqual(read_count(b"\x07", 0), (7, 1))
        self.assertEqual(read_count(b"\xff\xea\x01\x00\x00", 0), (490, 5))


class FacetTests(unittest.TestCase):
    def test_strip_with_normal(self):
        rec = facet(0x4A, [1], 4, (0, 0, 1), [4], [0, 1, 2, 3])
        h = parse_header(rec, 0, len(rec))
        self.assertEqual((h.strips, h.fans, h.triangles), ((4,), (), 0))
        r = read_record(rec, h, 1, len(rec), 4)
        self.assertEqual(triangulate(r), [(0, 1, 2), (2, 1, 3)])

    def test_mixed_triangles_and_strips(self):  # 0x4b: tri + strip + normal
        refs = [0, 1, 2, 3, 4, 5, 6]
        rec = facet(0x4B, [1, 1], 7, (0, 0, 1), [4], refs)
        h = parse_header(rec, 0, len(rec))
        self.assertEqual((h.triangles, h.strips), (1, (4,)))
        tris = triangulate(read_record(rec, h, 1, len(rec), 7))
        self.assertEqual(tris, [(0, 1, 2), (3, 4, 5), (5, 4, 6)])

    def test_strip_and_fan(self):  # 0x4e: strip + fan + normal; strips' refs come first
        rec = facet(0x4E, [1, 1], 9, (0, 0, -1), [4, 5], list(range(9)))
        tris = triangulate(read_record(rec, parse_header(rec, 0, len(rec)), 1, len(rec), 9))
        self.assertEqual(tris, [(0, 1, 2), (2, 1, 3), (4, 5, 6), (4, 6, 7), (4, 7, 8)])

    def test_escaped_header_count(self):
        rec = bytes([0x01, 0x41, 0xFF]) + struct.pack("<I", 2) + b"\xff" + struct.pack("<I", 6) + bytes(range(6))
        h = parse_header(rec, 0, len(rec))
        self.assertEqual(h.triangles, 2)

    def test_empty_face(self):
        rec = b"\x01\x40\xff\x00\x00\x00\x00"
        h = parse_header(rec, 0, len(rec))
        self.assertEqual(h.n, 0)

    def test_count_mismatch_rejected(self):
        rec = facet(0x4A, [1], 5, (0, 0, 1), [4], [0, 1, 2, 3, 4])
        self.assertIsNone(parse_header(rec, 0, len(rec)))

    def test_u16_refs(self):
        rec = facet(0x42, [1], 3, None, [3], [256, 300, 1], width=2)
        r = read_record(rec, parse_header(rec, 0, len(rec)), 2, len(rec), 301)
        self.assertEqual(r.refs, (256, 300, 1))


class FaceListTests(unittest.TestCase):
    def test_lod_before_list_is_not_a_face(self):
        lod = facet(0x43, [1, 0], 3, None, [], [0, 1, 2]) + struct.pack("<f", 3.3) + ATTR_GREY
        f1 = facet(0x4A, [1], 4, (0, 0, 1), [4], [0, 1, 2, 3])
        f2 = facet(0x40, [], 0)
        blob = b"\x00\x01" + struct.pack("<f", 0.2) + b"\x01" + lod + b"\x04\x02\x0a\x02" + f1 + f2 + ATTR_GREY + ATTR_SURF + b"\x01\x01\x00"
        fl = find_face_list(blob, 0, len(blob), 4)
        self.assertEqual(len(fl.faces), 2)
        self.assertEqual(fl.faces[0].record.header.marker, 0x4A)
        self.assertEqual(fl.faces[0].color, (0x5F, 0x5F, 0x5F, 0xFF))
        self.assertAlmostEqual(fl.sag, 0.2, places=6)
        self.assertEqual(len(fl.lods), 1)
        self.assertAlmostEqual(fl.lods[0].sag, 3.3, places=5)
        self.assertTrue(fl.faces[0].solid)
        self.assertFalse(fl.faces[1].solid)
        self.assertTrue(fl.faces[1].transparent)
        self.assertEqual(fl.faces[1].color, (0x83, 0xAA, 0xD6, 0xC6))

    def test_preamble_without_lods(self):
        f1 = facet(0x4A, [1], 4, (0, 0, 1), [4], [0, 1, 2, 3])
        blob = b"\x33" + struct.pack("<f", 0.05) + b"\x00\x03\x02\x0a\x01" + f1 + ATTR_GREY
        fl = find_face_list(blob, 0, len(blob), 4)
        self.assertAlmostEqual(fl.sag, 0.05, places=6)
        self.assertEqual((fl.lods, fl.list_kind), ((), 3))


class PointTests(unittest.TestCase):
    def node(self, body: bytes) -> bytes:
        return NODE_MAGIC + b"\x00" * 8 + LIST_TOKEN + body

    def test_table_encoding_prefix_law(self):
        # opcodes: 0 (xyz), 1 (repeat), 2 (z), 3 (y z) -> table byte 0b11_10_01_00
        words = struct.pack("<6f", 1, 2, 3, 4, 5, 6)
        body = b"\x04\x00" + b"\x00\x00\x00" + b"\x00" + bytes([0b11100100]) + b"\x06" + words
        data = self.node(body)
        pl = decode_point_list(data, 0, len(data))
        self.assertEqual(pl.points, ((1, 2, 3), (1, 2, 3), (1, 2, 4), (1, 5, 6)))

    def test_raw_encoding(self):
        body = b"\x02\x02" + b"\x00\x00\x00" + b"\x01" + struct.pack("<6f", 1, 2, 3, 4, 5, 6)
        data = self.node(body)
        self.assertEqual(decode_point_list(data, 0, len(data)).points, ((1, 2, 3), (4, 5, 6)))

    def test_texcoord_block(self):
        body = (b"\x01\x01" + b"\x01\x00\x01" + b"\x02" + struct.pack("<2f", 0.5, 0.25)
                + b"\x00" + b"\x00" + b"\x03" + struct.pack("<3f", 7, 8, 9))
        data = self.node(body)
        pl = decode_point_list(data, 0, len(data))
        self.assertEqual((pl.points, pl.texcoords), (((7, 8, 9),), (0.5, 0.25)))

    def test_word_count_mismatch_is_not_a_point_list(self):
        body = b"\x01\x01" + b"\x00\x00\x00" + b"\x00" + b"\x00" + b"\x02" + struct.pack("<3f", 1, 2, 3)
        data = self.node(body)
        with self.assertRaises(NotAPointList):
            decode_point_list(data, 0, len(data))


if __name__ == "__main__":
    unittest.main()
