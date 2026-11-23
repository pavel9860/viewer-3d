#!/usr/bin/env python3
"""Curate samples/github: keep only files whose content really is the claimed CAD format.

For every file collected by collect_samples.py:
  * git internals (.git) are removed - provenance is kept in SOURCES.tsv (repo + commit);
  * git-LFS pointers are replaced by the real object from media.githubusercontent.com
    when --fetch-lfs is given (the pointer's oid and size are verified);
  * the file's signature is checked against the format; files that fail are removed and
    listed in REJECTED.tsv with the reason, files that pass are listed in INVENTORY.tsv.

Signature checks are positive identifications only. A binary file of a format whose
signature is not published (CGR, CATIA V4 .model, NX .prt, SolidWorks 2015+) is kept as
"unverified" rather than accepted or dropped silently.
"""
from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import sys
import urllib.request
from pathlib import Path

CFB = b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1"
ZIP = b"PK\x03\x04"

FORMATS = {
    "acis": ("sat", "sab"),
    "catia-cgr": ("cgr",),
    "catia-v4": ("model",),
    "catia-v5": ("catpart", "catproduct", "catdrawing"),
    "catia-v6-3dxml": ("3dxml",),
    "catia-v4-exp": ("exp", "session"),
    "creo-nx-prt": ("prt", "asm"),
    "dwg": ("dwg",),
    "fusion": ("f3d", "f3z"),
    "ifc": ("ifc",),
    "iges": ("iges", "igs"),
    "inventor": ("ipt", "iam"),
    "jt": ("jt",),
    "parasolid": ("x_t", "x_b"),
    "rhino": ("3dm",),
    "sketchup": ("skp",),
    "solidedge": ("par", "psm", "asm"),
    "solidworks": ("sldprt", "sldasm", "slddrw"),
    "step": ("step", "stp"),
}


def is_text(b: bytes) -> bool:
    if not b:
        return True
    bad = sum(1 for x in b if x < 9 or (13 < x < 32 and x != 27))
    return bad == 0


def classify(ext: str, head: bytes) -> tuple[str, str]:
    """Returns (verdict, detail); verdict is verified | unverified | reject."""
    text = is_text(head[:4096])
    t = head.decode("latin-1")
    if ext == "sat":
        if re.match(r"\s*\d+ \d+ \d+ \d+\s", t):
            return "verified", "ACIS text header"
        return "reject", "no ACIS text header"
    if ext == "sab":
        return ("verified", "ACIS binary") if head.startswith(b"ACIS BinaryFile") else ("reject", "no 'ACIS BinaryFile' signature")
    if ext in ("ipt", "iam", "par", "psm"):
        return ("verified", "OLE compound file") if head.startswith(CFB) else ("reject", "not an OLE compound file")
    if ext == "asm":  # Creo / Solid Edge assemblies (both formats), never text sources
        if head.startswith(b"#UGC:"):
            return "verified", "Creo (#UGC header)"
        if head.startswith(CFB):
            return "verified", "OLE compound file (Solid Edge)"
        return "reject", "text file (assembly-language or other source)" if text else "unrecognised binary"
    if ext == "prt":
        if head.startswith(b"#UGC:"):
            return "verified", "Creo (#UGC header)"
        if text:
            return "reject", "text file, not a part"
        return "unverified", "binary .prt (NX?) - no published signature"
    if ext in ("sldprt", "sldasm", "slddrw"):
        if head.startswith(CFB):
            return "verified", "OLE compound file (SolidWorks <= 2014)"
        if text:
            return "reject", "text file"
        return "unverified", "binary SolidWorks 2015+ container - no published signature"
    if ext == "dwg":
        m = re.match(rb"AC10\d\d|AC1\.\d|AC2\.\d", head)
        return ("verified", head[:6].decode("latin-1")) if m else ("reject", "no AC10xx version signature")
    if ext in ("3dxml", "f3d", "f3z"):
        return ("verified", "zip container") if head.startswith(ZIP) else ("reject", "not a zip container")
    if ext == "ifc":
        if t.startswith("ISO-10303-21;") and re.search(r"FILE_SCHEMA\s*\(\s*\(\s*'IFC", t, re.I):
            return "verified", "STEP physical file, IFC schema"
        return "reject", "not ISO-10303-21 with an IFC schema"
    if ext in ("step", "stp"):
        if t.lstrip("﻿").startswith("ISO-10303-21;"):
            if re.search(r"FILE_SCHEMA\s*\(\s*\(\s*'IFC", t, re.I):
                return "reject", "IFC schema (belongs to ifc)"
            return "verified", "ISO-10303-21"
        return "reject", "no ISO-10303-21 header"
    if ext in ("iges", "igs"):
        first = t.splitlines()[0] if t else ""
        if len(first) >= 73 and first[72] in "SG":
            return "verified", "IGES fixed-format start section"
        return "reject", "no IGES 80-column start section"
    if ext == "jt":
        if head.startswith(b"Version ") and b"JT" in head[:80]:
            return "verified", head[:20].decode("latin-1").strip()
        return "reject", "no 'Version x.y JT' header"
    if ext == "x_t":
        if t.startswith("**ABCDEFGHIJKLMNOPQRSTUVWXYZ") or t.startswith("T51 : TRANSMIT FILE"):
            return "verified", "Parasolid text transmit"
        return "reject", "no Parasolid text header"
    if ext == "x_b":
        return ("verified", "Parasolid binary transmit") if head.startswith(b"PS") and not text else ("reject", "no Parasolid binary header")
    if ext == "3dm":
        return ("verified", "openNURBS 3DM") if head.startswith(b"3D Geometry File Format ") else ("reject", "no '3D Geometry File Format' header")
    if ext == "skp":
        return ("verified", "SketchUp model") if head.startswith(b"\xff\xfe\xff\x0eS\x00k\x00e\x00t\x00c\x00h\x00U\x00p") else ("reject", "no SketchUp signature")
    if ext in ("catpart", "catproduct", "catdrawing"):
        if head.startswith(b"V5_CFV2\x00"):
            m = re.search(rb"V5R\d+SP\d+|V5-6R20\d\d", head)
            return "verified", "V5_CFV2 " + (m.group(0).decode() if m else "release?")
        return "reject", "no V5_CFV2 signature"
    if ext in ("cgr", "model", "exp", "session"):
        if text:
            return "reject", "text file"
        return "unverified", "binary - no published signature"
    return "reject", "unknown extension"


def lfs_pointer(data: bytes) -> tuple[str, int] | None:
    if not data.startswith(b"version https://git-lfs.github.com/spec/v1"):
        return None
    m = re.search(rb"oid sha256:([0-9a-f]{64})\s+size (\d+)", data)
    return (m.group(1).decode(), int(m.group(2))) if m else None


def fetch_lfs(repo: str, commit: str, rel: str, oid: str, size: int, limit: int) -> bytes | None:
    if size > limit:
        return None
    url = f"https://media.githubusercontent.com/media/{repo}/{commit}/{urllib.request.quote(rel)}"
    req = urllib.request.Request(url, headers={"User-Agent": "cadcore-sample-curator"})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read(size + 1)
    except Exception:
        return None
    if len(data) != size or hashlib.sha256(data).hexdigest() != oid:
        return None
    return data


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--fetch-lfs", action="store_true")
    ap.add_argument("--lfs-limit", type=int, default=50 * 1024 * 1024)
    args = ap.parse_args()
    root: Path = args.root

    commits = {}
    for line in (root / "SOURCES.tsv").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fmt, repo, commit, *_ = line.split("\t")
        commits[(fmt, repo.replace("/", "__"))] = (repo, commit)

    inventory, rejected = [], []
    for fmt, exts in FORMATS.items():
        base = root / fmt
        if not base.is_dir():
            continue
        for git in list(base.glob("*/.git")):
            shutil.rmtree(git)
        for p in sorted(base.rglob("*")):
            if p.is_symlink():
                p.unlink()
                continue
            if not p.is_file():
                continue
            rel = p.relative_to(root).as_posix()
            ext = p.suffix.lower().lstrip(".")
            if ext not in exts:
                p.unlink()  # sparse checkouts can carry unrelated files (e.g. .gitattributes)
                continue
            data = p.read_bytes()
            ptr = lfs_pointer(data[:512])
            if ptr:
                repo_dir = p.relative_to(base).parts[0]
                repo, commit = commits.get((fmt, repo_dir), (None, None))
                got = None
                if args.fetch_lfs and repo:
                    got = fetch_lfs(repo, commit, "/".join(p.relative_to(base / repo_dir).parts), ptr[0], ptr[1], args.lfs_limit)
                if got is None:
                    rejected.append((fmt, rel, f"git-LFS pointer, object not fetched ({ptr[1]} bytes)"))
                    p.unlink()
                    continue
                p.write_bytes(got)
                data = got
            if not data:
                rejected.append((fmt, rel, "empty file"))
                p.unlink()
                continue
            verdict, detail = classify(ext, data[:8192])
            if verdict == "reject":
                rejected.append((fmt, rel, detail))
                p.unlink()
            else:
                inventory.append((fmt, rel, str(len(data)), verdict, detail))
        for d in sorted((x for x in base.rglob("*") if x.is_dir()), key=lambda x: -len(x.parts)):
            if not any(d.iterdir()):
                d.rmdir()

    with open(root / "INVENTORY.tsv", "w") as f:
        f.write("# format\tpath\tbytes\tverdict\tdetail\n")
        for row in inventory:
            f.write("\t".join(row) + "\n")
    with open(root / "REJECTED.tsv", "w") as f:
        f.write("# format\tpath (removed)\treason\n")
        for row in rejected:
            f.write("\t".join(row) + "\n")
    summary = {}
    for fmt, _, _, verdict, _ in inventory:
        summary.setdefault(fmt, {"verified": 0, "unverified": 0})[verdict] += 1
    for fmt in FORMATS:
        s = summary.get(fmt, {"verified": 0, "unverified": 0})
        nrej = sum(1 for r in rejected if r[0] == fmt)
        print(f"{fmt:16} verified {s['verified']:5}  unverified {s['unverified']:5}  rejected {nrej:5}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
