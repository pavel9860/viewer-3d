"""Collect public CAD sample files from GitHub for formats that lack test coverage.

Usage: python3 tools/collect_samples.py OUT_DIR [--per-format N] [--max-files M]

For every extension below, GitHub code search finds repositories containing such
files. The script keeps up to N repositories per format (distinct owners first,
for diversity of authoring software and versions), then makes blobless sparse
checkouts that fetch only the matching files (at most M per repository). It
writes OUT_DIR/SOURCES.tsv (format, repository, commit, licence, file count).
Files are data only; nothing downloaded is executed.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

FORMATS = {
    "solidworks": ["sldprt", "sldasm", "slddrw"],
    "inventor": ["ipt", "iam"],
    "creo-nx-prt": ["prt", "asm"],
    "catia-v5": ["catpart", "catproduct", "catdrawing"],
    "catia-v6-3dxml": ["3dxml"],
    "catia-v4-exp": ["exp", "session"],
    "catia-v4": ["model"],
    "catia-cgr": ["cgr"],
    "parasolid": ["x_t", "x_b"],
    "acis": ["sat", "sab"],
    "jt": ["jt"],
    "rhino": ["3dm"],
    "ifc": ["ifc"],
    "step": ["step", "stp"],
    "iges": ["iges", "igs"],
    "fusion": ["f3d", "f3z"],
    "dwg": ["dwg"],
    "sketchup": ["skp"],
    "solidedge": ["par", "psm"],
}


def gh_json(args: list[str]):
    for attempt in range(6):
        p = subprocess.run(["gh", "api", "-X", "GET", *args], capture_output=True, text=True)
        if p.returncode == 0:
            return json.loads(p.stdout)
        if "rate limit" in p.stderr.lower() or "secondary" in p.stderr.lower():
            time.sleep(20 * (attempt + 1))
            continue
        sys.stderr.write(p.stderr)
        return None
    return None


def search_repos(ext: str, pages: int = 3) -> dict[str, int]:
    repos: dict[str, int] = defaultdict(int)
    for page in range(1, pages + 1):
        res = gh_json(["search/code", "-f", f"q=extension:{ext}", "-f", "per_page=100", "-f", f"page={page}"])
        time.sleep(7)  # code search: about 10 requests per minute
        if not res or not res.get("items"):
            break
        for it in res["items"]:
            repos[it["repository"]["full_name"]] += 1
    return repos


def checkout(repo: str, dest: Path, exts: list[str], max_files: int) -> tuple[str, str, int]:
    if dest.exists():
        return "", "", sum(1 for _ in dest.rglob("*") if _.is_file() and _.suffix.lower().lstrip(".") in exts)
    url = f"https://github.com/{repo}.git"
    r = subprocess.run(["git", "clone", "-q", "--depth", "1", "--filter=blob:none", "--no-checkout", url, str(dest)], capture_output=True, text=True)
    if r.returncode:
        return "", "", 0
    names = subprocess.run(["git", "-C", str(dest), "ls-tree", "-r", "--name-only", "HEAD"], capture_output=True, text=True).stdout.splitlines()
    picked = [n for n in names if n.rsplit(".", 1)[-1].lower() in exts][:max_files]
    if not picked:
        subprocess.run(["rm", "-rf", str(dest)])
        return "", "", 0
    spec = "\n".join("/" + n.replace("[", "\\[").replace("*", "\\*").replace("?", "\\?") for n in picked) + "\n"
    subprocess.run(["git", "-C", str(dest), "sparse-checkout", "set", "--no-cone", "--stdin"], input=spec, text=True, capture_output=True)
    subprocess.run(["git", "-C", str(dest), "checkout", "-q", "HEAD"], capture_output=True)
    commit = subprocess.run(["git", "-C", str(dest), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    lic = gh_json([f"repos/{repo}"]) or {}
    licence = ((lic.get("license") or {}).get("spdx_id")) or "none"
    return commit, licence, len(picked)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("out", type=Path)
    ap.add_argument("--per-format", type=int, default=12)
    ap.add_argument("--max-files", type=int, default=60)
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    sources = args.out / "SOURCES.tsv"
    if not sources.exists():
        sources.write_text("# format\trepository\tcommit\tlicence\tfiles\n")
    for fmt, exts in FORMATS.items():
        if args.only and fmt not in args.only:
            continue
        counts: dict[str, int] = defaultdict(int)
        for ext in exts:
            for repo, n in search_repos(ext).items():
                counts[repo] += n
        # Diversity: at most one repository per owner first, then by file count.
        chosen, owners = [], set()
        for repo, _ in sorted(counts.items(), key=lambda kv: -kv[1]):
            owner = repo.split("/")[0]
            if owner in owners:
                continue
            owners.add(owner)
            chosen.append(repo)
            if len(chosen) == args.per_format:
                break
        for repo in chosen:
            dest = args.out / fmt / repo.replace("/", "__")
            commit, licence, n = checkout(repo, dest, exts, args.max_files)
            if n and commit:
                with open(sources, "a") as fh:
                    fh.write(f"{fmt}\t{repo}\t{commit}\t{licence}\t{n}\n")
            print(f"{fmt:16s} {repo:50s} files={n} licence={licence}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
