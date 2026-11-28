"""CGM corpus run: decode every CATPart's exact B-rep to IR (cadcli cgm), mesh it with the kernel's `cgm` binary and
compare with the tessellation CATIA stored in the part (cadcli info). Indicators only: closure and area are not
ground truth.

    python3 -I tools/cgm_corpus.py run  LIST|DIR... --cadcli CADCLI --cgm CGM --out RUN_DIR [--ir IR_DIR] [-j 4]
    python3 -I tools/cgm_corpus.py sum  RUN_DIR
    python3 -I tools/cgm_corpus.py diff OLD_RUN NEW_RUN

LIST is a file with one CATPart path per line; DIR is searched for *.CATPart (symlinks followed, duplicates by real
path removed). RUN_DIR/parts.jsonl gets one line per part: {"part", "decode": cadcli cgm summary, "stored": cadcli
info, "mesh": cgm JSON, "fail": {reason: faces}}. IRs are kept in IR_DIR (default RUN_DIR/ir) so a second kernel can
be run on the same IRs with `--reuse-ir`. Each child process gets an 8 GB address-space cap and a time limit.
"""
import argparse
import collections
import concurrent.futures as cf
import hashlib
import json
import os
import re
import resource
import subprocess
import sys

MEM = 8 << 30
TIMEOUT = 120

# failure messages of cgm_mesh::tessellate -> reason class (numbers removed)
REASONS = [
    (re.compile(r"off the shared edge points"), "boundary off the shared edge points"),
    (re.compile(r"boundary segments are not triangle edges"), "boundary segments lost"),
    (re.compile(r"triangulation covered"), "triangles do not cover the UV area"),
    (re.compile(r"outside the surface domain"), "boundary outside the surface domain"),
    (re.compile(r"weight"), "invalid NURBS weights"),
    (re.compile(r"finite"), "non-finite parameters"),
    (re.compile(r"knot"), "decreasing knots"),
    (re.compile(r"no loop with 3 points"), "no loop with 3 points"),
    (re.compile(r"chain"), "loop does not chain"),
]


def reason(msg):
    for rx, name in REASONS:
        if rx.search(msg):
            return name
    return "other: " + re.sub(r"-?\d[\d.e+-]*", "N", msg)[:80]


def limit():
    resource.setrlimit(resource.RLIMIT_AS, (MEM, MEM))


def sh(cmd):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, preexec_fn=limit, timeout=TIMEOUT)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "timeout"


def last_json(text):
    for line in reversed(text.splitlines()):
        line = line.strip()
        if line.startswith("{"):
            try:
                return json.loads(line)
            except json.JSONDecodeError:
                return {"error": "unparsable JSON"}
    return None


def parts(args):
    seen, out = set(), []
    for a in args:
        if os.path.isdir(a):
            found = sorted(os.path.join(d, f) for d, _, fs in os.walk(a, followlinks=True) for f in fs if f.lower().endswith(".catpart"))
        else:
            with open(a) as f:
                found = [x.strip() for x in f if x.strip()]
        for p in found:
            rp = os.path.realpath(p)
            if rp not in seen:
                seen.add(rp)
                out.append(p)
    return out


def one(p, a, ir_dir):
    ir = os.path.join(ir_dir, hashlib.sha1(os.path.realpath(p).encode()).hexdigest()[:16] + ".cgmir")
    row = {"part": p, "ir": ir}
    if not (a.reuse_ir and os.path.exists(ir)):
        rc, out, err = sh([a.cadcli, "cgm", p, ir])
        row["decode"] = last_json(out) or {"error": (err or "no output")[-200:]}
        if rc != 0 and os.path.exists(ir):
            os.remove(ir)
    rc, out, _ = sh([a.cadcli, "info", p])
    row["stored"] = last_json(out)
    if not os.path.exists(ir):
        return row
    rc, out, err = sh([a.cgm, ir])
    row["mesh"] = last_json(out) or {"error": "rc %d: %s" % (rc, err[-200:])}
    fail = collections.Counter()
    for line in err.splitlines():
        m = re.search(r": face \d+: (.*)$", line)
        if m:
            fail[reason(m.group(1))] += 1
    row["fail"] = dict(fail)
    return row


def run(a):
    ps = parts(a.inputs)
    ir_dir = a.ir or os.path.join(a.out, "ir")
    os.makedirs(ir_dir, exist_ok=True)
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "parts.jsonl"), "w") as f, cf.ThreadPoolExecutor(a.j) as ex:
        for i, row in enumerate(ex.map(lambda p: one(p, a, ir_dir), ps)):
            f.write(json.dumps(row) + "\n")
            if (i + 1) % 100 == 0:
                print(f"{i + 1}/{len(ps)}", file=sys.stderr)
    summary(a.out)


def load(d):
    with open(os.path.join(d, "parts.jsonl")) as f:
        return [json.loads(x) for x in f]


def area_ok(r, tol=0.01):
    m, s = r.get("mesh") or {}, r.get("stored") or {}
    a, b = m.get("area"), s.get("area")
    return a is not None and b and abs(a - b) <= tol * b


def complete(r):
    m = r.get("mesh") or {}
    return "faces" in m and m.get("failed") == 0 and m.get("skipped") == 0


def summary(d):
    rows = load(d)
    meshed = [r for r in rows if "faces" in (r.get("mesh") or {})]
    tot = lambda k: sum(r["mesh"][k] for r in meshed)
    fail = collections.Counter()
    for r in rows:
        fail.update(r.get("fail") or {})
    comp = [r for r in meshed if complete(r)]
    print(f"parts {len(rows)}; IR written {sum(os.path.exists(r['ir']) for r in rows)}; meshed {len(meshed)}; "
          f"mesh errors {sum('error' in (r.get('mesh') or {}) for r in rows)}")
    print(f"faces {tot('faces')}; skipped by the decoder {tot('skipped')}; failed in the kernel {tot('failed')}")
    print(f"open edges {tot('open_edges')}; non-manifold edges {tot('nonmanifold_edges')}")
    print(f"closed parts (indicator) {sum(r['mesh']['open_edges'] == 0 for r in meshed)} of {len(meshed)}")
    print(f"complete parts (no failed/skipped face) {len(comp)}; of those within 1% of the stored area (indicator) "
          f"{sum(area_ok(r) for r in comp)}")
    print("kernel failures by reason (faces):")
    for k, n in fail.most_common():
        print(f"  {n:7d}  {k}")
    no_ir = collections.Counter()
    for r in rows:
        if not os.path.exists(r["ir"]):
            dec = r.get("decode") or {}
            no_ir[(dec.get("error") or dec.get("see_body") or "?")[:70]] += 1
    if no_ir:
        print("parts without IR:")
        for k, n in no_ir.most_common(12):
            print(f"  {n:7d}  {k}")


def diff(a, b):
    old = {r["part"]: r for r in load(a)}
    new = {r["part"]: r for r in load(b)}
    keys = ["failed", "open_edges", "nonmanifold_edges"]
    t = collections.Counter()
    worse = []
    for p, n in new.items():
        o = old.get(p)
        mo, mn = (o or {}).get("mesh") or {}, n.get("mesh") or {}
        if "faces" not in mo or "faces" not in mn:
            continue
        for k in keys:
            t[k, 0] += mo[k]
            t[k, 1] += mn[k]
        if mn["failed"] > mo["failed"] or mn["open_edges"] > mo["open_edges"]:
            worse.append((mn["failed"] - mo["failed"], mn["open_edges"] - mo["open_edges"], p))
    for k in keys:
        print(f"{k}: {t[k, 0]} -> {t[k, 1]}")
    print(f"parts worse in failed faces or open edges: {len(worse)}")
    for w in sorted(worse, reverse=True)[:20]:
        print("  failed %+d open %+d  %s" % w)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("inputs", nargs="+")
    r.add_argument("--cadcli", required=True)
    r.add_argument("--cgm", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--ir")
    r.add_argument("--reuse-ir", action="store_true")
    r.add_argument("-j", type=int, default=4)
    s = sub.add_parser("sum")
    s.add_argument("run")
    d = sub.add_parser("diff")
    d.add_argument("old")
    d.add_argument("new")
    a = ap.parse_args()
    if a.j > 4 if a.cmd == "run" else False:
        sys.exit("at most 4 workers")
    {"run": run, "sum": lambda a: summary(a.run), "diff": lambda a: diff(a.old, a.new)}[a.cmd](a)


if __name__ == "__main__":
    main()
