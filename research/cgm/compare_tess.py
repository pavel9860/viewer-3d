"""CGM B-rep meshes vs CATIA's stored tessellation (ground truth) for parts that have both.

usage: python3 -I compare_tess.py PARTLIST OUTDIR [-j4]
Per part: cadcli cgm -> IR, kernel `cgm` -> mesh stats, cadcli info -> stored mesh stats.
Writes OUTDIR/compare.tsv and prints a summary. Each tool call has a 60 s limit."""
import json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
CADCLI = os.path.join(HERE, '..', '..', 'build', 'cadcli')
KERNEL = os.path.join(HERE, '..', '..', '..', 'kernel', 'target', 'release', 'cgm')

def run(cmd):
    try: r = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    except subprocess.TimeoutExpired: return {'error': 'timeout'}
    line = (r.stdout.strip().splitlines() or [''])[-1]
    try: return json.loads(line)
    except json.JSONDecodeError: return {'error': (r.stderr or r.stdout)[-200:]}

def one(args):
    i, part, out = args
    ir = os.path.join(out, f'{i}.cgmir')
    d = run([CADCLI, 'cgm', part, ir])
    if 'error' in d or not d.get('dialect'): return part, d, None, None
    k = run([KERNEL, ir])
    s = run([CADCLI, 'info', part])
    return part, d, k, s

def main():
    parts = [l for l in open(sys.argv[1]).read().splitlines() if l]
    out = sys.argv[2]; os.makedirs(out, exist_ok=True)
    with ThreadPoolExecutor(4) as ex: res = list(ex.map(one, [(i, p, out) for i, p in enumerate(parts)]))
    rows = []
    for part, d, k, s in res:
        if k is None or 'error' in (k or {}) or not s or 'error' in s: rows.append((part, 'nomesh', json.dumps(d)[:120] if k is None else json.dumps(k)[:120])); continue
        complete = d['faces_decoded'] == d['faces'] and not d['skipped']
        if not s.get('triangles'): rows.append((part, 'no stored tessellation', '', k['solids'], k['open_edges'], k['volume'])); continue
        rv = abs(k['volume']) / s['volume'] - 1 if s['volume'] else float('nan')
        ra = k['area'] / s['area'] - 1 if s['area'] else float('nan')
        rows.append((part, 'complete' if complete else 'partial', f'{rv:+.4f}', f'{ra:+.4f}', k['open_edges'], k['nonmanifold_edges'], 'neg' if k['volume'] < 0 else 'pos', k['solids'], k['tris']))
    with open(os.path.join(out, 'compare.tsv'), 'w') as f:
        for r in rows: f.write('\t'.join(map(str, r)) + '\n')
    print(len(rows), 'parts')

if __name__ == '__main__': main()
