#!/usr/bin/env python3
"""Ground truth for assembly decoders: the occurrence table of a STEP export.

Writes one TSV row per NEXT_ASSEMBLY_USAGE_OCCURRENCE:
  instance id, parent product-definition #, child product-definition #, 12 placement values
(x axis, y axis, z axis, translation of the child in the parent, from the
ITEM_DEFINED_TRANSFORMATION's target placement; the source placement must be the identity).
Only the STEP text is read; nothing is inferred.
"""
import math
import re
import sys


def main(src, dst):
    txt = open(src, encoding='latin-1').read()
    data = txt[txt.index('DATA;'):]
    ents = {int(m.group(1)): re.sub(r'\s+', ' ', m.group(2)) for m in re.finditer(r'#(\d+)\s*=\s*(.*?);\s*(?=#\d+\s*=|ENDSEC)', data, re.S)}

    def body(i):
        s = ents[i]
        return s[:s.index('(')].strip(), s[s.index('(') + 1:s.rindex(')')]

    def split(a):
        out, depth, cur, q = [], 0, '', False
        for ch in a:
            if ch == "'":
                q = not q
            if not q and ch == '(':
                depth += 1
            if not q and ch == ')':
                depth -= 1
            if ch == ',' and depth == 0 and not q:
                out.append(cur.strip())
                cur = ''
            else:
                cur += ch
        out.append(cur.strip())
        return out

    ref = lambda x: int(x.strip()[1:])

    def vec(i):
        return [float(v) for v in split(body(i)[1])[1].strip('()').split(',')]

    def axis2(i):
        f = split(body(i)[1])
        loc = vec(ref(f[1]))
        z = vec(ref(f[2])) if f[2] != '$' else [0.0, 0.0, 1.0]
        x = vec(ref(f[3])) if f[3] != '$' else [1.0, 0.0, 0.0]
        n = math.sqrt(sum(c * c for c in z)); z = [c / n for c in z]
        d = sum(a * b for a, b in zip(x, z)); x = [a - d * b for a, b in zip(x, z)]
        n = math.sqrt(sum(c * c for c in x)); x = [c / n for c in x]
        y = [z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0]]
        return x, y, z, loc

    nauo = {}
    for i, s in ents.items():
        if s.startswith('NEXT_ASSEMBLY_USAGE_OCCURRENCE'):
            f = split(body(i)[1])
            nauo[i] = (f[0].strip("'"), ref(f[3]), ref(f[4]))
    rows = {}
    for i, s in ents.items():
        if not s.startswith('CONTEXT_DEPENDENT_SHAPE_REPRESENTATION'):
            continue
        f = split(body(i)[1])
        rr, pds = ref(f[0]), ref(f[1])
        target = ref(split(body(pds)[1])[2])
        idt = [int(x) for x in re.findall(r'#(\d+)', ents[rr]) if ents[int(x)].startswith('ITEM_DEFINED_TRANSFORMATION')]
        if len(idt) != 1:
            raise SystemExit('representation relationship #%d without exactly one transformation' % rr)
        g = split(body(idt[0])[1])
        sx, sy, sz, sl = axis2(ref(g[2]))
        if max(abs(v) for v in sl) > 1e-12 or max(abs(a - b) for a, b in zip(sz, [0, 0, 1])) > 1e-12 or max(abs(a - b) for a, b in zip(sx, [1, 0, 0])) > 1e-12:
            raise SystemExit('non-identity source placement at #%d' % idt[0])
        rows[target] = axis2(ref(g[3]))
    with open(dst, 'w') as out:
        out.write('# instance\tparent_pd\tchild_pd\tx0 x1 x2 y0 y1 y2 z0 z1 z2 t0 t1 t2 (from %s)\n' % src.rsplit('/', 1)[-1])
        for i, (name, parent, child) in nauo.items():
            x, y, z, t = rows[i]
            out.write('%s\t%d\t%d\t%s\n' % (name, parent, child, ' '.join(repr(v) for v in x + y + z + t)))
    print(len(nauo), 'occurrences ->', dst)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
