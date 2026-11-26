"""Reader of the CGM flat IR (cadcore/src/formats/catia/cgm_model.hpp), for diagnostics."""
import struct
def read(path):
    b = open(path, 'rb').read(); assert b[:8] == b'CGMIR\0\0\x02'; k = 8
    def u(n=1):
        nonlocal k; v = struct.unpack_from('<%dI' % n, b, k); k += 4 * n; return v
    def f(n):
        nonlocal k; v = struct.unpack_from('<%dd' % n, b, k); k += 8 * n; return v
    c = u(10); m = dict(surfaces=[], pcurves=[], edges=[], coedges=[], loops=[], faces=[], shells=[], bodies=[], skipped=[], vertices=c[2])
    for _ in range(c[0]): kind, i, n = u(3); m['surfaces'].append((kind, i, f(n)))
    for _ in range(c[1]): kind, s, i, n = u(4); m['pcurves'].append((kind, s, i, f(n)))
    for _ in range(c[3]): m['edges'].append(u(3))
    for _ in range(c[4]):
        pc, e, sense, _ = struct.unpack_from('<IIiI', b, k); k += 16; t = f(2); m['coedges'].append((pc, e, sense, t))
    for _ in range(c[5]): m['loops'].append(u(2))
    for _ in range(c[6]):
        s, sense, l0, nl, i = struct.unpack_from('<IiIII', b, k); k += 20; m['faces'].append((s, sense, l0, nl, i))
    for _ in range(c[7]):
        i, solid, n = u(3); fs = [struct.unpack_from('<Ii', b, k + 8 * j) for j in range(n)]; k += 8 * n; m['shells'].append((i, solid, fs))
    for _ in range(c[8]):
        i, flags, n = u(3); m['bodies'].append((i, flags, u(n) if n else ()))
    for _ in range(c[9]):
        i, t, n = u(3); m['skipped'].append((i, t, b[k:k + n].decode(errors='replace'))); k += n
    return m
