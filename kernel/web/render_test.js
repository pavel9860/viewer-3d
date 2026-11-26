// Render tests for raster.js against an exact reference (ray casting through a BVH), over many views.
//
// Per view, every sampled pixel centre is classified:
//   ok       renderer and reference see the same surface (same triangle, or depth equal within tolerance, or the
//            reference triangle is drawn within one pixel: a sub-pixel triangle missing every pixel centre)
//   hole     reference hits a triangle, renderer shows background (crack, lost triangle, wrong clipping)
//   spill    renderer draws where the reference sees background
//   occluded renderer shows a different, farther or nearer surface (depth order or culling error)
// Plus: determinism (same frame twice, and 1 vs 4 workers, give identical buffers), timing per model and resolution
// against the budget, degenerate sizes.
import { prepare, camera, createRenderer, createLodRenderer, pickLevel, NEAR } from './raster.js';

// frame budget at 800x600: 60 fps x 5 (3.3 ms)
const BUDGET_MS = 1000 / 60 / 5;

const out = document.getElementById('o');
const lines = [];
const log = s => { lines.push(s); out.textContent += s + '\n'; };

// ---------------------------------------------------------------- reference ray caster (BVH over triangles)
function bvh(m) {
  const { nt, V, T } = m;
  const lo = new Float64Array(nt * 3), hi = new Float64Array(nt * 3), cen = new Float64Array(nt * 3);
  for (let t = 0; t < nt; t++) for (let k = 0; k < 3; k++) {
    const a = V[3 * T[3 * t] + k], b = V[3 * T[3 * t + 1] + k], c = V[3 * T[3 * t + 2] + k];
    lo[3 * t + k] = Math.min(a, b, c); hi[3 * t + k] = Math.max(a, b, c); cen[3 * t + k] = (a + b + c) / 3;
  }
  const idx = new Uint32Array(nt).map((_, i) => i);
  const nodes = []; // {lo, hi, l, r, s, e}
  const build = (s, e) => {
    const n = { lo: [Infinity, Infinity, Infinity], hi: [-Infinity, -Infinity, -Infinity], l: -1, r: -1, s, e };
    for (let i = s; i < e; i++) for (let k = 0; k < 3; k++) { n.lo[k] = Math.min(n.lo[k], lo[3 * idx[i] + k]); n.hi[k] = Math.max(n.hi[k], hi[3 * idx[i] + k]); }
    const me = nodes.length; nodes.push(n);
    if (e - s > 4) {
      const ax = [0, 1, 2].reduce((b, k) => (n.hi[k] - n.lo[k] > n.hi[b] - n.lo[b] ? k : b), 0);
      const part = Array.from(idx.subarray(s, e)).sort((a, b) => cen[3 * a + ax] - cen[3 * b + ax]);
      idx.set(part, s);
      const mid = (s + e) >> 1;
      n.l = build(s, mid); n.r = build(mid, e);
    }
    return me;
  };
  build(0, nt);
  return { nodes, idx };
}
// nearest triangle hit by the ray o + s d (s > 0) with camera depth beyond the near plane; returns [t, depth] or null
function cast(m, acc, o, d, cam) {
  const { V, T } = m, { nodes, idx } = acc;
  let best = Infinity, bt = -1;
  const inv = [1 / d[0], 1 / d[1], 1 / d[2]];
  const stack = [0];
  while (stack.length) {
    const n = nodes[stack.pop()];
    let t0 = 0, t1 = best;
    for (let k = 0; k < 3; k++) {
      let a = (n.lo[k] - o[k]) * inv[k], b = (n.hi[k] - o[k]) * inv[k];
      if (a > b) [a, b] = [b, a];
      t0 = Math.max(t0, a); t1 = Math.min(t1, b);
    }
    if (t0 > t1) continue;
    if (n.l >= 0) { stack.push(n.l, n.r); continue; }
    for (let i = n.s; i < n.e; i++) {
      const t = idx[i], a = 3 * T[3 * t], b = 3 * T[3 * t + 1], c = 3 * T[3 * t + 2];
      // Moller-Trumbore
      const e1 = [V[b] - V[a], V[b + 1] - V[a + 1], V[b + 2] - V[a + 2]], e2 = [V[c] - V[a], V[c + 1] - V[a + 1], V[c + 2] - V[a + 2]];
      const p = [d[1] * e2[2] - d[2] * e2[1], d[2] * e2[0] - d[0] * e2[2], d[0] * e2[1] - d[1] * e2[0]];
      const det = e1[0] * p[0] + e1[1] * p[1] + e1[2] * p[2];
      if (Math.abs(det) < 1e-300) continue;
      const iv = 1 / det, s = [o[0] - V[a], o[1] - V[a + 1], o[2] - V[a + 2]];
      const u = (s[0] * p[0] + s[1] * p[1] + s[2] * p[2]) * iv;
      if (u < -1e-9 || u > 1 + 1e-9) continue;
      const q = [s[1] * e1[2] - s[2] * e1[1], s[2] * e1[0] - s[0] * e1[2], s[0] * e1[1] - s[1] * e1[0]];
      const v = (d[0] * q[0] + d[1] * q[1] + d[2] * q[2]) * iv;
      if (v < -1e-9 || u + v > 1 + 1e-9) continue; // a ray exactly through a shared edge must hit one of its triangles
      const tt = (e2[0] * q[0] + e2[1] * q[1] + e2[2] * q[2]) * iv;
      if (tt <= 0 || tt >= best) continue;
      // beyond the near plane (camera depth along -z of the view)
      const x = o[0] + tt * d[0], y = o[1] + tt * d[1], z = o[2] + tt * d[2], r = cam.r;
      const cz = -(r[6] * x + r[7] * y + r[8] * z + cam.t[2]);
      if (cz <= NEAR) continue;
      best = tt; bt = t;
    }
  }
  return bt < 0 ? null : bt;
}
// eye position and ray direction (world) through pixel centre (px, py)
function ray(cam, w, h, px, py) {
  const r = cam.r, f = (h / 2) / cam.tanHalf;
  const dc = [(px + 0.5 - w / 2) / f, -(py + 0.5 - h / 2) / f, -1];
  const ec = [-cam.t[0], -cam.t[1], -cam.t[2]]; // eye in camera-translated frame
  // world = R^T (c - t)
  const tw = v => [r[0] * v[0] + r[3] * v[1] + r[6] * v[2], r[1] * v[0] + r[4] * v[1] + r[7] * v[2], r[2] * v[0] + r[5] * v[1] + r[8] * v[2]];
  return [tw(ec), tw(dc)];
}
// camera depth (1/z) of triangle t at pixel centre, from its plane
function invDepth(m, cam, w, h, t, px, py) {
  const [o, d] = ray(cam, w, h, px, py), { V, T } = m;
  const a = 3 * T[3 * t], n = [m.N[3 * t], m.N[3 * t + 1], m.N[3 * t + 2]];
  const den = n[0] * d[0] + n[1] * d[1] + n[2] * d[2];
  if (Math.abs(den) < 1e-300) return 0;
  const s = (n[0] * (V[a] - o[0]) + n[1] * (V[a + 1] - o[1]) + n[2] * (V[a + 2] - o[2])) / den;
  return 1 / s; // d has unit camera z, so s is the camera depth
}

// ---------------------------------------------------------------- one view
async function checkView(R, m, acc, cam, w, h, step) {
  const fr = await R.render(cam, w, h, { edges: false, twoSided: params.has('twosided') });
  const c = { ok: 0, hole: 0, spill: 0, occluded: 0, outline: 0 };
  const bad = [];
  for (let py = step >> 1; py < h; py += step) for (let px = step >> 1; px < w; px += step) {
    const [o, d] = ray(cam, w, h, px, py);
    const ref = cast(m, acc, o, d, cam), got = fr.id[py * w + px];
    let kind = 'ok';
    if (ref === null && got >= 0) kind = near(fr, w, h, px, py, -1) || snapHit(m, acc, cam, w, h, px, py, got) ? 'ok' : 'spill';
    else if (ref !== null && got < 0) {
      // background where the model is: a crack when the 8 neighbours are all drawn, else an outline difference
      // below the 1/16-pixel vertex snapping
      let filled = 0;
      for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) { const x = px + dx, y = py + dy; if ((dx || dy) && x >= 0 && y >= 0 && x < w && y < h && fr.id[y * w + x] >= 0) filled++; }
      kind = near(fr, w, h, px, py, ref) || nearBoundary(m, cam, w, h, px, py, ref) ? 'ok' : filled === 8 ? 'hole' : 'outline';
    }
    else if (ref !== null && got !== ref) {
      const zr = invDepth(m, cam, w, h, ref, px, py), zg = fr.depth[py * w + px];
      // coincident surfaces (both planes at the same depth here): either one is the right answer
      const zp = invDepth(m, cam, w, h, got, px, py);
      kind = Math.abs(zr - zp) <= 1e-6 * Math.abs(zr) || Math.abs(zr - zg) <= 1e-3 * Math.abs(zr) || snapHit(m, acc, cam, w, h, px, py, got) || near(fr, w, h, px, py, ref) || withinPixel(m, acc, cam, w, h, px, py, zg) ? 'ok' : 'occluded';
    }
    c[kind]++;
    if (kind !== 'ok' && bad.length < 3) {
      // exact depth of both triangles' planes at this pixel (relative difference shows near-coincident surfaces)
      const zr = ref === null ? 0 : invDepth(m, cam, w, h, ref, px, py), zg = got < 0 ? 0 : invDepth(m, cam, w, h, got, px, py);
      bad.push({ kind, px, py, ref, got, zRef: +zr.toPrecision(6), zGot: +zg.toPrecision(6), zDrawn: +fr.depth[py * w + px].toPrecision(6) });
    }
  }
  return { c, bad, fr };
}
// the drawn triangle is hit by a ray within the renderer's precision (1/16 pixel) of the pixel centre: vertices are
// snapped to that grid, so an edge-on wall or a sliver may cover the centre by that much
function snapHit(m, acc, cam, w, h, px, py, t) {
  // exact screen-space distance from the pixel centre to the projected triangle, clipped at the near plane
  const { V, T } = m, r = cam.r, tr = cam.t, f = (h / 2) / cam.tanHalf;
  const C = [0, 1, 2].map(k => { const i = 3 * T[3 * t + k], x = V[i], y = V[i + 1], z = V[i + 2]; return [r[0] * x + r[1] * y + r[2] * z + tr[0], r[3] * x + r[4] * y + r[5] * z + tr[1], -(r[6] * x + r[7] * y + r[8] * z + tr[2])]; });
  const Q = [];
  for (let j = 0; j < 3; j++) {
    const p = C[j], q = C[(j + 1) % 3], pin = p[2] > NEAR, qin = q[2] > NEAR;
    if (pin) Q.push(p);
    if (pin !== qin) { const u = (NEAR - p[2]) / (q[2] - p[2]); Q.push([p[0] + u * (q[0] - p[0]), p[1] + u * (q[1] - p[1]), NEAR]); }
  }
  if (Q.length < 3) return false;
  const P = Q.map(c => [w / 2 + f * c[0] / c[2], h / 2 - f * c[1] / c[2]]), n = P.length, q = [px + 0.5, py + 0.5];
  const cr = (a, b, c) => (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
  const sg = P.map((_, i) => cr(P[i], P[(i + 1) % n], q));
  if (sg.every(x => x >= 0) || sg.every(x => x <= 0)) return true;
  const seg = (a, b) => { const dx = b[0] - a[0], dy = b[1] - a[1], l = dx * dx + dy * dy; const u = l ? Math.max(0, Math.min(1, ((q[0] - a[0]) * dx + (q[1] - a[1]) * dy) / l)) : 0; return Math.hypot(a[0] + u * dx - q[0], a[1] + u * dy - q[1]); };
  return Math.min(...P.map((_, i) => seg(P[i], P[(i + 1) % n]))) <= 1 / 16;
}
// the pixel centre lies within the renderer's precision (1/16 px) of the boundary of triangle t: snapping can move
// it to either side (folds and outlines)
function nearBoundary(m, cam, w, h, px, py, t) {
  const { V, T } = m, r = cam.r, tr = cam.t, f = (h / 2) / cam.tanHalf;
  const P = [];
  for (let k = 0; k < 3; k++) {
    const i = 3 * T[3 * t + k], x = V[i], y = V[i + 1], z = V[i + 2];
    const cz = -(r[6] * x + r[7] * y + r[8] * z + tr[2]);
    if (cz <= NEAR) return true; // crosses the near plane: clipped edges, accept
    P.push([w / 2 + f * (r[0] * x + r[1] * y + r[2] * z + tr[0]) / cz, h / 2 - f * (r[3] * x + r[4] * y + r[5] * z + tr[1]) / cz]);
  }
  const q = [px + 0.5, py + 0.5];
  const seg = (a, b) => { const dx = b[0] - a[0], dy = b[1] - a[1], l = dx * dx + dy * dy; const u = l ? Math.max(0, Math.min(1, ((q[0] - a[0]) * dx + (q[1] - a[1]) * dy) / l)) : 0; return Math.hypot(a[0] + u * dx - q[0], a[1] + u * dy - q[1]); };
  return Math.min(seg(P[0], P[1]), seg(P[1], P[2]), seg(P[2], P[0])) <= 1 / 16;
}
// does the ray hit triangle t (any distance in front)
function castOne(m, o, d, t) {
  const { V, T } = m, a = 3 * T[3 * t], b = 3 * T[3 * t + 1], c = 3 * T[3 * t + 2];
  const e1 = [V[b] - V[a], V[b + 1] - V[a + 1], V[b + 2] - V[a + 2]], e2 = [V[c] - V[a], V[c + 1] - V[a + 1], V[c + 2] - V[a + 2]];
  const p = [d[1] * e2[2] - d[2] * e2[1], d[2] * e2[0] - d[0] * e2[2], d[0] * e2[1] - d[1] * e2[0]];
  const det = e1[0] * p[0] + e1[1] * p[1] + e1[2] * p[2];
  if (Math.abs(det) < 1e-300) return false;
  const iv = 1 / det, s = [o[0] - V[a], o[1] - V[a + 1], o[2] - V[a + 2]];
  const u = (s[0] * p[0] + s[1] * p[1] + s[2] * p[2]) * iv;
  const q = [s[1] * e1[2] - s[2] * e1[1], s[2] * e1[0] - s[0] * e1[2], s[0] * e1[1] - s[1] * e1[0]];
  const v = (d[0] * q[0] + d[1] * q[1] + d[2] * q[2]) * iv;
  return u >= -1e-9 && v >= -1e-9 && u + v <= 1 + 1e-9 && (e2[0] * q[0] + e2[1] * q[1] + e2[2] * q[2]) * iv > 0;
}
// the drawn depth lies within the range of true depths over the 3x3 pixel neighbourhood (steep, sub-pixel geometry)
function withinPixel(m, acc, cam, w, h, px, py, z) {
  let lo = Infinity, hi = -Infinity;
  for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
    const x = px + dx, y = py + dy;
    if (x < 0 || y < 0 || x >= w || y >= h) continue;
    const [o, d] = ray(cam, w, h, x, y), t = cast(m, acc, o, d, cam);
    const zt = t === null ? 0 : invDepth(m, cam, w, h, t, x, y);
    lo = Math.min(lo, zt); hi = Math.max(hi, zt);
  }
  return z >= lo * (1 - 1e-3) && z <= hi * (1 + 1e-3);
}
// the given triangle (or background, -1) appears within one pixel
function near(fr, w, h, px, py, t) {
  for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
    const x = px + dx, y = py + dy;
    if (x >= 0 && y >= 0 && x < w && y < h && fr.id[y * w + x] === t) return true;
  }
  return false;
}
const hash = a => { let x = 2166136261 >>> 0; const u = new Uint32Array(a.buffer, a.byteOffset, a.byteLength >> 2); for (let i = 0; i < u.length; i++) x = Math.imul(x ^ u[i], 16777619) >>> 0; return x.toString(16); };

// ---------------------------------------------------------------- suite
const params = new URLSearchParams(location.search);
const quick = params.has('quick');
const only = params.get('m');
const dirs = [[0, 0], [Math.PI / 2, 0], [-Math.PI / 2, 0], [0, Math.PI / 2], [0, Math.PI], [0, -Math.PI / 2],
  [0.6, -0.5], [-0.6, 0.5], [0.6, 2.6], [-0.6, -2.6], [1.2, 0.8], [-1.2, 2.2], [0.3, 1.7], [-0.9, -1.1]];
const dists = params.has('d') ? params.get('d').split(',').map(Number) : quick ? [0.6, 2.5] : [0.3, 0.6, 1.0, 1.5, 2.5, 6, 30];
const pans = [[0, 0], [0.5, -0.3]];

try {
  const list = only ? [only] : await (await fetch('meshes/index.json')).json();
  let failures = 0;
  const perf = [];
  for (const name of list) {
    const M = prepare(await (await fetch('meshes/' + name)).arrayBuffer());
    // every detail level is checked as a mesh of its own; frame times are of the level the camera picks
    const L = await createLodRenderer(M);
    for (let lv = 0; lv < M.levels.length; lv++) {
    if (params.has('lv') && +params.get('lv') !== lv) continue;
    const m = M.levels[lv];
    const R = await createRenderer(m), R1 = await createRenderer(m, 1); // reorders m's triangles: before the BVH
    const acc = bvh(m);
    const tot = { ok: 0, hole: 0, spill: 0, occluded: 0, outline: 0 };
    const worst = [];
    let views = 0;
    for (const dir of (quick ? dirs.slice(0, 6) : dirs)) for (const dist of dists) for (const pan of (quick ? pans.slice(0, 1) : pans)) {
      const cam = camera(dir, pan, dist);
      const { c, bad } = await checkView(R, m, acc, cam, 160, 120, 2);
      views++;
      for (const k in tot) tot[k] += c[k];
      if (c.hole + c.spill + c.occluded > 0) worst.push({ dir, dist, pan, ...c, bad });
    }
    // determinism and timing at viewer resolutions
    const cam = camera([0.6, -0.5], [0, 0], 2.5);
    const fh = async (r) => { const f = await r.render(cam, 640, 480, { edges: true }); return hash(f.color) + hash(f.id) + hash(f.depth); };
    const ha = await fh(R), hb = await fh(R), h1 = await fh(R1);
    // frame time over the 14 directions at the default distance: median and worst
    const time = async (w, h) => {
      const ts = [];
      for (const dir of dirs) {
        const c = camera(dir, [0, 0], 2.5), rs = [];
        await L.render(c, w, h, { edges: true });
        for (let i = 0; i < 5; i++) { const t0 = performance.now(); await L.render(c, w, h, { edges: true }); rs.push(performance.now() - t0); }
        ts.push(Math.min(...rs));
      }
      ts.sort((x, y) => x - y);
      return [ts[ts.length >> 1], ts[ts.length - 1]];
    };
    // frame times only once per model (they are of the picked level, not of this one)
    const [t800, w800] = lv ? [0, 0] : await time(800, 600), [t400] = lv ? [0] : await time(400, 300), [t1080] = quick || lv ? [0] : await time(1920, 1080);
    const picked = pickLevel(M, camera([0, 0], [0, 0], 2.5), 600);
    // degenerate sizes must not throw
    for (const [w, h] of [[1, 1], [7, 5], [3, 200]]) await R.render(cam, w, h, {});
    const errs = tot.hole + tot.spill + tot.occluded;
    const outlineRate = tot.outline / (tot.ok + errs + tot.outline);
    const det = ha === hb && ha === h1;
    const ok = errs === 0 && det && outlineRate < 1e-4 && w800 <= BUDGET_MS;
    failures += ok ? 0 : 1;
    if (!lv) perf.push({ name, nt: M.levels[picked].nt, t800: w800, t400, t1080 });
    log(`${ok ? 'PASS' : 'FAIL'} ${(lv ? `  level ${lv} (dev ${(m.err / 2).toExponential(1)} of diag)` : name).padEnd(24)} ${String(m.nt).padStart(7)} tris ${m.closed ? 'closed' : 'open  '} | ${views} views, ${tot.ok + errs + tot.outline} samples: crack ${tot.hole} spill ${tot.spill} occluded ${tot.occluded} (outline 1px ${tot.outline}) | deterministic ${det} (1 vs ${R.threads} threads) ${lv ? '' : ` | picked level ${picked} (${M.levels[picked].nt} tris) 800x600 ${t800.toFixed(2)} ms (worst ${w800.toFixed(2)}), 400x300 ${t400.toFixed(2)} ms${t1080 ? `, 1920x1080 ${t1080.toFixed(2)} ms` : ''}`}`);
    for (const wv of worst.slice(0, 3)) log(`     view dir ${wv.dir.map(x => x.toFixed(2))} dist ${wv.dist} pan ${wv.pan}: hole ${wv.hole} spill ${wv.spill} occluded ${wv.occluded} e.g. ${JSON.stringify(wv.bad)}`);
    }
  }
  const worstT = perf.reduce((a, b) => (b.t800 > a.t800 ? b : a), perf[0]);
  log(`\n${failures === 0 ? 'ALL PASS' : failures + ' FAILED'} | slowest 800x600 view: ${worstT.name} ${worstT.t800.toFixed(2)} ms (budget ${BUDGET_MS.toFixed(2)} ms) (${(worstT.nt / worstT.t800 / 1000).toFixed(2)} M tris/s)`);
} catch (e) { log('ERROR ' + e + '\n' + (e.stack || '')); }
log('DONE');
fetch('/result', { method: 'POST', body: lines.join('\n') }).catch(() => {});
