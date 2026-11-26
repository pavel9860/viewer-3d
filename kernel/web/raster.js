// CPU renderer for kernel meshes (used when WebGL is unavailable, and as the tested reference path).
//
// The rasterizer is kcore::raster (Rust, compiled into k.wasm): exact 1/16-pixel coverage with a top-left rule, 1/z
// depth, near-plane clipping, back faces culled on closed meshes. Sort-middle over workers in shared memory: each
// worker bins a contiguous range of triangles to screen bands, then draws its own interleaved bands from every
// range in triangle order. The frame is bit-identical for any number of workers; without cross-origin isolation it
// runs on the calling thread.

/**
 * parse a `brep mesh --out` file: u32 nv, u32 nt, f32 xyz*nv, u32 3*nt, u32 face*nt, normalized to the unit sphere;
 * then optional detail levels ("LODS", f32 chord of the mesh, u32 count, per level f32 error and the same block).
 * Returns the mesh with `chord` (its deviation from the exact model, normalized units) and `levels`: the mesh itself
 * (err 0), then each coarser level with `err`, its largest deviation from the mesh; all in the same frame.
 */
export function prepare(buffer) {
  const dv = new DataView(buffer);
  const block = off => {
    const nv = dv.getUint32(off, true), nt = dv.getUint32(off + 4, true), p = off + 8;
    return { nv, nt, P: new Float32Array(buffer.slice(p, p + nv * 12)), T: new Uint32Array(buffer.slice(p + nv * 12, p + nv * 12 + nt * 12)), F: new Uint32Array(buffer.slice(p + nv * 12 + nt * 12, p + nv * 12 + nt * 16)), end: p + nv * 12 + nt * 16 };
  };
  const b0 = block(0);
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < b0.nv; i++) for (let k = 0; k < 3; k++) { const x = b0.P[3 * i + k]; if (x < lo[k]) lo[k] = x; if (x > hi[k]) hi[k] = x; }
  const c = [0, 1, 2].map(k => (lo[k] + hi[k]) / 2), s = 2 / (Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1);
  const m = build(b0, c, s);
  m.err = 0; m.levels = [m];
  let off = b0.end;
  if (off + 12 <= buffer.byteLength && String.fromCharCode(...new Uint8Array(buffer, off, 4)) === 'LODS') {
    m.chord = dv.getFloat32(off + 4, true) * s;
    const n = dv.getUint32(off + 8, true);
    off += 12;
    for (let i = 0; i < n; i++) {
      const err = dv.getFloat32(off, true) * s, b = block(off + 4);
      const l = build(b, c, s);
      l.err = err; l.levels = null;
      m.levels.push(l);
      off = b.end;
    }
  }
  return m;
}

function build({ nv, nt, P, T, F }, c, s) {
  const V = new Float64Array(nv * 3);
  for (let i = 0; i < nv; i++) for (let k = 0; k < 3; k++) V[3 * i + k] = (P[3 * i + k] - c[k]) * s;
  // world normals and face colours per triangle
  const N = new Float64Array(nt * 3), C = new Float32Array(nt * 3);
  for (let t = 0; t < nt; t++) {
    const a = 3 * T[3 * t], b = 3 * T[3 * t + 1], d = 3 * T[3 * t + 2];
    const ux = V[b] - V[a], uy = V[b + 1] - V[a + 1], uz = V[b + 2] - V[a + 2], wx = V[d] - V[a], wy = V[d + 1] - V[a + 1], wz = V[d + 2] - V[a + 2];
    const nx = uy * wz - uz * wy, ny = uz * wx - ux * wz, nz = ux * wy - uy * wx, l = Math.hypot(nx, ny, nz) || 1;
    N[3 * t] = nx / l; N[3 * t + 1] = ny / l; N[3 * t + 2] = nz / l;
    C.set(hue(F[t]), 3 * t);
  }
  // closed (every directed edge has its reverse): back faces can be skipped
  const half = new Map();
  for (let t = 0; t < nt; t++) for (let k = 0; k < 3; k++) {
    const a = T[3 * t + k], b = T[3 * t + (k + 1) % 3], key = a * nv + b, rev = b * nv + a;
    if (half.has(rev)) { const n = half.get(rev) - 1; if (n) half.set(rev, n); else half.delete(rev); } else half.set(key, (half.get(key) || 0) + 1);
  }
  const blo = [Infinity, Infinity, Infinity], bhi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < nv; i++) for (let k = 0; k < 3; k++) { blo[k] = Math.min(blo[k], V[3 * i + k]); bhi[k] = Math.max(bhi[k], V[3 * i + k]); }
  return { nv, nt, V, T, F, N, C, closed: half.size === 0, lo: blo, hi: bhi };
}

/**
 * the detail level to draw: the coarsest whose deviation from the full mesh, projected at the nearest point of the
 * model's bounding sphere, stays within `px` pixels (the full mesh when the eye is inside the sphere)
 */
export function pickLevel(m, cam, h, px = 1) {
  const f = (h / 2) / cam.tanHalf, t = cam.t;
  // the model centre (origin) is at t in camera space; the sphere through the bounding box's farthest corner
  const near = Math.hypot(t[0], t[1], t[2]) - Math.hypot(...[0, 1, 2].map(k => Math.max(Math.abs(m.lo[k]), Math.abs(m.hi[k]))));
  if (!(near > NEAR)) return 0;
  let k = 0;
  for (let i = 1; i < m.levels.length; i++) if (m.levels[i].err * f / near <= px) k = i;
  return k;
}

export function hue(i) {
  const h = (i * 0.61803398875) % 1, s = 0.45, v = 0.9, k = n => (n + h * 6) % 6;
  const f = n => v - v * s * Math.max(0, Math.min(k(n), 4 - k(n), 1));
  return [f(5), f(3), f(1)];
}

/** camera: rotation about x then y, pan in view space, distance along -z; fov half-angle 0.4 rad */
export function camera(rot, pan, dist) {
  const [ax, ay] = rot, cx = Math.cos(ax), sx = Math.sin(ax), cy = Math.cos(ay), sy = Math.sin(ay);
  // rows of the world -> camera rotation, then translation
  return { r: [cy, 0, sy, sx * sy, cx, -sx * cy, -cx * sy, sx, cx * cy], t: [pan[0], pan[1], -dist], tanHalf: Math.tan(0.4) };
}

export const NEAR = 1e-3;
const MB = 1 << 20, HEAP = 16 * MB, CAM_BYTES = 128;

// Threads: the calling thread is thread 0, workers are threads 1..3. A frame is started through a control block in
// shared memory (no message round trips): CTL = [generation, argument block, bin done, draw done]. Each thread bins
// its triangle range, waits until every range is binned, then draws its bands.
const workerSrc = `
  onmessage = async ({ data: m }) => {
    const X = (await WebAssembly.instantiate(m.mod, { env: { memory: m.mem } })).exports;
    X.__stack_pointer.value = m.sp;
    const ctl = new Int32Array(m.mem.buffer, m.ctl, 8);
    postMessage(0);
    let seen = Atomics.load(ctl, 0);
    for (;;) {
      Atomics.wait(ctl, 0, seen);
      seen = Atomics.load(ctl, 0);
      const U = new Uint32Array(m.mem.buffer), a = Atomics.load(ctl, 1) >>> 2, nw = U[a + 9], i = m.index;
      if (i >= nw) continue;
      X.raster_bin(U[a], U[a + 1], U[a + 2], U[a + 3], U[a + 4], U[a + 5], U[a + 6], U[a + 7], U[a + 11 + i], U[a + 8] + 28 * i, U[a + 11 + nw + 2 * i], U[a + 12 + nw + 2 * i]);
      Atomics.add(ctl, 2, 1); Atomics.notify(ctl, 2);
      for (let v; (v = Atomics.load(ctl, 2)) < nw;) Atomics.wait(ctl, 2, v);
      // a band list pool overflowed: the calling thread grows it and starts the frame again
      let ok = 1;
      for (let k = 0; k < nw; k++) if (U[(U[a + 8] >>> 2) + 7 * k + 5] > U[(U[a + 8] >>> 2) + 7 * k + 4]) ok = 0;
      if (ok) X.raster_draw(U[a], U[a + 8], nw, U[a + 10], i, nw);
      Atomics.add(ctl, 3, 1); Atomics.notify(ctl, 3);
    }
  };`;

let shared = null; // { mod, mem, X, workers, ctl, top }
async function runtime(nthreads) {
  if (shared) return shared;
  const mod = await WebAssembly.compile(await (await fetch(new URL('k.wasm', import.meta.url))).arrayBuffer());
  const iso = self.crossOriginIsolated === true && typeof SharedArrayBuffer === 'function';
  const mem = new WebAssembly.Memory({ initial: HEAP / 65536, maximum: 65536, shared: iso });
  const X = (await WebAssembly.instantiate(mod, { env: { memory: mem } })).exports;
  const rt = { mod, mem, X, workers: [], top: HEAP };
  rt.ctl = alloc(rt, 64);
  if (iso && nthreads > 1) {
    const url = URL.createObjectURL(new Blob([workerSrc], { type: 'text/javascript' }));
    rt.workers = await Promise.all([...Array(nthreads - 1).keys()].map(i => new Promise(r => {
      const w = new Worker(url); w.onmessage = () => r(w); w.postMessage({ mod, mem, sp: 8 * MB + (i + 1) * MB, ctl: rt.ctl, index: i + 1 });
    })));
  }
  return (shared = rt);
}
function alloc(rt, bytes) {
  const p = rt.top; rt.top += Math.ceil(bytes / 64) * 64;
  if (rt.top > rt.mem.buffer.byteLength) rt.mem.grow(Math.ceil((rt.top - rt.mem.buffer.byteLength) / 65536));
  return p;
}

/**
 * a renderer for one mesh: `await createRenderer(m)`, then `await r.render(cam, w, h, opts)`.
 * Reorders m in place: triangles (T, N, C, F) into the culling-friendly draw order (m.perm: new -> old index) and
 * vertices (V) by first use (m.vperm: new -> old index).
 */
export async function createRenderer(m, nthreads = 4) {
  const rt = await runtime(4);
  const nw = Math.max(1, Math.min(nthreads, rt.workers.length + 1));
  const X = rt.X, { nv, nt } = m;
  const VTX = X.raster_sizes(0), REC = X.raster_sizes(1), CLB = X.raster_sizes(2), CL = X.raster_sizes(3), RPT = X.raster_sizes(4);
  const BH = 16, DESC = 7;
  const u8 = () => new Uint8Array(rt.mem.buffer), u32 = () => new Uint32Array(rt.mem.buffer);
  const put = arr => { const p = alloc(rt, arr.byteLength); u8().set(new Uint8Array(arr.buffer, arr.byteOffset, arr.byteLength), p); return p; };
  if (!m.perm) {
    const pv0 = put(m.V), pt0 = put(m.T), pf0 = put(m.F), pk = alloc(rt, 16 * nt), pp = alloc(rt, 4 * nt);
    X.raster_order(pv0, pt0, pf0, nv, nt, pk, pp);
    const perm = new Uint32Array(rt.mem.buffer, pp, nt).slice();
    const T = new Uint32Array(3 * nt), F = new Uint32Array(nt), N = new Float64Array(3 * nt), C = new Float32Array(3 * nt);
    for (let k = 0; k < nt; k++) {
      const o = perm[k];
      F[k] = m.F[o];
      for (let j = 0; j < 3; j++) { T[3 * k + j] = m.T[3 * o + j]; N[3 * k + j] = m.N[3 * o + j]; C[3 * k + j] = m.C[3 * o + j]; }
    }
    // vertices in first-use order: the per-frame vertex cache is then read nearly sequentially
    const ptn = put(T), pnew = alloc(rt, 4 * nv), pold = alloc(rt, 4 * nv);
    X.raster_renumber(ptn, nv, nt, pnew, pold);
    const old = new Uint32Array(rt.mem.buffer, pold, nv), V = new Float64Array(3 * nv);
    for (let i = 0; i < nv; i++) for (let j = 0; j < 3; j++) V[3 * i + j] = m.V[3 * old[i] + j];
    T.set(new Uint32Array(rt.mem.buffer, ptn, 3 * nt));
    Object.assign(m, { T, F, N, C, V, perm, vperm: old.slice() });
  }
  const pv = put(m.V);
  const pt = put(m.T), pn = put(m.N), pc = put(m.C);
  const pcl = alloc(rt, CLB * Math.ceil(nt / CL));
  X.raster_clusters(pv, pt, pn, pc, nv, nt, pcl);
  // contiguous triangle ranges (whole clusters), one per thread; their order is the draw order
  const ncl = Math.ceil(nt / CL);
  const ranges = [...Array(nw)].map((_, i) => [Math.min(nt, CL * Math.floor(ncl * i / nw)), Math.min(nt, CL * Math.floor(ncl * (i + 1) / nw))]);
  const caches = ranges.map(() => { const p = alloc(rt, VTX * nv); u8().fill(0, p, p + VTX * nv); return p; });
  const pcam = alloc(rt, CAM_BYTES), pdesc = alloc(rt, 4 * DESC * nw), parg = alloc(rt, 4 * (11 + 3 * nw));
  const recs = ranges.map(([a, b]) => RPT * (b - a) + 1);
  const chunk = recs.map(n => ({ rec: alloc(rt, REC * n), head: 0, tail: 0, pool: 0, np: 0 }));
  let nbCap = 0, frame = 0, cap = 0, frameNo = 0;
  const setDesc = () => { const U = u32(); chunk.forEach((c, i) => U.set([c.rec, c.head, c.tail, c.pool, c.np, 0, recs[i]], pdesc / 4 + DESC * i)); };
  const growPool = (c, n) => { c.np = n; c.pool = alloc(rt, 8 * n); };
  chunk.forEach((c, i) => growPool(c, recs[i] + 1024));
  const ctl = () => new Int32Array(rt.mem.buffer, rt.ctl, 8);

  // one frame on all threads; false when a band list pool overflowed (grown here: the caller renders again)
  // wait without spinning (a spinning thread competes with the workers for cores)
  async function until(c, k, n) {
    for (let v; (v = Atomics.load(c, k)) < n;) {
      const r = Atomics.waitAsync(c, k, v);
      if (r.async) await r.value;
    }
  }
  async function run() {
    setDesc();
    u32().set([pcam, pv, pt, pn, pc, nv, nt, pcl, pdesc, nw, frame, ...caches, ...ranges.flat()], parg / 4);
    const c = ctl();
    if (nw > 1) {
      Atomics.store(c, 1, parg); Atomics.store(c, 2, 0); Atomics.store(c, 3, 0);
      Atomics.add(c, 0, 1); Atomics.notify(c, 0);
    }
    X.raster_bin(pcam, pv, pt, pn, pc, nv, nt, pcl, caches[0], pdesc, ranges[0][0], ranges[0][1]);
    if (nw > 1) {
      Atomics.add(c, 2, 1); Atomics.notify(c, 2);
      await until(c, 2, nw);
    }
    const D = u32();
    const over = chunk.map((_, i) => D[pdesc / 4 + DESC * i + 5]).map((need, i) => need > chunk[i].np ? need : 0);
    if (!over.some(Boolean)) X.raster_draw(pcam, pdesc, nw, frame, 0, nw);
    if (nw > 1) await until(c, 3, nw - 1);
    over.forEach((need, i) => { if (need) growPool(chunk[i], need + (need >> 2)); });
    return !over.some(Boolean);
  }

  async function render(cam, w, h, opts = {}) {
    const nb = Math.ceil(h / BH);
    if (nb > nbCap) { nbCap = nb; for (const c of chunk) { c.head = alloc(rt, 4 * nb); c.tail = alloc(rt, 4 * nb); } }
    if (w * h > cap) { cap = w * h; frame = alloc(rt, 12 * cap); } // grows only
    const r = cam.r, tr = cam.t;
    // back faces are hidden by front faces only when the eye is outside the (closed) model's bounding box
    const e = [-(r[0] * tr[0] + r[3] * tr[1] + r[6] * tr[2]), -(r[1] * tr[0] + r[4] * tr[1] + r[7] * tr[2]), -(r[2] * tr[0] + r[5] * tr[1] + r[8] * tr[2])];
    const inside = [0, 1, 2].every(k => e[k] >= m.lo[k] && e[k] <= m.hi[k]);
    const cull = m.closed && !opts.twoSided && !inside;
    const t0 = performance.now();
    for (;;) {
      // a new frame number per attempt: every vertex cache entry is recomputed
      frameNo = (frameNo + 1) >>> 0 || 1;
      const dv = new DataView(rt.mem.buffer, pcam, CAM_BYTES);
      for (let i = 0; i < 9; i++) dv.setFloat64(8 * i, r[i], true);
      for (let i = 0; i < 3; i++) dv.setFloat64(72 + 8 * i, tr[i], true);
      dv.setFloat64(96, (h / 2) / cam.tanHalf, true);
      dv.setUint32(104, w, true); dv.setUint32(108, h, true);
      dv.setUint32(112, (cull ? 1 : 0) | (opts.edges ? 2 : 0) | (opts.faces !== false ? 4 : 0), true);
      dv.setUint32(116, frameNo, true); dv.setUint32(120, BH, true); dv.setUint32(124, nb, true);
      if (await run()) break;
    }
    const ms = performance.now() - t0;
    const n = w * h, buf = rt.mem.buffer;
    return { w, h, ms, color: new Uint32Array(buf, frame, n), depth: new Float32Array(buf, frame + 4 * n, n), id: new Int32Array(buf, frame + 8 * n, n) };
  }
  return { render, threads: nw };
}

/**
 * a renderer over a mesh's detail levels: each frame draws the level `pickLevel` chooses (renderers per level are
 * made on first use); the frame reports the level drawn
 */
export async function createLodRenderer(m, nthreads = 4, px = 1) {
  const rs = new Map();
  const get = async k => { if (!rs.has(k)) rs.set(k, await createRenderer(m.levels[k], nthreads)); return rs.get(k); };
  const r0 = await get(0);
  return {
    threads: r0.threads,
    async render(cam, w, h, opts = {}) {
      const k = opts.level ?? pickLevel(m, cam, h, opts.px ?? px);
      const fr = await (await get(k)).render(cam, w, h, opts);
      return { ...fr, level: k, nt: m.levels[k].nt };
    },
  };
}
