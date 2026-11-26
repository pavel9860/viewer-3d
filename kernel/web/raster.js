// CPU renderer for kernel meshes (used when WebGL is unavailable, and as the tested reference path).
//
// Exact coverage: vertices are snapped to 1/16 pixel and edge functions are evaluated in integers (exact in doubles),
// with a top-left fill rule, so the triangles of a watertight mesh cover every pixel centre exactly once: no cracks,
// no double coverage. Depth is 1/z interpolated in screen space (perspective-correct). Triangles crossing the near
// plane are clipped against it. Closed meshes skip back faces. Every vertex is transformed once per frame.

/** parse a `brep mesh --out` file: u32 nv, u32 nt, f32 xyz*nv, u32 3*nt, u32 face*nt; normalized to the unit sphere */
export function prepare(buffer) {
  const dv = new DataView(buffer);
  const nv = dv.getUint32(0, true), nt = dv.getUint32(4, true);
  const P = new Float32Array(buffer, 8, nv * 3), T = new Uint32Array(buffer.slice(8 + nv * 12, 8 + nv * 12 + nt * 12)), F = new Uint32Array(buffer.slice(8 + nv * 12 + nt * 12, 8 + nv * 12 + nt * 16));
  const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < nv; i++) for (let k = 0; k < 3; k++) { const x = P[3 * i + k]; if (x < lo[k]) lo[k] = x; if (x > hi[k]) hi[k] = x; }
  const c = [0, 1, 2].map(k => (lo[k] + hi[k]) / 2), s = 2 / (Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1);
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
const SUB = 16; // sub-pixel precision of snapped vertices

/**
 * render into colour (ABGR u32), depth (1/z, 0 = empty) and triangle id (-1 = empty) buffers of w x h.
 * opts: { edges: bool, faces: bool (face colours), twoSided: bool }
 */
export function render(m, cam, w, h, opts = {}, out = null) {
  const n = w * h;
  if (!out || out.w !== w || out.h !== h) out = { w, h, color: new Uint32Array(n), depth: new Float32Array(n), id: new Int32Array(n) };
  const { color, depth, id } = out;
  color.fill(0xff211d1b); depth.fill(0); id.fill(-1);
  const { nv, nt, V, T, N, C } = m, r = cam.r, tr = cam.t;
  const f = (h / 2) / cam.tanHalf, cxs = w / 2, cys = h / 2;
  // camera-space and screen positions of every vertex, once
  const X = new Float64Array(nv), Y = new Float64Array(nv), Z = new Float64Array(nv), SX = new Float64Array(nv), SY = new Float64Array(nv);
  for (let i = 0; i < nv; i++) {
    const x = V[3 * i], y = V[3 * i + 1], z = V[3 * i + 2];
    const cx = r[0] * x + r[1] * y + r[2] * z + tr[0], cy = r[3] * x + r[4] * y + r[5] * z + tr[1], cz = r[6] * x + r[7] * y + r[8] * z + tr[2];
    X[i] = cx; Y[i] = cy; Z[i] = -cz;
    if (-cz > NEAR) { SX[i] = Math.round((cxs + f * cx / -cz) * SUB); SY[i] = Math.round((cys - f * cy / -cz) * SUB); }
  }
  // back faces are hidden by front faces only when the eye is outside the (closed) model: test its bounding box
  const e = [-(r[0] * tr[0] + r[3] * tr[1] + r[6] * tr[2]), -(r[1] * tr[0] + r[4] * tr[1] + r[7] * tr[2]), -(r[2] * tr[0] + r[5] * tr[1] + r[8] * tr[2])];
  const inside = [0, 1, 2].every(k => e[k] >= m.lo[k] && e[k] <= m.hi[k]);
  const cull = m.closed && !opts.twoSided && !inside, edges = !!opts.edges, faces = opts.faces !== false;
  const L = [0.4 / 1.2329, 0.6 / 1.2329, 1 / 1.2329];
  const poly = new Float64Array(4 * 3); // clipped polygon in camera space
  const sx = new Float64Array(4), sy = new Float64Array(4), iz = new Float64Array(4);
  let drawn = 0;
  for (let t = 0; t < nt; t++) {
    const a = T[3 * t], b = T[3 * t + 1], c = T[3 * t + 2];
    const za = Z[a], zb = Z[b], zc = Z[c];
    if (za <= NEAR && zb <= NEAR && zc <= NEAR) continue;
    const nx = r[0] * N[3 * t] + r[1] * N[3 * t + 1] + r[2] * N[3 * t + 2], ny = r[3] * N[3 * t] + r[4] * N[3 * t + 1] + r[5] * N[3 * t + 2], nz = r[6] * N[3 * t] + r[7] * N[3 * t + 1] + r[8] * N[3 * t + 2];
    const d = 0.25 + 0.75 * Math.abs(nx * L[0] + ny * L[1] + nz * L[2]);
    const cr = faces ? C[3 * t] : 0.8, cg = faces ? C[3 * t + 1] : 0.8, cb = faces ? C[3 * t + 2] : 0.8;
    const rgb = (0xff000000 | (Math.min(255, (cb * d * 255) | 0) << 16) | (Math.min(255, (cg * d * 255) | 0) << 8) | Math.min(255, (cr * d * 255) | 0)) >>> 0;
    let k;
    if (za > NEAR && zb > NEAR && zc > NEAR) {
      sx[0] = SX[a]; sy[0] = SY[a]; iz[0] = 1 / za; sx[1] = SX[b]; sy[1] = SY[b]; iz[1] = 1 / zb; sx[2] = SX[c]; sy[2] = SY[c]; iz[2] = 1 / zc;
      k = 3;
    } else {
      // clip against the near plane (Sutherland-Hodgman), at most 4 vertices
      k = 0;
      const vs = [a, b, c];
      for (let j = 0; j < 3; j++) {
        const p = vs[j], q = vs[(j + 1) % 3], pin = Z[p] > NEAR, qin = Z[q] > NEAR;
        if (pin) { poly[3 * k] = X[p]; poly[3 * k + 1] = Y[p]; poly[3 * k + 2] = Z[p]; k++; }
        if (pin !== qin) {
          // computed from the lower vertex index, so both triangles of a shared edge get the identical point
          const [u, v] = p < q ? [p, q] : [q, p], s = (NEAR - Z[u]) / (Z[v] - Z[u]);
          poly[3 * k] = X[u] + s * (X[v] - X[u]); poly[3 * k + 1] = Y[u] + s * (Y[v] - Y[u]); poly[3 * k + 2] = NEAR; k++;
        }
      }
      for (let j = 0; j < k; j++) { const z = poly[3 * j + 2]; sx[j] = Math.round((cxs + f * poly[3 * j] / z) * SUB); sy[j] = Math.round((cys - f * poly[3 * j + 1] / z) * SUB); iz[j] = 1 / z; }
    }
    for (let j = 1; j + 1 < k; j++) drawn += tri(sx[0], sy[0], iz[0], sx[j], sy[j], iz[j], sx[j + 1], sy[j + 1], iz[j + 1], rgb, t, edges, cull, w, h, color, depth, id);
  }
  out.drawn = drawn;
  return out;
}

// one triangle in 1/16-pixel integer coordinates: span rasterization with exact edge functions and a top-left rule
function tri(x0, y0, z0, x1, y1, z1, x2, y2, z2, rgb, tid, edges, cull, w, h, color, depth, id) {
  // orient counter-clockwise in screen space (y down): area > 0
  let area = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
  if (area === 0) return 0;
  // back face, decided on the snapped screen triangle (the same data as coverage, so a closed mesh still covers
  // every pixel exactly once): a front face is counter-clockwise seen from the eye, i.e. negative area with y down
  if (cull && area > 0) return 0;
  if (area < 0) { let t = x1; x1 = x2; x2 = t; t = y1; y1 = y2; y2 = t; t = z1; z1 = z2; z2 = t; area = -area; }
  const minx = Math.max(0, Math.ceil((Math.min(x0, x1, x2) / SUB) - 0.5)), maxx = Math.min(w - 1, Math.floor((Math.max(x0, x1, x2) / SUB) - 0.5));
  const miny = Math.max(0, Math.ceil((Math.min(y0, y1, y2) / SUB) - 0.5)), maxy = Math.min(h - 1, Math.floor((Math.max(y0, y1, y2) / SUB) - 0.5));
  if (minx > maxx || miny > maxy) return 0;
  // edge e_k opposite vertex k: E(p) = (b - a) x (p - a) with (a, b) = next vertices; A = dE/dx, B = dE/dy (per pixel)
  const A0 = (y1 - y2) * SUB, B0 = (x2 - x1) * SUB, A1 = (y2 - y0) * SUB, B1 = (x0 - x2) * SUB, A2 = (y0 - y1) * SUB, B2 = (x1 - x0) * SUB;
  // top-left rule: on an edge that is not top or left, the boundary itself is outside (bias -1 on exact integers)
  const tl = (dx, dy) => (dy < 0 || (dy === 0 && dx > 0)) ? 0 : -1;
  const o0 = tl(x2 - x1, y2 - y1), o1 = tl(x0 - x2, y0 - y2), o2 = tl(x1 - x0, y1 - y0);
  const px = minx * SUB + SUB / 2, py = miny * SUB + SUB / 2;
  let r0 = (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1) + o0;
  let r1 = (x0 - x2) * (py - y2) - (y0 - y2) * (px - x2) + o1;
  let r2 = (x1 - x0) * (py - y0) - (y1 - y0) * (px - x0) + o2;
  // E_k increases by A_k per pixel in x and by B_k per row in y; inside: all three >= 0
  const inv = 1 / area;
  const w0 = 0.6 * SUB * Math.hypot(x2 - x1, y2 - y1), w1 = 0.6 * SUB * Math.hypot(x0 - x2, y0 - y2), w2 = 0.6 * SUB * Math.hypot(x1 - x0, y1 - y0); // 0.6 px edge lines
  const dzdx = (A0 * z0 + A1 * z1 + A2 * z2) * inv;
  const zlo = Math.min(z0, z1, z2), zhi = Math.max(z0, z1, z2);
  let covered = 0;
  for (let y = miny; y <= maxy; y++, r0 += B0, r1 += B1, r2 += B2) {
    // x range where r + A (x - minx) >= 0 for all three edges
    let xs = minx, xe = maxx;
    if (A0 > 0) xs = Math.max(xs, minx + Math.ceil(-r0 / A0)); else if (A0 < 0) xe = Math.min(xe, minx + Math.floor(r0 / -A0)); else if (r0 < 0) continue;
    if (A1 > 0) xs = Math.max(xs, minx + Math.ceil(-r1 / A1)); else if (A1 < 0) xe = Math.min(xe, minx + Math.floor(r1 / -A1)); else if (r1 < 0) continue;
    if (A2 > 0) xs = Math.max(xs, minx + Math.ceil(-r2 / A2)); else if (A2 < 0) xe = Math.min(xe, minx + Math.floor(r2 / -A2)); else if (r2 < 0) continue;
    // the divisions may round either way: widen by one pixel, then shrink with the exact integer test
    xs = Math.max(minx, xs - 1); xe = Math.min(maxx, xe + 1);
    while (xs <= xe && (r0 + A0 * (xs - minx) < 0 || r1 + A1 * (xs - minx) < 0 || r2 + A2 * (xs - minx) < 0)) xs++;
    while (xe >= xs && (r0 + A0 * (xe - minx) < 0 || r1 + A1 * (xe - minx) < 0 || r2 + A2 * (xe - minx) < 0)) xe--;
    if (xs > xe) continue;
    let e0 = r0 + A0 * (xs - minx), e1 = r1 + A1 * (xs - minx), e2 = r2 + A2 * (xs - minx);
    // depth from the unbiased edge values (the fill-rule bias of -1 is not geometry)
    let z = ((e0 - o0) * z0 + (e1 - o1) * z1 + (e2 - o2) * z2) * inv;
    let o = y * w + xs;
    for (let x = xs; x <= xe; x++, o++, e0 += A0, e1 += A1, e2 += A2, z += dzdx) {
      // a near-degenerate triangle's snapped barycentrics can extrapolate: depth stays within its vertex depths
      const zc = z < zlo ? zlo : z > zhi ? zhi : z;
      if (zc > depth[o]) {
        depth[o] = zc; id[o] = tid;
        color[o] = edges && (e0 < w0 || e1 < w1 || e2 < w2) ? 0xff100f0d : rgb;
        covered++;
      }
    }
  }
  return covered;
}
