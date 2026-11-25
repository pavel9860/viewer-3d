//! Robust orientation predicates: floating-point filter, then Shewchuk-style staged evaluation
//! (exact products of rounded differences, first-order tail correction) and an allocation-free exact
//! fallback on fixed stack expansions. Signs are always exact.

const EPS: f64 = 1.1102230246251565e-16; // 2^-53

#[inline] fn two_sum(a: f64, b: f64) -> (f64, f64) { let s = a + b; let bv = s - a; let av = s - bv; (s, (a - av) + (b - bv)) }
#[inline] fn fast_two_sum(a: f64, b: f64) -> (f64, f64) { let s = a + b; (s, b - (s - a)) }
#[inline] fn two_prod(a: f64, b: f64) -> (f64, f64) { let p = a * b; (p, a.mul_add(b, -p)) }
#[inline] fn two_diff(a: f64, b: f64) -> [f64; 2] { let (s, e) = two_sum(a, -b); [e, s] }
#[inline] fn diff_tail(a: f64, b: f64, x: f64) -> f64 { let bv = a - x; let av = x + bv; (a - av) + (bv - b) }
/// (a1 + a0) - (b1 + b0) as a 4-term expansion
#[inline] fn two_two_diff(a1: f64, a0: f64, b1: f64, b0: f64) -> [f64; 4] {
    let (i, x0) = two_sum(a0, -b0);
    let (j, r0) = two_sum(a1, i);
    let (i2, x1) = two_sum(r0, -b1);
    let (x3, x2) = two_sum(j, i2);
    [x0, x1, x2, x3]
}

/// h = e + f (nonoverlapping expansions, increasing magnitude, zeros eliminated); returns length
fn sum(e: &[f64], f: &[f64], h: &mut [f64]) -> usize {
    let (el, fl) = (e.len(), f.len());
    if el == 0 { h[..fl].copy_from_slice(f); return fl; }
    if fl == 0 { h[..el].copy_from_slice(e); return el; }
    let take_e = |ei: usize, fi: usize| fi >= fl || (ei < el && (f[fi] > e[ei]) == (f[fi] > -e[ei]));
    let (mut ei, mut fi, mut hi) = (0, 0, 0);
    let mut q = if take_e(ei, fi) { ei += 1; e[0] } else { fi += 1; f[0] };
    while ei < el || fi < fl {
        let x = if take_e(ei, fi) { ei += 1; e[ei - 1] } else { fi += 1; f[fi - 1] };
        let (nq, hh) = two_sum(q, x);
        q = nq;
        if hh != 0.0 { h[hi] = hh; hi += 1; }
    }
    if q != 0.0 || hi == 0 { h[hi] = q; hi += 1; }
    hi
}
/// h = e * b; h.len() >= 2 e.len()
fn scale(e: &[f64], b: f64, h: &mut [f64]) -> usize {
    let (mut q, hh) = two_prod(e[0], b);
    let mut hi = 0;
    if hh != 0.0 { h[hi] = hh; hi += 1; }
    for &x in &e[1..] {
        let (p1, p0) = two_prod(x, b);
        let (s, hh) = two_sum(q, p0);
        if hh != 0.0 { h[hi] = hh; hi += 1; }
        let (nq, hh) = fast_two_sum(p1, s);
        q = nq;
        if hh != 0.0 { h[hi] = hh; hi += 1; }
    }
    if q != 0.0 || hi == 0 { h[hi] = q; hi += 1; }
    hi
}
/// h = e * f for short f
fn mul(e: &[f64], f: &[f64], h: &mut [f64]) -> usize {
    let (mut acc, mut tmp, mut part) = ([0.0; 256], [0.0; 256], [0.0; 128]);
    let mut n = 0;
    for &b in f {
        let m = scale(e, b, &mut part);
        let k = sum(&acc[..n], &part[..m], &mut tmp);
        acc[..k].copy_from_slice(&tmp[..k]);
        n = k;
    }
    h[..n].copy_from_slice(&acc[..n]);
    n
}
fn neg(e: &mut [f64]) { e.iter_mut().for_each(|x| *x = -*x); }
fn sign_of(e: &[f64]) -> i8 { let t = e[e.len() - 1]; (t > 0.0) as i8 - (t < 0.0) as i8 }
pub fn sg(x: f64) -> i8 { (x > 0.0) as i8 - (x < 0.0) as i8 }

/// exact sign of (b-a) x (c-a)
pub fn orient2d_exact(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> i8 {
    let (mut l, mut r, mut h) = ([0.0; 8], [0.0; 8], [0.0; 16]);
    let nl = mul(&two_diff(b[0], a[0]), &two_diff(c[1], a[1]), &mut l);
    let nr = mul(&two_diff(b[1], a[1]), &two_diff(c[0], a[0]), &mut r);
    neg(&mut r[..nr]);
    let n = sum(&l[..nl], &r[..nr], &mut h);
    sign_of(&h[..n])
}
/// exact sign of det[a-p; b-p; c-p]
pub fn orient3d_exact(a: [f64; 3], b: [f64; 3], c: [f64; 3], p: [f64; 3]) -> i8 {
    let d = |u: [f64; 3]| [two_diff(u[0], p[0]), two_diff(u[1], p[1]), two_diff(u[2], p[2])];
    let (da, db, dc) = (d(a), d(b), d(c));
    let minor = |u: &[[f64; 2]; 3], v: &[[f64; 2]; 3], out: &mut [f64; 16]| {
        let (mut l, mut r) = ([0.0; 8], [0.0; 8]);
        let nl = mul(&u[1], &v[2], &mut l);
        let nr = mul(&u[2], &v[1], &mut r);
        neg(&mut r[..nr]);
        sum(&l[..nl], &r[..nr], out)
    };
    let mut m = [[0.0; 16]; 3];
    let n = [minor(&db, &dc, &mut m[0]), minor(&dc, &da, &mut m[1]), minor(&da, &db, &mut m[2])];
    let mut t = [[0.0; 64]; 3];
    let k = [mul(&m[0][..n[0]], &da[0], &mut t[0]), mul(&m[1][..n[1]], &db[0], &mut t[1]), mul(&m[2][..n[2]], &dc[0], &mut t[2])];
    let (mut s12, mut s) = ([0.0; 128], [0.0; 192]);
    let n12 = sum(&t[0][..k[0]], &t[1][..k[1]], &mut s12);
    let ns = sum(&s12[..n12], &t[2][..k[2]], &mut s);
    sign_of(&s[..ns])
}

/// adaptive orient2d; returns (sign, stage that decided: 0 filter, 1 B, 2 C, 3 exact)
pub fn orient2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> (i8, u8) {
    let (acx, bcx, acy, bcy) = (a[0] - c[0], b[0] - c[0], a[1] - c[1], b[1] - c[1]);
    let (dl, dr) = (acx * bcy, acy * bcx);
    let det = dl - dr;
    let detsum = dl.abs() + dr.abs();
    if det.abs() >= (3.0 + 16.0 * EPS) * EPS * detsum { return (sg(det), 0); }
    let ((l1, l0), (r1, r0)) = (two_prod(acx, bcy), two_prod(acy, bcx));
    let bexp = two_two_diff(l1, l0, r1, r0);
    let mut det: f64 = bexp.iter().sum();
    if det.abs() >= (2.0 + 12.0 * EPS) * EPS * detsum { return (sg(det), 1); }
    let (acxt, bcxt, acyt, bcyt) = (diff_tail(a[0], c[0], acx), diff_tail(b[0], c[0], bcx), diff_tail(a[1], c[1], acy), diff_tail(b[1], c[1], bcy));
    if acxt == 0.0 && acyt == 0.0 && bcxt == 0.0 && bcyt == 0.0 { return (sg(det), 1); }
    let errc = (9.0 + 64.0 * EPS) * EPS * EPS * detsum + (3.0 + 8.0 * EPS) * EPS * det.abs();
    det += (acx * bcyt + bcy * acxt) - (acy * bcxt + bcx * acyt);
    if det.abs() >= errc { return (sg(det), 2); }
    let (mut h1, mut h2, mut h3) = ([0.0; 8], [0.0; 12], [0.0; 16]);
    let tail = |x1: f64, y1: f64, x2: f64, y2: f64| { let ((s1, s0), (t1, t0)) = (two_prod(x1, y1), two_prod(x2, y2)); two_two_diff(s1, s0, t1, t0) };
    let n1 = sum(&bexp, &tail(acxt, bcy, acyt, bcx), &mut h1);
    let n2 = sum(&h1[..n1], &tail(acx, bcyt, acy, bcxt), &mut h2);
    let n3 = sum(&h2[..n2], &tail(acxt, bcyt, acyt, bcxt), &mut h3);
    (sign_of(&h3[..n3]), 3)
}

/// adaptive orient3d; returns (sign, stage)
pub fn orient3d(a: [f64; 3], b: [f64; 3], c: [f64; 3], p: [f64; 3]) -> (i8, u8) {
    let (ax, ay, az) = (a[0] - p[0], a[1] - p[1], a[2] - p[2]);
    let (bx, by, bz) = (b[0] - p[0], b[1] - p[1], b[2] - p[2]);
    let (cx, cy, cz) = (c[0] - p[0], c[1] - p[1], c[2] - p[2]);
    let det = ax * (by * cz - bz * cy) + bx * (cy * az - cz * ay) + cx * (ay * bz - az * by);
    let perm = ((by * cz).abs() + (bz * cy).abs()) * ax.abs() + ((cy * az).abs() + (cz * ay).abs()) * bx.abs() + ((ay * bz).abs() + (az * by).abs()) * cx.abs();
    if det.abs() >= (7.0 + 56.0 * EPS) * EPS * perm { return (sg(det), 0); }
    let mn = |u1: f64, v2: f64, u2: f64, v1: f64| { let ((p1, p0), (q1, q0)) = (two_prod(u1, v2), two_prod(u2, v1)); two_two_diff(p1, p0, q1, q0) };
    let (mut s1, mut s2, mut s3) = ([0.0; 8], [0.0; 8], [0.0; 8]);
    let (k1, k2, k3) = (scale(&mn(by, cz, bz, cy), ax, &mut s1), scale(&mn(cy, az, cz, ay), bx, &mut s2), scale(&mn(ay, bz, az, by), cx, &mut s3));
    let (mut t, mut bexp) = ([0.0; 16], [0.0; 24]);
    let nt = sum(&s1[..k1], &s2[..k2], &mut t);
    let nb = sum(&t[..nt], &s3[..k3], &mut bexp);
    let mut det: f64 = bexp[..nb].iter().sum();
    if det.abs() >= (3.0 + 28.0 * EPS) * EPS * perm { return (sg(det), 1); }
    let tl = [diff_tail(a[0], p[0], ax), diff_tail(a[1], p[1], ay), diff_tail(a[2], p[2], az), diff_tail(b[0], p[0], bx), diff_tail(b[1], p[1], by), diff_tail(b[2], p[2], bz), diff_tail(c[0], p[0], cx), diff_tail(c[1], p[1], cy), diff_tail(c[2], p[2], cz)];
    if tl.iter().all(|v| *v == 0.0) { return (sg(det), 1); }
    let [axt, ayt, azt, bxt, byt, bzt, cxt, cyt, czt] = tl;
    let errc = (26.0 + 288.0 * EPS) * EPS * EPS * perm + (3.0 + 8.0 * EPS) * EPS * det.abs();
    det += (axt * (by * cz - bz * cy) + ax * ((byt * cz + by * czt) - (bzt * cy + bz * cyt)))
         + (bxt * (cy * az - cz * ay) + bx * ((cyt * az + cy * azt) - (czt * ay + cz * ayt)))
         + (cxt * (ay * bz - az * by) + cx * ((ayt * bz + ay * bzt) - (azt * by + az * byt)));
    if det.abs() >= errc { return (sg(det), 2); }
    (orient3d_exact(a, b, c, p), 3)
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Rng(u64);
    impl Rng { fn f(&mut self) -> f64 { self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15); let mut z = self.0; z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB); ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64 } }
    /// independent reference: grow-expansion arithmetic on Vecs (simple, slow, obviously exact)
    fn grow(e: &[f64], b: f64) -> Vec<f64> { let mut q = b; let mut h = vec![]; for &x in e { let (s, er) = two_sum(q, x); if er != 0.0 { h.push(er); } q = s; } h.push(q); h }
    fn add(e: &[f64], f: &[f64]) -> Vec<f64> { f.iter().fold(e.to_vec(), |r, &x| grow(&r, x)) }
    fn mulv(e: &[f64], f: &[f64]) -> Vec<f64> { let mut r = vec![]; for &b in f { let mut s = vec![]; for &x in e { let (p, er) = two_prod(x, b); s = grow(&s, er); s = grow(&s, p); } r = add(&r, &s); } r }
    fn sgn(e: &[f64]) -> i8 { e.iter().rev().find(|x| **x != 0.0).map_or(0, |x| sg(*x)) }
    fn d(a: f64, b: f64) -> Vec<f64> { two_diff(a, b).to_vec() }
    fn ref2(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> i8 { let r: Vec<f64> = mulv(&d(b[1], a[1]), &d(c[0], a[0])).iter().map(|x| -x).collect(); sgn(&add(&mulv(&d(b[0], a[0]), &d(c[1], a[1])), &r)) }
    #[test]
    fn orient2d_matches_reference_on_near_degenerate_input() {
        let mut r = Rng(1);
        for _ in 0..200_000 {
            let (a, b, t) = ([r.f(), r.f()], [r.f(), r.f()], r.f());
            let c = [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])];
            let want = ref2(a, b, c);
            assert_eq!(orient2d_exact(a, b, c), want);
            assert_eq!(orient2d(a, b, c).0, want);
        }
    }
    #[test]
    fn orient3d_adaptive_matches_exact() {
        let mut r = Rng(2);
        for _ in 0..50_000 {
            let (a, b, c) = ([r.f(), r.f(), r.f()], [r.f(), r.f(), r.f()], [r.f(), r.f(), r.f()]);
            let (s, t) = (r.f(), r.f());
            let p = [0, 1, 2].map(|k| a[k] + s * (b[k] - a[k]) + t * (c[k] - a[k]));
            assert_eq!(orient3d(a, b, c, p).0, orient3d_exact(a, b, c, p));
        }
    }
}
