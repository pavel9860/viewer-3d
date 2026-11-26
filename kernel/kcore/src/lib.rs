//! Allocation-free, `no_std` kernel math shared by native and wasm32 builds, so determinism is tested on
//! identical code. Every function uses only IEEE add/sub/mul/div and own implementations of sin and sqrt.
#![no_std]

pub mod raster;

pub const SUM_N: usize = 4_000_000;
pub const SUM_CH: usize = 65_536;
pub const NT: usize = 4000;
pub const PN: u32 = 64;

/// splitmix64
pub struct Rng(pub u64);
impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    /// uniform in [-1, 1)
    pub fn f(&mut self) -> f64 { (self.next() >> 11) as f64 * (1.0 / (1u64 << 53) as f64) * 2.0 - 1.0 }
}

#[inline]
pub fn mix(h: u64, w: u64) -> u64 { let x = (h ^ w).wrapping_mul(0x9E37_79B9_7F4A_7C15); x ^ (x >> 29) }

/// 4-lane content hash (not cryptographic)
pub fn hash_bytes(b: &[u8]) -> u64 {
    let mut h = [0x243F_6A88_85A3_08D3u64, 0x1319_8A2E_0370_7344, 0xA409_3822_299F_31D0, 0x082E_FA98_EC4E_6C89];
    let word = |c: &[u8]| u64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]);
    let mut ch = b.chunks_exact(32);
    for c in &mut ch { for l in 0..4 { h[l] = mix(h[l], word(&c[l * 8..])); } }
    let mut t = [0u8; 32];
    let rem = ch.remainder();
    t[..rem.len()].copy_from_slice(rem);
    for l in 0..4 { h[l] = mix(h[l], word(&t[l * 8..])); }
    h.iter().fold(b.len() as u64, |o, &x| mix(o, x))
}
pub fn hash_f64(x: &[f64]) -> u64 { let mut h = x.len() as u64; for v in x { h = mix(h, v.to_bits()); } h }

/// sqrt: bit-trick seed + Newton; same operations on every target
#[inline]
pub fn own_sqrt(x: f64) -> f64 {
    if !(x > 0.0) { return 0.0; }
    let mut y = f64::from_bits((x.to_bits() >> 1) + 0x1FF8_0000_0000_0000);
    for _ in 0..5 { y = 0.5 * (y + x / y); }
    y
}
/// sin: fdlibm kernels, 2-part Cody-Waite reduction (|x| < ~1e5)
pub fn own_sin(x: f64) -> f64 {
    const INV: f64 = 6.366_197_723_675_813_4e-1;
    const P1: f64 = 1.570_796_326_734_125_614_17e0;
    const P2: f64 = 6.077_100_506_506_192_249_32e-11;
    let k = (x * INV + if x >= 0.0 { 0.5 } else { -0.5 }) as i64;
    let kf = k as f64;
    let r = (x - kf * P1) - kf * P2;
    let z = r * r;
    let s = r + r * z * (-1.666_666_666_666_663_243_48e-1 + z * (8.333_333_333_322_489_461_24e-3 + z * (-1.984_126_982_985_794_931_34e-4 + z * (2.755_731_370_707_006_767_89e-6 + z * (-2.505_076_025_340_686_341_95e-8 + z * 1.589_690_995_211_550_102_21e-10)))));
    let c = 1.0 - 0.5 * z + z * z * (4.166_666_666_666_660_190_37e-2 + z * (-1.388_888_888_887_410_957_49e-3 + z * (2.480_158_728_947_672_941_78e-5 + z * (-2.755_731_435_139_066_330_35e-7 + z * (2.087_572_321_298_174_827_90e-9 + z * -1.135_964_755_778_819_482_65e-11)))));
    match k & 3 { 0 => s, 1 => c, 2 => -s, _ => -c }
}

#[inline]
fn bern(t: f64) -> [f64; 4] { let s = 1.0 - t; [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t] }
/// synthetic "tessellate one face": random bicubic patch on an n x n grid, hashed
pub fn patch_work(seed: u64, n: u32) -> u64 {
    let mut r = Rng(seed);
    let mut cp = [[0.0; 3]; 16];
    for c in cp.iter_mut() { for x in c.iter_mut() { *x = r.f(); } }
    let d = 1.0 / (n - 1) as f64;
    let mut h = seed;
    for i in 0..n {
        let bu = bern(i as f64 * d);
        for j in 0..n {
            let bv = bern(j as f64 * d);
            let mut p = [0.0; 3];
            for a in 0..4 { for b in 0..4 { let w = bu[a] * bv[b]; for k in 0..3 { p[k] += w * cp[a * 4 + b][k]; } } }
            h = mix(h, p[0].to_bits() ^ p[1].to_bits().rotate_left(21) ^ p[2].to_bits().rotate_left(42));
        }
    }
    h
}

const MAG: [f64; 12] = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 1e1, 1e2, 1e3, 1e4, 1e5];
/// values of mixed magnitude (like real geometry sums)
pub fn fill(x: &mut [f64], seed: u64) { let mut r = Rng(seed); for v in x.iter_mut() { let a = r.f(); *v = a * MAG[(r.next() % 12) as usize]; } }
pub fn sum_seq(x: &[f64]) -> f64 { x.iter().fold(0.0, |s, v| s + v) }

/// hash of orient2d values on near-collinear triples (bit-exactness of a topology predicate)
pub fn orient_hash(seed: u64, n: u32) -> u64 {
    let mut r = Rng(seed);
    let mut h = 0;
    for _ in 0..n {
        let (ax, ay, bx, by, t) = (r.f(), r.f(), r.f(), r.f(), r.f());
        let (cx, cy) = (ax + t * (bx - ax), ay + t * (by - ay));
        h = mix(h, ((bx - ax) * (cy - ay) - (by - ay) * (cx - ax)).to_bits());
    }
    h
}
pub fn sin_hash(seed: u64, n: u32) -> u64 { let mut r = Rng(seed); let mut h = 0; for _ in 0..n { h = mix(h, own_sin(r.f() * 1000.0).to_bits()); } h }

/// scratch size (u64 words) for `mesh_work` on `nv` vertices
pub fn mesh_scratch_len(nv: usize) -> usize { (2 * nv).next_power_of_two() + 4 * nv }
/// per-part mesh processing on unindexed triangles (9 f32 each): exact-bit vertex welding, area, signed volume,
/// bbox, per-welded-vertex normals; returns a hash of all results
pub fn mesh_work(pos: &[f32], scratch: &mut [u64]) -> u64 {
    let nv = pos.len() / 3;
    let cap = (2 * nv).next_power_of_two();
    let (table, rest) = scratch.split_at_mut(cap);
    let (remap, nrm) = rest.split_at_mut(nv);
    table.fill(0);
    let key = |i: usize| (pos[3 * i].to_bits(), pos[3 * i + 1].to_bits(), pos[3 * i + 2].to_bits());
    let mut nu = 0usize;
    for i in 0..nv {
        let k = key(i);
        let mut h = (mix(mix(k.0 as u64, k.1 as u64), k.2 as u64) as usize) & (cap - 1);
        loop {
            let e = table[h];
            if e == 0 { table[h] = i as u64 + 1; remap[i] = nu as u64; nu += 1; break; }
            if key((e - 1) as usize) == k { remap[i] = remap[(e - 1) as usize]; break; }
            h = (h + 1) & (cap - 1);
        }
    }
    nrm[..3 * nu].fill(0f64.to_bits());
    let (mut area, mut vol) = (0.0f64, 0.0f64);
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for t in 0..nv / 3 {
        let p = |k: usize| [pos[9 * t + 3 * k] as f64, pos[9 * t + 3 * k + 1] as f64, pos[9 * t + 3 * k + 2] as f64];
        let (a, b, c) = (p(0), p(1), p(2));
        let (u, v) = ([b[0] - a[0], b[1] - a[1], b[2] - a[2]], [c[0] - a[0], c[1] - a[1], c[2] - a[2]]);
        let n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
        area += 0.5 * own_sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        vol += (a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0;
        for q in [a, b, c] { for k in 0..3 { lo[k] = lo[k].min(q[k]); hi[k] = hi[k].max(q[k]); } }
        for k in 0..3 { let w = remap[3 * t + k] as usize; for d in 0..3 { nrm[3 * w + d] = (f64::from_bits(nrm[3 * w + d]) + n[d]).to_bits(); } }
    }
    let mut h = mix(mix(mix(nv as u64, nu as u64), area.to_bits()), vol.to_bits());
    for k in 0..3 { h = mix(mix(h, lo[k].to_bits()), hi[k].to_bits()); }
    for w in 0..nu {
        let n = [f64::from_bits(nrm[3 * w]), f64::from_bits(nrm[3 * w + 1]), f64::from_bits(nrm[3 * w + 2])];
        let l = own_sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        let inv = if l > 0.0 { 1.0 / l } else { 0.0 };
        h = mix(h, (n[0] * inv).to_bits() ^ (n[1] * inv).to_bits().rotate_left(21) ^ (n[2] * inv).to_bits().rotate_left(42));
    }
    h
}
