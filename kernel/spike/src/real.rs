//! Real-data workloads on CATIA parts and an assembly exported by tools/export_parts.py:
//! per-part mesh processing (scaling + determinism vs the browser), hashing, transfer sizes, assembly picking.
//! Also packs the parts into DIR/parts.bin + parts.idx for the browser test.
use crate::bench::{as_bytes, par_for, slab, Bvh, WORKERS};
use crate::ms;
use kcore::{hash_bytes, mesh_scratch_len, mesh_work, Rng};
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::time::Instant;

/// vertex positions of a binary little-endian PLY written by cadcli (x,y,z f32, unindexed triangles)
fn ply_positions(path: &str) -> Vec<f32> {
    let b = std::fs::read(path).unwrap();
    let he = b.windows(11).position(|w| w == b"end_header\n").unwrap() + 11;
    let head = std::str::from_utf8(&b[..he]).unwrap();
    assert!(head.contains("property float x\nproperty float y\nproperty float z\nelement face"), "unexpected vertex layout in {path}");
    let nv: usize = head.lines().find(|l| l.starts_with("element vertex")).unwrap().split_whitespace().nth(2).unwrap().parse().unwrap();
    b[he..he + nv * 12].chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()
}

fn tri_hit(o: [f32; 3], d: [f32; 3], t: &[f32]) -> f32 {
    let s = |a: [f32; 3], b: [f32; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let c = |a: [f32; 3], b: [f32; 3]| [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
    let dt = |a: [f32; 3], b: [f32; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let p0 = [t[0], t[1], t[2]];
    let (e1, e2) = (s([t[3], t[4], t[5]], p0), s([t[6], t[7], t[8]], p0));
    let p = c(d, e2);
    let det = dt(e1, p);
    if det == 0.0 { return f32::MAX; }
    let tv = s(o, p0);
    let u = dt(tv, p) / det;
    if !(0.0..=1.0).contains(&u) { return f32::MAX; }
    let q = c(tv, e1);
    let v = dt(d, q) / det;
    if v < 0.0 || u + v > 1.0 { return f32::MAX; }
    let tt = dt(e2, q) / det;
    if tt > 0.0 { tt } else { f32::MAX }
}

pub fn run(dir: &str) {
    println!("== real CATIA data ({dir}) ==");
    let manifest = std::fs::read_to_string(format!("{dir}/parts/manifest.tsv")).unwrap();
    let parts: Vec<Vec<f32>> = manifest.lines().filter(|l| !l.is_empty()).map(|l| ply_positions(&format!("{dir}/parts/{}", l.split('\t').next().unwrap()))).filter(|p| !p.is_empty()).collect();
    let tris: usize = parts.iter().map(|p| p.len() / 9).sum();
    println!("  {} non-empty parts, {tris} triangles", parts.len());
    // pack for the browser test
    let (mut bin, mut idx, mut off) = (Vec::with_capacity(tris * 36), vec![], 0u64);
    for p in &parts { bin.extend_from_slice(as_bytes(p)); idx.extend_from_slice(&off.to_le_bytes()); idx.extend_from_slice(&(p.len() as u64).to_le_bytes()); off += p.len() as u64; }
    std::fs::write(format!("{dir}/parts.bin"), &bin).unwrap();
    std::fs::write(format!("{dir}/parts.idx"), &idx).unwrap();
    // per-part processing, largest first, 1..WORKERS threads
    let mut order: Vec<usize> = (0..parts.len()).collect();
    order.sort_by_key(|&i| (std::cmp::Reverse(parts[i].len()), i));
    let mut base = 0.0;
    for t in 1..=WORKERS {
        let res: Vec<AtomicU64> = (0..parts.len()).map(|_| AtomicU64::new(0)).collect();
        let t0 = Instant::now();
        par_for(order.len(), t, |k| { let i = order[k]; let mut scr = vec![0u64; mesh_scratch_len(parts[i].len() / 3)]; res[i].store(mesh_work(&parts[i], &mut scr), Relaxed); });
        let e = ms(t0);
        if t == 1 { base = e; }
        let v: Vec<u64> = res.iter().map(|x| x.load(Relaxed)).collect();
        println!("  per-part mesh processing, threads={t}: {e:6.1} ms speedup {:.2}x  results hash {:016x}", base / e, hash_bytes(as_bytes(&v)));
    }
    let t0 = Instant::now();
    let hx = parts.iter().fold(0u64, |h, p| h ^ hash_bytes(as_bytes(p)));
    let e = ms(t0);
    println!("  hashing all parts: {:.1} MB in {e:.1} ms = {:.2} GB/s (xor {hx:016x})", bin.len() as f64 / 1e6, bin.len() as f64 / e / 1e6);
    let mut sizes: Vec<usize> = parts.iter().map(|p| p.len() * 4).collect();
    sizes.sort_unstable();
    let pct = |q: f64| sizes[((sizes.len() - 1) as f64 * q) as usize] as f64 / 1e6;
    println!("  part mesh size p50 {:.3} MB, p95 {:.3} MB, max {:.3} MB", pct(0.5), pct(0.95), pct(1.0));
    // assembly picking on world-space triangles, checked against brute force
    let Ok(_) = std::fs::metadata(format!("{dir}/asm.ply")) else { return };
    let w = ply_positions(&format!("{dir}/asm.ply"));
    let nt = w.len() / 9;
    let bound = |t: usize| { let (mut l, mut h) = ([f32::MAX; 3], [f32::MIN; 3]); for v in 0..3 { for k in 0..3 { let x = w[9 * t + 3 * v + k]; l[k] = l[k].min(x); h[k] = h[k].max(x); } } (l, h) };
    let (lo, hi): (Vec<_>, Vec<_>) = (0..nt).map(bound).unzip();
    let t0 = Instant::now();
    let bvh = Bvh::build(&lo, &hi, 4);
    let tb = ms(t0);
    let (l0, h0) = (bvh.lo[0], bvh.hi[0]);
    let mut r = Rng(5);
    let rays: Vec<([f32; 3], [f32; 3])> = (0..2000).map(|_| ([l0[0] + (h0[0] - l0[0]) * (0.5 + 0.5 * r.f() as f32), l0[1] + (h0[1] - l0[1]) * (0.5 + 0.5 * r.f() as f32), h0[2] + 100.0], [0.05 * r.f() as f32, 0.05 * r.f() as f32, -1.0])).collect();
    let pick = |o: [f32; 3], d: [f32; 3]| {
        let inv = d.map(|x| 1.0 / x);
        let mut best = f32::MAX;
        let mut st = vec![0usize];
        while let Some(n) = st.pop() {
            if slab(bvh.lo[n], bvh.hi[n], o, inv) >= best { continue; }
            if bvh.count[n] > 0 { for &t in &bvh.items[bvh.first[n] as usize..(bvh.first[n] + bvh.count[n]) as usize] { best = best.min(tri_hit(o, d, &w[9 * t as usize..9 * t as usize + 9])); } }
            else { st.push(bvh.right[n] as usize); st.push(n + 1); }
        }
        best
    };
    let t0 = Instant::now();
    let hits: Vec<f32> = rays.iter().map(|&(o, d)| pick(o, d)).collect();
    let tp = ms(t0) * 1e3 / rays.len() as f64;
    let mismatch = rays.iter().take(100).zip(&hits).filter(|((o, d), h)| (0..nt).map(|t| tri_hit(*o, *d, &w[9 * t..9 * t + 9])).fold(f32::MAX, f32::min) != **h).count();
    println!("  assembly: {nt} triangles, BVH {tb:.0} ms, pick {tp:.2} us/ray ({}/{} hit), vs brute force on 100 rays: {mismatch} mismatches", hits.iter().filter(|h| **h < f32::MAX).count(), rays.len());
}
