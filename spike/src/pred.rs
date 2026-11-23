//! Exact predicate benchmark (agreement with the exact fallback on random and near-degenerate inputs).
use brep::exact::*;
/// benchmark + agreement check against the exact predicates on random and near-degenerate inputs
pub fn bench() -> usize {
    use kcore::Rng;
    use std::hint::black_box;
    use std::time::Instant;
    println!("== exact predicates (filter -> staged -> exact) ==");
    let mut bad = 0;
    for degen in [false, true] {
        let n = 1_000_000;
        let mut r = Rng(11);
        let pts: Vec<_> = (0..n).map(|_| { let a = [r.f(), r.f()]; let b = [r.f(), r.f()]; let c = if degen { let t = r.f(); [a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])] } else { [r.f(), r.f()] }; (a, b, c) }).collect();
        let t0 = Instant::now(); let plain: Vec<i8> = pts.iter().map(|&(a, b, c)| sg((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))).collect(); let tp = t0.elapsed().as_secs_f64();
        let mut st = [0u64; 4];
        let t0 = Instant::now(); let ad: Vec<i8> = pts.iter().map(|&(a, b, c)| { let (s, k) = orient2d(a, b, c); st[k as usize] += 1; s }).collect(); let ta = t0.elapsed().as_secs_f64();
        let ex: Vec<i8> = pts.iter().map(|&(a, b, c)| orient2d_exact(a, b, c)).collect();
        let wrong = ad.iter().zip(&ex).filter(|(x, y)| x != y).count();
        bad += wrong;
        println!("  orient2d {:<15} plain {:4.1} ns  adaptive {:5.1} ns  stages {:?}  plain wrong {:4.1}%  adaptive wrong {wrong}", if degen { "near-degenerate" } else { "random" }, tp * 1e9 / n as f64, ta * 1e9 / n as f64, st, 100.0 * plain.iter().zip(&ex).filter(|(x, y)| x != y).count() as f64 / n as f64);
        black_box(&plain);
    }
    for degen in [false, true] {
        let n = 200_000;
        let mut r = Rng(12);
        let pts: Vec<_> = (0..n).map(|_| { let a = [r.f(), r.f(), r.f()]; let b = [r.f(), r.f(), r.f()]; let c = [r.f(), r.f(), r.f()]; let p = if degen { let (s, t) = (r.f(), r.f()); [0, 1, 2].map(|k| a[k] + s * (b[k] - a[k]) + t * (c[k] - a[k])) } else { [r.f(), r.f(), r.f()] }; (a, b, c, p) }).collect();
        let mut st = [0u64; 4];
        let t0 = Instant::now(); let ad: Vec<i8> = pts.iter().map(|&(a, b, c, p)| { let (s, k) = orient3d(a, b, c, p); st[k as usize] += 1; s }).collect(); let ta = t0.elapsed().as_secs_f64();
        let t0 = Instant::now(); let ex: Vec<i8> = pts.iter().map(|&(a, b, c, p)| orient3d_exact(a, b, c, p)).collect(); let te = t0.elapsed().as_secs_f64();
        let wrong = ad.iter().zip(&ex).filter(|(x, y)| x != y).count();
        bad += wrong;
        println!("  orient3d {:<15} adaptive {:6.1} ns  exact {:6.1} ns  stages {:?}  adaptive wrong {wrong}", if degen { "near-degenerate" } else { "random" }, ta * 1e9 / n as f64, te * 1e9 / n as f64, st);
    }
    bad
}

