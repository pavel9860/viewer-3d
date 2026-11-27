//! Triangulation of real face polygons (UV loops dumped from CATIA STEP faces with BREP_DUMP).
use brep::tri::triangulate;
use brep::trim::area;

fn load(name: &str) -> Vec<Vec<[f64; 2]>> {
    std::fs::read_to_string(format!("{}/tests/data/{name}", env!("CARGO_MANIFEST_DIR"))).unwrap().lines()
        .map(|l| l.split(';').map(|p| { let v: Vec<f64> = p.split(' ').map(|x| x.parse().unwrap()).collect(); [v[0], v[1]] }).collect()).collect()
}
fn check(name: &str) {
    let r = load(name);
    let (pts, tris) = triangulate(&r);
    let want = area(&r[0]) + r[1..].iter().map(|h| area(h)).sum::<f64>();
    let got: f64 = tris.iter().map(|t| area(&[pts[t[0]], pts[t[1]], pts[t[2]]])).sum();
    let neg = tris.iter().filter(|t| area(&[pts[t[0]], pts[t[1]], pts[t[2]]]) < 0.0).count();
    let bad: Vec<_> = tris.iter().filter(|t| { let c = [(pts[t[0]][0] + pts[t[1]][0] + pts[t[2]][0]) / 3.0, (pts[t[0]][1] + pts[t[1]][1] + pts[t[2]][1]) / 3.0]; !brep::trim::in_poly(c, &r[0]) || r[1..].iter().any(|h| brep::trim::in_poly(c, h)) }).map(|t| [pts[t[0]], pts[t[1]], pts[t[2]]]).collect();
    for t in bad.iter().take(4) { eprintln!("triangle outside the region: {t:?}"); }
    // no triangle edge may cross a ring edge
    let rings: Vec<([f64; 2], [f64; 2])> = r.iter().flat_map(|q| (0..q.len()).map(move |i| (q[i], q[(i + 1) % q.len()]))).collect();
    let pr = &pts;
    let crossing = tris.iter().flat_map(|t| (0..3).map(move |k| (pr[t[k]], pr[t[(k + 1) % 3]]))).filter(|(a, b)| rings.iter().any(|(c, d)| {
        let o = |p: [f64; 2], q: [f64; 2], r: [f64; 2]| brep::exact::orient2d(p, q, r).0;
        o(*a, *b, *c) * o(*a, *b, *d) < 0 && o(*c, *d, *a) * o(*c, *d, *b) < 0
    })).count();
    if brep::debug("tri") {
        for t in &tris { for k in 0..3 { let (a, b) = (pts[t[k]], pts[t[(k + 1) % 3]]); for (c, d) in &rings {
            let o = |p: [f64; 2], q: [f64; 2], r: [f64; 2]| brep::exact::orient2d(p, q, r).0;
            if o(a, b, *c) * o(a, b, *d) < 0 && o(*c, *d, a) * o(*c, *d, b) < 0 { eprintln!("triangle {:?} edge {a:?}-{b:?} crosses ring edge {c:?}-{d:?}", t); }
        } } }
    }
    let flipped = tris.iter().filter(|t| brep::exact::orient2d(pts[t[0]], pts[t[1]], pts[t[2]]).0 < 0).count();
    assert_eq!((crossing, flipped), (0, 0), "{name}: triangle edges crossing the boundary, non-positive triangles");
    assert!((got - want).abs() <= 1e-9 * want, "{name}: area {got} vs {want}, {} triangles, {neg} inverted", tris.len());
}
#[test] fn face_with_two_holes() { check("face_2holes.txt") }
#[test] fn square_with_square_hole() {
    let r = vec![vec![[0.0, 0.0], [4.0, 0.0], [4.0, 4.0], [0.0, 4.0]], vec![[1.0, 1.0], [1.0, 3.0], [3.0, 3.0], [3.0, 1.0]]];
    let (pts, t) = triangulate(&r);
    assert_eq!(t.len(), 8);
    assert!((t.iter().map(|t| area(&[pts[t[0]], pts[t[1]], pts[t[2]]])).sum::<f64>() - 12.0).abs() < 1e-12);
}
#[test] fn rectangle_with_noisy_straight_edges() { check("rect_noisy_edges.txt") }

/// conforming: every ring point is a vertex, every ring segment a triangle edge, no edge in more than two triangles
/// (exactly one on the boundary), no degenerate triangles, area kept
fn conforming(name: &str, r: &[Vec<[f64; 2]>]) {
    let (pts, tris) = triangulate(r);
    let key = |p: [f64; 2]| (p[0].to_bits(), p[1].to_bits());
    let id: std::collections::HashMap<_, usize> = pts.iter().enumerate().map(|(i, &p)| (key(p), i)).collect();
    let mut uses: std::collections::HashMap<(usize, usize), usize> = Default::default();
    for t in &tris { for k in 0..3 { let (a, b) = (t[k], t[(k + 1) % 3]); *uses.entry((a.min(b), a.max(b))).or_default() += 1; } }
    let mut ring = std::collections::HashSet::new();
    for q in r { for i in 0..q.len() { let (a, b) = (id[&key(q[i])], id[&key(q[(i + 1) % q.len()])]); ring.insert((a.min(b), a.max(b))); } }
    for e in &ring {
        if uses.get(e) != Some(&1) {
            for t in tris.iter().filter(|t| t.contains(&e.0) || t.contains(&e.1)) { eprintln!("  triangle {:?}", t.map(|i| pts[i])); }
        }
        assert_eq!(uses.get(e), Some(&1), "{name}: ring segment {:?}-{:?}", pts[e.0], pts[e.1]);
    }
    for (e, n) in &uses { assert!(ring.contains(e) || *n == 2, "{name}: inner edge {:?}-{:?} in {n} triangles", pts[e.0], pts[e.1]); }
    assert!(tris.iter().all(|t| brep::exact::orient2d(pts[t[0]], pts[t[1]], pts[t[2]]).0 > 0), "{name}: degenerate or inverted triangle");
    let want: f64 = r.iter().map(|q| area(q)).sum();
    let got: f64 = tris.iter().map(|t| area(&[pts[t[0]], pts[t[1]], pts[t[2]]])).sum();
    assert!((got - want).abs() <= 1e-9 * want, "{name}: area {got} vs {want}");
}

#[test]
fn collinear_runs_stay_conforming() {
    // edges sampled into many collinear points (straight UV lines of cylinders and planes)
    // points o + t d as a straight p-curve gives them (not exactly collinear in floating point)
    let run = |a: [f64; 2], b: [f64; 2], n: usize| (0..n).map(move |k| { let t = 0.1 + 0.7 * k as f64 / n as f64; let d = [(b[0] - a[0]) / 0.7, (b[1] - a[1]) / 0.7]; let o = [a[0] - 0.1 * d[0], a[1] - 0.1 * d[1]]; [o[0] + t * d[0], o[1] + t * d[1]] });
    let rect: Vec<[f64; 2]> = run([0.0, 0.0], [12.5, 0.0], 8).chain(run([12.5, 0.0], [12.5, 1.8], 5)).chain(run([12.5, 1.8], [0.0, 1.8], 8)).chain(run([0.0, 1.8], [0.0, 0.0], 5)).collect();
    conforming("rectangle", &[rect.clone()]);
    // with a square hole whose sides are also runs
    let hole: Vec<[f64; 2]> = run([2.0, 0.5], [2.0, 1.2], 3).chain(run([2.0, 1.2], [4.0, 1.2], 4)).chain(run([4.0, 1.2], [4.0, 0.5], 3)).chain(run([4.0, 0.5], [2.0, 0.5], 4)).collect();
    conforming("rectangle with hole", &[rect, hole]);
    // L shape: a long side with a run meeting a reflex corner
    let l: Vec<[f64; 2]> = run([0.0, 0.0], [10.0, 0.0], 10).chain(run([10.0, 0.0], [10.0, 3.0], 3)).chain(run([10.0, 3.0], [3.0, 3.0], 7)).chain(run([3.0, 3.0], [3.0, 10.0], 7)).chain(run([3.0, 10.0], [0.0, 10.0], 3)).chain(run([0.0, 10.0], [0.0, 0.0], 10)).collect();
    conforming("L", &[l]);
}

#[test]
fn plate_with_two_holes_and_collinear_runs() { conforming("face_collinear", &load("face_collinear.txt")) }

#[test]
fn face_with_sampled_straight_edges() { conforming("face_flat", &load("face_flat.txt")) }
