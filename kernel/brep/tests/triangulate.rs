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
