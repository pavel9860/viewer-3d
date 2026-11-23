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
