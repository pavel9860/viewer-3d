//! Debug one imported face: loops in UV (areas, extents, items), triangulation area.
//! Usage: face_dbg FILE.stp STEP_FACE_ID [SAG] [DUMP_FILE]
use brep::{step, trim::area};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let im = step::import(&a[1]).unwrap();
    let id: u32 = a[2].parse().unwrap();
    let (fid, fc) = im.body.faces.iter().find(|(_, f)| f.origin.src == id).unwrap();
    let sag: f64 = a.get(3).and_then(|x| x.parse().ok()).unwrap_or(brep::trim::SAG);
    let loops: Vec<Vec<_>> = im.body.face_loops(fid);
    let t = brep::trim::Trim::with_sag(&im.body, fc.surface.as_ref(), fc.sense, &loops, sag);
    let rings: Vec<Vec<[f64; 2]>> = t.polygons().to_vec();
    let (pts, tris) = brep::tri::triangulate(&rings);
    let want: f64 = rings.iter().map(|r| area(r)).sum();
    let got: f64 = tris.iter().map(|x| area(&[pts[x[0]], pts[x[1]], pts[x[2]]])).sum();
    println!("triangulation: {} triangles, area {got:.6e} of {want:.6e}", tris.len());
    if let Some(d) = a.get(4) { std::fs::write(d, rings.iter().map(|r| r.iter().map(|q| format!("{} {}", q[0], q[1])).collect::<Vec<_>>().join(";")).collect::<Vec<_>>().join("\n")).unwrap(); }
    println!("face #{id} sense {} dom {:?} per {:?} sg {}", fc.sense, fc.surface.dom(), t.per, t.sg);
    for (i, (l, p)) in t.loops.iter().zip(t.polygons()).enumerate() {
        println!(" loop {i}: area {:.4e}, {} pts, items {}", area(p), p.len(), l.len());
        for it in l { println!("   edge #{} sense {} {} pts uv {:?} -> {:?}", im.body.edges.get(it.e).origin.src, it.sense, it.uv.len(), it.uv[0], it.uv[it.uv.len() - 1]); }
    }
}
