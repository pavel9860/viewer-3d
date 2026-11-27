//! cgm IR... [--sag S] [--sheets] [--all-solids] [--bodies] [--stl OUT.stl]
//!
//! Meshes CATIA CGM B-reps (IR from `cadcli cgm PART OUT`) and prints one JSON line per part:
//! faces, left-out faces, triangles, area, volume, closure. By default the shells the part shows
//! (`cgm_mesh::shown_shells`); `--all-solids` meshes every solid; `--bodies` prints one line per body
//! with the volume of its solids and its flags instead.
use brep::cgm::Model;
use brep::cgm_mesh::{shown_shells, stats, tessellate};
use std::io::Write;
use std::time::Instant;

/// a JSON number with 6 decimals, null when not finite
fn num(x: f64) -> String { if x.is_finite() { format!("{x:.6}") } else { "null".into() } }

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut files = vec![];
    let (mut sag_rel, mut stl, mut sheets, mut all, mut bodies) = (1e-3, None, false, false, false);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--sag" => { sag_rel = args[i + 1].parse().expect("--sag S (fraction of the part's diagonal)"); i += 1; }
            "--stl" => { stl = Some(args[i + 1].clone()); i += 1; }
            "--sheets" => sheets = true,
            "--all-solids" => all = true,
            "--bodies" => bodies = true,
            f => files.push(f.to_string()),
        }
        i += 1;
    }
    if files.is_empty() { eprintln!("usage: cgm IR... [--sag S] [--sheets] [--all-solids] [--bodies] [--stl OUT.stl]"); std::process::exit(2); }
    let mut bad = 0;
    for f in &files {
        let t0 = Instant::now();
        let m = match std::fs::read(f).map_err(|e| e.to_string()).and_then(|b| Model::read(&b)) {
            Ok(m) => m,
            Err(e) => { println!("{{\"ir\":\"{f}\",\"error\":\"{e}\"}}"); bad += 1; continue; }
        };
        let solids: Vec<usize> = (0..m.shells.len()).filter(|&i| m.shells[i].solid).collect();
        // sag relative to the size of the solids' edge end points
        let b = stats(&tessellate(&m, f64::MAX, &solids)).bbox;
        let diag = ((b[3] - b[0]).powi(2) + (b[4] - b[1]).powi(2) + (b[5] - b[2]).powi(2)).sqrt();
        let sag = sag_rel * if diag.is_finite() && diag > 0.0 { diag } else { 1.0 };
        if bodies {
            for bd in &m.bodies {
                let sh: Vec<usize> = bd.shells.iter().copied().filter(|&i| m.shells[i].solid).collect();
                let st = stats(&tessellate(&m, sag, &sh));
                println!("{{\"ir\":\"{f}\",\"body\":{},\"flags\":{},\"solid_shells\":{},\"volume\":{:.6},\"open_edges\":{}}}", bd.id, bd.flags, sh.len(), st.volume, st.open_edges);
            }
            continue;
        }
        if brep::debug("gaps") { for (k, (n, g)) in brep::cgm_mesh::edge_gaps(&m) { eprintln!("{f}: gap (surface, p-curve) {:?} x {:?}: {n} edges, max {:.1e} of diag", k[0], k[1], g / sag * sag_rel); } }
        let shells = if all { solids.clone() } else { shown_shells(&m, sheets) };
        let mesh = tessellate(&m, sag, &shells);
        let s = stats(&mesh);
        let bb = s.bbox;
        let bb = bb.map(num);
        println!("{{\"ir\":\"{f}\",\"solids\":{},\"shells\":{},\"faces\":{},\"skipped\":{},\"failed\":{},\"tris\":{},\"area\":{},\"volume\":{},\"open_edges\":{},\"nonmanifold_edges\":{},\"bbox\":[{}],\"ms\":{:.1}}}",
            solids.len(), shells.len(), m.faces.len(), m.skipped.len(), mesh.failed.len(), s.tris, num(s.area), num(s.volume), s.open_edges, s.nonmanifold_edges, bb.join(","), t0.elapsed().as_secs_f64() * 1e3);
        if brep::debug("orient") {
            match brep::cgm_mesh::orientation_fix(&mesh) {
                Some(fix) => for (f, x) in mesh.faces.iter().zip(fix) { eprintln!("orient {} {} {} {} {} {}", f.0, f.1, f.2, f.3, f.4, x as i32); },
                None => eprintln!("orient none"),
            }
        }
        if brep::debug("nonman") {
            let mut u: std::collections::HashMap<(u32, u32), Vec<usize>> = Default::default();
            for (i, t) in mesh.tris.iter().enumerate() { for k in 0..3 { let (a, b) = (t[k], t[(k + 1) % 3]); u.entry((a.min(b), a.max(b))).or_default().push(i); } }
            for (e, ts) in u.iter().filter(|x| x.1.len() > 2) { eprintln!("non-manifold edge {:?} -> {:?}: faces {:?}", mesh.p[e.0 as usize], mesh.p[e.1 as usize], ts.iter().map(|&t| mesh.faces[mesh.tri_face[t] as usize].0).collect::<Vec<_>>()); }
        }
        if brep::debug("open") {
            let c = brep::cgm_mesh::open_causes(&mesh);
            eprintln!("open edges: {} on edges whose other face is not meshed ({} edges), {} on edges of two faces ({} edges: {:?}), {} inside faces", c.lone, c.lone_edges.len(), c.joined, c.joined_edges.len(), &c.joined_edges[..c.joined_edges.len().min(20)], c.inside);
        }
        if brep::debug("faces") {
            // largest faces by mesh area, with surface kind and p-curve kinds
            let mut area = vec![0.0f64; mesh.faces.len()];
            for (t, &f) in mesh.tris.iter().zip(&mesh.tri_face) {
                let [a, b, c] = t.map(|i| mesh.p[i as usize]);
                area[f as usize] += 0.5 * brep::geom::norm(brep::geom::cross(brep::geom::sub(b, a), brep::geom::sub(c, a)));
            }
            let by_id: std::collections::HashMap<u32, &brep::cgm::Face> = m.faces.iter().map(|f| (f.id, f)).collect();
            let mut order: Vec<usize> = (0..area.len()).collect();
            order.sort_by(|&x, &y| area[y].total_cmp(&area[x]));
            for &i in order.iter().take(5) {
                let fc = by_id[&mesh.faces[i].0];
                let mut pk: Vec<u32> = fc.loops.iter().flat_map(|l| l.clone()).map(|c| m.pcurve_kinds[m.coedges[c].pcurve]).collect();
                pk.sort_unstable();
                pk.dedup();
                eprintln!("face {} surface kind {} p-curves {pk:?} area {:.3e}", fc.id, m.surface_kinds[fc.surface], area[i]);
            }
        }
        for (id, why) in &mesh.failed { eprintln!("{f}: face {id}: {why}"); }
        if let Some(out) = &stl {
            let mut w = std::io::BufWriter::new(std::fs::File::create(out).expect("cannot create the STL"));
            w.write_all(&[0u8; 80]).unwrap();
            w.write_all(&(mesh.tris.len() as u32).to_le_bytes()).unwrap();
            for t in &mesh.tris {
                w.write_all(&[0u8; 12]).unwrap();
                for &k in t { for c in mesh.p[k as usize] { w.write_all(&(c as f32).to_le_bytes()).unwrap(); } }
                w.write_all(&[0u8; 2]).unwrap();
            }
        }
    }
    std::process::exit((bad > 0) as i32);
}
