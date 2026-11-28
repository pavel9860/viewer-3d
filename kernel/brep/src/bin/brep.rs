//! brep CLI.
//!   brep check FILE.stp|LIST.tsv          import + validate
//!   brep bool FILE.stp|LIST.tsv [SHELLS]  booleans of the largest closed shells with a cylinder through them
//!   brep mesh FILE.stp|LIST.tsv [CHORD] [--out DIR]
//!       sew open shells, tessellate (chord relative to the model diagonal, default 1e-3); --out writes
//!       DIR/<name>.mesh with detail levels and DIR/index.json for the web viewer
use brep::boolean::{boolean, Op};
use brep::geom::*;
use brep::topo::*;
use brep::{check, make, mesh, step};
use std::time::Instant;

const WORKERS: usize = 4;
fn ms(t: Instant) -> f64 { t.elapsed().as_secs_f64() * 1e3 }
fn base(f: &str) -> &str { f.rsplit('/').next().unwrap() }
fn import(f: &str) -> Option<step::Import> { step::import(f).map_err(|e| println!("{f}: {e}")).ok() }

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 { eprintln!("usage: brep check|bool|mesh FILE.stp|LIST.tsv [SHELLS|CHORD] [--out DIR]"); std::process::exit(2); }
    let files: Vec<String> = if args[2].ends_with(".tsv") { std::fs::read_to_string(&args[2]).unwrap().lines().filter(|l| !l.is_empty()).map(|l| l.rsplit('\t').next().unwrap().to_string()).collect() } else { vec![args[2].clone()] };
    let bad = match args[1].as_str() {
        "check" => files.iter().map(|f| check_file(f)).sum::<usize>(),
        "bool" => { let n = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(2); files.iter().map(|f| bool_file(f, n)).sum() }
        "mesh" => {
            let out = args.iter().position(|a| a == "--out").and_then(|i| args.get(i + 1)).cloned();
            let c = args.get(3).filter(|a| !a.starts_with("--")).and_then(|s| s.parse().ok()).unwrap_or(1e-3);
            files.iter().map(|f| mesh_file(f, c, out.as_deref())).sum()
        }
        _ => { eprintln!("unknown command {}", args[1]); 2 }
    };
    std::process::exit((bad > 0) as i32);
}

fn check_file(f: &str) -> usize {
    let t0 = Instant::now();
    let Some(im) = import(f) else { return 1 };
    let ti = ms(t0);
    let t1 = Instant::now();
    let r = check::check(&im.body, 7, WORKERS);
    let mm = |x: f64| x * im.diag;
    println!("{:<26} {} | import {ti:.0} ms, check {:.0} ms | violations {} (open loops {}, edge use {}, orientation {}, vertex off curve {}, edge off surface {}) | non-manifold {} | max dev vertex {:.1e} mm edge {:.1e} mm | no curve {} | unsupported {:?}",
        base(f), im.body.stats(), ms(t1), r.violations(), r.loops_open, r.edge_use, r.orientation, r.vertex_off_curve, r.edge_off_surface, r.nonmanifold, mm(r.max_vertex_dev), mm(r.max_edge_dev), r.edges_without_curve, im.unsupported);
    r.violations()
}

/// A (op) B for the `n` largest closed shells A of the file and a z-cylinder B through A's box centre.
/// No ground truth: checks are the validator, arrangement health and the volume identities
/// V(A-B) + V(A^B) = V(A) and V(AuB) = V(A) + V(B) - V(A^B).
fn bool_file(f: &str, n: usize) -> usize {
    let Some(im) = import(f) else { return 1 };
    let b = &im.body;
    let mut shells: Vec<&Shell> = b.shells.iter().map(|x| x.1).filter(|s| s.closed).collect();
    shells.sort_by_key(|s| std::cmp::Reverse(s.faces.len()));
    let mut bad = 0;
    for (si, sh) in shells.iter().take(n).enumerate() {
        let a = b.extract(&sh.faces);
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for (_, v) in a.vertices.iter() { for k in 0..3 { lo[k] = lo[k].min(v.p[k]); hi[k] = hi[k].max(v.p[k]); } }
        let c = mulv(addv(lo, hi), 0.5);
        let r = 0.2 * (hi[0] - lo[0]).min(hi[1] - lo[1]).max(1e-3);
        let cyl = make::cylinder(c, r, (hi[2] - lo[2]).max(r), [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        let vol = |x: &Body| mesh::volume(x, &x.faces.ids().collect::<Vec<_>>());
        let (va, ba) = vol(&a);
        let vb = vol(&cyl).0;
        let mut v = [0.0; 3];
        let mut line = format!("{:<26} shell {si} ({} faces) V(A) {va:.4e}{}", base(f), sh.faces.len(), if ba > 0 { format!(" [{ba} faces lost area]") } else { String::new() });
        for (k, op) in [Op::Difference, Op::Intersection, Op::Union].into_iter().enumerate() {
            let t0 = Instant::now();
            let (res, st) = boolean(&a, &cyl, op, WORKERS);
            let t = ms(t0);
            let rep = check::check(&res, 5, WORKERS);
            let (vr, br) = vol(&res);
            v[k] = vr;
            if brep::debug("check") && rep.edge_use > 0 { unmatched(&res, &a, &sh.faces, b); }
            let fail = st.dangling + st.on + rep.violations() + br;
            bad += fail;
            line += &format!(" | {op:?} {} faces {t:.0} ms{}", res.faces.len(), if fail > 0 { format!(" FAIL dangling {} on {} violations {} (open {}, use {}, orient {}, vtx {}, edge {}) lost-area {br}", st.dangling, st.on, rep.violations(), rep.loops_open, rep.edge_use, rep.orientation, rep.vertex_off_curve, rep.edge_off_surface) } else { String::new() });
        }
        let id1 = (v[0] + v[1] - va).abs() / va.abs().max(1e-30);
        let id2 = (v[2] - (va + vb - v[1])).abs() / v[2].abs().max(1e-30);
        println!("{line} | identities {id1:.1e} {id2:.1e}");
    }
    bad
}

/// debug: closed-shell edges with unbalanced use, with the STEP ids of the faces they came from
fn unmatched(r: &Body, a: &Body, shell: &[Id<Face>], src: &Body) {
    let mut uses: std::collections::HashMap<Id<Edge>, Vec<(u32, bool)>> = Default::default();
    for (fid, _) in r.faces.iter() { for ce in r.face_coedges(fid) { uses.entry(ce.edge).or_default().push((fid.0, ce.sense)); } }
    let step_of = |wf: u32| if (wf as usize) < a.faces.len() { format!("A#{}", src.faces.get(shell[wf as usize]).origin.src) } else { format!("B{}", wf as usize - a.faces.len()) };
    for (e, u) in uses { let f = u.iter().filter(|x| x.1).count(); if 2 * f != u.len() {
        let ed = r.edges.get(e);
        let fs: Vec<String> = u.iter().map(|(fi, s)| format!("{}{}", step_of(r.faces.get(Id::new(*fi as usize)).origin.src - 1), if *s { "+" } else { "-" })).collect();
        eprintln!("  unmatched edge {e:?} origin {:?} {:?} -> {:?} faces {fs:?}", ed.origin, r.vertices.get(ed.v[0]).p, r.vertices.get(ed.v[1]).p);
    } }
}

/// tessellate every shell; closed shells must give closed meshes. Volume and area vs the exact face integrals are
/// indicators (area per face catches facets that cut across their surface, which volume does not).
fn mesh_file(f: &str, chord: f64, out: Option<&str>) -> usize {
    use brep::tess::{mesh_volume, open_edges, tessellate, Opts};
    let Some(mut im) = import(f) else { return 1 };
    let tol = im.body.tol;
    let sewn = brep::sew::sew(&mut im.body, tol);
    let sew_note = if sewn.free_before > 0 {
        format!(" | sewn: {} free edges -> {}, {} shells -> {} ({} closed), {} orientation conflicts", sewn.free_before, sewn.free_after, sewn.shells_before, sewn.shells_after, sewn.closed_after, sewn.conflicts)
    } else { String::new() };
    let b = &im.body;
    let t0 = Instant::now();
    let o = Opts { chord, ..Default::default() };
    let (mut tris, mut open, mut nonman, mut empty, mut coarse, mut dev, mut vrel) = (0, 0, 0, 0, 0, 0.0f64, 0.0f64);
    let (mut mv, mut ev) = (0.0, 0.0);
    let (mut ma, mut ea, mut arel) = (0.0, 0.0, 0.0f64);
    let mut all = brep::tess::Mesh::default();
    for (_, sh) in b.shells.iter() {
        let (m, st) = tessellate(b, &sh.faces, &o, WORKERS);
        if out.is_some() {
            let (base, fbase) = (all.pos.len() as u32, all.face.iter().max().map_or(0, |x| x + 1));
            all.pos.extend(&m.pos);
            all.tri.extend(m.tri.iter().map(|t| t.map(|i| i + base)));
            all.face.extend(m.face.iter().map(|x| x + fbase));
        }
        tris += m.tri.len();
        // area per face: mesh vs exact
        let mut fa = vec![0.0f64; sh.faces.len()];
        for (t, &fi) in m.tri.iter().zip(&m.face) { let [p, q, r] = t.map(|i| m.pos[i as usize]); fa[fi as usize] += 0.5 * norm(cross(sub(q, p), sub(r, p))); }
        let ex = brep::par::map(sh.faces.len(), WORKERS, |i| mesh::area_of(b, &sh.faces[i..i + 1]).0);
        for (a, e) in fa.iter().zip(&ex) { ma += a; ea += e; arel = arel.max((a - e).abs() / e.abs().max(1e-30)); }
        empty += st.empty;
        coarse += st.coarse;
        dev = dev.max(st.max_dev);
        if sh.closed {
            let (op, nm) = open_edges(&m);
            if op + nm > 0 && brep::debug("check") { mesh_defects(b, &sh.faces, &m); }
            open += op;
            nonman += nm;
            let (v, _) = mesh::volume(b, &sh.faces);
            mv += mesh_volume(&m);
            ev += v;
            vrel = vrel.max((mesh_volume(&m) - v).abs() / v.abs().max(1e-30));
        }
    }
    let t = ms(t0);
    let mut lods = String::new();
    if let Some(dir) = out {
        // tolerances from a quarter of the chord, doubling
        let t1 = Instant::now();
        let tols: Vec<f64> = (0..10).map(|k| 0.25 * chord * (1u64 << k) as f64).collect();
        let lv = brep::simplify::levels(&all.pos, &all.tri, &all.face, &tols, WORKERS);
        lods = format!(" | levels {} in {:.0} ms", lv.iter().map(|l| format!("{}@{:.1e}", l.tri.len(), l.err)).collect::<Vec<_>>().join(" "), ms(t1));
        write_mesh(dir, base(f), &all, chord, &lv, im.diag);
    }
    let bad = open + empty + coarse;
    println!("{:<26} {} faces -> {tris} triangles in {t:.0} ms | open edges {open} non-manifold {nonman} | empty faces {empty} coarse {coarse} | max dev {:.1e} mm | volume mesh/exact {:.6} (worst shell rel {vrel:.1e}) | area mesh/exact {:.6} (worst face rel {arel:.1e}){sew_note}{lods}{}",
        base(f), b.faces.len(), dev * im.diag, mv / ev, ma / ea, if bad > 0 { " FAIL" } else { "" });
    bad
}

/// debug: non-manifold mesh edges and open mesh edges by face (STEP ids)
fn mesh_defects(b: &Body, faces: &[Id<Face>], m: &brep::tess::Mesh) {
    let src = |i: usize| b.faces.get(faces[m.face[i] as usize]).origin.src;
    let mut d: std::collections::HashMap<(u32, u32), Vec<usize>> = Default::default();
    for (i, t) in m.tri.iter().enumerate() { for k in 0..3 { d.entry((t[k], t[(k + 1) % 3])).or_default().push(i); } }
    let mut open: std::collections::BTreeMap<u32, usize> = Default::default();
    for (e, ts) in &d {
        let back = d.get(&(e.1, e.0)).map_or(0, |x| x.len());
        if back == 0 { for &i in ts { *open.entry(src(i)).or_default() += 1; } }
        if e.0 < e.1 && ts.len() + back > 2 { eprintln!("  non-manifold mesh edge {:?} -> {:?} faces STEP {:?}", m.pos[e.0 as usize], m.pos[e.1 as usize], ts.iter().chain(d.get(&(e.1, e.0)).into_iter().flatten()).map(|&i| src(i)).collect::<Vec<_>>()); }
    }
    if !open.is_empty() { eprintln!("  open edges by face (STEP id: count): {open:?}"); }
}

/// little-endian block: u32 vertex count, u32 triangle count, f32 xyz per vertex, u32 x3 per triangle, u32 face id
/// per triangle; then "LODS", f32 chord (deviation of the mesh from the exact model), u32 level count, and per level
/// f32 error (measured deviation from the mesh) and a block. Lengths in mm. DIR/index.json lists the files.
fn write_mesh(dir: &str, name: &str, m: &brep::tess::Mesh, err: f64, lv: &[brep::simplify::Level], scale: f64) {
    std::fs::create_dir_all(dir).unwrap();
    let mut b: Vec<u8> = vec![];
    let block = |b: &mut Vec<u8>, pos: &[[f64; 3]], tri: &[[u32; 3]], face: &[u32]| {
        b.extend((pos.len() as u32).to_le_bytes());
        b.extend((tri.len() as u32).to_le_bytes());
        for p in pos { for k in 0..3 { b.extend(((p[k] * scale) as f32).to_le_bytes()); } }
        for t in tri { for k in 0..3 { b.extend(t[k].to_le_bytes()); } }
        for f in face { b.extend(f.to_le_bytes()); }
    };
    block(&mut b, &m.pos, &m.tri, &m.face);
    b.extend(b"LODS");
    b.extend(((err * scale) as f32).to_le_bytes());
    b.extend((lv.len() as u32).to_le_bytes());
    for l in lv { b.extend(((l.err * scale) as f32).to_le_bytes()); block(&mut b, &l.pos, &l.tri, &l.tag); }
    let file = format!("{}.mesh", name.trim_end_matches(".stp"));
    std::fs::write(format!("{dir}/{file}"), b).unwrap();
    let idx = format!("{dir}/index.json");
    let mut list: Vec<String> = std::fs::read_to_string(&idx).ok().map(|s| s.trim().trim_start_matches('[').trim_end_matches(']').split(',').map(|x| x.trim().trim_matches('"').to_string()).filter(|x| !x.is_empty()).collect()).unwrap_or_default();
    if !list.contains(&file) { list.push(file); }
    list.sort();
    std::fs::write(idx, format!("[{}]", list.iter().map(|x| format!("\"{x}\"")).collect::<Vec<_>>().join(","))).unwrap();
}
