//! brep CLI.
//!   brep check FILE.stp|LIST.tsv          import + validate
//!   brep bool FILE.stp|LIST.tsv [SHELLS]  booleans of the largest closed shells with a cylinder through them
//!   brep mesh FILE.stp|LIST.tsv [CHORD]   watertight tessellation (chord relative to the model diagonal, default 1e-3)
use brep::boolean::{boolean, Op};
use brep::geom::*;
use brep::topo::*;
use brep::{check, make, mesh, step};
use std::time::Instant;

const WORKERS: usize = 4;
fn ms(t: Instant) -> f64 { t.elapsed().as_secs_f64() * 1e3 }

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 { eprintln!("usage: brep check|bool FILE.stp|LIST.tsv [SHELLS]"); std::process::exit(2); }
    let files: Vec<String> = if args[2].ends_with(".tsv") { std::fs::read_to_string(&args[2]).unwrap().lines().filter(|l| !l.is_empty()).map(|l| l.rsplit('\t').next().unwrap().to_string()).collect() } else { vec![args[2].clone()] };
    let bad = match args[1].as_str() {
        "check" => files.iter().map(|f| check_file(f)).sum::<usize>(),
        "bool" => { let n = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(2); files.iter().map(|f| bool_file(f, n)).sum() }
        "mesh" => { let c = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1e-3); files.iter().map(|f| mesh_file(f, c)).sum() }
        _ => { eprintln!("unknown command {}", args[1]); 2 }
    };
    std::process::exit((bad > 0) as i32);
}

fn check_file(f: &str) -> usize {
    let t0 = Instant::now();
    let im = match step::import(f) { Ok(i) => i, Err(e) => { println!("{f}: {e}"); return 1 } };
    let ti = ms(t0);
    let t1 = Instant::now();
    let r = check::check(&im.body, 7, WORKERS);
    let mm = |x: f64| x * im.diag;
    println!("{:<26} {} | import {ti:.0} ms, check {:.0} ms | violations {} (open loops {}, edge use {}, orientation {}, vertex off curve {}, edge off surface {}) | non-manifold {} | max dev vertex {:.1e} mm edge {:.1e} mm | no curve {} | unsupported {:?}",
        f.rsplit('/').next().unwrap(), im.body.stats(), ms(t1), r.violations(), r.loops_open, r.edge_use, r.orientation, r.vertex_off_curve, r.edge_off_surface, r.nonmanifold, mm(r.max_vertex_dev), mm(r.max_edge_dev), r.edges_without_curve, im.unsupported);
    r.violations()
}

/// A (op) B for the `n` largest closed shells A of the file and a z-cylinder B through A's box centre.
/// No ground truth: checks are the validator, arrangement health and the volume identities
/// V(A-B) + V(A^B) = V(A) and V(AuB) = V(A) + V(B) - V(A^B).
fn bool_file(f: &str, n: usize) -> usize {
    let im = match step::import(f) { Ok(i) => i, Err(e) => { println!("{f}: {e}"); return 1 } };
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
        let mut line = format!("{:<26} shell {si} ({} faces) V(A) {va:.4e}{}", f.rsplit('/').next().unwrap(), sh.faces.len(), if ba > 0 { format!(" [{ba} faces lost area]") } else { String::new() });
        for (k, op) in [Op::Difference, Op::Intersection, Op::Union].into_iter().enumerate() {
            let t0 = Instant::now();
            let (res, st) = boolean(&a, &cyl, op, WORKERS);
            let t = ms(t0);
            let rep = check::check(&res, 5, WORKERS);
            let (vr, br) = vol(&res);
            v[k] = vr;
            if std::env::var("BREP_DEBUG").is_ok() && rep.edge_use > 0 { unmatched(&res, &a, &sh.faces, b); }
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
    for (fid, f) in r.faces.iter() { for &l in &f.loops { for &c in &r.loops.get(l).coedges { let ce = r.coedges.get(c); uses.entry(ce.edge).or_default().push((fid.0, ce.sense)); } } }
    let step_of = |wf: u32| if (wf as usize) < a.faces.len() { format!("A#{}", src.faces.get(shell[wf as usize]).origin.src) } else { format!("B{}", wf as usize - a.faces.len()) };
    for (e, u) in uses { let f = u.iter().filter(|x| x.1).count(); if 2 * f != u.len() {
        let ed = r.edges.get(e);
        let fs: Vec<String> = u.iter().map(|(fi, s)| format!("{}{}", step_of(r.faces.get(Id::new(*fi as usize)).origin.src - 1), if *s { "+" } else { "-" })).collect();
        eprintln!("  unmatched edge {e:?} origin {:?} {:?} -> {:?} faces {fs:?}", ed.origin, r.vertices.get(ed.v[0]).p, r.vertices.get(ed.v[1]).p);
    } }
}

/// tessellate every shell; closed shells must give closed meshes. Volume vs the exact face integral is an indicator.
fn mesh_file(f: &str, chord: f64) -> usize {
    use brep::tess::{mesh_volume, open_edges, tessellate, Opts};
    let im = match step::import(f) { Ok(i) => i, Err(e) => { println!("{f}: {e}"); return 1 } };
    let b = &im.body;
    let t0 = Instant::now();
    let o = Opts { chord, ..Default::default() };
    let (mut tris, mut open, mut nonman, mut empty, mut coarse, mut dev, mut vrel) = (0, 0, 0, 0, 0, 0.0f64, 0.0f64);
    let (mut mv, mut ev) = (0.0, 0.0);
    for (_, sh) in b.shells.iter() {
        let (m, st) = tessellate(b, &sh.faces, &o, WORKERS);
        tris += m.tri.len();
        empty += st.empty;
        coarse += st.coarse;
        dev = dev.max(st.max_dev);
        if sh.closed {
            let (op, nm) = open_edges(&m);
            if op > 0 && std::env::var("BREP_DEBUG").is_ok() {
                let mut d: std::collections::HashMap<(u32, u32), Vec<u32>> = Default::default();
                for (i, t) in m.tri.iter().enumerate() { for k in 0..3 { d.entry((t[k], t[(k + 1) % 3])).or_default().push(m.face[i]); } }
                let mut fs: std::collections::BTreeMap<u32, usize> = Default::default();
                for (e, f) in &d { if !d.contains_key(&(e.1, e.0)) { for x in f { *fs.entry(b.faces.get(sh.faces[*x as usize]).origin.src).or_default() += 1; } } }
                eprintln!("  open edges by face (STEP id: count): {fs:?}");
            }
            open += op;
            nonman += nm;
            let (v, _) = mesh::volume(b, &sh.faces);
            mv += mesh_volume(&m);
            ev += v;
            vrel = vrel.max((mesh_volume(&m) - v).abs() / v.abs().max(1e-30));
        }
    }
    let t = ms(t0);
    let bad = open + empty + coarse;
    println!("{:<26} {} faces -> {tris} triangles in {t:.0} ms | open edges {open} non-manifold {nonman} | empty faces {empty} coarse {coarse} | max dev {:.1e} mm | volume mesh/exact {:.6} (worst shell rel {vrel:.1e}){}",
        f.rsplit('/').next().unwrap(), b.faces.len(), dev * im.diag, mv / ev, if bad > 0 { " FAIL" } else { "" });
    bad
}
