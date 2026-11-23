//! STEP exact B-rep edge check (test harness, not a decoder).
//!
//! Ground truth: the edges CATIA wrote. Every edge between two faces is classified as
//! - Coincident:   one surface split into two faces (no intersection exists)
//! - PatchJoin:    shared parameter boundary of both patches
//! - BoundaryOn:   one patch's boundary laid onto the other surface within the model tolerance
//! - Reproduced:   the exact intersection A ∩ B lies within the model tolerance of the stored edge
//! - Tolerant:     surfaces tangent within tolerance; the tolerant contact reproduces the edge
//! - Band:         tangent within tolerance across a band: edge valid, position set by topology
//! - GapWithinTol: transversal faces that never touch, gap below the model tolerance (tolerant G0 join)
//! - Miss:         none of the above (a kernel or data problem)
//!
//! Local mode (default, fast) projects each edge sample onto A ∩ B; `--full` runs a complete SSI of the
//! two surfaces (tests seeding and tracing too).
//! Usage: spike step LIST [--full] [--pair SUBSTR] [--file SUBSTR] [--debug]
use crate::geom::*;
use crate::ms;
use crate::nurbs::*;
use crate::p21::*;
use crate::ssi::{gauss_newton, intersect, intersect_tol};
use std::collections::{BTreeMap, HashMap};
use std::f64::consts::PI;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

/// CATIA's G1 angular tolerance is about 0.5 degrees
const ANG_TOL: f64 = 1e-2;
const WORKERS: usize = 4;

use brep::step::{curve, surface, M};

/// an edge between two distinct faces, sampled from its stored curve
struct Edge { fa: usize, fb: usize, samples: Vec<V>, cat: String }
struct Model { name: String, faces: usize, surfs: Vec<Option<(Box<dyn Surf>, &'static str)>>, edges: Vec<Edge>, diag: f64, tol: f64 }

fn load(path: &str, pair: Option<&str>) -> Model {
    let raw = std::fs::read(path).unwrap();
    let e = parse(&raw);
    let txt = String::from_utf8_lossy(&raw);
    let model_tol = txt.find("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(").and_then(|i| txt[i + 45..].split(')').next()?.trim().parse::<f64>().ok()).unwrap_or(1e-6);
    // normalize by the vertex bounding box
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for ent in e.values() {
        if let Ent::Simple(n, a) = ent { if n == "VERTEX_POINT" { if let Some(Ent::Simple(_, pa)) = e.get(&rf(&a[1])) { let c = lst(&pa[1]); for k in 0..3 { lo[k] = lo[k].min(num(&c[k])); hi[k] = hi[k].max(num(&c[k])); } } } }
    }
    let diag = norm(sub(hi, lo)).max(1e-9);
    let m = M { e, sc: 1.0 / diag, off: mulv(addv(lo, hi), 0.5) };
    let mut faces: Vec<u32> = m.e.iter().filter(|(_, v)| matches!(v, Ent::Simple(n, _) if n == "ADVANCED_FACE")).map(|(k, _)| *k).collect();
    faces.sort_unstable();
    // face -> edge curves (via bounds, loops, oriented edges)
    let mut edge_faces: HashMap<u32, Vec<usize>> = HashMap::new();
    let face_edges: Vec<Vec<u32>> = faces.iter().enumerate().map(|(fi, &fid)| {
        let mut es = vec![];
        for b in lst(&m.simple(fid).1[1]) {
            let (_, la) = m.simple(rf(&m.simple(rf(b)).1[1]));
            if let Val::List(oes) = &la[1] { for oe in oes { let ec = rf(&m.simple(rf(oe)).1[3]); es.push(ec); edge_faces.entry(ec).or_default().push(fi); } }
        }
        es
    }).collect();
    // 7 interior samples of each stored edge curve, in edge direction
    let samples: HashMap<u32, Vec<V>> = edge_faces.keys().filter_map(|&ec| {
        let (_, a) = m.simple(ec);
        let (p1, p2) = (m.vertex(rf(&a[1])), m.vertex(rf(&a[2])));
        let c = curve(&m, rf(&a[3]))?;
        let same = matches!(&a[4], Val::Enum(s) if s == "T");
        let (ta, mut tb) = if same { (t_of(c.as_ref(), p1), t_of(c.as_ref(), p2)) } else { (t_of(c.as_ref(), p2), t_of(c.as_ref(), p1)) };
        if c.periodic() { while tb <= ta + 1e-12 { tb += 2.0 * PI; } }
        Some((ec, (1..8).map(|i| c.ev(ta + (tb - ta) * i as f64 / 8.0).0).collect()))
    }).collect();
    let surfs: Vec<_> = faces.iter().enumerate().map(|(fi, &fid)| {
        let mut pts: Vec<V> = face_edges[fi].iter().filter_map(|e| samples.get(e)).flatten().cloned().collect();
        pts.extend(face_edges[fi].iter().map(|e| m.vertex(rf(&m.simple(*e).1[1]))));
        surface(&m, rf(&m.simple(fid).1[2]), &pts)
    }).collect();
    let mut ids: Vec<&u32> = edge_faces.keys().collect();
    ids.sort_unstable();
    let edges = ids.into_iter().filter_map(|ec| {
        let mut d = edge_faces[ec].clone();
        d.sort_unstable(); d.dedup();
        if d.len() != 2 { return None; } // seam or non-manifold edge
        let (s, Some(a), Some(b)) = (samples.get(ec)?, &surfs[d[0]], &surfs[d[1]]) else { return None };
        let mut k = [a.1, b.1]; k.sort_unstable();
        let cat = format!("{}-{}", k[0], k[1]);
        if pair.is_some_and(|p| !cat.contains(p)) { return None; }
        Some(Edge { fa: d[0], fb: d[1], samples: s.clone(), cat })
    }).collect();
    Model { name: path.rsplit('/').next().unwrap().to_string(), faces: faces.len(), surfs, edges, diag, tol: model_tol / diag }
}

// ------------------------------------------------------------------ checks
#[derive(Clone, Copy, PartialEq, Debug)]
enum Verdict { Reproduced, Tolerant, Band, Coincident, PatchJoin, BoundaryOn, GapWithinTol, Miss }
impl Verdict {
    fn ok(self) -> bool { !matches!(self, Verdict::Miss) }
    fn label(self) -> &'static str {
        match self { Verdict::Reproduced => "", Verdict::Tolerant | Verdict::Band => " [tangent within tol]", Verdict::Coincident => " [coincident: face split]", Verdict::PatchJoin => " [patch join]", Verdict::BoundaryOn => " [patch boundary on surface]", Verdict::GapWithinTol => " [transversal, gap within tol]", Verdict::Miss => "" }
    }
}
struct Check { v: Verdict, on: f64, sin: f64, dist: f64, gap: f64 }

/// local: project each sample onto A ∩ B in the plane through it orthogonal to the stored edge.
/// Returns (max distance sample -> projection, max gap A-B at the projections).
fn project(sa: &dyn Surf, sb: &dyn Surf, s: &[V]) -> (f64, f64) {
    let (mut ga, mut gb) = (sa.uv_of(s[0]), sb.uv_of(s[0]));
    let (mut worst, mut worst_gap) = (0.0f64, 0.0f64);
    for i in 0..s.len() {
        let p = s[i];
        let t = unit(sub(s[(i + 1).min(s.len() - 1)], s[i.saturating_sub(1)]));
        let (ua, ub) = (uv_refine(sa, p, ga), uv_refine(sb, p, gb));
        let res = |x: &[f64; 4], f: &mut [f64]| { let (pa, pb) = (sa.ev(x[0], x[1]).p, sb.ev(x[2], x[3]).p); f[..3].copy_from_slice(&sub(pa, pb)); f[3] = dot(t, sub(pa, p)); };
        let x = gauss_newton([ua[0], ua[1], ub[0], ub[1]], 4, &res, 60);
        let (pa, pb) = (sa.ev(x[0], x[1]).p, sb.ev(x[2], x[3]).p);
        worst = worst.max(norm(sub(pa, p)));
        worst_gap = worst_gap.max(norm(sub(pa, pb)));
        (ga, gb) = ([x[0], x[1]], [x[2], x[3]]);
    }
    (worst, worst_gap)
}

/// full: SSI of the two surfaces, then each sample projected onto the nearest traced curve.
/// Returns (max distance, max gap).
fn trace(sa: &dyn Surf, sb: &dyn Surf, s: &[V], tolerant: Option<(f64, f64)>) -> (f64, f64) {
    let cv = match tolerant { Some((tol, ang)) => intersect_tol(sa, sb, tol, ang).0, None => intersect(sa, sb).0 };
    let pts: Vec<V> = cv.iter().flat_map(|c| c.pts.iter().copied()).filter(|q| q.iter().all(|v| v.is_finite())).collect();
    s.iter().fold((0.0f64, 0.0f64), |(wd, wg), p| {
        let Some(q0) = pts.iter().copied().min_by(|a, b| norm(sub(*a, *p)).total_cmp(&norm(sub(*b, *p)))) else { return (f64::MAX, f64::MAX) };
        let (d, g) = project(sa, sb, &[q0, *p]); // project from the traced point toward the sample
        (wd.max(d.max(norm(sub(q0, *p)).min(d + norm(sub(q0, *p))))), wg.max(g))
    })
}

fn check(m: &Model, e: &Edge, full: bool) -> Check {
    let (sa, sb) = (m.surfs[e.fa].as_ref().unwrap().0.as_ref(), m.surfs[e.fb].as_ref().unwrap().0.as_ref());
    let s = &e.samples;
    let on = s.iter().map(|p| dist(sa, *p).max(dist(sb, *p))).fold(0.0, f64::max);
    let sin = s.iter().map(|p| {
        let (ua, ub) = (sa.uv_of(*p), sb.uv_of(*p));
        let (ea, eb) = (sa.ev(ua[0], ua[1]), sb.ev(ub[0], ub[1]));
        let (na, nb) = (cross(ea.u, ea.v), cross(eb.u, eb.v));
        norm(cross(na, nb)) / (norm(na) * norm(nb)).max(1e-300)
    }).fold(f64::MAX, f64::min);
    let mk = |v, dist, gap| Check { v, on, sin, dist, gap };
    let da = sa.dom();
    let coincident = (0..25).all(|k| dist(sb, sa.ev(da[0] + (da[1] - da[0]) * (0.1 + 0.2 * (k / 5) as f64), da[2] + (da[3] - da[2]) * (0.1 + 0.2 * (k % 5) as f64)).p) < 1e-9);
    if coincident { return mk(Verdict::Coincident, 0.0, 0.0); }
    let on_bound = |srf: &dyn Surf| {
        let d = srf.dom();
        let (eu, ev) = (1e-6 * (d[1] - d[0]), 1e-6 * (d[3] - d[2]));
        s.iter().all(|p| { let uv = srf.uv_of(*p); (!srf.per_u() && ((uv[0] - d[0]).abs() < eu || (uv[0] - d[1]).abs() < eu)) || (!srf.per_v() && ((uv[1] - d[2]).abs() < ev || (uv[1] - d[3]).abs() < ev)) })
    };
    let (ba, bb) = (on_bound(sa), on_bound(sb));
    if ba && bb { return mk(Verdict::PatchJoin, 0.0, 0.0); }
    let tolerant = sin < ANG_TOL;
    let (dist, gap) = if full { trace(sa, sb, s, tolerant.then_some((m.tol, ANG_TOL))) } else { project(sa, sb, s) };
    let v = match () {
        // the exact intersection lies within the model tolerance of the stored edge
        _ if !tolerant && gap <= 1e-10 && dist <= m.tol => Verdict::Reproduced,
        _ if tolerant && gap <= m.tol && dist <= m.tol => Verdict::Tolerant,
        _ if tolerant && on <= m.tol => Verdict::Band,
        // one patch's boundary laid onto the other surface: they do not cross inside the patch domain
        _ if (ba || bb) && on <= m.tol => Verdict::BoundaryOn,
        // transversal faces that meet the edge but never touch: a tolerant G0 join (gap below the model tolerance)
        _ if gap <= m.tol && dist <= m.tol && on <= m.tol => Verdict::GapWithinTol,
        _ => Verdict::Miss,
    };
    mk(v, dist, gap)
}

/// run `f` over 0..n on WORKERS threads, results in index order
fn par_map<T: Send>(n: usize, f: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let out: Vec<Mutex<Option<T>>> = (0..n).map(|_| Mutex::new(None)).collect();
    let next = AtomicUsize::new(0);
    std::thread::scope(|sc| for _ in 0..WORKERS { sc.spawn(|| loop {
        let i = next.fetch_add(1, Ordering::Relaxed);
        if i >= n { break; }
        *out[i].lock().unwrap() = Some(f(i));
    }); });
    out.into_iter().map(|m| m.into_inner().unwrap().unwrap()).collect()
}

// ------------------------------------------------------------------ report
pub fn run(args: &[String]) -> usize {
    let flag = |k: &str| args.iter().any(|a| a == k);
    let opt = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).map(|s| s.as_str());
    let (full, debug, pair, file) = (flag("--full"), flag("--debug"), opt("--pair"), opt("--file"));
    println!("== STEP exact B-rep edge check ({} mode) ==", if full { "full SSI" } else { "local" });
    let files: Vec<String> = std::fs::read_to_string(&args[0]).unwrap().lines().filter(|l| !l.is_empty()).map(|l| l.rsplit('\t').next().unwrap().to_string()).filter(|f| file.map_or(true, |s| f.contains(s))).collect();
    let mut table: BTreeMap<String, [usize; 3]> = BTreeMap::new(); // edges, ok, stored within 1e-3 mm
    let (mut n_edges, mut n_miss) = (0, 0);
    let mut misses = vec![];
    let t_all = Instant::now();
    for path in &files {
        let t0 = Instant::now();
        let m = load(path, pair);
        let checks = par_map(m.edges.len(), |i| check(&m, &m.edges[i], full));
        for (e, c) in m.edges.iter().zip(&checks) {
            let row = table.entry(format!("{}{}", e.cat, c.v.label())).or_default();
            row[0] += 1;
            row[1] += c.v.ok() as usize;
            row[2] += (c.on * m.diag <= 1e-3) as usize;
            if !c.v.ok() {
                n_miss += 1;
                if debug || misses.len() < 20 { misses.push(format!("{} {}: on-surface {:.1e} mm, sin {:.1e}, distance {:.1e} mm, gap {:.1e} mm (model tol {:.1e} mm)", m.name, e.cat, c.on * m.diag, c.sin, c.dist * m.diag, c.gap * m.diag, m.tol * m.diag)); }
            }
        }
        n_edges += m.edges.len();
        println!("  {:<26} faces {:5} edges {:6} ok {:6} | {:6.0} ms", m.name, m.faces, m.edges.len(), checks.iter().filter(|c| c.v.ok()).count(), ms(t0));
    }
    println!("  {:<48} {:>7} {:>7} {:>9}", "category", "edges", "ok", "<=1e-3mm");
    for (k, v) in &table { println!("  {k:<48} {:7} {:7} {:9}", v[0], v[1], v[2]); }
    println!("  TOTAL edges {n_edges} | MISS {n_miss} | {:.1} s", ms(t_all) / 1e3);
    for s in &misses { println!("    MISS {s}"); }
    n_miss
}
