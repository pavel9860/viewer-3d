//! STEP exact B-rep edge check (test harness, not a decoder).
//! Ground truth: the edges CATIA wrote. For every edge between two faces:
//!   (a) do the stored edge samples lie on both adjacent surfaces?
//!   (b) local mode (default, fast): does each sample project onto the intersection A ∩ B (exact, or the
//!       tolerant contact for surfaces tangent within the model tolerance)?
//!   (c) --full: does a complete SSI of the two surfaces reproduce the edge (tests seeding/tracing too)?
//! Usage: spike step LIST [--full] [--pair SUBSTR] [--file SUBSTR] [--debug]
use crate::geom::*;
use crate::ms;
use crate::nurbs::*;
use crate::p21::*;
use crate::ssi::{gauss_newton, intersect, intersect_tol};
use std::collections::{BTreeMap, HashMap};
use std::f64::consts::PI;
use std::time::Instant;

struct M { e: HashMap<u32, Ent>, sc: f64, off: V }
impl M {
    fn simple(&self, id: u32) -> (&str, &Vec<Val>) { match &self.e[&id] { Ent::Simple(n, a) => (n, a), Ent::Complex(p) => (&p[0].0, &p[0].1) } }
    fn part(&self, id: u32, name: &str) -> Option<&Vec<Val>> { match &self.e[&id] { Ent::Simple(n, a) => if n == name { Some(a) } else { None }, Ent::Complex(p) => p.iter().find(|x| x.0 == name).map(|x| &x.1) } }
    fn kind(&self, id: u32) -> String {
        match &self.e[&id] {
            Ent::Simple(n, _) => n.clone(),
            Ent::Complex(p) => p.iter().map(|x| x.0.as_str()).find(|n| n.starts_with("B_SPLINE_SURFACE_WITH_KNOTS") || n.starts_with("B_SPLINE_CURVE_WITH_KNOTS")).map(|n| format!("RATIONAL_{n}")).unwrap_or_else(|| p[0].0.clone()),
        }
    }
    fn pt(&self, id: u32) -> V { let (_, a) = self.simple(id); let c = lst(&a[1]); [0, 1, 2].map(|k| (num(&c[k]) - self.off[k]) * self.sc) }
    fn dir(&self, id: u32) -> V { let (_, a) = self.simple(id); let c = lst(&a[1]); unit([num(&c[0]), num(&c[1]), num(&c[2])]) }
    /// AXIS2_PLACEMENT_3D -> (origin, x, y, z)
    fn ax2(&self, id: u32) -> (V, V, V, V) {
        let (_, a) = self.simple(id);
        let o = self.pt(rf(&a[1]));
        let z = if let Val::Ref(r) = a[2] { self.dir(r) } else { [0.0, 0.0, 1.0] };
        let xr = if let Val::Ref(r) = a[3] { self.dir(r) } else if z[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let x = unit(sub(xr, mulv(z, dot(xr, z))));
        (o, x, cross(z, x), z)
    }
}
fn expand(mults: &[Val], knots: &[Val]) -> Vec<f64> { let mut k = vec![]; for (m, v) in mults.iter().zip(knots) { for _ in 0..num(m) as usize { k.push(num(v)); } } k }

fn curve(m: &M, id: u32) -> Option<Box<dyn Crv>> {
    let k = m.kind(id);
    let (_, a) = m.simple(id);
    match k.as_str() {
        "LINE" => { let o = m.pt(rf(&a[1])); let (_, va) = m.simple(rf(&a[2])); let d = m.dir(rf(&va[1])); Some(Box::new(LineC { o, d: mulv(d, num(&va[2]) * m.sc) })) }
        "CIRCLE" => { let (c, x, y, _) = m.ax2(rf(&a[1])); let r = num(&a[2]) * m.sc; Some(Box::new(ConicC { c, x, y, a: r, b: r })) }
        "ELLIPSE" => { let (c, x, y, _) = m.ax2(rf(&a[1])); Some(Box::new(ConicC { c, x, y, a: num(&a[2]) * m.sc, b: num(&a[3]) * m.sc })) }
        "SURFACE_CURVE" | "SEAM_CURVE" => curve(m, rf(&a[1])),
        "B_SPLINE_CURVE_WITH_KNOTS" | "RATIONAL_B_SPLINE_CURVE_WITH_KNOTS" => {
            let (deg, pts, mults, knots) = if let Some(b) = m.part(id, "B_SPLINE_CURVE") { let kk = m.part(id, "B_SPLINE_CURVE_WITH_KNOTS").unwrap(); (num(&b[0]), lst(&b[1]).clone(), lst(&kk[0]).clone(), lst(&kk[1]).clone()) } else { (num(&a[1]), lst(&a[2]).clone(), lst(&a[6]).clone(), lst(&a[7]).clone()) };
            let w: Vec<f64> = m.part(id, "RATIONAL_B_SPLINE_CURVE").map(|r| lst(&r[0]).iter().map(num).collect()).unwrap_or(vec![1.0; pts.len()]);
            let cp = pts.iter().zip(&w).map(|(p, w)| { let q = m.pt(rf(p)); [q[0] * w, q[1] * w, q[2] * w, *w] }).collect();
            Some(Box::new(BspC { p: deg as usize, k: expand(&mults, &knots), cp }))
        }
        _ => None,
    }
}

/// face surface; `pts` (boundary samples) bound the infinite analytic surfaces
fn surface(m: &M, id: u32, pts: &[V]) -> Option<(Box<dyn Surf>, &'static str)> {
    let k = m.kind(id);
    let (_, a) = m.simple(id);
    let margin = 0.05;
    let proj = |o: V, z: V| { let t: Vec<f64> = pts.iter().map(|p| dot(sub(*p, o), z)).collect(); (t.iter().cloned().fold(f64::MAX, f64::min), t.iter().cloned().fold(f64::MIN, f64::max)) };
    let unitz = |r: f64, h: f64| Box::new(Cyl { r, h, x: [1.0, 0.0, 0.0], y: [0.0, 1.0, 0.0], z: [0.0, 0.0, 1.0] }) as Box<dyn Surf>;
    match k.as_str() {
        "PLANE" => {
            let (o, x, y, _) = m.ax2(rf(&a[1]));
            let ((ux, vx), (uy, vy)) = (proj(o, x), proj(o, y));
            let c = addv(o, addv(mulv(x, (ux + vx) / 2.0), mulv(y, (uy + vy) / 2.0)));
            Some((Box::new(Plane { o: c, a: x, b: y, ext: 0.5 * (vx - ux).max(vy - uy) + margin }), "plane"))
        }
        "CYLINDRICAL_SURFACE" => {
            let (o, x, y, z) = m.ax2(rf(&a[1]));
            let (lo, hi) = proj(o, z);
            Some((Box::new(Moved { s: unitz(num(&a[2]) * m.sc, (hi - lo) / 2.0 + margin), r: frame(x, y, z), t: addv(o, mulv(z, (lo + hi) / 2.0)) }), "cylinder"))
        }
        "CONICAL_SURFACE" => {
            let (o, x, y, z) = m.ax2(rf(&a[1]));
            let (r, t) = (num(&a[2]) * m.sc, num(&a[3]).tan());
            if t.abs() < 1e-12 { return None; }
            let apex = sub(o, mulv(z, r / t));
            let (lo, hi) = proj(apex, z);
            let (v0, v1) = if t > 0.0 { ((lo - margin).max(1e-6), hi + margin) } else { (lo - margin, (hi + margin).min(-1e-6)) };
            Some((Box::new(Moved { s: Box::new(Cone { t, v0, v1 }), r: frame(x, y, z), t: apex }), "cone"))
        }
        "SPHERICAL_SURFACE" => { let (o, x, y, z) = m.ax2(rf(&a[1])); Some((Box::new(Moved { s: Box::new(Sphere { r: num(&a[2]) * m.sc }), r: frame(x, y, z), t: o }), "sphere")) }
        "TOROIDAL_SURFACE" => { let (o, x, y, z) = m.ax2(rf(&a[1])); Some((Box::new(Moved { s: Box::new(Torus { big: num(&a[2]) * m.sc, r: num(&a[3]) * m.sc }), r: frame(x, y, z), t: o }), "torus")) }
        "B_SPLINE_SURFACE_WITH_KNOTS" | "RATIONAL_B_SPLINE_SURFACE_WITH_KNOTS" => {
            let (pu, pv, rows, ku, kv) = if let Some(b) = m.part(id, "B_SPLINE_SURFACE") {
                let kk = m.part(id, "B_SPLINE_SURFACE_WITH_KNOTS").unwrap();
                (num(&b[0]), num(&b[1]), lst(&b[2]).clone(), expand(lst(&kk[0]), lst(&kk[2])), expand(lst(&kk[1]), lst(&kk[3])))
            } else { (num(&a[1]), num(&a[2]), lst(&a[3]).clone(), expand(lst(&a[8]), lst(&a[10])), expand(lst(&a[9]), lst(&a[11]))) };
            let w: Option<Vec<Vec<f64>>> = m.part(id, "RATIONAL_B_SPLINE_SURFACE").map(|r| lst(&r[0]).iter().map(|row| lst(row).iter().map(num).collect()).collect());
            let cp = rows.iter().enumerate().map(|(i, row)| lst(row).iter().enumerate().map(|(j, p)| { let q = m.pt(rf(p)); let ww = w.as_ref().map_or(1.0, |w| w[i][j]); [q[0] * ww, q[1] * ww, q[2] * ww, ww] }).collect()).collect();
            Some((Box::new(BspS { pu: pu as usize, pv: pv as usize, ku, kv, cp }), "bspline"))
        }
        "SURFACE_OF_LINEAR_EXTRUSION" => {
            let c = curve(m, rf(&a[1]))?;
            let (_, va) = m.simple(rf(&a[2]));
            let d = mulv(m.dir(rf(&va[1])), num(&va[2]) * m.sc);
            let o = c.ev(c.range().0).0;
            let (lo, hi) = proj(o, unit(d));
            let l = norm(d);
            // u range: the profile parameters of the face's boundary samples, unless the profile is periodic
            let (u0, u1) = if c.periodic() { c.range() } else {
                let ts: Vec<f64> = pts.iter().map(|p| t_of(c.as_ref(), sub(*p, mulv(unit(d), dot(sub(*p, o), unit(d)))))).collect();
                let (lo_t, hi_t) = (ts.iter().cloned().fold(f64::MAX, f64::min), ts.iter().cloned().fold(f64::MIN, f64::max));
                let (r0, r1) = c.range();
                let pad = 0.05 * (hi_t - lo_t).max(1e-6);
                ((lo_t - pad).max(r0), (hi_t + pad).min(r1))
            };
            Some((Box::new(ExtrS { c, d, u0, u1, v0: (lo - margin) / l, v1: (hi + margin) / l }), "extrusion"))
        }
        "SURFACE_OF_REVOLUTION" => {
            let c = curve(m, rf(&a[1]))?;
            let (_, ax) = m.simple(rf(&a[2]));
            Some((Box::new(RevS { c, a: m.pt(rf(&ax[1])), z: m.dir(rf(&ax[2])) }), "revolution"))
        }
        _ => None,
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Verdict { Reproduced, Tolerant, Band, Coincident, PatchJoin, Miss }
struct Opts { full: bool, pair: Option<String>, file: Option<String>, debug: bool }

/// local check: project each sample onto A ∩ B in the plane through the sample orthogonal to the stored edge.
/// Returns (max distance sample->intersection, max gap A-B at the projections).
fn local(sa: &dyn Surf, sb: &dyn Surf, s: &[V]) -> (f64, f64) {
    let (mut worst, mut worst_gap) = (0.0f64, 0.0f64);
    let (mut ga, mut gb) = (sa.uv_of(s[0]), sb.uv_of(s[0]));
    for i in 0..s.len() {
        let p = s[i];
        let t = unit(sub(s[(i + 1).min(s.len() - 1)], s[i.saturating_sub(1)]));
        let (ua, ub) = (uv_refine(sa, p, ga), uv_refine(sb, p, gb));
        let res = |x: &[f64; 4], f: &mut [f64]| { let (pa, pb) = (sa.ev(x[0], x[1]).p, sb.ev(x[2], x[3]).p); f[..3].copy_from_slice(&sub(pa, pb)); f[3] = dot(t, sub(pa, p)); };
        let x = gauss_newton([ua[0], ua[1], ub[0], ub[1]], 4, &res, 60);
        let (pa, pb) = (sa.ev(x[0], x[1]).p, sb.ev(x[2], x[3]).p);
        worst = worst.max(norm(sub(pa, p)));
        worst_gap = worst_gap.max(norm(sub(pa, pb)));
        ga = [x[0], x[1]];
        gb = [x[2], x[3]];
    }
    (worst, worst_gap)
}

/// full check: SSI of the two surfaces, then exact projection of each sample onto the traced curve
fn full(sa: &dyn Surf, sb: &dyn Surf, s: &[V], tolerant: Option<(f64, f64)>) -> f64 {
    let cv = match tolerant { Some((tol, ang)) => intersect_tol(sa, sb, tol, ang).0, None => intersect(sa, sb).0 };
    s.iter().map(|p| {
        let Some(&q0) = cv.iter().flat_map(|c| c.pts.iter()).filter(|q| q.iter().all(|v| v.is_finite())).min_by(|a, b| norm(sub(**a, *p)).total_cmp(&norm(sub(**b, *p)))) else { return f64::MAX };
        if norm(sub(q0, *p)) > 0.05 { return norm(sub(q0, *p)); }
        let (ua, ub) = (sa.uv_of(q0), sb.uv_of(q0));
        let (ea, eb) = (sa.ev(ua[0], ua[1]), sb.ev(ub[0], ub[1]));
        let tv = cross(cross(ea.u, ea.v), cross(eb.u, eb.v));
        let t = if norm(tv) > 1e-9 * norm(cross(ea.u, ea.v)) * norm(cross(eb.u, eb.v)) { unit(tv) } else { unit(sub(q0, *p)) };
        let res = |x: &[f64; 4], f: &mut [f64]| { let (pa, pb) = (sa.ev(x[0], x[1]).p, sb.ev(x[2], x[3]).p); f[..3].copy_from_slice(&sub(pa, pb)); f[3] = dot(t, sub(pa, *p)); };
        let x = gauss_newton([ua[0], ua[1], ub[0], ub[1]], 4, &res, 80);
        let (pa, pb) = (sa.ev(x[0], x[1]).p, sb.ev(x[2], x[3]).p);
        if norm(sub(pa, pb)) > tolerant.map_or(1e-9, |t| t.0) { return norm(sub(q0, *p)); }
        norm(sub(pa, *p))
    }).fold(0.0, f64::max)
}

pub fn run(args: &[String]) -> usize {
    let list = &args[0];
    let get = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).cloned();
    let o = Opts { full: args.iter().any(|a| a == "--full"), pair: get("--pair"), file: get("--file"), debug: args.iter().any(|a| a == "--debug") };
    println!("== STEP exact B-rep edge check ({} mode) ==", if o.full { "full SSI" } else { "local" });
    let files: Vec<String> = std::fs::read_to_string(list).unwrap().lines().filter(|l| !l.is_empty()).map(|l| l.split('\t').last().unwrap().to_string()).filter(|f| o.file.as_ref().map_or(true, |s| f.contains(s.as_str()))).collect();
    let mut table: BTreeMap<String, [usize; 6]> = BTreeMap::new(); // per category: edges, reproduced, tolerant, band, miss, on<=1e-3mm
    let mut totals = [0usize; 6];
    let mut misses: Vec<String> = vec![];
    let t_all = Instant::now();
    for f in &files {
        let t0 = Instant::now();
        let raw = std::fs::read(f).unwrap();
        let e = parse(&raw);
        let txt = String::from_utf8_lossy(&raw);
        let model_tol = txt.find("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(").and_then(|i| txt[i + 45..].split(')').next().and_then(|v| v.trim().parse::<f64>().ok())).unwrap_or(1e-6);
        let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
        for ent in e.values() { if let Ent::Simple(n, a) = ent { if n == "VERTEX_POINT" { if let Some(Ent::Simple(_, pa)) = e.get(&rf(&a[1])) { let c = lst(&pa[1]); for k in 0..3 { let x = num(&c[k]); lo[k] = lo[k].min(x); hi[k] = hi[k].max(x); } } } } }
        let diag = norm(sub(hi, lo)).max(1e-9);
        let (tol_n, ang) = (model_tol / diag, 1e-2); // CATIA's G1 angular tolerance is ~0.5 deg
        let m = M { e, sc: 1.0 / diag, off: mulv(addv(lo, hi), 0.5) };
        // faces and the edges they use
        let mut faces: Vec<u32> = m.e.iter().filter(|(_, v)| matches!(v, Ent::Simple(n, _) if n == "ADVANCED_FACE")).map(|(k, _)| *k).collect();
        faces.sort();
        let mut edge_faces: HashMap<u32, Vec<usize>> = HashMap::new();
        let mut face_edges: Vec<Vec<u32>> = vec![];
        for (fi, &fid) in faces.iter().enumerate() {
            let (_, a) = m.simple(fid);
            let mut es = vec![];
            for b in lst(&a[1]) { let (_, ba) = m.simple(rf(b)); let (_, la) = m.simple(rf(&ba[1])); if let Val::List(oes) = &la[1] { for oe in oes { let (_, oa) = m.simple(rf(oe)); let ec = rf(&oa[3]); es.push(ec); edge_faces.entry(ec).or_default().push(fi); } } }
            face_edges.push(es);
        }
        // edge samples from the stored curve between its vertices
        let mut samples: HashMap<u32, Vec<V>> = HashMap::new();
        for &ec in edge_faces.keys() {
            let (_, a) = m.simple(ec);
            let (p1, p2) = (m.pt(rf(&m.simple(rf(&a[1])).1[1])), m.pt(rf(&m.simple(rf(&a[2])).1[1])));
            let same = matches!(&a[4], Val::Enum(s) if s == "T");
            let Some(c) = curve(&m, rf(&a[3])) else { continue };
            let (ta, mut tb) = if same { (t_of(c.as_ref(), p1), t_of(c.as_ref(), p2)) } else { (t_of(c.as_ref(), p2), t_of(c.as_ref(), p1)) };
            if c.periodic() { while tb <= ta + 1e-12 { tb += 2.0 * PI; } }
            samples.insert(ec, (1..8).map(|i| c.ev(ta + (tb - ta) * i as f64 / 8.0).0).collect());
        }
        let surfs: Vec<Option<(Box<dyn Surf>, &'static str)>> = faces.iter().enumerate().map(|(fi, &fid)| {
            let mut pts: Vec<V> = face_edges[fi].iter().filter_map(|e| samples.get(e)).flatten().cloned().collect();
            for e in &face_edges[fi] { let (_, a) = m.simple(*e); pts.push(m.pt(rf(&m.simple(rf(&a[1])).1[1]))); }
            surface(&m, rf(&m.simple(fid).1[2]), &pts)
        }).collect();
        // edges between two distinct faces with supported surfaces
        let mut recs: Vec<(usize, usize, Vec<V>, String)> = vec![];
        let mut keys: Vec<&u32> = edge_faces.keys().collect();
        keys.sort();
        for ec in keys {
            let mut d = edge_faces[ec].clone(); d.sort(); d.dedup();
            if d.len() != 2 { continue; }
            let (Some(s), Some(a), Some(b)) = (samples.get(ec), &surfs[d[0]], &surfs[d[1]]) else { continue };
            let mut k = [a.1, b.1]; k.sort();
            let cat = format!("{}-{}", k[0], k[1]);
            if o.pair.as_ref().map_or(false, |p| !cat.contains(p.as_str())) { continue; }
            recs.push((d[0], d[1], s.clone(), cat));
        }
        // check edges on 6 workers
        let out: Vec<std::sync::Mutex<Option<(Verdict, f64, f64, f64)>>> = (0..recs.len()).map(|_| std::sync::Mutex::new(None)).collect();
        let ctr = std::sync::atomic::AtomicUsize::new(0);
        std::thread::scope(|sc| { for _ in 0..6 { sc.spawn(|| loop {
            let i = ctr.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            if i >= recs.len() { break; }
            let (fa, fb, s, _) = &recs[i];
            let (sa, sb) = (surfs[*fa].as_ref().unwrap().0.as_ref(), surfs[*fb].as_ref().unwrap().0.as_ref());
            let on = s.iter().map(|p| dist(sa, *p).max(dist(sb, *p))).fold(0.0, f64::max);
            let tang = s.iter().map(|p| { let (ua, ub) = (sa.uv_of(*p), sb.uv_of(*p)); let (ea, eb) = (sa.ev(ua[0], ua[1]), sb.ev(ub[0], ub[1])); let (na, nb) = (cross(ea.u, ea.v), cross(eb.u, eb.v)); norm(cross(na, nb)) / (norm(na) * norm(nb)).max(1e-300) }).fold(f64::MAX, f64::min);
            // one surface split into two faces
            let da = sa.dom();
            let coincident = (0..5).all(|ii| (0..5).all(|jj| dist(sb, sa.ev(da[0] + (da[1] - da[0]) * (0.1 + 0.2 * ii as f64), da[2] + (da[3] - da[2]) * (0.1 + 0.2 * jj as f64)).p) < 1e-9));
            // shared boundary of two patches
            let on_bound = |srf: &dyn Surf| { let d = srf.dom(); s.iter().all(|p| { let uv = srf.uv_of(*p); let (eu, ev) = (1e-6 * (d[1] - d[0]), 1e-6 * (d[3] - d[2])); (!srf.per_u() && ((uv[0] - d[0]).abs() < eu || (uv[0] - d[1]).abs() < eu)) || (!srf.per_v() && ((uv[1] - d[2]).abs() < ev || (uv[1] - d[3]).abs() < ev)) }) };
            let tolerant = tang < ang;
            let v = if coincident { (Verdict::Coincident, 0.0, 0.0) }
                else if on_bound(sa) && on_bound(sb) { (Verdict::PatchJoin, 0.0, 0.0) }
                else if o.full {
                    let cov = full(sa, sb, s, if tolerant { Some((tol_n, ang)) } else { None });
                    (if !tolerant && cov <= 2.0 * on + 1e-7 { Verdict::Reproduced } else if tolerant && cov <= tol_n { Verdict::Tolerant } else if tolerant && on <= tol_n { Verdict::Band } else { Verdict::Miss }, cov, 0.0)
                } else {
                    let (d, gap) = local(sa, sb, s);
                    let v = if !tolerant && gap <= 1e-10 && d <= 2.0 * on + 1e-7 { Verdict::Reproduced }
                        else if tolerant && gap <= tol_n && d <= tol_n { Verdict::Tolerant }
                        else if tolerant && on <= tol_n { Verdict::Band }
                        else { Verdict::Miss };
                    (v, d, gap)
                };
            *out[i].lock().unwrap() = Some((v.0, v.1, v.2, on));
            if o.debug && v.0 == Verdict::Miss { eprintln!("MISS {} edge#{i}: on-surface {:.1e} sin {:.1e} distance {:.1e} gap {:.1e} (model units: diag {diag:.1})", recs[i].3, on, tang, v.1, v.2); }
        }); } });
        let mut n_ok = 0;
        for (i, (_, _, _, cat)) in recs.iter().enumerate() {
            let (v, d, _, on) = out[i].lock().unwrap().unwrap();
            let label = match v { Verdict::Coincident => format!("{cat} [coincident: face split]"), Verdict::PatchJoin => format!("{cat} [patch join]"), Verdict::Band | Verdict::Tolerant => format!("{cat} [tangent within tolerance]"), _ => cat.clone() };
            let t = table.entry(label).or_default();
            t[0] += 1;
            if on * diag <= 1e-3 { t[5] += 1; }
            let slot = match v { Verdict::Reproduced => 1, Verdict::Tolerant => 2, Verdict::Band => 3, Verdict::Miss => 4, _ => 0 };
            if slot > 0 { t[slot] += 1; totals[slot] += 1; }
            if matches!(v, Verdict::Reproduced | Verdict::Tolerant) { n_ok += 1; }
            if v == Verdict::Miss && misses.len() < 30 { misses.push(format!("{}: {cat} edge#{i} on-surface {:.1e} mm, distance {:.1e} mm", f.rsplit('/').next().unwrap(), on * diag, d * diag)); }
            totals[0] += 1;
        }
        println!("  {:<26} faces {:5} edges {:6} ok {:6} | {:6.0} ms", f.rsplit('/').next().unwrap(), faces.len(), recs.len(), n_ok, ms(t0));
    }
    println!("  categories: edges | reproduced exact | tolerant reproduced | tolerant band (valid, position set by topology) | MISS | stored within 1e-3 mm");
    for (k, v) in &table { println!("    {k:<44} {:6} {:6} {:6} {:6} {:6} {:6}", v[0], v[1], v[2], v[3], v[4], v[5]); }
    println!("  TOTAL edges {} | reproduced {} | tolerant {} | band {} | MISS {} | {:.1} s", totals[0], totals[1], totals[2], totals[3], totals[4], ms(t_all) / 1e3);
    for s in &misses { println!("    MISS {s}"); }
    totals[4]
}
