//! STEP (ISO 10303-21) exact B-rep import: geometry builders and a `Body` with shared edges and vertices.
use crate::geom::*;
use crate::nurbs::*;
use crate::p21::*;
use crate::topo::*;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::Arc;

// ------------------------------------------------------------------ model
/// STEP model access; points are mapped by (p - off) * sc
pub struct M { pub e: HashMap<u32, Ent>, pub sc: f64, pub off: V }
impl M {
    pub fn simple(&self, id: u32) -> (&str, &Vec<Val>) { match &self.e[&id] { Ent::Simple(n, a) => (n, a), Ent::Complex(p) => (&p[0].0, &p[0].1) } }
    pub fn part(&self, id: u32, name: &str) -> Option<&Vec<Val>> { match &self.e[&id] { Ent::Simple(n, a) => (n == name).then_some(a), Ent::Complex(p) => p.iter().find(|x| x.0 == name).map(|x| &x.1) } }
    pub fn kind(&self, id: u32) -> String {
        match &self.e[&id] {
            Ent::Simple(n, _) => n.clone(),
            Ent::Complex(p) => p.iter().map(|x| x.0.as_str()).find(|n| n.starts_with("B_SPLINE_SURFACE_WITH_KNOTS") || n.starts_with("B_SPLINE_CURVE_WITH_KNOTS")).map(|n| format!("RATIONAL_{n}")).unwrap_or_else(|| p[0].0.clone()),
        }
    }
    /// point in normalized model units (centred, scaled to unit diagonal)
    pub fn pt(&self, id: u32) -> V { let c = lst(&self.simple(id).1[1]); [0, 1, 2].map(|k| (num(&c[k]) - self.off[k]) * self.sc) }
    pub fn dir(&self, id: u32) -> V { let c = lst(&self.simple(id).1[1]); unit([num(&c[0]), num(&c[1]), num(&c[2])]) }
    pub fn vertex(&self, id: u32) -> V { self.pt(rf(&self.simple(id).1[1])) }
    /// AXIS2_PLACEMENT_3D -> (origin, x, y, z)
    pub fn ax2(&self, id: u32) -> (V, V, V, V) {
        let (_, a) = self.simple(id);
        let o = self.pt(rf(&a[1]));
        let z = if let Val::Ref(r) = a[2] { self.dir(r) } else { [0.0, 0.0, 1.0] };
        let xr = if let Val::Ref(r) = a[3] { self.dir(r) } else if z[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
        let x = unit(sub(xr, mulv(z, dot(xr, z))));
        (o, x, cross(z, x), z)
    }
}
pub fn expand(mults: &[Val], knots: &[Val]) -> Vec<f64> { mults.iter().zip(knots).flat_map(|(m, v)| std::iter::repeat(num(v)).take(num(m) as usize)).collect() }

pub fn curve(m: &M, id: u32) -> Option<Box<dyn Crv>> {
    let (_, a) = m.simple(id);
    match m.kind(id).as_str() {
        "LINE" => { let (_, va) = m.simple(rf(&a[2])); Some(Box::new(LineC { o: m.pt(rf(&a[1])), d: mulv(m.dir(rf(&va[1])), num(&va[2]) * m.sc) })) }
        "CIRCLE" => { let (c, x, y, _) = m.ax2(rf(&a[1])); let r = num(&a[2]) * m.sc; Some(Box::new(ConicC { c, x, y, a: r, b: r })) }
        "ELLIPSE" => { let (c, x, y, _) = m.ax2(rf(&a[1])); Some(Box::new(ConicC { c, x, y, a: num(&a[2]) * m.sc, b: num(&a[3]) * m.sc })) }
        "SURFACE_CURVE" | "SEAM_CURVE" => curve(m, rf(&a[1])),
        "B_SPLINE_CURVE_WITH_KNOTS" | "RATIONAL_B_SPLINE_CURVE_WITH_KNOTS" => {
            let (deg, pts, mults, knots) = match m.part(id, "B_SPLINE_CURVE") {
                Some(b) => { let kk = m.part(id, "B_SPLINE_CURVE_WITH_KNOTS")?; (num(&b[0]), lst(&b[1]), lst(&kk[0]), lst(&kk[1])) }
                None => (num(&a[1]), lst(&a[2]), lst(&a[6]), lst(&a[7])),
            };
            let w: Vec<f64> = m.part(id, "RATIONAL_B_SPLINE_CURVE").map(|r| lst(&r[0]).iter().map(num).collect()).unwrap_or(vec![1.0; pts.len()]);
            let cp = pts.iter().zip(&w).map(|(p, w)| { let q = m.pt(rf(p)); [q[0] * w, q[1] * w, q[2] * w, *w] }).collect();
            Some(Box::new(BspC { p: deg as usize, k: expand(mults, knots), cp }))
        }
        _ => None,
    }
}

/// face surface; `pts` (boundary samples) bound the infinite analytic surfaces
pub fn surface(m: &M, id: u32, pts: &[V]) -> Option<(Box<dyn Surf>, &'static str)> {
    const MARGIN: f64 = 0.05;
    let (_, a) = m.simple(id);
    let range = |o: V, z: V| pts.iter().map(|p| dot(sub(*p, o), z)).fold((f64::MAX, f64::MIN), |(l, h), t| (l.min(t), h.max(t)));
    let unitz = |r: f64, h: f64| Box::new(Cyl { r, h, x: [1.0, 0.0, 0.0], y: [0.0, 1.0, 0.0], z: [0.0, 0.0, 1.0] }) as Box<dyn Surf>;
    let placed = |s: Box<dyn Surf>, x: V, y: V, z: V, t: V| Box::new(Moved { s, r: frame(x, y, z), t }) as Box<dyn Surf>;
    Some(match m.kind(id).as_str() {
        "PLANE" => {
            let (o, x, y, _) = m.ax2(rf(&a[1]));
            let ((ux, vx), (uy, vy)) = (range(o, x), range(o, y));
            let c = addv(o, addv(mulv(x, (ux + vx) / 2.0), mulv(y, (uy + vy) / 2.0)));
            (Box::new(Plane { o: c, a: x, b: y, ext: 0.5 * (vx - ux).max(vy - uy) + MARGIN }), "plane")
        }
        "CYLINDRICAL_SURFACE" => {
            let (o, x, y, z) = m.ax2(rf(&a[1]));
            let (lo, hi) = range(o, z);
            (placed(unitz(num(&a[2]) * m.sc, (hi - lo) / 2.0 + MARGIN), x, y, z, addv(o, mulv(z, (lo + hi) / 2.0))), "cylinder")
        }
        "CONICAL_SURFACE" => {
            let (o, x, y, z) = m.ax2(rf(&a[1]));
            let (r, t) = (num(&a[2]) * m.sc, num(&a[3]).tan());
            if t.abs() < 1e-12 { return None; }
            let apex = sub(o, mulv(z, r / t));
            let (lo, hi) = range(apex, z);
            let (v0, v1) = if t > 0.0 { ((lo - MARGIN).max(1e-6), hi + MARGIN) } else { (lo - MARGIN, (hi + MARGIN).min(-1e-6)) };
            (placed(Box::new(Cone { t, v0, v1 }), x, y, z, apex), "cone")
        }
        "SPHERICAL_SURFACE" => { let (o, x, y, z) = m.ax2(rf(&a[1])); (placed(Box::new(Sphere { r: num(&a[2]) * m.sc }), x, y, z, o), "sphere") }
        "TOROIDAL_SURFACE" => { let (o, x, y, z) = m.ax2(rf(&a[1])); (placed(Box::new(Torus { big: num(&a[2]) * m.sc, r: num(&a[3]) * m.sc }), x, y, z, o), "torus") }
        "B_SPLINE_SURFACE_WITH_KNOTS" | "RATIONAL_B_SPLINE_SURFACE_WITH_KNOTS" => {
            let (pu, pv, rows, ku, kv) = match m.part(id, "B_SPLINE_SURFACE") {
                Some(b) => { let kk = m.part(id, "B_SPLINE_SURFACE_WITH_KNOTS")?; (num(&b[0]), num(&b[1]), lst(&b[2]), expand(lst(&kk[0]), lst(&kk[2])), expand(lst(&kk[1]), lst(&kk[3]))) }
                None => (num(&a[1]), num(&a[2]), lst(&a[3]), expand(lst(&a[8]), lst(&a[10])), expand(lst(&a[9]), lst(&a[11]))),
            };
            let w: Option<Vec<Vec<f64>>> = m.part(id, "RATIONAL_B_SPLINE_SURFACE").map(|r| lst(&r[0]).iter().map(|row| lst(row).iter().map(num).collect()).collect());
            let cp = rows.iter().enumerate().map(|(i, row)| lst(row).iter().enumerate().map(|(j, p)| { let q = m.pt(rf(p)); let ww = w.as_ref().map_or(1.0, |w| w[i][j]); [q[0] * ww, q[1] * ww, q[2] * ww, ww] }).collect()).collect();
            (Box::new(BspS { pu: pu as usize, pv: pv as usize, ku, kv, cp }), "bspline")
        }
        "SURFACE_OF_LINEAR_EXTRUSION" => {
            let c = curve(m, rf(&a[1]))?;
            let (_, va) = m.simple(rf(&a[2]));
            let d = mulv(m.dir(rf(&va[1])), num(&va[2]) * m.sc);
            let (o, l, du) = (c.ev(c.range().0).0, norm(d), unit(d));
            let (lo, hi) = range(o, du);
            // u: profile parameters of the boundary samples (projected along d), unless the profile is periodic
            let (u0, u1) = if c.periodic() { c.range() } else {
                let (tl, th) = pts.iter().map(|p| t_of(c.as_ref(), sub(*p, mulv(du, dot(sub(*p, o), du))))).fold((f64::MAX, f64::MIN), |(l, h), t| (l.min(t), h.max(t)));
                let (r0, r1) = c.range();
                let pad = 0.05 * (th - tl).max(1e-6);
                ((tl - pad).max(r0), (th + pad).min(r1))
            };
            (Box::new(ExtrS { c, d, u0, u1, v0: (lo - MARGIN) / l, v1: (hi + MARGIN) / l }), "extrusion")
        }
        "SURFACE_OF_REVOLUTION" => {
            let c = curve(m, rf(&a[1]))?;
            let (_, ax) = m.simple(rf(&a[2]));
            (Box::new(RevS { c, a: m.pt(rf(&ax[1])), z: m.dir(rf(&ax[2])) }), "revolution")
        }
        _ => return None,
    })
}


/// result of importing one STEP file
pub struct Import { pub body: Body, pub diag: f64, pub unsupported: HashMap<String, usize> }

/// Import all B-rep faces of a STEP file into one `Body`. Coordinates are normalized to a unit model diagonal
/// (`diag` gives the scale back to model units); edges and vertices shared by faces are shared in the body.
pub fn import(path: &str) -> std::io::Result<Import> {
    let raw = std::fs::read(path)?;
    let e = parse(&raw);
    let txt = String::from_utf8_lossy(&raw);
    let model_tol = txt.find("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(").and_then(|i| txt[i + 45..].split(')').next()?.trim().parse::<f64>().ok()).unwrap_or(1e-6);
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for ent in e.values() {
        if let Ent::Simple(n, a) = ent { if n == "VERTEX_POINT" { if let Some(Ent::Simple(_, pa)) = e.get(&rf(&a[1])) { let c = lst(&pa[1]); for k in 0..3 { lo[k] = lo[k].min(num(&c[k])); hi[k] = hi[k].max(num(&c[k])); } } } }
    }
    let diag = norm(sub(hi, lo)).max(1e-9);
    let m = M { e, sc: 1.0 / diag, off: mulv(addv(lo, hi), 0.5) };
    let mut b = Body { tol: model_tol / diag, ..Default::default() };
    let mut unsupported: HashMap<String, usize> = HashMap::new();
    let is_t = |v: &Val| matches!(v, Val::Enum(s) if s == "T");
    // vertices and edges (shared by every face that uses them)
    let mut vmap: HashMap<u32, Id<Vertex>> = HashMap::new();
    let mut emap: HashMap<u32, Id<Edge>> = HashMap::new();
    let mut ids: Vec<u32> = m.e.keys().copied().collect();
    ids.sort_unstable();
    for &id in &ids {
        if !matches!(&m.e[&id], Ent::Simple(n, _) if n == "EDGE_CURVE") { continue; }
        let (_, a) = m.simple(id);
        let mut vid = |r: u32, b: &mut Body| *vmap.entry(r).or_insert_with(|| b.vertices.push(Vertex { p: m.vertex(r), tol: b.tol }));
        let (v0, v1) = (vid(rf(&a[1]), &mut b), vid(rf(&a[2]), &mut b));
        let crv = curve(&m, rf(&a[3]));
        if crv.is_none() { *unsupported.entry(format!("curve {}", m.kind(rf(&a[3])))).or_default() += 1; }
        // parameters at the edge's start and end vertex, in edge direction (v0 -> v1)
        let t = crv.as_ref().map_or([0.0, 0.0], |c| {
            let (p0, p1) = (b.vertices.get(v0).p, b.vertices.get(v1).p);
            let (mut t0, mut t1) = (t_of(c.as_ref(), p0), t_of(c.as_ref(), p1));
            if c.periodic() {
                // the edge runs along the curve when same_sense, else against it; full circles when v0 == v1
                if is_t(&a[4]) { while t1 <= t0 + 1e-12 { t1 += 2.0 * PI; } } else { while t0 <= t1 + 1e-12 { t0 += 2.0 * PI; } }
            }
            [t0, t1]
        });
        let eid = b.edges.push(Edge { curve: crv.map(Arc::from), v: [v0, v1], t, tol: b.tol, origin: Origin { op: 0, src: id } });
        emap.insert(id, eid);
    }
    // faces with loops of coedges
    let mut fmap: HashMap<u32, Id<Face>> = HashMap::new();
    for &id in &ids {
        if !matches!(&m.e[&id], Ent::Simple(n, _) if n == "ADVANCED_FACE") { continue; }
        let (_, a) = m.simple(id);
        let mut loops = vec![];
        let mut samples = vec![];
        for bound in lst(&a[1]) {
            let (_, ba) = m.simple(rf(bound));
            let (_, la) = m.simple(rf(&ba[1]));
            let Val::List(oes) = &la[1] else { continue }; // VERTEX_LOOP: point bound (apex), no coedges
            let mut ces: Vec<(Id<Edge>, bool)> = oes.iter().map(|oe| { let (_, oa) = m.simple(rf(oe)); (emap[&rf(&oa[3])], is_t(&oa[4])) }).collect();
            if !is_t(&ba[2]) { ces.reverse(); for c in ces.iter_mut() { c.1 = !c.1; } }
            for &(e, _) in &ces {
                let ed = b.edges.get(e);
                samples.push(b.vertices.get(ed.v[0]).p);
                if let Some(c) = &ed.curve { samples.extend((1..8).map(|k| c.ev(ed.t[0] + (ed.t[1] - ed.t[0]) * k as f64 / 8.0).0)); }
            }
            let coedges = ces.into_iter().map(|(edge, sense)| b.coedges.push(Coedge { edge, sense })).collect();
            loops.push(b.loops.push(Loop { coedges }));
        }
        let sid = rf(&a[2]);
        let Some((s, _)) = surface(&m, sid, &samples) else { *unsupported.entry(format!("surface {}", m.kind(sid))).or_default() += 1; continue };
        fmap.insert(id, b.faces.push(Face { surface: Arc::from(s), sense: is_t(&a[3]), loops, tol: b.tol, origin: Origin { op: 0, src: id } }));
    }
    // shells
    for &id in &ids {
        let (n, a) = m.simple(id);
        let closed = match n { "CLOSED_SHELL" => true, "OPEN_SHELL" => false, _ => continue };
        let faces = lst(&a[1]).iter().filter_map(|f| fmap.get(&rf(f)).copied()).collect();
        b.shells.push(Shell { faces, closed });
    }
    Ok(Import { body: b, diag, unsupported })
}
