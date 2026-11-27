//! Edges between a surface of the given IR kind and another face: how far the kind's coedge points lie off the
//! other face's surface (ends and middle), against the end gap between the two coedges.
//! Usage: cgm_edges IR KIND
use brep::cgm::*;
use brep::geom::*;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let m = Model::read(&std::fs::read(&a[1]).unwrap()).unwrap();
    let kind: u32 = a[2].parse().unwrap();
    let mut by_edge: std::collections::HashMap<usize, Vec<(usize, usize)>> = Default::default();
    for f in &m.faces { for l in &f.loops { for c in l.clone() { by_edge.entry(m.coedges[c].edge).or_default().push((c, f.surface)); } } }
    let (mut n, mut off, mut mid_off, mut gap) = (0, 0.0f64, 0.0f64, 0.0f64);
    let mut worst = vec![];
    for cs in by_edge.values().filter(|cs| cs.len() == 2) {
        let Some(&(cf, _)) = cs.iter().find(|x| m.surface_kinds[x.1] == kind) else { continue };
        let &(co, so) = cs.iter().find(|x| x.0 != cf).unwrap();
        if m.surface_kinds[so] == kind { continue; }
        let (a, b) = (&m.coedges[cf], &m.coedges[co]);
        let s = m.surfaces[so].as_ref();
        let on = |p: V| { let uv = s.uv_of(p); norm(sub(s.ev(uv[0], uv[1]).p, p)) };
        let e = [on(m.point(a, a.t[0])), on(m.point(a, a.t[1]))];
        let mo = on(m.point(a, (a.t[0] + a.t[1]) / 2.0));
        let (pa, pb) = ([m.point(a, a.t[0]), m.point(a, a.t[1])], [m.point(b, b.t[0]), m.point(b, b.t[1])]);
        let d = |x: V, y: V| norm(sub(x, y));
        let g = d(pa[0], pb[0]).max(d(pa[1], pb[1])).min(d(pa[0], pb[1]).max(d(pa[1], pb[0])));
        n += 1;
        off = off.max(e[0].max(e[1]));
        mid_off = mid_off.max(mo);
        gap = gap.max(g);
        worst.push((mo, m.surface_kinds[so], m.pcurve_kinds[b.pcurve], g));
    }
    worst.sort_by(|x, y| y.0.total_cmp(&x.0));
    println!("{n} edges: off the other surface ends {off:.1e} middle {mid_off:.1e}; end gap {gap:.1e}; worst (middle off, other kind, other p-curve kind, gap) {:?}", &worst[..worst.len().min(4)]);
}
