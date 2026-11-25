//! SSI tests: hand-made hard cases and a randomized stress test against analytic ground truth.
use crate::geom::*;
use crate::k::Rng;
use crate::ms;
use crate::ssi::{intersect, End, Kind};
use std::f64::consts::PI;
use std::time::Instant;

fn bump(z: impl Fn(usize, usize) -> f64) -> Bez {
    let mut cp = [[0.0; 3]; 16];
    for i in 0..4 { for j in 0..4 { cp[i * 4 + j] = [i as f64 / 3.0, j as f64 / 3.0, z(i, j)]; } }
    Bez { cp }
}
fn implicit_case(name: &str, p: V) -> f64 {
    let cz = (p[0] * p[0] + p[1] * p[1]).sqrt() - 1.0;
    match name.as_bytes()[0] {
        b'a' => cz.abs().max((p[2] - 0.3 * p[0]).abs()),
        b'b' => cz.abs().max(((p[1] * p[1] + p[2] * p[2]).sqrt() - 0.7).abs()),
        b'c' => cz.abs().max(((p[1] * p[1] + p[2] * p[2]).sqrt() - 1.0).abs()),
        b'd' => cz.abs().max((p[0] - 1.0).abs()),
        b'e' => (p[2] - 0.2).abs(),
        b'h' => 0.0,
        _ => (p[0] - 0.1).abs(),
    }
}

/// hand-made cases: transversal, branch points, tangential curves, near-tangency, gaps. Returns failures.
pub fn cases() -> usize {
    println!("== SSI cases ==");
    let (x, y, z) = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
    let bz = |i: usize, j: usize| if (1..=2).contains(&i) && (1..=2).contains(&j) { 0.6 } else { 0.0 };
    // (name, expected (curves, branch points), a, b)
    let list: Vec<(&str, (usize, usize), Box<dyn Surf>, Box<dyn Surf>)> = vec![
        ("a plane tilted x cylinder", (1, 0), Box::new(Plane { o: [0.0; 3], a: [1.0, 0.0, 0.3], b: y, ext: 2.0 }), Box::new(Cyl { r: 1.0, h: 2.0, x, y, z })),
        ("b cyl r1 x cyl r0.7 perpendicular", (2, 0), Box::new(Cyl { r: 1.0, h: 2.0, x, y, z }), Box::new(Cyl { r: 0.7, h: 2.0, x: y, y: z, z: x })),
        ("c Steinmetz: 4 arcs, 2 branch points", (4, 2), Box::new(Cyl { r: 1.0, h: 2.0, x, y, z }), Box::new(Cyl { r: 1.0, h: 2.0, x: y, y: z, z: x })),
        ("d plane tangent to cylinder: 1 line", (1, 0), Box::new(Plane { o: [1.0, 0.0, 0.0], a: y, b: z, ext: 2.0 }), Box::new(Cyl { r: 1.0, h: 2.0, x, y, z })),
        ("e bicubic bump x plane z=0.2", (1, 0), Box::new(bump(bz)), Box::new(Plane { o: [0.0, 0.0, 0.2], a: x, b: y, ext: 2.0 })),
        ("f bicubic x near-tangent copy (1e-4 tilt)", (1, 0), Box::new(bump(bz)), Box::new(bump(move |i, j| bz(i, j) + 1e-4 * (i as f64 / 3.0) - 1e-5))),
        ("h plane 1e-9 off tangent: no contact", (0, 0), Box::new(Plane { o: [1.0 + 1e-9, 0.0, 0.0], a: y, b: z, ext: 2.0 }), Box::new(Cyl { r: 1.0, h: 2.0, x, y, z })),
    ];
    let mut fails = 0;
    for (name, exp, a, b) in &list {
        let t0 = Instant::now();
        let (cv, nb) = intersect(a.as_ref(), b.as_ref());
        let res = cv.iter().flat_map(|c| c.pts.iter()).map(|p| implicit_case(name, *p)).fold(0.0f64, f64::max);
        let ok = (cv.len(), nb) == *exp && res < 1e-9;
        if !ok { fails += 1; }
        println!("  {name:<42} curves {} branch {} (expect {:?}) residual {res:.1e} {:.1} ms {}", cv.len(), nb, exp, ms(t0), if ok { "PASS" } else { "FAIL" });
    }
    fails
}

fn rand_rot(r: &mut Rng) -> [V; 3] {
    let mut q = [r.f(), r.f(), r.f(), r.f()];
    let n = (q.iter().map(|x| x * x).sum::<f64>()).sqrt();
    for x in q.iter_mut() { *x /= n; }
    let [w, x, y, z] = q;
    [[1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - w * z), 2.0 * (x * z + w * y)],
     [2.0 * (x * y + w * z), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - w * x)],
     [2.0 * (x * z - w * y), 2.0 * (y * z + w * x), 1.0 - 2.0 * (x * x + y * y)]]
}
fn mv(s: Box<dyn Surf>, r: [V; 3], t: V) -> Box<dyn Surf> { Box::new(Moved { s, r, t }) }
fn rot_to(z: V) -> [V; 3] {
    let z = unit(z);
    let h = if z[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    let e1 = unit(cross(h, z));
    frame(e1, cross(z, e1), z)
}
/// random rotation about z (moves the seam of an axially symmetric surface)
fn about_z(r: &mut Rng) -> [V; 3] { let a = PI * r.f(); let (s, c) = a.sin_cos(); [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]] }
fn numeric_len(f: impl Fn(f64) -> V, a: f64, b: f64) -> f64 {
    let n = 20000;
    let (mut l, mut p) = (0.0, f(a));
    for i in 1..=n { let q = f(a + (b - a) * i as f64 / n as f64); l += norm(sub(q, p)); p = q; }
    l
}
#[derive(Default, Debug, Clone, Copy, PartialEq)]
struct Expect { closed: usize, open: usize, tangential: usize, points: usize, branch_arcs: usize }

/// 10 analytic families x `trials`, random parameters (every 5th at the critical value), random rigid motion.
pub fn stress(trials: usize) -> usize {
    println!("== SSI stress: analytic families, random rigid motion of both surfaces ==");
    let mut rng = Rng(2026);
    let unitz = |r: f64, h: f64| Box::new(Cyl { r, h, x: [1.0, 0.0, 0.0], y: [0.0, 1.0, 0.0], z: [0.0, 0.0, 1.0] }) as Box<dyn Surf>;
    let pl = |o: V, a: V, b: V, ext: f64| Box::new(Plane { o, a, b, ext }) as Box<dyn Surf>;
    let fams = ["plane x sphere", "sphere x sphere", "cyl x cyl (perp, offset)", "plane x cyl (tilted)", "plane x cyl (parallel)", "sphere x cyl (coaxial)", "torus x plane (perp)", "torus x plane (axial)", "cone x plane (perp)", "plane x plane (coincident)"];
    let (mut fails, mut tot_ms, mut n) = (0, 0.0, 0);
    let mut msgs = vec![];
    for (fi, fam) in fams.iter().enumerate() {
        let (mut ok, mut worst_res, mut worst_len) = (0, 0.0f64, 0.0f64);
        for k in 0..trials {
            let crit = k % 5 == 4;
            let (a, b, exp, len): (Box<dyn Surf>, Box<dyn Surf>, Expect, f64);
            let e0 = Expect::default();
            match fi {
                0 => {
                    let d = if crit { 1.0 } else { 1.5 * rng.f() };
                    a = pl([0.0, 0.0, d], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 3.0);
                    b = mv(Box::new(Sphere { r: 1.0 }), rand_rot(&mut rng), [0.0; 3]);
                    (exp, len) = if d.abs() < 1.0 { (Expect { closed: 1, ..e0 }, 2.0 * PI * (1.0 - d * d).sqrt()) } else if d.abs() == 1.0 { (Expect { points: 1, ..e0 }, 0.0) } else { (e0, 0.0) };
                }
                1 => {
                    let r2 = 0.3 + 1.2 * (0.5 + 0.5 * rng.f());
                    let e = if crit { 1.0 + r2 } else { 3.0 * (0.5 + 0.5 * rng.f()) };
                    a = mv(Box::new(Sphere { r: 1.0 }), rand_rot(&mut rng), [0.0; 3]);
                    b = mv(Box::new(Sphere { r: r2 }), rand_rot(&mut rng), [e, 0.0, 0.0]);
                    (exp, len) = if e < 1.0 + r2 && e > (1.0 - r2).abs() { let x = (e * e + 1.0 - r2 * r2) / (2.0 * e); (Expect { closed: 1, ..e0 }, 2.0 * PI * (1.0 - x * x).sqrt()) } else if e == 1.0 + r2 { (Expect { points: 1, ..e0 }, 0.0) } else { (e0, 0.0) };
                }
                2 => {
                    let r2 = 0.2 + 0.7 * (0.5 + 0.5 * rng.f());
                    let e = if crit { 1.0 - r2 } else { 2.0 * (0.5 + 0.5 * rng.f()) };
                    a = unitz(1.0, 3.0);
                    b = mv(unitz(r2, 3.0), rot_to([1.0, 0.0, 0.0]), [0.0, e, 0.0]);
                    let br = move |f: f64| { let y = e + r2 * f.cos(); [(1.0 - y * y).max(0.0).sqrt(), y, r2 * f.sin()] };
                    (exp, len) = if e + r2 < 1.0 { (Expect { closed: 2, ..e0 }, 2.0 * numeric_len(br, 0.0, 2.0 * PI)) }
                        else if e + r2 == 1.0 { (Expect { branch_arcs: 2, ..e0 }, 2.0 * numeric_len(br, 0.0, 2.0 * PI)) }
                        else if e - r2 < 1.0 { let f0 = ((1.0 - e) / r2).clamp(-1.0, 1.0).acos(); (Expect { closed: 1, ..e0 }, 2.0 * numeric_len(br, f0, 2.0 * PI - f0)) }
                        else if e - r2 == 1.0 { (Expect { points: 1, ..e0 }, 0.0) } else { (e0, 0.0) };
                }
                3 => {
                    let (th, d) = (0.6 * rng.f(), 0.5 * rng.f());
                    a = mv(pl([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 4.0), rot_to([th.sin(), 0.0, th.cos()]), [0.0, 0.0, d]);
                    b = mv(unitz(1.0, 3.0), about_z(&mut rng), [0.0; 3]);
                    exp = Expect { closed: 1, ..e0 };
                    len = numeric_len(|f| { let (x, y) = (f.cos(), f.sin()); [x, y, d / th.cos() - x * th.tan()] }, 0.0, 2.0 * PI);
                }
                4 => {
                    let d = if crit { 1.0 } else { 1.5 * (0.5 + 0.5 * rng.f()) };
                    a = pl([d, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 4.0);
                    b = mv(unitz(1.0, 2.0), about_z(&mut rng), [0.0; 3]);
                    (exp, len) = if d < 1.0 { (Expect { open: 2, ..e0 }, 8.0) } else if d == 1.0 { (Expect { tangential: 1, ..e0 }, 4.0) } else { (e0, 0.0) };
                }
                5 => {
                    let r = if crit { 1.0 } else { 0.3 + (0.5 + 0.5 * rng.f()) };
                    a = mv(Box::new(Sphere { r: 1.0 }), rand_rot(&mut rng), [0.0; 3]);
                    b = mv(unitz(r, 2.0), about_z(&mut rng), [0.0; 3]);
                    (exp, len) = if r < 1.0 { (Expect { closed: 2, ..e0 }, 4.0 * PI * r) } else if r == 1.0 { (Expect { tangential: 1, ..e0 }, 2.0 * PI) } else { (e0, 0.0) };
                }
                6 => {
                    let h = if crit { 0.5 } else { 0.7 * rng.f() };
                    a = Box::new(Torus { big: 2.0, r: 0.5 });
                    b = pl([0.0, 0.0, h], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 4.0);
                    (exp, len) = if h.abs() < 0.5 { let w = (0.25 - h * h).sqrt(); (Expect { closed: 2, ..e0 }, 2.0 * PI * (2.0 + w) + 2.0 * PI * (2.0 - w)) } else if h.abs() == 0.5 { (Expect { tangential: 1, ..e0 }, 4.0 * PI) } else { (e0, 0.0) };
                }
                7 => {
                    let ang = PI * rng.f();
                    a = Box::new(Torus { big: 2.0, r: 0.5 });
                    b = pl([0.0; 3], [ang.cos(), ang.sin(), 0.0], [0.0, 0.0, 1.0], 4.0);
                    (exp, len) = (Expect { closed: 2, ..e0 }, 2.0 * PI);
                }
                8 => {
                    let h = 0.3 + 1.6 * (0.5 + 0.5 * rng.f());
                    let t = (PI / 6.0).tan();
                    a = mv(Box::new(Cone { t, v0: 0.2, v1: 2.0 }), about_z(&mut rng), [0.0; 3]);
                    b = pl([0.0, 0.0, h], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 4.0);
                    (exp, len) = (Expect { closed: 1, ..e0 }, 2.0 * PI * h * t);
                }
                _ => {
                    a = pl([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0);
                    let (c, s) = (0.3f64.cos(), 0.3f64.sin());
                    b = pl([0.1, 0.2, 0.0], [c, s, 0.0], [-s, c, 0.0], 1.0);
                    (exp, len) = (e0, 0.0);
                }
            }
            let (g, gt) = (rand_rot(&mut rng), [3.0 * rng.f(), 3.0 * rng.f(), 3.0 * rng.f()]);
            let (a, b) = (mv(a, g, gt), mv(b, g, gt));
            let t0 = Instant::now();
            let (cv, _) = intersect(a.as_ref(), b.as_ref());
            tot_ms += ms(t0);
            n += 1;
            let (mut got, mut coincident, mut glen) = (Expect::default(), false, 0.0);
            for c in &cv {
                match (c.kind, c.ends) {
                    (Kind::Isolated, (End::Fail, _)) => coincident = true,
                    (Kind::Isolated, _) => got.points += 1,
                    (Kind::Tangential, _) => got.tangential += 1,
                    (_, (End::Closed, _)) => got.closed += 1,
                    (_, (End::Branch, _)) | (_, (_, End::Branch)) => got.branch_arcs += 1,
                    _ => got.open += 1,
                }
                glen += c.pts.windows(2).map(|w| norm(sub(w[1], w[0]))).sum::<f64>();
            }
            let res = cv.iter().flat_map(|c| c.pts.iter()).map(|p| dist(a.as_ref(), *p).max(dist(b.as_ref(), *p))).fold(0.0, f64::max);
            let len_err = if len > 0.0 { (glen - len).abs() / len } else { glen };
            let pass = if fi == 9 { coincident } else { got == exp && len_err < 2e-3 && res < 1e-9 };
            if pass { ok += 1; } else { msgs.push(format!("{fam} #{k}{}: expect {exp:?} len {len:.4}, got {got:?} len {glen:.4} res {res:.1e}", if crit { " (critical)" } else { "" })); }
            worst_res = worst_res.max(res);
            if len > 0.0 { worst_len = worst_len.max(len_err); }
        }
        fails += trials - ok;
        println!("  {fam:<28} {ok:3}/{trials} pass   worst residual {worst_res:.1e}   worst length error {worst_len:.1e}");
    }
    println!("  mean time per pair {:.2} ms", tot_ms / n as f64);
    for m in msgs.iter().take(20) { println!("    FAIL {m}"); }
    fails
}
