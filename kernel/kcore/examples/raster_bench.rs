//! native timing of kcore::raster on a `brep mesh --out` file: raster_bench FILE.mesh [W H]
//! Reports, per phase, the sum over 6 views of the fastest of 9 repetitions (robust to a busy machine).
use kcore::raster::*;
use std::time::Instant;

// ---- self-sampling profiler (PROF=1): SIGPROF every 100 us of CPU time records the interrupted instruction
// pointer; addresses are printed relative to the executable's load base, for `addr2line -f -i -e`.
mod prof {
    use std::sync::atomic::{AtomicUsize, Ordering};
    pub static N: AtomicUsize = AtomicUsize::new(0);
    pub static mut PC: [usize; 1 << 20] = [0; 1 << 20];
    #[repr(C)]
    struct SigAction { handler: usize, mask: [u64; 16], flags: i32, restorer: usize }
    #[repr(C)]
    struct ITimer { iv_s: i64, iv_us: i64, v_s: i64, v_us: i64 }
    extern "C" {
        fn sigaction(sig: i32, act: *const SigAction, old: *mut SigAction) -> i32;
        fn setitimer(which: i32, new: *const ITimer, old: *mut ITimer) -> i32;
    }
    extern "C" fn on_prof(_sig: i32, _info: usize, uc: *const u8) {
        // x86_64 glibc ucontext_t: uc_mcontext.gregs at byte 40, REG_RIP = 16
        let rip = unsafe { *(uc.add(40 + 16 * 8) as *const usize) };
        let i = N.fetch_add(1, Ordering::Relaxed);
        if i < 1 << 20 { unsafe { PC[i] = rip; } }
    }
    pub fn start() {
        const SA_SIGINFO: i32 = 4;
        const SA_RESTART: i32 = 0x1000_0000;
        let act = SigAction { handler: on_prof as usize, mask: [0; 16], flags: SA_SIGINFO | SA_RESTART, restorer: 0 };
        unsafe {
            sigaction(27, &act, std::ptr::null_mut());
            setitimer(2, &ITimer { iv_s: 0, iv_us: 100, v_s: 0, v_us: 100 }, std::ptr::null_mut());
        }
    }
    pub fn dump(path: &str) {
        unsafe { setitimer(2, &ITimer { iv_s: 0, iv_us: 0, v_s: 0, v_us: 0 }, std::ptr::null_mut()); }
        let maps = std::fs::read_to_string("/proc/self/maps").unwrap();
        let exe = std::fs::read_link("/proc/self/exe").unwrap();
        let base = maps.lines().find(|l| l.ends_with(exe.to_str().unwrap())).map(|l| usize::from_str_radix(l.split('-').next().unwrap(), 16).unwrap()).unwrap();
        let n = N.load(Ordering::Relaxed).min(1 << 20);
        let mut out = String::new();
        for i in 0..n { out += &format!("{:x}\n", unsafe { PC[i] }.wrapping_sub(base)); }
        std::fs::write(path, out).unwrap();
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let b = std::fs::read(&a[1]).expect("mesh file");
    let (w, h): (u32, u32) = (a.get(2).map_or(800, |s| s.parse().unwrap()), a.get(3).map_or(600, |s| s.parse().unwrap()));
    let u = |o: usize| u32::from_le_bytes(b[o..o + 4].try_into().unwrap());
    let (nv, nt) = (u(0) as usize, u(4) as usize);
    let p: Vec<f64> = (0..3 * nv).map(|i| f32::from_le_bytes(b[8 + 4 * i..12 + 4 * i].try_into().unwrap()) as f64).collect();
    let t0: Vec<u32> = (0..3 * nt).map(|i| u(8 + 12 * nv + 4 * i)).collect();
    let face: Vec<u32> = (0..nt).map(|i| u(8 + 12 * nv + 12 * nt + 4 * i)).collect();
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for i in 0..nv { for k in 0..3 { lo[k] = lo[k].min(p[3 * i + k]); hi[k] = hi[k].max(p[3 * i + k]); } }
    let s = 2.0 / ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2)).sqrt();
    let v: Vec<f64> = (0..3 * nv).map(|i| (p[i] - (lo[i % 3] + hi[i % 3]) / 2.0) * s).collect();
    let (mut keys, mut perm) = (vec![0u128; nt], vec![0u32; nt]);
    if std::env::var("NOORDER").is_ok() { for k in 0..nt { perm[k] = k as u32; } } else { order(&v, &t0, &face, &mut keys, &mut perm); }
    let mut t: Vec<u32> = (0..3 * nt).map(|i| t0[3 * perm[i / 3] as usize + i % 3]).collect();
    let (mut vnew, mut vold) = (vec![0u32; nv], vec![0u32; nv]);
    renumber(&mut t, &mut vnew, &mut vold);
    let v: Vec<f64> = (0..3 * nv).map(|i| v[3 * vold[i / 3] as usize + i % 3]).collect();
    let mut n = vec![0.0; 3 * nt];
    for k in 0..nt {
        let (a, b, c) = (3 * t[3 * k] as usize, 3 * t[3 * k + 1] as usize, 3 * t[3 * k + 2] as usize);
        let (u, q) = ([v[b] - v[a], v[b + 1] - v[a + 1], v[b + 2] - v[a + 2]], [v[c] - v[a], v[c + 1] - v[a + 1], v[c + 2] - v[a + 2]]);
        let x = [u[1] * q[2] - u[2] * q[1], u[2] * q[0] - u[0] * q[2], u[0] * q[1] - u[1] * q[0]];
        let l = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
        if l > 0.0 { n[3 * k..3 * k + 3].copy_from_slice(&[x[0] / l, x[1] / l, x[2] / l]); }
    }
    let c = vec![0.7f32; 3 * nt];
    let m = Mesh { v: &v, t: &t, n: &n, c: &c };
    let mut vc = vec![Vtx::default(); nv];
    let mut cls = vec![Cluster::default(); nt.div_ceil(CL)];
    clusters(&m, &mut cls);
    let bh = 16u32;
    let nb = h.div_ceil(bh) as usize;
    let mut rec = vec![Rec::default(); RECS_PER_TRI * nt];
    let (mut head, mut tail, mut pool, mut need) = (vec![0u32; nb], vec![0u32; nb], vec![0u32; 2 * (4 * nt + 1024)], 0u32);
    let npx = (w * h) as usize;
    let (mut col, mut dep, mut id) = (vec![0u32; npx], vec![0f32; npx], vec![0i32; npx]);
    let (mut tb, mut td, mut px) = (0.0, 0.0, 0u64);
    let mut frame = 0;
    if std::env::var("PROF").is_ok() { prof::start(); }
    for &(ax, ay) in &[(0.0f64, 0.0f64), (0.6, -0.5), (-0.6, 0.5), (1.2, 0.8), (0.3, 1.7), (-0.9, -1.1)] {
        let (mut bb, mut dd) = (f64::MAX, f64::MAX);
        for _ in 0..9 {
            frame += 1;
            let (cx, sx, cy, sy) = (ax.cos(), ax.sin(), ay.cos(), ay.sin());
            let cam = Cam { r: [cy, 0.0, sy, sx * sy, cx, -sx * cy, -cx * sy, sx, cx * cy], t: [0.0, 0.0, -2.5], f: (h as f64 / 2.0) / 0.4f64.tan(), w, h, flags: 1 | 2 | 4, frame, bh, nb: nb as u32 };
            let t0 = Instant::now();
            let mut bi = Bins { rec: &mut rec, head: &mut head, tail: &mut tail, pool: &mut pool, need: &mut need };
            bin(&m, &cls, &cam, &mut vc, &mut bi, 0, nt);
            let t1 = Instant::now();
            let ch = [bi];
            let mut fr = Frame { color: &mut col, depth: &mut dep, id: &mut id };
            let mut p = 0;
            for band in 0..nb { p += raster(&cam, &ch, &mut fr, band); }
            let t2 = Instant::now();
            bb = bb.min((t1 - t0).as_secs_f64());
            dd = dd.min((t2 - t1).as_secs_f64());
            px = p;
        }
        tb += bb / 6.0;
        td += dd / 6.0;
    }
    if std::env::var("ALL").is_ok() {
        // every band drawn twice in a row, the second timed: draw cost with the band hot in cache
        let cam = Cam { r: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], t: [0.0; 3], f: 1.0, w, h, flags: 6, frame: 1, bh, nb: nb as u32 };
        let ch = [Bins { rec: &mut rec, head: &mut head, tail: &mut tail, pool: &mut pool, need: &mut need }];
        let mut fr = Frame { color: &mut col, depth: &mut dep, id: &mut id };
        let (mut cold, mut hot) = (0.0, 0.0);
        for band in 0..nb {
            let (mut c, mut hh) = (f64::MAX, f64::MAX);
            for _ in 0..20 {
                let t0 = Instant::now(); raster(&cam, &ch, &mut fr, band); let t1 = Instant::now(); raster(&cam, &ch, &mut fr, band);
                c = c.min((t1 - t0).as_secs_f64()); hh = hh.min(t1.elapsed().as_secs_f64());
            }
            cold += c; hot += hh;
        }
        println!("last view, per band best of 20: first draw {:.2} ms, repeated {:.2} ms", cold * 1e3, hot * 1e3);
        return;
    }
    if let Ok(ri) = std::env::var("RI") {
        // one record drawn alone into its band, fastest of many: tri cost without noise (minus the band clear)
        let rr = rec[ri.parse::<usize>().unwrap()];
        let cam = Cam { r: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], t: [0.0; 3], f: 1.0, w, h, flags: 6, frame: 1, bh, nb: nb as u32 };
        let band = (rr.y.iter().cloned().fold(f64::MAX, f64::min) / 16.0 - 0.5).ceil().max(0.0) as usize / bh as usize;
        let (mut best, mut best0) = (f64::MAX, f64::MAX);
        for with in [true, false] {
            let mut rec1 = vec![rr; 1];
            let (mut hd, mut tl, mut pl, mut nd) = (vec![u32::MAX; nb], vec![0u32; nb], vec![0u32, u32::MAX], 0u32);
            if with { hd[band] = 0; }
            let ch = [Bins { rec: &mut rec1, head: &mut hd, tail: &mut tl, pool: &mut pl, need: &mut nd }];
            let mut fr = Frame { color: &mut col, depth: &mut dep, id: &mut id };
            for _ in 0..3000 {
                let t0 = Instant::now();
                let p = raster(&cam, &ch, &mut fr, band);
                let e = t0.elapsed().as_secs_f64();
                if with { px = p; }
                if with { best = best.min(e) } else { best0 = best0.min(e) }
            }
        }
        println!("record {ri}: x {:?} y {:?}: {:.0} ns, {px} px", rr.x, rr.y, (best - best0) * 1e9);
    }
    if let Ok(p) = std::env::var("PROF") { prof::dump(&p); }
    println!("{nt} tris {w}x{h}: bin {:.2} ms, draw {:.2} ms per frame (one thread, best of 9), {px} px written (last view)", tb * 1e3, td * 1e3);
}
