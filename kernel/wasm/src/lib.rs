//! wasm32 exports of kcore for the browser tests (shared memory, threads via workers).
#![no_std]
use core::slice::{from_raw_parts as sl, from_raw_parts_mut as slm};

#[panic_handler]
fn panic(_: &core::panic::PanicInfo) -> ! { core::arch::wasm32::unreachable() }

#[no_mangle] pub extern "C" fn patch_work(seed: u64, n: u32) -> u64 { kcore::patch_work(seed, n) }
#[no_mangle] pub unsafe extern "C" fn fill(p: *mut f64, n: usize, seed: u64) { kcore::fill(slm(p, n), seed) }
#[no_mangle] pub unsafe extern "C" fn sum(p: *const f64, n: usize) -> f64 { kcore::sum_seq(sl(p, n)) }
#[no_mangle] pub unsafe extern "C" fn hash(p: *const u8, n: usize) -> u64 { kcore::hash_bytes(sl(p, n)) }
#[no_mangle] pub extern "C" fn orient_hash(seed: u64, n: u32) -> u64 { kcore::orient_hash(seed, n) }
#[no_mangle] pub extern "C" fn sin_hash(seed: u64, n: u32) -> u64 { kcore::sin_hash(seed, n) }
#[no_mangle] pub extern "C" fn mesh_scratch_len(nv: usize) -> usize { kcore::mesh_scratch_len(nv) }
#[no_mangle] pub unsafe extern "C" fn mesh_work(pos: *const f32, nf: usize, scr: *mut u64, ns: usize) -> u64 { kcore::mesh_work(sl(pos, nf), slm(scr, ns)) }

// ---- CPU rasterizer (kcore::raster). Blocks laid out by raster.js:
// mesh: v f64[3nv], t u32[3nt], n f64[3nt], c f32[3nt], clusters Cluster[nt / CL]; cache: Vtx[nv];
// chunk descriptor: u32 [rec ptr, head ptr, tail ptr, pool ptr, pool entries, need, record capacity];
// frame: color u32[w h], depth f32[w h], id i32[w h].
use kcore::raster::{Bins, Cam, Cluster, Frame, Mesh, Rec, Vtx, CL};
const MAXCH: usize = 64;
const DESC: usize = 7;
unsafe fn frame<'a>(p: *mut u8, n: usize) -> Frame<'a> {
    Frame { color: slm(p as *mut u32, n), depth: slm(p.add(4 * n) as *mut f32, n), id: slm(p.add(8 * n) as *mut i32, n) }
}
unsafe fn bins<'a>(d: *mut u32, nb: usize) -> Bins<'a> {
    let (rec, head, tail, pool, np, nrec) = (*d as *mut Rec, *d.add(1) as *mut u32, *d.add(2) as *mut u32, *d.add(3) as *mut u32, *d.add(4) as usize, *d.add(6) as usize);
    Bins { rec: slm(rec, nrec), head: slm(head, nb), tail: slm(tail, nb), pool: slm(pool, 2 * np), need: &mut *d.add(5) }
}
unsafe fn mesh<'a>(v: *const f64, t: *const u32, n: *const f64, c: *const f32, nv: usize, nt: usize) -> Mesh<'a> {
    Mesh { v: sl(v, 3 * nv), t: sl(t, 3 * nt), n: sl(n, 3 * nt), c: sl(c, 3 * nt) }
}
#[no_mangle] pub extern "C" fn raster_sizes(what: u32) -> usize {
    match what { 0 => core::mem::size_of::<Vtx>(), 1 => core::mem::size_of::<Rec>(), 2 => core::mem::size_of::<Cluster>(), 3 => CL, _ => kcore::raster::RECS_PER_TRI }
}
/// draw order (new -> old triangle index) into `perm`; `keys` is scratch of 16 nt bytes
#[no_mangle] pub unsafe extern "C" fn raster_order(v: *const f64, t: *const u32, face: *const u32, nv: usize, nt: usize, keys: *mut u128, perm: *mut u32) {
    kcore::raster::order(sl(v, 3 * nv), sl(t, 3 * nt), sl(face, nt), slm(keys, nt), slm(perm, nt))
}
/// renumber vertices by first use in `t` (in place); `old` receives new -> old, `new` is scratch (old -> new)
#[no_mangle] pub unsafe extern "C" fn raster_renumber(t: *mut u32, nv: usize, nt: usize, new: *mut u32, old: *mut u32) {
    kcore::raster::renumber(slm(t, 3 * nt), slm(new, nv), slm(old, nv))
}
#[no_mangle] pub unsafe extern "C" fn raster_clusters(v: *const f64, t: *const u32, n: *const f64, c: *const f32, nv: usize, nt: usize, out: *mut Cluster) {
    kcore::raster::clusters(&mesh(v, t, n, c, nv, nt), slm(out, nt.div_ceil(CL)))
}
/// phase 1 for triangles [t0, t1) into the chunk described at `desc`
#[no_mangle] pub unsafe extern "C" fn raster_bin(cam: *const Cam, v: *const f64, t: *const u32, n: *const f64, c: *const f32, nv: usize, nt: usize,
    cl: *const Cluster, cache: *mut Vtx, desc: *mut u32, t0: usize, t1: usize) {
    let cam = &*cam;
    kcore::raster::bin(&mesh(v, t, n, c, nv, nt), sl(cl, nt.div_ceil(CL)), cam, slm(cache, nv), &mut bins(desc, cam.nb as usize), t0, t1);
}
/// phase 2: bands first, first + step, ... of the frame at `fr`, from `k` chunks (descriptors at `descs`)
#[no_mangle] pub unsafe extern "C" fn raster_draw(cam: *const Cam, descs: *mut u32, k: usize, fr: *mut u8, first: usize, step: usize) -> f64 {
    let cam = &*cam;
    let nb = cam.nb as usize;
    let k = k.min(MAXCH);
    let mut ch: [core::mem::MaybeUninit<Bins>; MAXCH] = core::mem::MaybeUninit::uninit().assume_init();
    for i in 0..k { ch[i].write(bins(descs.add(DESC * i), nb)); }
    let ch = core::slice::from_raw_parts(ch.as_ptr() as *const Bins, k);
    let mut fr = frame(fr, cam.w as usize * cam.h as usize);
    let mut drawn = 0u64;
    let mut b = first;
    while b < nb { drawn += kcore::raster::raster(cam, ch, &mut fr, b); b += step.max(1); }
    drawn as f64
}
