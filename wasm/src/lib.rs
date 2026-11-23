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
