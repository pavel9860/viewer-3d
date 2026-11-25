#!/bin/sh
# Build kwasm with real atomics (nightly -Zbuild-std): shared imported memory, exported stack pointer for workers.
set -e
cd "$(dirname "$0")"
RUSTFLAGS="-C target-feature=+atomics,+bulk-memory,+mutable-globals -C link-arg=--shared-memory -C link-arg=--import-memory -C link-arg=--max-memory=4294967296 -C link-arg=--export=__stack_pointer" \
  cargo +nightly build -j2 --release --target wasm32-unknown-unknown -Zbuild-std=core,panic_abort
cp target/wasm32-unknown-unknown/release/kwasm.wasm ../web/k.wasm
