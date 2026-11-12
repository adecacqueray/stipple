#!/bin/sh
# Compile the Rust crate to WebAssembly for the site: ./build.sh
# Needs Rust and the wasm target:  rustup target add wasm32-unknown-unknown
set -e
cd "$(dirname "$0")"
RUSTFLAGS="-C target-feature=+simd128" cargo build --release --lib --target wasm32-unknown-unknown
cp target/wasm32-unknown-unknown/release/stipple.wasm web/stipple.wasm
ls -l web/stipple.wasm
