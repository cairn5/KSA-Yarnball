# Advanced Flight Planner

Advanced flight planning for KSA and beyond.

Rust trajectory engine compiled to WebAssembly, called from a TypeScript web interface.

```
engine/   Rust library (cdylib + rlib). Headless maths, exposed via wasm-bindgen.
web/      Vite + TypeScript interface. Imports the wasm build from web/src/engine/ (generated).
          src/plot/  full-screen canvas plan view: log radial scale, pan/zoom, layers of lines and markers.
          src/routeWorker.ts  runs the flyby-sequence search (engine/src/mga.rs) in a Web Worker.
          src/refineWorker.ts one SADE island (engine/src/mga1dsm.rs, sade.rs) per Web Worker.
tools/    pykep_reference.py regenerates engine/tests/data/pykep_reference.json, which the
          engine's tests compare against (needs pykep 2 in a conda environment).
```

## Prerequisites

- Rust with the wasm target: `rustup target add wasm32-unknown-unknown`
- wasm-pack: `cargo install wasm-pack`
- Node.js LTS

## Run

```bash
cd web
npm install
npm run dev
```

`npm run dev` rebuilds the engine with wasm-pack, then starts Vite on http://localhost:5173.
After changing Rust code, run `npm run wasm` (Vite reloads automatically).

## Test

```bash
cargo test
```
