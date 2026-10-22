#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
cargo build -p guisdap-web --target wasm32-unknown-unknown --release
wasm-bindgen --target web --out-dir web/pkg target/wasm32-unknown-unknown/release/guisdap_web.wasm

# These public library classes are imported by the shared plotting module.
for plot_class in FieldView SeriesView PlotViewport; do
  grep -q "export class ${plot_class}" web/pkg/guisdap_web.js || {
    echo "Missing Rustplotlib export: ${plot_class}" >&2
    exit 1
  }
done
