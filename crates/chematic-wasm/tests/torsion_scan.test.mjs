// WASM public-contract test for the playground's torsion scan. This verifies
// the wasm-bindgen export itself, not merely the Rust helper behind it.
import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "..", "..", "..");
const wasm = await import(path.join(root, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));

assert.equal(typeof wasm.torsion_scan_json, "function");

const butane = wasm.parse_smiles("CCCC");
try {
  const points = JSON.parse(wasm.torsion_scan_json(butane, 0, 1, 2, 3, 12));
  assert.equal(points.length, 12);
  assert.equal(points[0].angle, 0);
  assert.ok(points.every(({ angle, energy }) => Number.isFinite(angle) && Number.isFinite(energy)));
  assert.ok(new Set(points.map(({ energy }) => energy)).size > 1, "scan must vary with torsion angle");

  const outOfRange = JSON.parse(wasm.torsion_scan_json(butane, 0, 1, 2, 4, 12));
  assert.match(outOfRange.error, /atom index out of range/);
} finally {
  butane.free();
}

console.log("torsion_scan.test.mjs: all assertions passed");
