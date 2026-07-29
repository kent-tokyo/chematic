// Tests for retro_disconnect_json (issue #91): single-step retrosynthetic
// disconnection exposed to WASM/browser consumers, mirroring the Rust
// chematic_rxn::retro::retro_disconnect API and the Python
// Mol.retro_disconnect() binding's JSON field shape.
//
// Not wired into any CI workflow (this repo has no WASM-test CI job). Run
// manually after building the Node target:
//
//   wasm-pack build crates/chematic-wasm --target nodejs --out-dir pkg-node --release
//   node crates/chematic-wasm/tests/retro_disconnect.test.mjs

import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");

const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const { parse_smiles, retro_disconnect_json } = wasm;

// acetanilide: has an amide bond, the standard chematic-rxn::retro doc example.
const acetanilide = parse_smiles("CC(=O)Nc1ccccc1");

// Empty result (no disconnectable bond) is a valid `[]`, not an error.
{
  const methane = parse_smiles("C");
  const json = retro_disconnect_json(methane, 20, "");
  assert.deepEqual(JSON.parse(json), []);
}

// Basic shape + field names, matching the Python binding's dict keys.
{
  const json = retro_disconnect_json(acetanilide, 20, "");
  const results = JSON.parse(json);
  assert.ok(results.length >= 1, "acetanilide should have >=1 disconnection");
  for (const r of results) {
    assert.equal(typeof r.template, "string");
    assert.equal(typeof r.reaction_class, "string");
    assert.ok(Array.isArray(r.precursors));
    assert.ok(Array.isArray(r.sa_scores));
    assert.equal(r.sa_scores.length, r.precursors.length);
    assert.equal(typeof r.max_sa_score, "number");
  }
}

// reaction_class filter narrows results; an unrelated class returns [].
{
  const filtered = JSON.parse(retro_disconnect_json(acetanilide, 0, "AmideBond"));
  assert.ok(filtered.length >= 1);
  for (const r of filtered) {
    assert.equal(r.reaction_class, "AmideBond");
  }
  const unrelated = JSON.parse(retro_disconnect_json(acetanilide, 0, "Ether"));
  assert.deepEqual(unrelated, []);
}

// max_results caps output.
{
  const capped = JSON.parse(retro_disconnect_json(acetanilide, 1, ""));
  assert.ok(capped.length <= 1);
}

// Unknown reaction_class is a real thrown error, not a silent empty result.
{
  assert.throws(
    () => retro_disconnect_json(acetanilide, 20, "NotARealClass"),
    "unrecognized reaction_class must throw"
  );
}

console.log("retro_disconnect.test.mjs: all assertions passed");
