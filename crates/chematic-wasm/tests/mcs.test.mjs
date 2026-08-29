// Tests for mcs_smiles_json_with_config (Track A/99-point directive Phase 3):
// full McsConfig + McsOutcome exposed to WASM, mirroring the Rust
// chematic_smarts::find_mcs_with_config_checked API and the Python
// find_mcs_checked() binding's (mol, was_timed_out) shape.
//
// Run manually after building the Node target:
//
//   wasm-pack build crates/chematic-wasm --target nodejs --out-dir pkg-node --release
//   node crates/chematic-wasm/tests/mcs.test.mjs

import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");

const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const { mcs_smiles_json, mcs_smiles_json_with_config } = wasm;

// Default config (`{}`) matches the bare mcs_smiles_json's own default.
// (mcs_smiles_json returns a plain SMILES string, or the literal string
// "null" -- not JSON -- so it's compared directly, not JSON.parse'd.)
{
  const smiles = JSON.stringify(["CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1"]);
  const full = JSON.parse(mcs_smiles_json_with_config(smiles, "{}"));
  const bare = mcs_smiles_json(smiles);
  assert.equal(full.smiles, bare);
  assert.equal(full.wasTimedOut, false);
}

// No common substructure -> smiles: null, not an error.
{
  const smiles = JSON.stringify(["[He]", "[Ne]"]);
  const result = JSON.parse(mcs_smiles_json_with_config(smiles, "{}"));
  assert.equal(result.smiles, null);
  assert.equal(result.wasTimedOut, false);
}

// matchCharge:true changes the result vs the false default.
{
  const smiles = JSON.stringify(["CC(=O)[O-]", "CC(=O)O"]);
  const without = mcs_smiles_json_with_config(smiles, JSON.stringify({ matchCharge: false }));
  const withCharge = mcs_smiles_json_with_config(smiles, JSON.stringify({ matchCharge: true }));
  assert.notEqual(without, withCharge, "matchCharge:true must change the MCS result");
}

// atomCompare:"any_heavy_atom" widens the match vs the "elements" default.
{
  const smiles = JSON.stringify(["c1ccccc1", "c1ccncc1"]);
  const elements = mcs_smiles_json_with_config(smiles, JSON.stringify({ atomCompare: "elements" }));
  const anyHeavy = mcs_smiles_json_with_config(
    smiles,
    JSON.stringify({ atomCompare: "any_heavy_atom" })
  );
  assert.notEqual(elements, anyHeavy, "any_heavy_atom compare must change the MCS result");
}

// timeoutMs:0 is reported as timed out.
{
  const smiles = JSON.stringify(["CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1"]);
  const result = JSON.parse(mcs_smiles_json_with_config(smiles, JSON.stringify({ timeoutMs: 0 })));
  assert.equal(result.wasTimedOut, true);
}

// Invalid atomCompare string is a real thrown error, not a silent default.
{
  const smiles = JSON.stringify(["CCO", "CCO"]);
  assert.throws(
    () => mcs_smiles_json_with_config(smiles, JSON.stringify({ atomCompare: "bogus" })),
    "invalid atomCompare must throw"
  );
}

// Unknown config field is rejected (deny_unknown_fields), not silently ignored.
{
  const smiles = JSON.stringify(["CCO", "CCO"]);
  assert.throws(
    () => mcs_smiles_json_with_config(smiles, JSON.stringify({ notAField: true })),
    "unknown config field must throw"
  );
}

// Fewer than 2 SMILES is rejected.
{
  assert.throws(
    () => mcs_smiles_json_with_config(JSON.stringify(["CCO"]), "{}"),
    "fewer than 2 SMILES must throw"
  );
}

console.log("mcs.test.mjs: all assertions passed");
