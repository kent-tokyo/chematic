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
const { mcs_smiles_json, mcs_smiles_json_with_config, mcs_smarts_json_with_config } = wasm;

{
  const result = JSON.parse(mcs_smarts_json_with_config(
    JSON.stringify(["NO", "CO"]), JSON.stringify({ atomCompare: "any" })
  ));
  assert.match(result.smarts, /,/);
  assert.equal(result.wasTimedOut, false);
  assert.equal(result.wasResourceLimited, false);
  const constrained = JSON.parse(mcs_smarts_json_with_config(
    JSON.stringify(["[13CH3][NH3+]", "[NH3+][13CH3]"]),
    JSON.stringify({ matchIsotope: true, matchCharge: true })
  ));
  assert.match(constrained.smarts, /;13\*/);
  assert.match(constrained.smarts, /;\+1/);
  assert.equal(constrained.wasTimedOut, false);
  assert.equal(constrained.wasResourceLimited, false);
  assert.throws(() => mcs_smarts_json_with_config(
    JSON.stringify(["C[C@H](O)Cl", "C[C@H](O)Cl"]),
    JSON.stringify({ matchChiralTag: true })
  ), /matchChiralTag/);
}

{
  const inputs = JSON.stringify(["CCCCC.C12C3C1C23", "CCCCC.C12C3C1C23"]);
  const atomFirst = JSON.parse(mcs_smarts_json_with_config(inputs, "{}"));
  const bondFirst = JSON.parse(mcs_smarts_json_with_config(
    inputs, JSON.stringify({ bondCountFirst: true })
  ));
  assert.equal(atomFirst.smarts.match(/#6/g).length, 5);
  assert.equal(bondFirst.smarts.match(/#6/g).length, 4);
  assert.equal(bondFirst.wasTimedOut, false);
  assert.equal(bondFirst.wasResourceLimited, false);
}

// Default config (`{}`) matches the bare mcs_smiles_json's own default.
// (mcs_smiles_json returns a plain SMILES string, or the literal string
// "null" -- not JSON -- so it's compared directly, not JSON.parse'd.)
{
  const smiles = JSON.stringify(["CC(=O)Oc1ccccc1C(=O)O", "CC(=O)Nc1ccc(O)cc1"]);
  const full = JSON.parse(mcs_smiles_json_with_config(smiles, "{}"));
  const bare = mcs_smiles_json(smiles);
  assert.equal(full.smiles, bare);
  assert.equal(full.wasTimedOut, false);
  assert.equal(full.wasResourceLimited, false);
  // Regression test: qmol_to_molecule once matched only the bare
  // AtomQuery::Primitive(AtomicNum(n)) case, but at the time find_mcs's query
  // atoms were always the compound And(AtomicNum(n), Aromatic(bool)) -- so
  // every atom silently fell through to carbon. Aspirin/paracetamol's shared
  // hydroxyphenyl MCS must keep its oxygen atom.
  assert.match(full.smiles, /[Oo]/, `MCS lost its oxygen atom: ${full.smiles}`);
}

// Both result formats retain an input dummy atom.
{
  const input = JSON.stringify(["[*]O", "[*]O"]);
  const query = JSON.parse(mcs_smarts_json_with_config(input, "{}"));
  const projected = JSON.parse(mcs_smiles_json_with_config(input, "{}"));
  assert.match(query.smarts, /\[#0\]/);
  assert.match(projected.smiles, /\*/);
}

// Resource limits are distinct from timeouts and cannot be mistaken for an
// exhaustive result.
{
  const smiles = JSON.stringify(["CCC", "CCC"]);
  const result = JSON.parse(
    mcs_smiles_json_with_config(smiles, JSON.stringify({ maxSearchNodes: 0 }))
  );
  assert.equal(result.wasTimedOut, false);
  assert.equal(result.wasResourceLimited, true);
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
