// Tests for 2D E/Z (double-bond cis/trans) direction perception wired into
// the WASM MOL/SDF bindings (mol_from_sdf_block, mol_from_v3000_block,
// sdf_to_records_json). See docs/stereo2d_reader_integration_rfc.md and
// crates/chematic-mol/tests/stereo2d_ez_reader_integration.rs for the
// Rust-side equivalents of these same fixtures (the (Z)-but-2-ene MOL block
// below is the exact same coordinates/atom order as that file's `but2ene_v2000`
// helper with Z_COORDS, so the expected canonical SMILES is a value already
// independently confirmed there and against a live RDKit oracle, not a fresh
// assumption made here).
//
// Not wired into any CI workflow (this repo has no WASM-test CI job). Run
// manually after building the Node target:
//
//   wasm-pack build crates/chematic-wasm --target nodejs --out-dir pkg-node --release
//   node crates/chematic-wasm/tests/stereo2d_ez_directions.test.mjs

import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");

const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const { mol_from_sdf_block, mol_from_v3000_block, sdf_to_records_json, mol_block_from_smiles } = wasm;

// (Z)-but-2-ene, V2000: Me0-C1=C2-Me3, zigzag "cis" layout -- same
// coordinates as stereo2d_ez_reader_integration.rs's Z_COORDS constant.
const Z_2BUTENE_V2000 = `but2ene
  chematic

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    2.3660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
M  END
`;

const Z_2BUTENE_V3000 = `but2ene
  chematic

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -0.8660 0.5000 0.0 0
M  V30 2 C 0.0 0.0 0.0 0
M  V30 3 C 1.5 0.0 0.0 0
M  V30 4 C 2.366 0.5 0.0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 END BOND
M  V30 END CTAB
M  END
`;

// Terminal alkene negative control (propene): no substituent at the =CH2
// end, so no E/Z is stereogenic at all -- matches
// stereo2d_ez_reader_integration.rs's negative_17_terminal_alkene fixture.
const PROPENE_V2000 = `propene
  chematic

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    2.3660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
M  END
`;

const EXPECTED_Z_2BUTENE_SMILES = "C(=C/C)/C";

// 1. V2000 positive: exact-match canonical SMILES (not a substring check --
// a substring check would still pass if the marker were on the wrong bond
// or backwards).
{
  const mol = mol_from_sdf_block(Z_2BUTENE_V2000); // V2000-single-block parser, despite its name
  const smi = mol.canonical_smiles();
  assert.equal(smi, EXPECTED_Z_2BUTENE_SMILES, `V2000 canonical SMILES mismatch: ${smi}`);
}

// 2. V3000 positive: same fixture, same expected canonical SMILES --
// cross-dialect semantic consistency (already verified on the Rust side;
// confirming it holds through WASM too).
{
  const mol = mol_from_v3000_block(Z_2BUTENE_V3000);
  const smi = mol.canonical_smiles();
  assert.equal(smi, EXPECTED_Z_2BUTENE_SMILES, `V3000 canonical SMILES mismatch: ${smi}`);
}

// 3. SDF positive: same fixture wrapped as one SDF record.
{
  const sdf = Z_2BUTENE_V2000 + "$$$$\n";
  const records = JSON.parse(sdf_to_records_json(sdf));
  assert.equal(records.length, 1);
  assert.equal(records[0].smiles, EXPECTED_Z_2BUTENE_SMILES, `SDF record smiles mismatch: ${records[0].smiles}`);
}

// 4. Write -> parse -> write stability: bridge through mol_block_from_smiles
// (SMILES -> MOL block) since this WASM build has no direct SMILES ->
// MolHandle entry point, then re-parse via mol_from_sdf_block and confirm
// the canonical SMILES is stable across the round trip. Not full RDKit
// oracle validation (RDKit isn't available inside a WASM/Node test) -- a
// basic self-consistency check only.
{
  const roundTripBlock = mol_block_from_smiles(EXPECTED_Z_2BUTENE_SMILES);
  const roundTripMol = mol_from_sdf_block(roundTripBlock);
  const roundTripSmiles = roundTripMol.canonical_smiles();
  assert.equal(
    roundTripSmiles,
    EXPECTED_Z_2BUTENE_SMILES,
    `write->parse->write round trip unstable: ${roundTripSmiles}`
  );
}

// 5. Negative control: terminal alkene (propene) must produce NO E/Z
// directional token at all.
{
  const mol = mol_from_sdf_block(PROPENE_V2000);
  const smi = mol.canonical_smiles();
  assert.ok(!smi.includes("/") && !smi.includes("\\"), `expected no E/Z token in terminal alkene SMILES: ${smi}`);
}

console.log("stereo2d_ez_directions.test.mjs: all assertions passed");
