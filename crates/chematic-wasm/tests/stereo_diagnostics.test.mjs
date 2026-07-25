// Tests for 2D wedge/hash stereo perception + diagnostics wired into the
// WASM MOL/SDF bindings (mol_from_sdf_block, mol_from_v3000_block,
// mol_block_stereo_diagnostics_json, mol_v3000_stereo_diagnostics_json,
// sdf_to_records_json). See docs/stereo2d_reader_integration_rfc.md and
// crates/chematic-mol/tests/stereo_reader_integration.rs for the Rust-side
// equivalents of these same fixtures.
//
// Not wired into any CI workflow (this repo has no WASM-test CI job). Run
// manually after building the Node target:
//
//   wasm-pack build crates/chematic-wasm --target nodejs --out-dir pkg-node --release
//   node crates/chematic-wasm/tests/stereo_diagnostics.test.mjs

import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");

const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const {
  mol_from_sdf_block,
  mol_from_v3000_block,
  mol_block_stereo_diagnostics_json,
  mol_v3000_stereo_diagnostics_json,
  sdf_to_records_json,
} = wasm;

function atomLine(x, y, sym) {
  return (
    x.toFixed(4).padStart(10) +
    y.toFixed(4).padStart(10) +
    (0.0).toFixed(4).padStart(10) +
    " " +
    sym.padEnd(3) +
    " 0  0  0  0  0  0  0  0  0  0  0  0"
  );
}

function bondLine(a1, a2, btype, stereo) {
  return String(a1).padStart(3) + String(a2).padStart(3) + String(btype).padStart(3) + String(stereo).padStart(3);
}

function chfclbrV2000(wedgeBonds) {
  const lines = [
    "test",
    "  chematic",
    "",
    "  5  4  0  0  0  0  0  0  0  0999 V2000",
    atomLine(0.0, 0.0, "C"),
    atomLine(-1.0, 0.4, "F"),
    atomLine(0.9, 0.7, "Cl"),
    atomLine(-0.5, -1.1, "Br"),
    atomLine(0.8, -0.6, "I"),
  ];
  for (const sub of [2, 3, 4, 5]) {
    lines.push(bondLine(1, sub, 1, wedgeBonds[sub] ?? 0));
  }
  lines.push("M  END");
  lines.push("");
  return lines.join("\n");
}

const VALID_WEDGE_BLOCK = chfclbrV2000({ 2: 1 });
const CONTRADICTORY_WEDGE_BLOCK = chfclbrV2000({ 2: 1, 3: 1 });

const V3000_WEDGE_BLOCK = `wedge
  chematic

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 5 4 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0 0.0 0 0
M  V30 2 F -1.0 0.4 0 0
M  V30 3 Cl 0.9 0.7 0 0
M  V30 4 Br -0.5 -1.1 0 0
M  V30 5 I 0.8 -0.6 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2 CFG=1
M  V30 2 1 1 3
M  V30 3 1 1 4
M  V30 4 1 1 5
M  V30 END BOND
M  V30 END CTAB
M  END
`;

// mol_from_sdf_block automatically perceives a valid wedge -> chirality
// survives to canonical SMILES as `@`, with no diagnostics.
{
  const mol = mol_from_sdf_block(VALID_WEDGE_BLOCK);
  const smi = mol.canonical_smiles();
  assert.ok(smi.includes("@"), `expected chirality in SMILES: ${smi}`);
  const diag = JSON.parse(mol_block_stereo_diagnostics_json(VALID_WEDGE_BLOCK));
  assert.deepEqual(diag, []);
}

// Contradictory wedges: no chirality, one diagnostic.
{
  const mol = mol_from_sdf_block(CONTRADICTORY_WEDGE_BLOCK);
  const smi = mol.canonical_smiles();
  assert.ok(!smi.includes("@"), `expected no chirality in SMILES: ${smi}`);
  const diag = JSON.parse(mol_block_stereo_diagnostics_json(CONTRADICTORY_WEDGE_BLOCK));
  assert.deepEqual(diag, [{ atom_idx: 0, reason: "contradictory_wedges" }]);
}

// V3000 wedge matches the V2000 reading of the same drawing.
{
  const molV3000 = mol_from_v3000_block(V3000_WEDGE_BLOCK);
  const molV2000 = mol_from_sdf_block(VALID_WEDGE_BLOCK);
  assert.equal(molV3000.canonical_smiles(), molV2000.canonical_smiles());
  const diag = JSON.parse(mol_v3000_stereo_diagnostics_json(V3000_WEDGE_BLOCK));
  assert.deepEqual(diag, []);
}

// SDF records JSON carries the same diagnostics as the direct parse.
{
  const sdf = CONTRADICTORY_WEDGE_BLOCK + "$$$$\n";
  const records = JSON.parse(sdf_to_records_json(sdf));
  assert.equal(records.length, 1);
  assert.deepEqual(records[0].stereo_diagnostics, [
    { atom_idx: 0, reason: "contradictory_wedges" },
  ]);
}

console.log("stereo_diagnostics.test.mjs: all assertions passed");
