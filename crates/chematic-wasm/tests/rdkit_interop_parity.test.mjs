// RDKit-compatible WASM binding parity (issue #784).
//
// Build and run with:
//   wasm-pack build crates/chematic-wasm --target nodejs --out-dir pkg-node --release
//   node crates/chematic-wasm/tests/rdkit_interop_parity.test.mjs

import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");
const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));

const mol = wasm.parse_smiles("O=C1CC[C@H](C)N1CCc1ccccc1");

// Writers and structured chemistry use zero-based MolHandle atom indices.
assert.match(mol.rdkit_smarts(true), /^\[/);
assert.match(mol.rdkit_cx_smarts(), /^\[/);
assert.match(mol.rdkit_pdb_block(), /HETATM/);
assert.equal(mol.rdkit_murcko_scaffold(), "O=C1CCCN1CCc1ccccc1");
assert.deepEqual(JSON.parse(mol.rdkit_chiral_centers_json(true)), [[4, "S"]]);

const stereo = wasm.parse_smiles("CC(F)C(Cl)Br");
assert.equal(stereo.rdkit_stereoisomer_count(), "4");
assert.equal(JSON.parse(stereo.rdkit_stereoisomer_smiles_json(16)).length, 4);
assert.equal(
  wasm.parse_smiles("Cc1ccccc1CC(=O)O").rdkit_mol_hash("ExtendedMurcko", false),
  "*c1ccccc1*",
);
const bitInfo = JSON.parse(mol.rdkit_morgan_bit_info_json(2, 2048, false));
assert.ok(Object.keys(bitInfo).length > 0);

// PDB, XYZ and MOL2 readers return a handle and coordinates in the same order.
const pdbResult = wasm.rdkit_from_pdb_block(mol.rdkit_pdb_block(), true, true, 0, true);
assert.equal(pdbResult.molecule().atom_count(), mol.atom_count());
assert.equal(JSON.parse(pdbResult.coords_json()).length, mol.atom_count());

const xyzResult = wasm.rdkit_from_xyz_block("2\nethane\nC 0 0 0\nC 1.5 0 0\n");
assert.equal(xyzResult.molecule().atom_count(), 2);
assert.deepEqual(JSON.parse(xyzResult.coords_json())[1], [1.5, 0, 0]);

const mol2 = wasm.smiles_to_mol2("CCO");
const mol2Result = wasm.rdkit_from_mol2_block(mol2, true, true, true);
assert.equal(mol2Result.molecule().atom_count(), 3);
assert.equal(JSON.parse(mol2Result.coords_json()).length, 3);

// Alignment calls are explicitly bounded and return transforms/maps as JSON.
const ethanol = wasm.parse_smiles("CCO");
const reference = JSON.stringify([[0, 0, 0], [1, 0, 0], [2, 0, 0]]);
const probe = JSON.stringify([[1, 0, 0], [2, 0, 0], [3, 0, 0]]);
const aligned = JSON.parse(ethanol.rdkit_align_json(probe, reference, undefined, undefined, false, 50));
assert.ok(aligned.rmsd < 1e-12);
assert.equal(aligned.atom_map.length, 3);
assert.ok(ethanol.rdkit_best_rms(probe, reference, 100, true) < 1e-12);
assert.ok(ethanol.rdkit_calc_rms(reference, reference, 100, true) < 1e-12);
const best = JSON.parse(ethanol.rdkit_best_alignment_json(probe, reference, 100, true));
assert.equal(best.transform.length, 4);
assert.equal(best.atom_map.length, 3);

// Seeded embedding and bounds are separate from the generic 3D browser API.
const ethanolH = wasm.add_hydrogens(ethanol);
const embedded = JSON.parse(ethanolH.rdkit_embed_json(42, 0));
const bounds = JSON.parse(ethanolH.rdkit_bounds_matrix_json(true, true, false));
assert.equal(embedded.length, ethanolH.atom_count());
assert.equal(bounds.length, ethanolH.atom_count());
assert.ok(bounds.every((row) => row.length === ethanolH.atom_count()));

// Every argument-bearing operation rejects malformed or unbounded inputs by
// throwing a stable JS error, never by panicking or silently coercing values.
assert.throws(() => mol.rdkit_smarts(true, 999), /rdkit_compatibility_error/);
assert.throws(() => mol.rdkit_pdb_block("[]"), /rdkit_coordinate_count/);
assert.throws(() => stereo.rdkit_stereoisomer_smiles_json(0), /rdkit_stereoisomer_limit/);
assert.throws(() => mol.rdkit_mol_hash("NoSuchHash", false), /rdkit_hash_function/);
assert.throws(() => mol.rdkit_morgan_bit_info_json(9, 2048, false), /rdkit_morgan_config/);
assert.throws(
  () => ethanol.rdkit_align_json("[]", reference, undefined, undefined, false, 50),
  /rdkit_coordinate_count/,
);
assert.throws(() => ethanol.rdkit_best_rms(probe, reference, 0, true), /rdkit_alignment_limit/);
assert.throws(
  () => ethanol.rdkit_best_alignment_json(probe, reference, 1_000_001, true),
  /rdkit_alignment_limit/,
);
assert.throws(() => ethanol.rdkit_calc_rms(probe, reference, 0, true), /rdkit_alignment_limit/);
assert.throws(() => ethanolH.rdkit_embed_json(-1, 0), /rdkit_compatibility_error/);
assert.throws(() => wasm.rdkit_from_pdb_block("not pdb", true, true, 0, true), /rdkit_parse_failed|rdkit_compatibility_error/);
assert.throws(() => wasm.rdkit_from_xyz_block("not xyz"), /rdkit_compatibility_error/);
assert.throws(() => wasm.rdkit_from_mol2_block("not mol2", true, true, true), /rdkit_compatibility_error/);

// Bounds output is O(n^2), so the WASM contract refuses more than 512 atoms.
const oversized = wasm.parse_smiles("C".repeat(513));
assert.throws(() => oversized.rdkit_bounds_matrix_json(true, true, false), /rdkit_3d_limit/);

oversized.free();
ethanolH.free();
ethanol.free();
mol2Result.molecule().free();
mol2Result.free();
xyzResult.molecule().free();
xyzResult.free();
pdbResult.molecule().free();
pdbResult.free();
stereo.free();
mol.free();

console.log("RDKit-compatible WASM parity bindings passed");
