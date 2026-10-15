import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";
const url = process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href;
const wasm = await import(url);
const mol = wasm.parse_smiles("CCO");
const coords = [[0,0,0],[3,0,0],[3,4,0]];
assert.deepEqual(JSON.parse(mol.rdkit_distance_matrix_3d_checked_json(JSON.stringify(coords), false)), [[0,3,5],[3,0,4],[5,4,0]]);
assert.equal(JSON.parse(mol.rdkit_distance_matrix_3d_checked_json(JSON.stringify(coords), true))[2][2], 0.75);
const pair = wasm.parse_smiles("CC");
for (const length of [1e308, 1e-160, 1e-300, Number.MIN_VALUE]) {
  const matrix = JSON.parse(pair.rdkit_distance_matrix_3d_checked_json(JSON.stringify([[0,0,0],[length,0,0]]),false));
  assert.deepEqual(matrix, [[0,length],[length,0]]);
}
for (const bad of ["[[NaN,0,0],[0,0,0]]", "[[1e309,0,0],[0,0,0]]", "[]", "[[0,0],[0,0,0]]", " ".repeat(1_000_001), JSON.stringify([[1e308,0,0],[-1e308,0,0]])]) {
  assert.throws(() => pair.rdkit_distance_matrix_3d_checked_json(bad, false));
}
const dummy = wasm.parse_smiles("*");
assert.throws(() => dummy.rdkit_distance_matrix_3d_checked_json("[[0,0,0]]", true));
assert.deepEqual(JSON.parse(dummy.rdkit_distance_matrix_3d_checked_json("[[0,0,0]]", false)), [[0]]);
const explicit = wasm.parse_smiles("[H]C([H])([H])[H]");
assert.deepEqual(JSON.parse(explicit.rdkit_distance_matrix_3d_checked_json("[[0,0,0]]", false)), [[0]]);
const large = wasm.parse_smiles(Array(1025).fill("C").join("."));
assert.throws(() => large.rdkit_distance_matrix_3d_checked_json("[]", false));
for (const m of [mol,pair,dummy,explicit,large]) m.free();
const packet = JSON.parse(readFileSync(new URL("../../../validation/results/rdkit-checked-3d-matrix-2026-10-16.packet.json", import.meta.url), "utf8"));
for (const test of packet.ordinary) {
  const m = wasm.parse_smiles(test.smiles);
  const matrix = JSON.parse(m.rdkit_distance_matrix_3d_checked_json(JSON.stringify(test.coords), test.weights));
  assert.deepEqual(matrix, test.checked);
  assert.deepEqual(matrix, test.rdkit);
  m.free();
}
console.log("Checked 3D matrix: normal, extreme, input, allocation and 128-case cross-runtime contracts passed");
