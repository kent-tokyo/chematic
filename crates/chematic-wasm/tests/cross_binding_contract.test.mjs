// Consume the same fixture as the Rust and Python binding contract tests.

import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");
const fixture = JSON.parse(
  readFileSync(path.join(repoRoot, "validation/cross_binding_contract.json"), "utf8"),
);
const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));

assert.equal(fixture.schema_version, 1);
assert.equal(fixture.fixtures.length, 4);
assert.equal(fixture.descriptor_contract.schema_version, 1);
assert.equal(fixture.descriptor_contract.fields.tpsa.unit, "A2");

for (const expected of fixture.fixtures) {
  const mol = wasm.parse_smiles(expected.smiles);
  assert.equal(mol.canonical_smiles(), expected.canonical_smiles, expected.id);
  assert.equal(mol.atom_count(), expected.heavy_atoms, expected.id);
  mol.free();
}

for (const expected of fixture.descriptor_contract.fixtures) {
  const mol = wasm.parse_smiles(expected.smiles);
  assert.ok(Math.abs(mol.molecular_weight() - expected.molecular_weight) <= 1e-6, expected.id);
  assert.ok(Math.abs(mol.tpsa() - expected.tpsa) <= 1e-6, expected.id);
  assert.equal(mol.hbd_count(), expected.hbd, expected.id);
  assert.equal(mol.hba_count(), expected.hba, expected.id);
  assert.equal(mol.heavy_atom_count(), expected.heavy_atoms, expected.id);
  mol.free();
}

for (const expected of fixture.adversarial) {
  assert.throws(
    () => wasm.convert_common_format(expected.input, expected.format, "smiles"),
    undefined,
    expected.id,
  );
}
