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
assert.equal(fixture.standardization_contract.schema_version, 1);
assert.equal(fixture.descriptor_contract.fields.tpsa.unit, "A2");
assert.equal(fixture.fingerprint_contract.schema_version, 1);
assert.equal(fixture.fingerprint_contract.operations.ecfp4.bytes, 256);
assert.equal(fixture.fingerprint_contract.operations.maccs.bytes, 21);
assert.equal(fixture.fingerprint_detail_contract.schema_version, 1);
assert.equal(
  fixture.fingerprint_detail_contract.operations.rdkit_ecfp4_detail.configuration.radius,
  2,
);
assert.equal(fixture.batch_canonicalization_contract.schema_version, 1);
assert.equal(fixture.extxyz_contract.schema_version, 1);

const extxyz = fixture.extxyz_contract;
const extxyzActual = JSON.parse(wasm.extxyz_frame_json(extxyz.input));
assert.deepEqual(extxyzActual.coords, extxyz.expected.coords);
assert.deepEqual(extxyzActual.lattice, extxyz.expected.lattice);
assert.deepEqual(extxyzActual.properties, extxyz.expected.properties);
assert.deepEqual(extxyzActual.info, extxyz.expected.info);

const batchInputs = fixture.batch_canonicalization_contract.inputs.join("\n");
const batchManifest = JSON.parse(wasm.canonicalize_smiles_batch_json(batchInputs, "\n"));
assert.equal(batchManifest.schema_version, 1);
assert.equal(batchManifest.operation, "canonicalize_smiles");
assert.equal(batchManifest.status, "complete");
assert.equal(batchManifest.record_count, fixture.batch_canonicalization_contract.expected.length);
const batchActual = batchManifest.records;
assert.deepEqual(
  batchActual.map(({ error, ...record }) => record),
  fixture.batch_canonicalization_contract.expected,
);

// The batch boundary is deliberately a stable JSON manifest: malformed
// records stay inline, later records are processed, and input/shape limits
// reject before parsing rather than returning an ambiguous partial result.
const malformedBatch = JSON.parse(
  wasm.canonicalize_smiles_batch_json("CC\nC1CC\nCCO", "\n"),
);
assert.equal(malformedBatch.status, "complete");
assert.equal(malformedBatch.record_count, 3);
assert.deepEqual(
  malformedBatch.records.map(({ input_index: inputIndex, status }) => ({ inputIndex, status })),
  [
    { inputIndex: 0, status: "accepted" },
    { inputIndex: 1, status: "rejected" },
    { inputIndex: 2, status: "accepted" },
  ],
);
assert.equal(typeof malformedBatch.records[1].error, "string");
assert.throws(() => wasm.canonicalize_smiles_batch_json("CC", ""));
assert.throws(() => wasm.canonicalize_smiles_batch_json("CC\n".repeat(1024), "\n"));
assert.throws(() => wasm.canonicalize_smiles_batch_json("C".repeat(1_000_001), "\n"));

const screeningBatch = JSON.parse(wasm.screen_smiles_json("CC\nC1CC\nCCO", "\n"));
assert.equal(screeningBatch.records.length, 3);
assert.equal(screeningBatch.records[0].error, null);
assert.equal(screeningBatch.records[1].report, null);
assert.match(screeningBatch.records[1].error, /SMILES parse failed/);
assert.equal(screeningBatch.records[2].error, null);
assert.deepEqual(JSON.parse(wasm.screen_smiles_json("CC", "")), {
  error: "delimiter must not be empty",
});
assert.deepEqual(JSON.parse(wasm.screen_smiles_json("CC\n".repeat(1024), "\n")), {
  error: "smiles_batch exceeds maximum item count (1025 > 1024)",
});

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

for (const expected of fixture.standardization_contract.fixtures) {
  const actual = JSON.parse(
    wasm.standardize_smiles_report_json(expected.smiles, true, true, true, true),
  );
  assert.equal(actual.smiles, expected.output_smiles, expected.id);
}

for (const expected of fixture.fingerprint_contract.fixtures) {
  const mol = wasm.parse_smiles(expected.smiles);
  const ecfp4 = wasm.ecfp4_bitvec(mol);
  const topoPath = wasm.topo_path_bitvec(mol);
  const maccs = wasm.maccs_bitvec(mol);
  assert.equal(ecfp4.length, 256, `${expected.id} ECFP4 shape`);
  assert.equal(topoPath.length, 256, `${expected.id} topo_path shape`);
  assert.equal(maccs.length, 21, `${expected.id} MACCS shape`);
  assert.deepEqual(
    Array.from({ length: 2048 }, (_, bit) => (ecfp4[bit >> 3] >> (bit & 7)) & 1 ? bit : null).filter((bit) => bit !== null),
    expected.ecfp4_bits,
    `${expected.id} ECFP4 bits`,
  );
  assert.deepEqual(
    Array.from({ length: 2048 }, (_, bit) => (topoPath[bit >> 3] >> (bit & 7)) & 1 ? bit : null).filter((bit) => bit !== null),
    expected.topo_path_bits,
    `${expected.id} topo_path bits`,
  );
  assert.equal(Buffer.from(maccs).toString("hex"), expected.maccs_hex, `${expected.id} MACCS bytes`);
  assert.ok(ecfp4.some((byte) => byte !== 0), `${expected.id} ECFP4 non-empty`);
  assert.ok(maccs.some((byte) => byte !== 0), `${expected.id} MACCS non-empty`);
  mol.free();
}

for (const expected of fixture.fingerprint_detail_contract.fixtures) {
  const mol = wasm.parse_smiles(expected.smiles);
  const detail = JSON.parse(wasm.rdkit_ecfp4_detail_json(mol));
  assert.equal(detail.fingerprint.length, 256, `${expected.id} detail shape`);
  assert.ok(Object.keys(detail.sparseCounts).length > 0, `${expected.id} sparse counts`);
  assert.ok(Object.keys(detail.rawBitInfo).length > 0, `${expected.id} raw explanation`);
  assert.ok(Object.keys(detail.foldedBitInfo).length > 0, `${expected.id} folded explanation`);
  for (const provenance of [detail.rawBitInfo, detail.foldedBitInfo]) {
    for (const pairs of Object.values(provenance)) {
      for (const [atom, radius] of pairs) {
        assert.ok(atom >= 0 && atom < mol.atom_count(), `${expected.id} atom provenance`);
        assert.ok(radius >= 0 && radius <= 2, `${expected.id} radius provenance`);
      }
    }
  }
  mol.free();
}

for (const expected of fixture.adversarial) {
  assert.throws(
    () => wasm.convert_common_format(expected.input, expected.format, "smiles"),
    undefined,
    expected.id,
  );
}
