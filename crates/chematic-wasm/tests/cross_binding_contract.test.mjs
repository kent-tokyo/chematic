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
assert.equal(fixture.operation_manifest.operations.length, 56);
assert.equal(
  new Set(fixture.operation_manifest.operations.map(({ id }) => id)).size,
  56,
);

const orcaOutput = fixture.orca_output_contract;
const orcaOutputJson = JSON.parse(wasm.orca_output_to_json(orcaOutput.input));
assert.equal(orcaOutputJson.charge, orcaOutput.expected.charge);
assert.equal(orcaOutputJson.multiplicity, orcaOutput.expected.multiplicity);
assert.equal(orcaOutputJson.final_energy_hartree, orcaOutput.expected.final_energy_hartree);
assert.equal(orcaOutputJson.termination.kind, orcaOutput.expected.termination);
assert.equal(orcaOutputJson.optimization_convergence, orcaOutput.expected.optimization_convergence);
assert.equal(orcaOutputJson.trajectory.length, orcaOutput.expected.trajectory_count);
assert.throws(() => wasm.orca_output_to_json(orcaOutput.malformed_input));
for (const operation of fixture.operation_manifest.operations) {
  assert.equal(operation.bindings.length, 4, operation.id);
  assert.ok(operation.test_anchors.length > 0, operation.id);
}
assert.equal(fixture.descriptor_contract.schema_version, 1);
assert.equal(fixture.standardization_contract.schema_version, 1);
assert.equal(fixture.descriptor_contract.fields.tpsa.unit, "A2");
assert.equal(fixture.descriptor_contract.fields.exact_mass.unit, "Da");
assert.equal(fixture.descriptor_contract.fields.formula.unit, "Hill notation");
assert.equal(fixture.fingerprint_contract.schema_version, 1);
assert.equal(fixture.fingerprint_contract.operations.ecfp4.bytes, 256);
assert.equal(fixture.fingerprint_contract.operations.maccs.bytes, 21);
assert.equal(fixture.fingerprint_contract.operations.rdkit_rdk.bytes, 256);
assert.equal(fixture.fingerprint_contract.operations.rdkit_path.bytes, 256);
assert.equal(fixture.fingerprint_detail_contract.schema_version, 1);
assert.equal(
  fixture.fingerprint_detail_contract.operations.rdkit_ecfp4_detail.configuration.radius,
  2,
);
assert.equal(fixture.batch_canonicalization_contract.schema_version, 1);
assert.equal(fixture.smiles_validity_contract.schema_version, 1);
assert.equal(fixture.inchi_contract.schema_version, 1);
for (const testCase of fixture.inchi_contract.cases) {
  assert.equal(wasm.inchi_from_smiles(testCase.smiles), testCase.inchi, testCase.smiles);
  assert.equal(wasm.inchikey_from_smiles(testCase.smiles), testCase.inchikey, testCase.smiles);
}
assert.equal(fixture.fragment_parent_contract.schema_version, 1);
for (const testCase of fixture.fragment_parent_contract.cases) {
  const mol = wasm.parse_smiles(testCase.smiles);
  const result = JSON.parse(wasm.fragment_parent_json(mol));
  assert.equal(result.status, "completed", testCase.smiles);
  assert.equal(result.smiles, testCase.canonical_smiles, testCase.smiles);
  mol.free();
}
for (const [key, transform] of [
  ["charge_parent_contract", wasm.charge_parent_json],
  ["isotope_parent_contract", wasm.isotope_parent_json],
  ["stereo_parent_contract", wasm.stereo_parent_json],
]) {
  assert.equal(fixture[key].schema_version, 1);
  for (const testCase of fixture[key].cases) {
    const mol = wasm.parse_smiles(testCase.smiles);
    const result = JSON.parse(transform(mol));
    assert.equal(result.status, "completed", testCase.smiles);
    assert.equal(result.smiles, testCase.canonical_smiles, `${key}: ${testCase.smiles}`);
    mol.free();
  }
}
assert.equal(fixture.hydrogen_contract.schema_version, 1);
for (const testCase of fixture.hydrogen_contract.cases) {
  const mol = wasm.parse_smiles(testCase.smiles);
  const added = wasm.add_hydrogens(mol);
  assert.equal(added.canonical_smiles(), testCase.added_canonical_smiles, `add: ${testCase.smiles}`);
  const removed = wasm.remove_hydrogens(added);
  assert.equal(removed.canonical_smiles(), testCase.removed_canonical_smiles, `remove: ${testCase.smiles}`);
  mol.free();
  added.free();
  removed.free();
}
assert.equal(fixture.super_parent_report_contract.schema_version, 1);
for (const testCase of fixture.super_parent_report_contract.cases) {
  const mol = wasm.parse_smiles(testCase.smiles);
  const limits = fixture.super_parent_report_contract.limits;
  const report = JSON.parse(wasm.super_parent_report_json(
    mol, limits.max_transforms, limits.max_tautomers, null,
  ));
  assert.equal(report.status, "completed", testCase.smiles);
  assert.equal(report.smiles, testCase.canonical_smiles, testCase.smiles);
  assert.deepEqual(
    report.stages.map(({ name }) => name),
    fixture.super_parent_report_contract.stage_names,
  );
  assert.ok(report.stages.every(({ smiles }) => smiles.length > 0));
  mol.free();
}
assert.equal(fixture.canonical_tautomer_contract.schema_version, 1);
for (const testCase of fixture.canonical_tautomer_contract.cases) {
  const mol = wasm.parse_smiles(testCase.smiles);
  const tautomer = wasm.canonical_tautomer(mol);
  assert.equal(tautomer.canonical_smiles(), testCase.canonical_smiles, testCase.smiles);
  mol.free();
  tautomer.free();
}
assert.equal(fixture.neutralize_charges_contract.schema_version, 1);
for (const testCase of fixture.neutralize_charges_contract.cases) {
  const mol = wasm.parse_smiles(testCase.smiles);
  const neutral = wasm.neutralize_charges(mol);
  assert.equal(neutral.canonical_smiles(), testCase.canonical_smiles, testCase.smiles);
  mol.free();
  neutral.free();
}
assert.equal(fixture.largest_fragment_contract.schema_version, 1);
for (const testCase of fixture.largest_fragment_contract.cases) {
  const mol = wasm.parse_smiles(testCase.smiles);
  const fragment = wasm.largest_fragment(mol);
  assert.equal(fragment.canonical_smiles(), testCase.canonical_smiles, testCase.smiles);
  mol.free();
  fragment.free();
}
for (const [key, transform] of [
  ["tautomer_parent_contract", wasm.tautomer_parent_json],
  ["super_parent_contract", wasm.super_parent_json],
]) {
  assert.equal(fixture[key].schema_version, 1);
  const limits = fixture[key].limits;
  for (const testCase of fixture[key].cases) {
    const mol = wasm.parse_smiles(testCase.smiles);
    const result = JSON.parse(transform(mol, limits.max_transforms, limits.max_tautomers, null));
    assert.equal(result.status, "completed", `${key}: ${testCase.smiles}`);
    assert.equal(result.smiles, testCase.canonical_smiles, `${key}: ${testCase.smiles}`);
    mol.free();
  }
}
for (const smiles of fixture.smiles_validity_contract.accepted) {
  assert.equal(wasm.is_valid_smiles(smiles), true, `accepted: ${smiles}`);
}
for (const smiles of fixture.smiles_validity_contract.rejected) {
  assert.equal(wasm.is_valid_smiles(smiles), false, `rejected: ${smiles}`);
}
assert.equal(fixture.mol_v2000_contract.schema_version, 1);
assert.equal(fixture.mol_v3000_contract.schema_version, 1);
assert.equal(fixture.xyz_contract.schema_version, 1);
assert.equal(fixture.extxyz_contract.schema_version, 1);
assert.equal(fixture.pdb_contract.schema_version, 1);
assert.equal(fixture.cml_contract.schema_version, 1);
assert.equal(fixture.cdxml_contract.schema_version, 1);
assert.equal(fixture.mmcif_contract.schema_version, 1);
assert.equal(fixture.moljson_contract.schema_version, 1);
assert.equal(fixture.xyz_batch_contract.schema_version, 1);

for (const [key, parser] of [
  ["mol_v2000_contract", wasm.mol_from_sdf_block],
  ["mol_v3000_contract", wasm.mol_from_v3000_block],
  ["xyz_contract", wasm.mol_from_xyz],
]) {
  const contract = fixture[key];
  const mol = parser(contract.input);
  assert.equal(mol.atom_count(), contract.expected.atom_count, key);
  assert.equal(mol.canonical_smiles(), contract.expected.canonical_smiles, key);
  if (key === "mol_v2000_contract") {
    const roundtripped = wasm.mol_from_sdf_block(wasm.to_mol_block(mol));
    assert.equal(roundtripped.atom_count(), contract.expected.atom_count, `${key} roundtrip`);
    assert.equal(roundtripped.canonical_smiles(), contract.expected.canonical_smiles, `${key} roundtrip`);
    roundtripped.free();
  }
  if (key === "mol_v3000_contract") {
    const roundtripped = wasm.mol_from_v3000_block(wasm.to_mol_v3000_block(mol));
    assert.equal(roundtripped.atom_count(), contract.expected.atom_count, `${key} roundtrip`);
    assert.equal(roundtripped.canonical_smiles(), contract.expected.canonical_smiles, `${key} roundtrip`);
    roundtripped.free();
  }
  mol.free();
}

const xyzRoundtrip = fixture.xyz_contract;
const xyzRoundtripSource = wasm.mol_from_xyz(xyzRoundtrip.input);
const xyzRoundtripMol = wasm.mol_from_xyz(wasm.to_xyz(xyzRoundtripSource));
assert.equal(xyzRoundtripMol.atom_count(), xyzRoundtrip.expected.atom_count);
assert.equal(xyzRoundtripMol.canonical_smiles(), xyzRoundtrip.expected.canonical_smiles);
xyzRoundtripMol.free();
xyzRoundtripSource.free();

const extxyz = fixture.extxyz_contract;
const extxyzActual = JSON.parse(wasm.extxyz_frame_json(extxyz.input));
assert.deepEqual(extxyzActual.coords, extxyz.expected.coords);
assert.deepEqual(extxyzActual.lattice, extxyz.expected.lattice);
assert.deepEqual(extxyzActual.properties, extxyz.expected.properties);
assert.deepEqual(extxyzActual.info, extxyz.expected.info);

const extxyzWriter = extxyz.writer;
const extxyzWriterMol = wasm.parse_smiles(extxyzWriter.smiles);
const extxyzWriterOptions = JSON.stringify({
  lattice: extxyzWriter.lattice,
  properties: extxyzWriter.properties,
  info: extxyzWriter.info,
});
const extxyzWritten = wasm.to_extxyz_json(
  extxyzWriterMol,
  JSON.stringify(extxyzWriter.coords),
  extxyzWriterOptions,
);
const extxyzWrittenActual = JSON.parse(wasm.extxyz_frame_json(extxyzWritten));
assert.deepEqual(extxyzWrittenActual.coords, extxyzWriter.expected.coords);
assert.deepEqual(extxyzWrittenActual.lattice, extxyzWriter.expected.lattice);
assert.deepEqual(extxyzWrittenActual.properties, extxyzWriter.expected.properties);
assert.deepEqual(extxyzWrittenActual.info, extxyzWriter.expected.info);
extxyzWriterMol.free();

const pdb = fixture.pdb_contract;
const pdbMol = wasm.mol_from_pdb(pdb.input);
assert.equal(pdbMol.atom_count(), pdb.expected.atom_count);
assert.deepEqual(JSON.parse(wasm.pdb_coords_json(pdb.input)), pdb.expected.coords);
pdbMol.free();

const qcschema = fixture.qcschema_contract;
const qcschemaMol = wasm.mol_from_qcschema_molecule(qcschema.input);
assert.equal(qcschemaMol.atom_count(), qcschema.expected.atom_count);
assert.deepEqual(
  JSON.parse(wasm.qcschema_molecule_coords_json(qcschema.input)).coords,
  qcschema.expected.coords_angstrom,
);
qcschemaMol.free();

const atomicInput = fixture.qcschema_atomic_input_contract;
const normalizedAtomicInput = JSON.parse(
  wasm.qcschema_validate_atomic_input(atomicInput.input),
);
assert.equal(normalizedAtomicInput.driver, atomicInput.expected.driver);
assert.equal(normalizedAtomicInput.model.method, atomicInput.expected.model_method);
assert.equal(normalizedAtomicInput.molecule.symbols.length, atomicInput.expected.atom_count);
assert.equal(normalizedAtomicInput.extras.contract_marker, atomicInput.expected.extra_marker);
assert.equal(normalizedAtomicInput.vendor_extension, atomicInput.expected.vendor_extension);

const atomicResult = fixture.qcschema_atomic_result_contract;
const normalizedAtomicResult = JSON.parse(
  wasm.qcschema_validate_atomic_result(atomicResult.input),
);
assert.equal(normalizedAtomicResult.driver, atomicResult.expected.driver);
assert.equal(normalizedAtomicResult.return_result, atomicResult.expected.return_result);
assert.equal(normalizedAtomicResult.success, atomicResult.expected.success);
assert.equal(normalizedAtomicResult.molecule.symbols.length, atomicResult.expected.atom_count);
assert.equal(normalizedAtomicResult.vendor_result_marker, atomicResult.expected.vendor_result_marker);

const orca = fixture.orca_input_contract;
const orcaJson = JSON.parse(wasm.orca_input_to_json(orca.input));
assert.deepEqual(orcaJson.keywords, orca.expected.keywords);
assert.equal(orcaJson.coords.type, orca.expected.kind);
assert.equal(orcaJson.coords.charge, orca.expected.charge);
assert.equal(orcaJson.coords.multiplicity, orca.expected.multiplicity);
assert.equal(orcaJson.coords.atoms.length, orca.expected.atom_count);
assert.equal(orcaJson.coords.atoms[0].element, orca.expected.element);
assert.deepEqual(
  JSON.parse(wasm.orca_input_coords_json(orca.input)).coords,
  orca.expected.coords,
);

const mol2 = fixture.mol2_contract;
const mol2Smiles = wasm.mol2_to_smiles(mol2.input);
assert.equal(typeof mol2Smiles, "string");
assert.notEqual(mol2Smiles, "");
const mol2Roundtrip = wasm.mol2_to_smiles(wasm.smiles_to_mol2(mol2Smiles));
assert.equal(mol2Roundtrip, mol2Smiles);

const cml = fixture.cml_contract;
const cmlMol = wasm.mol_from_cml(cml.input);
assert.equal(cmlMol.atom_count(), cml.expected.atom_count);
const cmlRoundtrip = wasm.mol_from_cml(wasm.to_cml(cmlMol));
assert.equal(cmlRoundtrip.atom_count(), cml.expected.atom_count);
cmlRoundtrip.free();
cmlMol.free();

const cdxml = fixture.cdxml_contract;
const cdxmlMol = wasm.mol_from_cdxml(cdxml.input);
assert.equal(cdxmlMol.atom_count(), cdxml.expected.atom_count);
const cdxmlRoundtrip = wasm.convert_common_format(cdxml.input, "cdxml", "cdxml");
const cdxmlRoundtripMol = wasm.mol_from_cdxml(cdxmlRoundtrip);
assert.equal(cdxmlRoundtripMol.atom_count(), cdxml.expected.atom_count);
cdxmlRoundtripMol.free();
cdxmlMol.free();

const mmcif = fixture.mmcif_contract;
const mmcifMol = wasm.mol_from_mmcif(mmcif.input);
assert.equal(mmcifMol.atom_count(), mmcif.expected.atom_count);
assert.deepEqual(JSON.parse(wasm.mmcif_coords_json(mmcif.input)), mmcif.expected.coords);
mmcifMol.free();

const moljson = fixture.moljson_contract;
const moljsonMol = wasm.mol_from_moljson(moljson.input);
assert.equal(moljsonMol.atom_count(), moljson.expected.atom_count);
assert.equal(moljsonMol.canonical_smiles(), moljson.expected.canonical_smiles);
const moljsonRoundtrip = wasm.mol_from_moljson(wasm.to_moljson(moljsonMol));
assert.equal(moljsonRoundtrip.atom_count(), moljson.expected.atom_count);
assert.equal(moljsonRoundtrip.canonical_smiles(), moljson.expected.canonical_smiles);
moljsonRoundtrip.free();
moljsonMol.free();

for (const format of ["xyz", "extxyz"]) {
  const contract = fixture.xyz_batch_contract[format];
  const batch = format === "xyz"
    ? wasm.xyz_frames_batch_json(contract.input, 0, fixture.xyz_batch_contract.batch_size)
    : wasm.extxyz_frames_batch_json(contract.input, 0, fixture.xyz_batch_contract.batch_size);
  const actual = JSON.parse(batch);
  assert.equal(actual.format, format);
  assert.equal(actual.status, contract.expected.status);
  assert.equal(actual.record_count, contract.expected.record_count);
  assert.equal(actual.rejected_count, contract.expected.rejected_count);
  assert.equal(actual.records[0].status, "rejected");
  assert.equal(actual.records[1].status, "accepted");
  assert.deepEqual(actual.records[1].frame.coords, contract.expected.records[1].coords);
}

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

const sdfInput = readFileSync(path.join(repoRoot, "benchmarks/fixtures/streaming.sdf"), "utf8");
const sdfPage = JSON.parse(wasm.sdf_records_batch_json(sdfInput, 0, 1));
assert.equal(sdfPage.status, "partial");
assert.equal(sdfPage.record_count, 1);
assert.equal(sdfPage.next_offset, 1);
const sdfNextPage = JSON.parse(wasm.sdf_records_batch_json(sdfInput, sdfPage.next_offset, 1));
assert.equal(sdfNextPage.offset, 1);
assert.equal(sdfNextPage.records[0].input_index, 1);
assert.throws(() => wasm.sdf_records_batch_json(sdfInput, 0, 0));
assert.throws(() => wasm.sdf_records_batch_json("x".repeat(1_000_001), 0, 1));

const xyzInput = "1\nfirst\nC 0 0 0\n1\nsecond\nO 1 0 0\n";
const xyzPage = JSON.parse(wasm.xyz_frames_batch_json(xyzInput, 0, 1));
assert.equal(xyzPage.status, "partial");
assert.equal(xyzPage.records[0].input_index, 0);
const xyzNextPage = JSON.parse(wasm.xyz_frames_batch_json(xyzInput, xyzPage.next_offset, 1));
assert.equal(xyzNextPage.records[0].input_index, 1);
assert.throws(() => wasm.xyz_frames_batch_json(xyzInput, 0, 0));

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
  assert.ok(Math.abs(mol.exact_mass() - expected.exact_mass) <= 1e-6, expected.id);
  assert.equal(mol.formula(), expected.formula, expected.id);
  assert.ok(Math.abs(mol.tpsa() - expected.tpsa) <= 1e-6, expected.id);
  assert.equal(mol.hbd_count(), expected.hbd, expected.id);
  assert.equal(mol.hba_count(), expected.hba, expected.id);
  assert.equal(mol.aromatic_ring_count(), expected.aromatic_ring_count, expected.id);
  assert.equal(mol.rotatable_bond_count(), expected.rotatable_bonds, expected.id);
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
  const torsion = wasm.torsion_bitvec(mol);
  const rdkitTorsion = wasm.rdkit_torsion_bitvec(mol);
  const maccs = wasm.maccs_bitvec(mol);
  assert.equal(ecfp4.length, 256, `${expected.id} ECFP4 shape`);
  assert.equal(topoPath.length, 256, `${expected.id} topo_path shape`);
  assert.equal(torsion.length, 256, `${expected.id} torsion shape`);
  assert.equal(rdkitTorsion.length, 256, `${expected.id} RDKit torsion shape`);
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
  assert.deepEqual(
    Array.from({ length: 2048 }, (_, bit) => (torsion[bit >> 3] >> (bit & 7)) & 1 ? bit : null).filter((bit) => bit !== null),
    expected.torsion_bits,
    `${expected.id} torsion bits`,
  );
  assert.deepEqual(
    Array.from({ length: 2048 }, (_, bit) => (rdkitTorsion[bit >> 3] >> (bit & 7)) & 1 ? bit : null).filter((bit) => bit !== null),
    expected.rdkit_torsion_bits,
    `${expected.id} RDKit torsion bits`,
  );
  assert.equal(Buffer.from(maccs).toString("hex"), expected.maccs_hex, `${expected.id} MACCS bytes`);
  assert.ok(ecfp4.some((byte) => byte !== 0), `${expected.id} ECFP4 non-empty`);
  assert.ok(maccs.some((byte) => byte !== 0), `${expected.id} MACCS non-empty`);
  mol.free();
}

for (const expected of fixture.fingerprint_contract.rdkit_rdk_fixtures) {
  const mol = wasm.parse_smiles(expected.smiles);
  const actual = wasm.rdkit_rdk_bitvec(mol);
  assert.equal(actual.length, 256, `${expected.id} RDKit RDK shape`);
  assert.deepEqual(
    Array.from({ length: 2048 }, (_, bit) => (actual[bit >> 3] >> (bit & 7)) & 1 ? bit : null).filter((bit) => bit !== null),
    expected.rdkit_rdk_bits,
    `${expected.id} RDKit RDK bits`,
  );
  mol.free();
}

for (const expected of fixture.fingerprint_contract.rdkit_path_fixtures) {
  const mol = wasm.parse_smiles(expected.smiles);
  const actual = wasm.rdkit_path_bitvec(mol);
  assert.equal(actual.length, 256, `${expected.id} RDKit path shape`);
  assert.deepEqual(
    Array.from({ length: 2048 }, (_, bit) => (actual[bit >> 3] >> (bit & 7)) & 1 ? bit : null).filter((bit) => bit !== null),
    expected.rdkit_path_bits,
    `${expected.id} RDKit path bits`,
  );
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
