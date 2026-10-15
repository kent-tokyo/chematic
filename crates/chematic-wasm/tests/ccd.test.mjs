import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { createHash } from "node:crypto";
import { pathToFileURL } from "node:url";
const url = process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href;
const wasm = await import(url);
const text = readFileSync(new URL("../../../validation/openbabel_issues/ccd/water.cif", import.meta.url), "utf8");
const record = JSON.parse(wasm.parse_ccd_json(text));
assert.equal(record.component_id, "SYN");
assert.equal(record.atoms.length, 3);
assert.equal(record.bonds.length, 2);
assert.equal(record.atoms[0].atom_id, "O#1");
assert.deepEqual(record.atoms[0].model_coordinates, [1, 2, 3]);
assert.deepEqual(record.atoms[0].ideal_coordinates, [0, 0, 0]);
const mol = wasm.mol_from_ccd(text, false);
assert.equal(mol.atom_count(), 3);
assert.equal(mol.bond_count(), 2);
mol.free();
for (const bad of ["", text + "\ndata_OTHER\n", text.replaceAll("SING N N", "DELO N N"), text.replace("N N 1 2 3", "N N NaN 2 3"), " ".repeat(1_000_001)]) {
  assert.throws(() => wasm.parse_ccd_json(bad));
}
const chiral = "data_T\n_chem_comp.id T\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_stereo_config\nT C C 0 R\nT F F 0 N\nT CL Cl 0 N\nT BR Br 0 N\nT H H 0 N\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nT C F SING\nT C CL SING\nT C BR SING\nT C H SING\n";
assert.equal(JSON.parse(wasm.parse_ccd_json(chiral)).atoms[0].stereo_config, "R");
assert.throws(() => wasm.mol_from_ccd(chiral, false));
const graph = wasm.mol_from_ccd(chiral, true);
assert.equal(graph.atom_count(), 5);
graph.free();
console.log("CCD record, explicit graph, malformed input, bounds and stereo-loss contracts passed");

const evidence = JSON.parse(readFileSync(new URL("../../../validation/results/openbabel-ccd-2026-10-16.json", import.meta.url), "utf8"));
assert.deepEqual(record, evidence.synthetic_record);
if (process.env.CHEMATIC_CCD_SOURCES) {
  for (const test of evidence.official_components) {
    const source = readFileSync(`${process.env.CHEMATIC_CCD_SOURCES}/${test.component_id}.cif`);
    assert.equal(createHash("sha256").update(source).digest("hex"), test.input_sha256);
    assert.deepEqual(JSON.parse(wasm.parse_ccd_json(source.toString("utf8"))), test.record);
    if (test.graph_refusal) {
      assert.throws(() => wasm.mol_from_ccd(source.toString("utf8"), true), /formal charge is unknown/);
    } else {
      const graph = wasm.mol_from_ccd(source.toString("utf8"), test.stereo_loss_consent_required);
      assert.equal(graph.atom_count(), test.atoms);
      assert.equal(graph.bond_count(), test.bonds);
      graph.free();
    }
  }
  console.log("Eight official CCD records match the Python/Gemmi evidence in real WASM");
}

const unknownCharge = text.replace("'O#1' O 0", "'O#1' O ?");
assert.equal(JSON.parse(wasm.parse_ccd_json(unknownCharge)).atoms[0].charge, null);
assert.throws(() => wasm.mol_from_ccd(unknownCharge, true));

const longTag = "_chem_comp_atom.extra_" + "x".repeat(20_000);
const expanded = "data_A\n_chem_comp.id A\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n" + longTag + "\n" + Array.from({length: 2000}, (_,i) => `A C${i} C v\n`).join("");
assert.throws(() => wasm.parse_ccd_json(expanded), /expanded loop bytes/);
