import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";
const moduleUrl = process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href;
const wasm = await import(moduleUrl);
const source = readFileSync(new URL("../../../validation/mol2_preservation/amide.mol2", import.meta.url), "utf8");
for (const policy of ["strict", "allow_loss"]) {
  const report = JSON.parse(wasm.convert_with_report_json(source, "mol2", "smiles", policy));
  assert.equal(report.schema, "chematic.conversion-report.v1");
  assert.equal(report.status, policy === "strict" ? "refused" : "completed_with_loss");
  assert.equal(report.candidate_graph_verified, true);
  assert.equal(report.candidate_atom_mapping.length, 4);
  assert.equal(report.candidate_bond_mapping.length, 3);
  const fields = Object.fromEntries(report.observations.map(o => [o.field, o.status]));
  for (const field of ["coordinates", "partial_charges", "opaque_sections", "atom_types"]) assert.equal(fields[field], "dropped");
  for (const field of ["formal_charge_interpretation", "stereochemistry"]) assert.equal(fields[field], "unverified");
  if (policy === "strict") assert.equal(report.output, null);
  else assert.equal(typeof report.output, "string");
}
for (const policy of ["strict", "allow_loss", "preserve_source"]) {
  const input = source.replaceAll("\n", "\r\n");
  const report = JSON.parse(wasm.convert_with_report_json(input, " .MOL2 ", "mol2", policy));
  assert.equal(report.output, input);
  assert.equal(report.verification_scope, "exact_source_bytes");
  assert.equal(report.candidate_graph_verified, false);
  assert.ok(report.observations.some(o => o.status === "opaque_retained"));
}
for (const args of [
  ["invalid", "mol2", "mol2", "strict"],
  [source, "sdf", "smiles", "strict"],
  [source, "mol2", "cml", "allow_loss"],
  [source, "mol2", "smiles", "preserve_source"],
  [source, "mol2", "smiles", "typo"],
  ["C".repeat(1_000_001), "mol2", "mol2", "strict"],
]) assert.throws(() => wasm.convert_with_report_json(...args));
console.log("Checked conversion policy, losses, mapping and refusal contract passed");

for (const policy of ["strict", "allow_loss"]) {
  const input = source.replace("0.123456789 0 0 C.3", "0.123456789 2.5 -3.75 N.4");
  const report = JSON.parse(wasm.convert_with_report_json(input, "mol2", "cjson", policy));
  const fields = Object.fromEntries(report.observations.map(o => [o.field, o.status]));
  assert.equal(fields.coordinates, "preserved");
  assert.equal(fields.projected_formal_charges, "preserved");
  assert.equal(fields.partial_charges, "preserved");
  assert.equal(fields.charge_method, "preserved");
  assert.equal(fields.molecule_name, "preserved");
  assert.equal(fields.formal_charge_interpretation, "unverified");
  assert.equal(fields.implicit_hydrogens, undefined);
  assert.equal(report.candidate_graph_verified, true);
  if (policy === "strict") {
    assert.equal(report.status, "refused");
    assert.equal(report.output, null);
  } else {
    const candidate = JSON.parse(report.output);
    assert.deepEqual(candidate.atoms.coords["3d"].slice(0, 3), [0.123456789, 2.5, -3.75]);
    assert.deepEqual(candidate.atoms.formalCharges, [1, 0, 0, 0]);
  }
}
for (const value of ["-1.4175225513503917", "-0.0", "5e-324", "1.7976931348623157e308"]) {
  const input = source.replace("0.123456789", value);
  const report = JSON.parse(wasm.convert_with_report_json(input, "mol2", "cjson", "allow_loss"));
  const candidate = JSON.parse(report.output);
  assert.ok(Object.is(candidate.atoms.coords["3d"][0], Number(value)));
  assert.ok(report.observations.some(o => o.field === "coordinates" && o.status === "preserved"));
}
for (const token of ["nc", "du", "un"]) {
  const input = source.replace("12 7 3 am", `12 7 3 ${token}`);
  const report = JSON.parse(wasm.convert_with_report_json(input, "mol2", "cjson", "allow_loss"));
  assert.equal(report.candidate_graph_verified, false);
  assert.deepEqual(report.candidate_atom_mapping, []);
}
assert.throws(() => wasm.convert_with_report_json(source, "mol2", "cjson", "preserve_source"));
console.log("CJSON coordinates, projected charges, precision and refusal contract passed");

for (const literal of ["-0.0", "-1.4175225513503917", "5e-324", "1.7976931348623157e308"]) {
  const input = source.replace("ACM 0.70", `ACM ${literal}`);
  const report = JSON.parse(wasm.convert_with_report_json(input, "mol2", "cjson", "allow_loss"));
  const candidate = JSON.parse(report.output);
  assert.ok(Object.is(candidate.partialCharges.USER_CHARGES[0], Number(literal)));
  assert.equal(candidate.properties.name, "acetamide");
  assert.ok(report.observations.some(o => o.field === "partial_charges" && o.status === "preserved"));
}
for (const input of [source.replace("ACM 0.70", "ACM"), source.replace("USER_CHARGES", "NO_CHARGES")]) {
  const report = JSON.parse(wasm.convert_with_report_json(input, "mol2", "cjson", "allow_loss"));
  assert.equal(JSON.parse(report.output).partialCharges, undefined);
  assert.ok(report.observations.some(o => o.field === "partial_charges" && o.status === "dropped"));
}
console.log("CJSON declared partial-charge values, method, name and missing-data contract passed");
