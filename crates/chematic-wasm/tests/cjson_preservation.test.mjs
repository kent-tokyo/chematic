import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";
const wasm = await import(process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href);
const source = readFileSync(new URL("../../../validation/cjson_preservation/charged.cjson", import.meta.url), "utf8");
for (const input of [source, source.replaceAll("\n", "\r\n")]) {
  assert.equal(wasm.roundtrip_cjson_record(input), input);
  const record = JSON.parse(wasm.parse_cjson_record_json(input));
  assert.equal(record.schema, "chematic.cjson-record.v1");
  assert.equal(record.source, input);
  assert.deepEqual(record.formal_charges, [0, 1, -1, 0]);
  assert.deepEqual(record.isotopes, [13, 15, null, 2]);
  assert.deepEqual(record.opaque_fields, ["/atoms/coords/2d", "/atoms/elements/symbol", "/atoms/partialCharges", "/name", "/vendor~1a~0b"]);
  assert.ok(Object.is(record.coordinates[0][0], -1.4175225513503917));
  assert.ok(Object.is(record.coordinates[0][1], -0));
  assert.equal(record.coordinates[0][2], Number.MIN_VALUE);
  for (const policy of ["strict", "allow_loss", "preserve_source"]) {
    const report = JSON.parse(wasm.convert_with_report_json(input, "cjson", "cjson", policy));
    assert.equal(report.output, input);
    assert.equal(report.status, "completed");
    assert.equal(report.candidate_graph_verified, false);
  }
}
for (const policy of ["strict", "allow_loss"]) {
  const report = JSON.parse(wasm.convert_with_report_json(source, "cjson", "smiles", policy));
  assert.equal(report.candidate_graph_verified, true);
  assert.equal(report.status, policy === "strict" ? "refused" : "completed_with_loss");
  assert.equal(report.candidate_atom_mapping.length, 4);
  if (policy === "strict") assert.equal(report.output, null);
  else for (const isotope of ["13C", "15N", "2H"]) assert.ok(report.output.includes(isotope));
}
for (const input of [
  source.replace("[0, 1, -1, 0]", "[]"),
  source.replace("[13, 15, 0, 2]", "null"),
  source.replace('"order": [1, 1]', '"order": [2.2, 1]'),
  source.replace("[0, 1, 1, 3]", "[0, 1, 1, 0]"),
  source.replace('"stereo": "opaque"', '"stereo": "opaque", "stereo": "other"'),
  source.replace('"chemicalJson": 1', '"chemicalJson": 0'),
]) {
  assert.throws(() => wasm.roundtrip_cjson_record(input));
  assert.throws(() => wasm.parse_cjson_record_json(input));
  assert.throws(() => wasm.convert_with_report_json(input, "cjson", "cjson", "strict"));
}
assert.throws(() => wasm.roundtrip_cjson_record(" ".repeat(1_000_001)));
const partial = JSON.parse(readFileSync(new URL("../../../validation/cjson_preservation/partial-charges.cjson", import.meta.url), "utf8"));
for (const n of [-0, -1.4175225513503917, Number.MIN_VALUE, Number.MAX_VALUE]) {
  // JSON.stringify normalizes negative zero, so use explicit text for that case.
  const input = structuredClone(partial);
  input.partialCharges.USER_CHARGES[0] = n;
  input.atoms.coords["3d"][0] = n;
  let text = JSON.stringify(input);
  if (Object.is(n, -0)) text = text.replace('"USER_CHARGES":[0,', '"USER_CHARGES":[-0.0,').replace('"3d":[0,', '"3d":[-0.0,');
  const record = JSON.parse(wasm.parse_cjson_record_json(text));
  assert.ok(Object.is(record.partial_charge_datasets.USER_CHARGES[0], n));
  assert.deepEqual(record.opaque_fields, ["/properties/vendor"]);
  const report = JSON.parse(wasm.convert_with_report_json(text, "cjson", "mol2", "allow_loss"));
  assert.equal(report.status, "completed_with_loss");
  assert.equal(report.candidate_graph_verified, true);
  const after = JSON.parse(wasm.parse_mol2_record_json(report.output));
  assert.ok(Object.is(after.atoms[0].partial_charge, n));
  assert.ok(Object.is(after.coordinates[0][0], n));
  assert.equal(after.name, "acetamide charges");
  assert.equal(after.charge_type, "USER_CHARGES");
  const fields = Object.fromEntries(report.observations.map(o => [o.field, o.status]));
  for (const field of ["coordinates", "partial_charges", "charge_method", "molecule_name"]) assert.equal(fields[field], "preserved");
  assert.equal(JSON.parse(wasm.convert_with_report_json(text, "cjson", "mol2", "strict")).output, null);
}
for (const mode of ["multiple", "missing-coordinates", "name-newline", "method-whitespace"]) {
  const input = structuredClone(partial);
  if (mode === "multiple") input.partialCharges.SECOND = [0,0,0,0];
  else if (mode === "missing-coordinates") delete input.atoms.coords;
  else if (mode === "name-newline") input.properties.name = "a\n@<TRIPOS>ATOM";
  else input.partialCharges = {"two methods": [0,0,0,0]};
  const text = JSON.stringify(input);
  assert.equal(wasm.roundtrip_cjson_record(text), text);
  const report = JSON.parse(wasm.convert_with_report_json(text, "cjson", "mol2", "allow_loss"));
  assert.equal(report.status, "refused");
  assert.equal(report.output, null);
  assert.deepEqual(report.candidate_atom_mapping, []);
}
for (const bad of [null, [], {method:[]}, {method:[0,0,0]}, {method:[0,0,0,"NaN"]}]) {
  const input = structuredClone(partial); input.partialCharges = bad;
  const text = JSON.stringify(input);
  assert.throws(() => wasm.roundtrip_cjson_record(text));
  assert.throws(() => wasm.convert_with_report_json(text, "cjson", "mol2", "allow_loss"));
}
console.log("CJSON source retention, strict validation, isotopes and checked conversion passed");
