// The same hand-derived source is exercised by Rust, Python, and CLI tests.
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";

const moduleUrl = process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href;
const wasm = await import(moduleUrl);
const source = readFileSync(new URL("../../../validation/mol2_preservation/amide.mol2", import.meta.url), "utf8");
const record = JSON.parse(wasm.parse_mol2_record_json(source));
assert.equal(record.schema, "chematic.mol2-record.v1");
assert.equal(record.name, "acetamide");
assert.equal(record.charge_type, "USER_CHARGES");
assert.equal(record.atoms[2].source_id, 90);
assert.equal(record.atoms[2].atom_type, "O.2");
assert.equal(record.atoms[2].partial_charge, -0.70);
assert.deepEqual(record.formal_charges, [0, 0, 0, 0]);
assert.equal(record.bonds[2].bond_type, "am");
assert.deepEqual(record.coordinates[0], [0.123456789, 0.0, 0.0]);
assert.deepEqual(record.opaque_sections, ["SUBSTRUCTURE", "VENDOR_DATA"]);
assert.equal(record.graph_projection.opaque_sections, "not_interpreted");
assert.equal(record.source, source);

for (const input of [source, source.replaceAll("\n", "\r\n")]) {
  assert.equal(wasm.roundtrip_mol2_record(input), input);
  assert.equal(wasm.convert_common_format(input, ".mol2", "mol2"), input);
}
for (const input of [
  source.replace("ACM -0.70", "ACM NaN"),
  source.replace("12 7 3 am", "12 7 3 vendor"),
  source.replace("4 3 1 0 0", "4 4 1 0 0"),
  source + source,
  "C".repeat(1_000_001),
]) {
  assert.throws(() => wasm.parse_mol2_record_json(input));
  assert.throws(() => wasm.roundtrip_mol2_record(input));
  assert.throws(() => wasm.convert_common_format(input, "mol2", "mol2"));
}
console.log("MOL2 source preservation and rejection contract passed");
