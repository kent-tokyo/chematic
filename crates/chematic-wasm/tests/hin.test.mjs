import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";
const url = process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href;
const wasm = await import(url);
const source = readFileSync(new URL("../../../validation/openbabel_issues/hin/aromatic-ring.hin", import.meta.url), "utf8");
const record = JSON.parse(wasm.parse_hin_record_json(source));
assert.equal(record.schema, "chematic.hin-record.v1");
assert.deepEqual(record.atoms.map(a => a.source_id), [1, 2, 3, 4, 5, 6]);
assert.deepEqual(record.atoms[0].coordinates, [1.2345678901234567, 0, 0.125]);
assert.equal(record.atoms[0].partial_charge, -0.75);
assert.ok(record.bonds.every(b => b.bond_type === "a"));
assert.equal(record.graph_projection.formal_charges, "default_zero_not_inferred_from_partial_charges");
assert.equal(record.source, source);
for (const text of [source, source.replaceAll("\n", "\r\n"), source.trimEnd()]) {
  assert.equal(wasm.roundtrip_hin_record(text), text);
}
for (const bad of ["", source + source, source.replace("endmol 7", "endmol 8"),
  source.replace("-0.75", "NaN"), source.replace("6 a 2 a", "6 s 2 a"),
  source.replace("6 a 2 a", "6 a 99 a"), source.replace("6 a 2 a", "6 q 2 a"),
  source.replace("endmol 7", "res 1 ALA 1\nendmol 7"), " ".repeat(1_000_001)]) {
  assert.throws(() => wasm.parse_hin_record_json(bad));
  assert.throws(() => wasm.roundtrip_hin_record(bad));
}
console.log("HIN source, aromatic bond, precision and rejection contract passed");
