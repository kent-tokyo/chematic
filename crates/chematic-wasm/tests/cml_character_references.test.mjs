import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import {pathToFileURL} from "node:url";
const wasm = await import(process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href);
const packet = JSON.parse(readFileSync(new URL("../../../validation/cml_character_references/cases.json", import.meta.url), "utf8"));
for (const reader of [wasm.mol_from_cml, wasm.mol_from_cml_strict]) {
  const literal = reader(packet.cases[0].source);
  const expected = wasm.to_cml(literal);
  literal.free();
  for (const row of packet.cases) {
    if (row.expected_refusal) { assert.throws(() => reader(row.source), undefined, row.id); continue; }
    const mol = reader(row.source);
    assert.equal(mol.atom_count(), 2, row.id);
    assert.equal(wasm.to_cml(mol), expected, row.id);
    mol.free();
  }
}
console.log("CML character-reference/typed-field packet passed: 28 public reader cases");
