import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import {pathToFileURL} from "node:url";
const wasm = await import(process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href);
const fixture = name => readFileSync(new URL(`../../../validation/rdkit_issues/cdxml/${name}.cdxml`, import.meta.url), "utf8");
const text = fixture("nondirectional-stereo");
const codes = mol => JSON.parse(wasm.cip_assignments_json(mol)).map(v=>v.cipCode);
for (const [display, expected] of [["Bold","S"],["Hash","R"],["Dash","R"]]) {
  const input = text.replaceAll("Bold",display);
  for (const read of [wasm.mol_from_cdxml, t=>wasm.mol_from_cdxml_with_options(t,false)]) {
    const m = read(input);
    assert.deepEqual(codes(m),[]);
    assert.equal(m.atom_count(),6);
    m.free();
  }
  for (const input2 of [input,input.replace('B="1" E="2" Display=', 'B="2" E="1" Display=')]) {
    const m = wasm.mol_from_cdxml_with_options(input2,true);
    assert.deepEqual(codes(m),[expected]);
    assert.equal(m.atom_count(),6);
    m.free();
  }
}
for (const input of [fixture("ambiguous-nondirectional"),text.replace(/ p="[^"]*"/g,"")]) {
  const m = wasm.mol_from_cdxml_with_options(input,true);
  assert.deepEqual(codes(m),[]);
  m.free();
}
const multi = text.replace("</page>", `<fragment><n id="99" Element="8"/></fragment></page>`);
const first = wasm.mol_from_cdxml_with_options(multi,true);
assert.equal(first.atom_count(),6);
first.free();
const defaults = JSON.parse(wasm.cdxml_to_smiles_json(multi));
assert.deepEqual(JSON.parse(wasm.cdxml_to_smiles_json_with_options(multi,false)),defaults);
const enabled = JSON.parse(wasm.cdxml_to_smiles_json_with_options(multi,true));
assert.equal(enabled.length,2);
assert.equal(enabled[1],defaults[1]);
assert(!defaults[0].includes("@"));
assert(enabled[0].includes("@"));
for (const read of [wasm.mol_from_cdxml_with_options,wasm.cdxml_to_smiles_json_with_options]) {
  assert.throws(()=>read('<CDXML><fragment><n id="1" Element="999"/></fragment></CDXML>',true));
  assert.throws(()=>read(" ".repeat(1000001),true));
}
assert.throws(()=>wasm.cdxml_to_smiles_json_with_options(`<CDXML>${'<fragment><n id="1"/></fragment>'.repeat(1025)}</CDXML>`,true));
console.log("CDXML stereo options: default compatibility, CIP, endpoint reversal, abstention, batch and limits passed");
