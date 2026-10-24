import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const w = await import(process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL('../pkg-node/chematic_wasm.js', import.meta.url).href);
const cases=JSON.parse(readFileSync(new URL('../../../validation/cdxml_bond_meaning/cases.json',import.meta.url),'utf8')).cases;
function smiles(m) {try{return m.canonical_smiles()}finally{m.free()}}
for(const c of cases) {
  if(c.group==='refusal') {
    for(const infer of [false,true]) {
      assert.throws(()=>w.mol_from_cdxml_with_options(c.input,infer),/unsupported CDXML bond order/);
      assert.throws(()=>w.cdxml_to_smiles_json_with_options(c.input,infer),/unsupported CDXML bond order/);
    }
    for(const format of ['smiles','cdxml'])assert.throws(()=>w.convert_common_format(c.input,'cdxml',format),/unsupported CDXML bond order/);
    continue;
  }
  const expected=smiles(w.parse_smiles(c.expected.smiles));
  for(const infer of [false,true]) {
    assert.equal(smiles(w.mol_from_cdxml_with_options(c.input,infer)),expected,c.id);
    assert.deepEqual(JSON.parse(w.cdxml_to_smiles_json_with_options(c.input,infer)),[expected],c.id);
  }
  assert.equal(w.convert_common_format(c.input,'cdxml','smiles'),expected,c.id);
  assert.equal(smiles(w.mol_from_cdxml(w.convert_common_format(c.input,'cdxml','cdxml'))),expected,c.id);
  assert.equal(smiles(w.mol_from_cdxml(w.convert_common_format(c.expected.smiles,'smiles','cdxml'))),expected,c.id);
}
console.log('CDXML bond meaning: 12 graph/class/direction cases and 12 refusals passed');
