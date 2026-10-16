import assert from 'node:assert/strict';
import { readFileSync, writeFileSync } from 'node:fs';
import { pathToFileURL } from 'node:url';
const moduleUrl = process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL('../pkg-node/chematic_wasm.js', import.meta.url).href;
const wasm = await import(moduleUrl);
const folder = new URL('../../../validation/openbabel_issues/mol2-section-order/', import.meta.url);
const fixture = JSON.parse(readFileSync(new URL('cases.json',folder),'utf8'));
const rows = [];
let checks = 0;
for (const item of fixture.cases.filter(c => c.smiles)) {
  const output = wasm.smiles_to_mol2(item.smiles,true);
  const original = wasm.smiles_to_mol2(item.smiles);
  assert.equal(original,wasm.smiles_to_mol2(item.smiles,false));
  if (item.formal_charges.some(c => c !== 0)) {
    assert.ok(original.indexOf('@<TRIPOS>UNITY_ATOM_ATTR') < original.indexOf('@<TRIPOS>BOND'));
  }
  const markers = output.split(/\r?\n/).filter(line => line.startsWith('@<TRIPOS>'));
  const expected = ['@<TRIPOS>MOLECULE','@<TRIPOS>ATOM','@<TRIPOS>BOND'];
  if (item.formal_charges.some(charge => charge !== 0)) expected.push('@<TRIPOS>UNITY_ATOM_ATTR');
  assert.deepEqual(markers,expected);
  const parsed = JSON.parse(wasm.parse_mol2_record_json(output));
  assert.deepEqual(parsed.formal_charges,item.formal_charges);
  assert.ok(parsed.atoms.every(atom => atom.partial_charge === 0));
  rows.push({id:item.id,output,default_output:original});
  checks++;
}
const source = readFileSync(new URL('charged-source.mol2',folder),'utf8');
for (const text of [source,source.replaceAll('\n','\r\n')]) {
  assert.equal(wasm.roundtrip_mol2_record(text),text);
  assert.equal(wasm.convert_common_format(text,'mol2','mol2'),text);
  checks++;
}
assert.equal(wasm.smiles_to_mol2('CO'),readFileSync(new URL('neutral-expected.mol2',folder),'utf8'));
checks++;
if (process.env.CHEMATIC_MOL2_ORDER_DUMP) {
  writeFileSync(process.env.CHEMATIC_MOL2_ORDER_DUMP,JSON.stringify({schema:'chematic.mol2-section-order-node.v1',rows},null,2)+'\n');
}
console.log(`MOL2 section-order checks passed: ${checks}`);
