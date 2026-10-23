import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(process.env.CHEMATIC_WASM_MODULE?pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href:new URL('../pkg-node/chematic_wasm.js',import.meta.url).href);
const packet=JSON.parse(readFileSync(new URL('../../../validation/orca-checked-2026-10-24/cases.json',import.meta.url),'utf8'));
const bits=v=>{const b=Buffer.alloc(8);b.writeDoubleBE(v);return b.toString('hex');};
for(const c of packet.cases){
  if(c.expected==='refuse') { assert.throws(()=>wasm.write_orca_input_json_checked(c.input_json),c.id); continue; }
  // Feed raw JSON: JS JSON.stringify would normalize negative zero.
  const text=wasm.write_orca_input_json_checked(c.input_json);
  const got=JSON.parse(wasm.orca_input_to_json(text)),expected=JSON.parse(JSON.stringify(c.fields));
  expected.coords.type=expected.coords.kind;delete expected.coords.kind;
  // Check metadata separately from floating values and sign bits.
  assert.deepEqual(got.comments,expected.comments,c.id);
  assert.deepEqual(got.keywords,expected.keywords,c.id);
  assert.deepEqual(got.blocks,expected.blocks,c.id);
  assert.equal(got.coords.type,'xyz');
  assert.equal(got.coords.charge,expected.coords.charge);
  assert.equal(got.coords.multiplicity,expected.coords.multiplicity);
  assert.equal(got.coords.atoms.length,expected.coords.atoms.length);
  for(const [i,a] of got.coords.atoms.entries()){
    for(const key of ['element','frozen','extra'])assert.deepEqual(a[key],expected.coords.atoms[i][key]);
    assert.deepEqual([a.x,a.y,a.z].map(bits),c.coordinate_bits[i],c.id);
  }
  assert.equal(wasm.write_orca_input_json_checked(wasm.orca_input_to_json(text)),text,c.id);
}
const base=packet.cases[0].input_json;
const exact=base+' '.repeat(1_000_000-Buffer.byteLength(base));
assert.equal(wasm.write_orca_input_json_checked(exact),wasm.write_orca_input_json_checked(base));
assert.throws(()=>wasm.write_orca_input_json_checked(exact+' '));
assert.throws(()=>wasm.write_orca_input_json_checked('é'.repeat(500_001)));
const doc=JSON.parse(base);doc.comments=[];doc.keywords=[];
const atom=doc.coords.atoms[1];doc.coords.atoms=Array(1022).fill(atom);
const boundary=wasm.write_orca_input_json_checked(JSON.stringify(doc));
assert.equal(JSON.parse(wasm.orca_input_to_json(boundary)).coords.atoms.length,1022);
doc.coords.atoms.push(atom);
assert.throws(()=>wasm.write_orca_input_json_checked(JSON.stringify(doc)),/lines/);
console.log('ORCA checked output: 24 field/bit-exact controls, 38 refusals, raw JSON and UTF-8 input boundaries passed');
