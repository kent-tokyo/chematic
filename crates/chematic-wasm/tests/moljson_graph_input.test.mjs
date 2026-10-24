import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
import test from 'node:test';
const require=createRequire(import.meta.url);
const wasm=require(process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/moljson_graph_input/cases.json',import.meta.url)));
function snapshot(text) {
 const v=JSON.parse(text), ids=new Map(v.atoms.map((a,i)=>[a.id,i]));
 return {atoms:v.atoms.map(a=>({element:a.element,charge:a.charge,isotope:a.isotope,hydrogens:a.hydrogens,aromatic:a.aromatic})),
   bonds:v.bonds.map(b=>({source:ids.get(b.source_id),target:ids.get(b.target_id),order:b.order,aromatic:b.aromatic}))};
}
test('local MolJSON field graph survives direct/common read-write cycles',()=>{
 for(const c of packet.valid){
  let text=c.input;
  for(let cycle=0;cycle<3;cycle++){
   const mol=wasm.mol_from_moljson(text);
   try{
    const direct=wasm.to_moljson(mol);
    assert.deepEqual(snapshot(direct),c.expected,c.id);
    assert.deepEqual(snapshot(wasm.convert_common_format(text,'moljson','moljson')),c.expected,c.id);
    text=direct;
   }finally{mol.free();}
  }
 }
});
test('malformed fields and ambiguous graph references throw',()=>{
 for(const c of packet.invalid){
  assert.throws(()=>wasm.mol_from_moljson(c.input),c.id);
  assert.throws(()=>wasm.convert_common_format(c.input,'moljson','moljson'),c.id);
 }
});
