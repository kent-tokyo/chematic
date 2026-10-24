import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/moljson_hydrogens/cases.json',import.meta.url)));
for(const c of packet.cases){
 const reference=wasm.parse_smiles(c.expected_smiles);
 const expected=reference.canonical_smiles();
 const roundtrip=wasm.mol_from_moljson(wasm.to_moljson(reference));
 try{assert.equal(roundtrip.canonical_smiles(),expected,c.id);}finally{roundtrip.free();reference.free();}
 let mol=wasm.mol_from_moljson(c.input);
 try{
  for(let i=0;i<3;i++){
   assert.equal(mol.canonical_smiles(),expected,c.id);
   const encoded=wasm.to_moljson(mol);
   assert.deepEqual(JSON.parse(encoded).atoms.map(a=>a.hydrogens),c.expected_hydrogens,c.id);
   assert.equal(wasm.convert_common_format(encoded,'moljson','smiles'),expected,c.id);
   const next=wasm.mol_from_moljson(encoded);mol.free();mol=next;
  }
 }finally{mol.free();}
}
for(const hydrogens of packet.invalid_hydrogens){
 const text=JSON.stringify({atoms:[{id:'a1',element:'C',hydrogens}],bonds:[]});
 assert.throws(()=>wasm.mol_from_moljson(text),/hydrogens/);
 assert.throws(()=>wasm.convert_common_format(text,'moljson','smiles'),/hydrogens/);
}
console.log('MolJSON non-node H counts: 13 graphs, repeated direct/common routes and 9 invalid values passed');
