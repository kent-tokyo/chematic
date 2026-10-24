import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/v3000_atom_scalars/cases.json',import.meta.url)));
for(const c of packet.cases){
 const mol=wasm.mol_from_v3000_block(c.input);
 try{
  const a=JSON.parse(wasm.to_moljson(mol)).atoms[0],e=c.expected;
  assert.equal(a.element,({1:'H',6:'C',7:'N',8:'O'})[e.number],c.id);
  assert.equal(a.charge,e.charge,c.id);assert.equal(a.isotope??0,e.isotope,c.id);
  const smiles=mol.canonical_smiles();
  assert.equal(e.atom_map?smiles.includes(`:${e.atom_map}]`):!smiles.includes(':'),true,c.id);
  assert.equal(wasm.convert_common_format(c.input,'v3000','smiles'),smiles,c.id);
  assert.equal(wasm.roundtrip_mol_v3000_block(c.input),c.input,c.id);
  assert.equal(wasm.convert_common_format(c.input,'v3000','v3000'),c.input,c.id);
 }finally{mol.free();}
}
for(const c of packet.invalid)for(const call of [wasm.mol_from_v3000_block,wasm.roundtrip_mol_v3000_block,
 s=>wasm.convert_common_format(s,'v3000','smiles'),s=>wasm.convert_common_format(s,'v3000','v3000')]){
 assert.throws(()=>call(c.input),new RegExp(c.field),c.id);
}
console.log('V3000 atom scalars: 8 ordinary + 2 quoted controls and 25 refusals passed');
