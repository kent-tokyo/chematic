import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/v2000_isotope_output/cases.json',import.meta.url)));
const numbers={H:1,C:6,N:7,O:8,Cl:17,Br:35,U:92};
for(const c of packet.cases){
 const mol=wasm.parse_smiles(c.smiles);
 try{
  const outputs=[wasm.to_mol_block(mol),wasm.to_mol_block_strict(mol),
   wasm.convert_common_format(c.smiles,'smiles','mol'),
   wasm.convert_common_format(c.smiles,'smiles','sdf'),
   wasm.smiles_array_to_sdf(JSON.stringify([c.smiles]))];
  for(const output of outputs){
   const after=wasm.mol_from_sdf_block(output);
   try{
    const value=JSON.parse(wasm.to_moljson(after)),atoms=value.atoms,e=c.expected;
    assert.deepEqual(atoms.map(a=>numbers[a.element]),e.numbers,c.id);
    assert.deepEqual(atoms.map(a=>a.charge),e.charges,c.id);
    assert.deepEqual(atoms.map(a=>a.isotope??0),e.isotopes,c.id);
    const ids=new Map(atoms.map((a,i)=>[a.id,i]));
    const bonds=value.bonds.map(b=>[...[ids.get(b.source_id),ids.get(b.target_id)].sort((a,b)=>a-b),b.order]).sort((a,b)=>a[0]-b[0]||a[1]-b[1]||a[2]-b[2]);
    assert.deepEqual(bonds,e.bonds,c.id);
    const maps=[...after.canonical_smiles().matchAll(/:(\d+)\]/g)].map(m=>Number(m[1])).sort((a,b)=>a-b);
    assert.deepEqual(maps,e.atom_maps.filter(Boolean).sort((a,b)=>a-b),c.id);
   }finally{after.free();}
  }
 }finally{mol.free();}
}
console.log('V2000 isotope output: 15 ordinary + 4 storage controls passed on five routes');
