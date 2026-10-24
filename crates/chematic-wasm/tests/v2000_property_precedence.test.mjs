import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/v2000_property_precedence/cases.json',import.meta.url)));
const numbers={C:6,N:7,O:8};
const sortRows=rows=>rows.sort((a,b)=>a[0]-b[0]||a[1]-b[1]||a[2]-b[2]);
function check(mol,e,id){
 const value=JSON.parse(wasm.to_moljson(mol)),atoms=value.atoms;
 assert.deepEqual(sortRows(atoms.map(a=>[numbers[a.element],a.charge,a.isotope??0])),sortRows(e.numbers.map((n,i)=>[n,e.charges[i],e.isotopes[i]])),id);
 const ids=new Map(atoms.map(a=>[a.id,numbers[a.element]]));
 const bonds=value.bonds.map(b=>[...[ids.get(b.source_id),ids.get(b.target_id)].sort((a,b)=>a-b),b.order]);
 const wanted=e.bonds.map(([a,b,o])=>[...[e.numbers[a],e.numbers[b]].sort((a,b)=>a-b),o]);
 assert.deepEqual(sortRows(bonds),sortRows(wanted),id);
 assert.deepEqual([...mol.canonical_smiles().matchAll(/:(\d+)\]/g)].map(m=>Number(m[1])).sort((a,b)=>a-b),e.atom_maps.filter(Boolean).sort((a,b)=>a-b),id);
}
for(const c of packet.cases){
 const direct=wasm.mol_from_sdf_block(c.input);
 try{check(direct,c.expected,c.id);}finally{direct.free();}
 const smiles=[wasm.convert_common_format(c.input,'mol','smiles'),...JSON.parse(wasm.sdf_to_smiles_json(c.input+'$$$$\n'))];
 for(const text of smiles){const m=wasm.parse_smiles(text);try{check(m,c.expected,c.id);}finally{m.free();}}
}
const cases=[...packet.cases,...packet.cases.toReversed()];
const text=cases.map(c=>c.input+'> <note>\nM  CHG  1   1   3\nM  ISO  1   1  99\n\n$$$$\n').join('');
const smiles=JSON.parse(wasm.sdf_to_smiles_json(text));assert.equal(smiles.length,cases.length);
for(let i=0;i<cases.length;i++){const m=wasm.parse_smiles(smiles[i]);try{check(m,cases[i].expected,cases[i].id);}finally{m.free();}}
console.log('V2000 property precedence: 14 cases on direct/common/SDF routes and 28 alternating records passed');
