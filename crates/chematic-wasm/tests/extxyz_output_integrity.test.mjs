import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url),w=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const p=JSON.parse(readFileSync(new URL('../../../validation/extxyz_output_integrity/cases.json',import.meta.url)));
// JSON.stringify normalizes -0. Use a replacer token restored to a JSON number.
function encode(x){return JSON.stringify(x,(_,v)=>Object.is(v,-0)?'__negative_zero__':v).replaceAll('"__negative_zero__"','-0.0');}
function bits(x){return x.flat().map(v=>{const b=Buffer.alloc(8);b.writeDoubleBE(v);return b.toString('hex');});}
const mol=w.parse_smiles('CO');
for(const c of p.valid){
 const out=w.to_extxyz_json(mol,encode(c.coords),encode({lattice:c.lattice,properties:c.properties,info:c.info})),r=JSON.parse(w.extxyz_frame_json(out));
 assert.deepEqual(bits(r.coords),bits(c.coords),c.id);assert.deepEqual(r.info,c.info,c.id);
 assert.deepEqual(r.lattice===null?null:bits(r.lattice),c.lattice===null?null:bits(c.lattice),c.id);
 assert.deepEqual(Object.keys(r.properties).sort(),Object.keys(c.properties).sort(),c.id);
 for(const [name,values] of Object.entries(c.properties))assert.deepEqual(bits(r.properties[name]),bits(values),c.id);
}
for(const c of p.invalid)assert.throws(()=>w.to_extxyz_json(mol,encode(c.coords),encode({lattice:c.lattice,properties:c.properties,info:c.info})),undefined,c.id);
for(const options of ['null','[]','[null,null,null]','true','{"properties":[]}','{"info":3}','{"unknown":1}',
 '{"info":{"note":"a","note":"b"}}','{"properties":{"x":[[1],[2]],"x":[[3],[4]]}}','{"lattice":null,"lattice":null}',
 '{"properties":{"x":[[],[]]}}','{"lattice":[1e999,0,0,0,1,0,0,0,1]}']){
 assert.throws(()=>w.to_extxyz_json(mol,'[[0,0,0],[1,0,0]]',options),undefined,options);
}
mol.free();console.log('ExtXYZ output preflight, metadata and strict option decoding passed');
