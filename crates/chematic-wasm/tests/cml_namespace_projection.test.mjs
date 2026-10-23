import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const w=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const p=JSON.parse(readFileSync(new URL('../../../validation/cml_namespace_projection/cases.json',import.meta.url)));
const variants=s=>[s,s.replaceAll('><','>\n<'),s.replaceAll('><','>\r\n<')];
const bits=values=>values.flat().map(v=>{const b=Buffer.alloc(8);b.writeDoubleBE(v);return b.toString('hex');});
for(const c of p.valid)for(const source of variants(c.input)){
 const d=JSON.parse(w.convert_common_format(source,'cml','cjson')),a=d.atoms,b=d.bonds,e=c.expected,n=e.numbers.length;
 assert.deepEqual(a.elements.number,e.numbers,c.id);
 assert.deepEqual(bits(a.coords['3d']),bits(e.coordinates),c.id);
 assert.deepEqual(a.formalCharges??Array(n).fill(0),e.formal_charges,c.id);
 assert.deepEqual(a.isotopes??Array(n).fill(0),e.isotopes,c.id);
 assert.deepEqual(b.order.map((o,i)=>[b.connections.index[2*i],b.connections.index[2*i+1],o]),e.bonds,c.id);
}
for(const c of p.invalid)for(const source of variants(c.input))assert.throws(()=>w.convert_common_format(source,'cml','cjson'),undefined,c.id);
console.log(`CML namespace projection passed: ${p.valid.length*3} valid and ${p.invalid.length*3} refused cases`);
