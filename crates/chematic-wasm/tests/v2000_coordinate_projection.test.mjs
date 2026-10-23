import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/v2000_coordinate_projection/cases.json',import.meta.url)));
const bits=values=>values.map(v=>{const b=Buffer.alloc(8);b.writeDoubleBE(Number(v));return b.toString('hex');});
for(const c of packet.cases) {
 for(const text of [c.input,c.input.replaceAll('\n','\r\n'),c.input.trimEnd()]) {
  for(const format of ['.MOL','sdf']) {
   const source=format==='sdf'?text+'\n$$$$\n':text;
   const v=JSON.parse(wasm.convert_common_format(source,format,'cjson'));
   const a=v.atoms,b=v.bonds,e=c.expected,n=e.numbers.length;
   assert.deepEqual(bits(a.coords['3d']),bits(e.coordinates.flat()),c.id);
   assert.deepEqual(a.elements.number,e.numbers);
   assert.deepEqual(a.formalCharges??Array(n).fill(0),e.formal_charges);
   assert.deepEqual(a.isotopes??Array(n).fill(0),e.isotopes);
   assert.deepEqual(b.order.map((o,i)=>[b.connections.index[2*i],b.connections.index[2*i+1],o]),e.bonds);
  }
 }
}
for(let axis=0;axis<3;axis++) for(const raw of ['', 'NaN','inf','-inf','1e999','not-a-num']) {
 const rows=packet.cases[0].input.split('\n'),r=rows[4];
 rows[4]=r.slice(0,axis*10)+raw.padStart(10)+r.slice(axis*10+10);
 assert.throws(()=>wasm.convert_common_format(rows.join('\n'),'mol','cjson'),/coordinate/);
}
const source=packet.cases[0].input;
assert.throws(()=>wasm.convert_common_format(source+'$$$$\n'+source+'$$$$\n','sdf','cjson'),/records exceeds limit 1/);
assert.throws(()=>wasm.convert_common_format(' '.repeat(1000001),'mol','cjson'),/input size/);
console.log(`V2000 source geometry contract passed: ${packet.cases.length*6} projections, 18 malformed fields`);
