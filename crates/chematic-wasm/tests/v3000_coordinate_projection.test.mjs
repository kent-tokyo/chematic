import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/v3000_coordinate_projection/cases.json',import.meta.url)));
const bits=values=>values.map(value=>{const b=Buffer.alloc(8);b.writeDoubleBE(Number(value));return b.toString('hex');});
for(const c of packet.cases){
 for(const source of [c.input,c.input.replaceAll('\n','\r\n'),c.input.trimEnd()]){
  const out=JSON.parse(wasm.convert_common_format(source,'.V3000','cjson'));
  const a=out.atoms,b=out.bonds,e=c.expected,n=e.numbers.length;
  assert.deepEqual(bits(a.coords['3d']),bits(e.coordinates.flat()),c.id);
  assert.deepEqual(a.elements.number,e.numbers,c.id);
  assert.deepEqual(a.formalCharges??Array(n).fill(0),e.formal_charges,c.id);
  assert.deepEqual(a.isotopes??Array(n).fill(0),e.isotopes,c.id);
  assert.deepEqual(b.order.map((o,i)=>[b.connections.index[2*i],b.connections.index[2*i+1],o]),e.bonds,c.id);
 }
}
const source=packet.cases[0].input;
for(const bad of [source.replace('M  END',''),source.replace('0.123456789012345','NaN'),source.replace('31 1 1 2','31 1 1 99')]){
 assert.throws(()=>wasm.convert_common_format(bad,'v3000','cjson'));
}
assert.throws(()=>wasm.convert_common_format(' '.repeat(1000001),'v3000','cjson'),/input size/);
const many=n=>['count test','','',source.split('\n')[3],'M  V30 BEGIN CTAB',`M  V30 COUNTS ${n} 0 0 0 0`,'M  V30 BEGIN ATOM',...Array.from({length:n},(_,i)=>`M  V30 ${i+1} C 0 0 0 0`),'M  V30 END ATOM','M  V30 BEGIN BOND','M  V30 END BOND','M  V30 END CTAB','M  END',''].join('\n');
assert.equal(JSON.parse(wasm.convert_common_format(many(10000),'v3000','cjson')).atoms.coords['3d'].length,30000);
assert.throws(()=>wasm.convert_common_format(many(10001),'v3000','cjson'),/atom count/);
console.log('V3000 coordinate projection: 24 inputs preserve xyz bits and indexed graph, limits/invalid input passed');
