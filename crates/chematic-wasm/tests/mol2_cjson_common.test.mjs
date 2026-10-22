import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/mol2_cjson_common/cases.json',import.meta.url)));
const bits=values=>values.map(value=>{const b=Buffer.alloc(8);b.writeDoubleBE(Number(value));return b.toString('hex');});
for(const c of packet.cases){
 for(const source of [c.input,c.input.replaceAll('\n','\r\n'),c.input.trimEnd()]){
  if(c.outcome!=='preserved'){
   assert.throws(()=>wasm.convert_common_format(source,'.MOL2','cjson'),undefined,c.id);
   continue;
  }
  const value=JSON.parse(wasm.convert_common_format(source,'.MOL2','cjson'));
  const a=value.atoms,e=c.expected;
  assert.deepEqual(a.elements.number,e.numbers,c.id);
  assert.deepEqual(a.formalCharges??Array(e.numbers.length).fill(0),e.formal_charges,c.id);
  assert.deepEqual(value.bonds.connections.index,[0,1],c.id);
  assert.deepEqual(value.bonds.order,[1],c.id);
  assert.deepEqual(bits(a.coords['3d']),bits(e.coordinates.flat()),c.id);
  assert.equal(value.properties.name,e.name,c.id);
  if(e.method===null)assert.equal(value.partialCharges,undefined,c.id);
  else{
   assert.deepEqual(Object.keys(value.partialCharges),[e.method],c.id);
   assert.deepEqual(bits(value.partialCharges[e.method]),bits(e.partial_charges),c.id);
  }
 }
}
const generated=wasm.convert_common_format('CCO','smiles','mol2');
assert.equal(JSON.parse(wasm.convert_common_format(generated,'mol2','cjson')).partialCharges,undefined);
assert.throws(()=>wasm.convert_common_format(' '.repeat(1000001),'mol2','cjson'),/input size/);
const many=n=>`@<TRIPOS>MOLECULE\ncount test\n${n} 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n`+Array.from({length:n},(_,i)=>`${i+1} C${i+1} 0 0 0 C.3`).join('\n')+'\n';
assert.equal(JSON.parse(wasm.convert_common_format(many(10000),'mol2','cjson')).atoms.elements.number.length,10000);
assert.throws(()=>wasm.convert_common_format(many(10001),'mol2','cjson'),/limit|maximum|atoms/);
console.log('MOL2 -> CJSON: 27 selected-field projections, 12 ambiguous-charge and 9 invalid refusals; limits passed');
