import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.argv[2]??process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/cml_3d/cases.json',import.meta.url)));
const variants=c=>[c.input,c.input.replaceAll('\n','\r\n'),c.input.split('\n').map(s=>s.trim()).join('')];
const bits=values=>values.flat().map(v=>{const b=Buffer.alloc(8);b.writeDoubleBE(v);return b.toString('hex')});
for(const c of packet.valid)for(const source of variants(c)){
 const data=JSON.parse(wasm.convert_common_format(source,'.CML','cjson'));
 assert.deepEqual(bits(data.atoms.coords['3d']),bits(c.coords),c.id);
 assert.deepEqual(data.atoms.elements.number,c.numbers,c.id);
 assert.deepEqual(data.bonds.connections.index,c.bond_indices,c.id);
 assert.deepEqual(data.bonds.order,c.bond_orders,c.id);
}
for(const c of packet.absent)for(const source of variants(c)){
 assert.equal(JSON.parse(wasm.convert_common_format(source,'cml','cjson')).atoms.coords,undefined,c.id);
}
for(const c of packet.invalid)for(const source of variants(c))assert.throws(()=>wasm.convert_common_format(source,'cml','cjson'),undefined,c.id);
assert.throws(()=>wasm.convert_common_format(' '.repeat(1000001),'cml','cjson'),/input size/);
const data=JSON.parse(wasm.convert_common_format(packet.valid[0].input,'cml','cjson'));
assert.deepEqual(data.atoms.formalCharges,[0,0,1]);assert.deepEqual(data.atoms.isotopes,[0,13,0]);
const many=n=>'<molecule><atomArray>'+Array.from({length:n},(_,i)=>`<atom id="a${i}" elementType="C" x3="${i}" y3="0" z3="0"/>`).join('')+'</atomArray></molecule>';
assert.equal(JSON.parse(wasm.convert_common_format(many(10000),'cml','cjson')).atoms.elements.number.length,10000);
assert.throws(()=>wasm.convert_common_format(many(10001),'cml','cjson'),/atom elements/);
console.log('CML 3D: 24 source-ordered coordinate/graph cases, 6 absent, 42 invalid, charge/isotope and WASM byte/atom limits passed');
