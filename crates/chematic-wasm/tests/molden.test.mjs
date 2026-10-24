import {createRequire} from 'node:module';
import fs from 'node:fs';
import assert from 'node:assert/strict';
const require=createRequire(import.meta.url);
const wasm=require(process.env.CHEMATIC_WASM_MODULE || '../pkg-node/chematic_wasm.js');
const packet=JSON.parse(fs.readFileSync(new URL('../../../validation/molden_vibrations/cases.json',import.meta.url),'utf8'));
for(const c of packet.cases){
    const record=JSON.parse(wasm.molden_record_json(c.input));assert.deepEqual(record,c.expected,c.id);
    assert.deepEqual(JSON.parse(wasm.molden_record_json(wasm.molden_from_record_json(JSON.stringify(record)))),record,c.id);
}
for(const c of packet.refusals)assert.throws(()=>wasm.molden_record_json(c.input),undefined,c.id);
const edited=structuredClone(packet.cases[2].expected);
edited.frequencies[0]=-123.45678901234567;edited.intensities[0].raman=0;
edited.normal_modes[0].displacements_bohr[0]=[1e-130,.12345678901234567,-3];
assert.deepEqual(JSON.parse(wasm.molden_record_json(wasm.molden_from_record_json(JSON.stringify(edited)))),edited);
edited.orbitals=[];assert.throws(()=>wasm.molden_from_record_json(JSON.stringify(edited)));
assert.throws(()=>wasm.molden_record_json(' '.repeat(1000001)));
const large='[Molden Format]\n[Atoms] AU\n'+Array.from({length:10000},(_,i)=>`H ${i+1} 1 .12345678901234567 .23456789012345678 .34567890123456789\n`).join('');
const json=wasm.molden_record_json(large);assert.ok(json.length>1000000);
assert.deepEqual(JSON.parse(wasm.molden_record_json(wasm.molden_from_record_json(json))),JSON.parse(json));
console.log(`Molden: ${packet.cases.length} roundtrips, ${packet.refusals.length} refusals, edits and large JSON passed`);
