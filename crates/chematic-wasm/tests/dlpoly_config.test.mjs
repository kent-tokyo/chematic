import {createRequire} from 'node:module';
import fs from 'node:fs';
import assert from 'node:assert/strict';
const require=createRequire(import.meta.url);
const wasm=require(process.env.CHEMATIC_WASM_MODULE || '../pkg-node/chematic_wasm.js');
const packet=JSON.parse(fs.readFileSync(new URL('../../../validation/dlpoly_config/cases.json',import.meta.url),'utf8'));
for(const c of packet.cases){
    const record=JSON.parse(wasm.dlpoly_config_record_json(c.input));
    assert.deepEqual(record,c.expected,c.id);
    assert.deepEqual(JSON.parse(wasm.dlpoly_config_record_json(wasm.dlpoly_config_from_record_json(JSON.stringify(record)))),record,c.id);
}
for(const c of packet.refusals)assert.throws(()=>wasm.dlpoly_config_record_json(c.input),undefined,c.id);
const edited=structuredClone(packet.cases[2].expected);
edited.atoms[0].velocity=[1.2345678901234567,0,0];
edited.atoms[0].force=[-123.45678901234567,0,0];
edited.atoms[0].header_extra=['opaque','-0.1250'];
assert.deepEqual(JSON.parse(wasm.dlpoly_config_record_json(wasm.dlpoly_config_from_record_json(JSON.stringify(edited)))),edited);
edited.units='angstrom';assert.throws(()=>wasm.dlpoly_config_from_record_json(JSON.stringify(edited)));
assert.throws(()=>wasm.dlpoly_config_record_json(' '.repeat(1000001)));
const large='large\n2 0 10000\n'+Array.from({length:10000},(_,i)=>`DPD-A ${i+1}\n0 0 0\n1 2 3\n4 5 6\n`).join('');
const json=wasm.dlpoly_config_record_json(large);assert.ok(json.length>1000000);
assert.deepEqual(JSON.parse(wasm.dlpoly_config_record_json(wasm.dlpoly_config_from_record_json(json))),JSON.parse(json));
console.log(`DL-POLY CONFIG: ${packet.cases.length} roundtrips, ${packet.refusals.length} refusals, edited vectors and large JSON passed`);
