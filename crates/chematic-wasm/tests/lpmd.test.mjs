import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
const require=createRequire(import.meta.url);
const wasm=require(process.env.CHEMATIC_WASM_MODULE || '../pkg-node/chematic_wasm.js');
const packet=JSON.parse(fs.readFileSync(fileURLToPath(new URL('../../../validation/lpmd_preservation/cases.json',import.meta.url)),'utf8'));
test('LPMD all-frame decimal and column preservation, refusals and explicit projection',()=>{
 for(const c of packet.valid){
  let record=wasm.parse_lpmd_json(c.text);
  for(let i=0;i<3;i++){
   assert.deepEqual(JSON.parse(record),c.record,c.id);
   assert.equal(JSON.parse(wasm.lpmd_positions_angstrom_json(record)).length,c.record.frames.length);
   record=wasm.parse_lpmd_json(wasm.write_lpmd_from_json(record));
  }
  assert.deepEqual(JSON.parse(record),c.record);
 }
 for(const c of packet.refused) assert.throws(()=>wasm.parse_lpmd_json(c.text),e=>String(e).startsWith('LPMD line '),c.id);
 for(const c of packet.projection_refused){
  const record=wasm.parse_lpmd_json(c.text);
  assert.deepEqual(JSON.parse(wasm.parse_lpmd_json(wasm.write_lpmd_from_json(record))),c.record);
  assert.throws(()=>wasm.lpmd_positions_angstrom_json(record),e=>String(e).startsWith('LPMD line '));
 }
});
test('LPMD writer uses edits and rejects ambiguous JSON and resource excess',()=>{
 const r=structuredClone(packet.valid[0].record);r.frames[0].cell[0][0]='10';r.frames[0].atoms[0][1]='.5';
 assert.deepEqual(JSON.parse(wasm.lpmd_positions_angstrom_json(JSON.stringify(r)))[0][0],[5,2,8]);
 const s=JSON.stringify(r);
 for(const bad of [s.replace('"columns":','"columns":[],"columns":'),s.replace('"cell":','"c\\u0065ll":[],"cell":'),JSON.stringify({...r,unknown:1})]) assert.throws(()=>wasm.write_lpmd_from_json(bad));
 assert.throws(()=>wasm.parse_lpmd_json(' '.repeat(1000001)));
 const huge=structuredClone(r);huge.frames=Array.from({length:1025},()=>r.frames[0]);
 assert.throws(()=>wasm.write_lpmd_from_json(JSON.stringify(huge)));
});
