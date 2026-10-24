import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
import test from 'node:test';
const require=createRequire(import.meta.url);
const wasm=require(process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const packet=JSON.parse(readFileSync(new URL('../../../validation/fract_preservation/cases.json',import.meta.url)));
test('fract tokens survive three cycles',()=>{
 for(const c of packet.valid){
  let r=JSON.parse(wasm.parse_fract_json(c.text));
  assert.deepEqual(r,c.expected,c.id);
  for(let i=0;i<3;i++){
   r=JSON.parse(wasm.parse_fract_json(wasm.write_fract_from_json(JSON.stringify(r))));
   assert.deepEqual(r,c.expected,c.id);
  }
 }
});
test('fract refuses malformed records, unknown/duplicate keys and over-limit input',()=>{
 for(const c of packet.refused) assert.throws(()=>wasm.parse_fract_json(c.text),/^fract line /,c.id);
 const r=structuredClone(packet.valid[0].expected);
 r.title='edited';r.cell[0]='7.25(2)';r.atoms[0].fractional[0]='-1.2345678901234567(003)';
 const s=JSON.stringify(r);
 assert.deepEqual(JSON.parse(wasm.parse_fract_json(wasm.write_fract_from_json(s))),r);
 for(const bad of [s.replace('"title":','"title":"a","title":'),s.replace('"title":','"t\\u0069tle":"a","title":'),s.replace('"element":','"charge":1,"element":'),s.replace('"element":','"element":"O","element":'),s.replace('"7.25(2)"','7.25')]) assert.throws(()=>wasm.write_fract_from_json(bad),/^fract line /);
 assert.throws(()=>wasm.parse_fract_json('x'.repeat(1000001)),/^fract line /);
 const many='large\n4 5 6 90 90 90\n'+'C 0 0 0\n'.repeat(10001);
 assert.throws(()=>wasm.parse_fract_json(many),/^fract line /);
 r.atoms=Array(10001).fill(r.atoms[0]);
 assert.throws(()=>wasm.write_fract_from_json(JSON.stringify(r)),/^fract line /);
});
