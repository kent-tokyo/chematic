import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
import test from 'node:test';
const require=createRequire(import.meta.url);
const w=require(process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const p=JSON.parse(readFileSync(new URL('../../../validation/unixyz_preservation/cases.json',import.meta.url)));
test('UniChem XYZ title and coordinate strings survive three cycles',()=>{
 for(const c of p.valid){let r=JSON.parse(w.parse_unixyz_json(c.text));assert.deepEqual(r,c.expected,c.id);
  for(let i=0;i<3;i++){r=JSON.parse(w.parse_unixyz_json(w.write_unixyz_from_json(JSON.stringify(r))));assert.deepEqual(r,c.expected,c.id);}
 }
});
test('UniChem XYZ refuses malformed/trailing data, strict JSON violations and limits',()=>{
 for(const c of p.refused)assert.throws(()=>w.parse_unixyz_json(c.text),/^UniChem XYZ line /,c.id);
 const r=structuredClone(p.valid[0].expected);r.title='edited';r.atoms[0].atomic_number=7;r.atoms[0].position_angstrom[0]='-.12345678901234567';
 const s=JSON.stringify(r);assert.deepEqual(JSON.parse(w.parse_unixyz_json(w.write_unixyz_from_json(s))),r);
 for(const bad of [s.replace('"title":','"title":"x","title":'),s.replace('"title":','"t\\u0069tle":"x","title":'),s.replace('"atomic_number":','"charge":1,"atomic_number":'),s.replace('"atomic_number":','"atomic_number":6,"atomic_number":'),s.replace('"-.12345678901234567"','-.125')])assert.throws(()=>w.write_unixyz_from_json(bad),/^UniChem XYZ line /);
 assert.throws(()=>w.parse_unixyz_json('x'.repeat(1000001)),/^UniChem XYZ line /);
 assert.throws(()=>w.parse_unixyz_json('title\n10001\n'),/^UniChem XYZ line /);
 r.atoms=Array(10001).fill(r.atoms[0]);assert.throws(()=>w.write_unixyz_from_json(JSON.stringify(r)),/^UniChem XYZ line /);
});
