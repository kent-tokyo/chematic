import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
import test from 'node:test';
const require=createRequire(import.meta.url);
const w=require(process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const p=JSON.parse(readFileSync(new URL('../../../validation/gpr_preservation/cases.json',import.meta.url)));
test('GPR records retain all fields across three cycles',()=>{
 for(const c of p.valid){let r=JSON.parse(w.parse_gpr_json(c.text));assert.deepEqual(r,c.expected,c.id);
  for(let i=0;i<3;i++){r=JSON.parse(w.parse_gpr_json(w.write_gpr_from_json(JSON.stringify(r))));assert.deepEqual(r,c.expected,c.id);}
 }
});
test('GPR errors, edited output, strict JSON and fixed WASM limits',()=>{
 for(const c of p.refused)assert.throws(()=>w.parse_gpr_json(c.text),/^GPR line /,c.id);
 const r=structuredClone(p.valid[2].expected);r.coordinates[0].position_nm[0]='.12345678901234567';r.partial_charges[0].value='-.125';r.bonds[0].order='D';
 const s=JSON.stringify(r);assert.deepEqual(JSON.parse(w.parse_gpr_json(w.write_gpr_from_json(s))),r);
 for(const bad of [s.replace('"atoms":','"atoms":[],"atoms":'),s.replace('"atoms":','"at\\u006fms":[],"atoms":'),s.replace('"id":','"id":0,"id":'),s.replace('"id":','"formal_charge":1,"id":'),s.replace('".12345678901234567"','0.125')])assert.throws(()=>w.write_gpr_from_json(bad),/^GPR line /);
 r.partial_charges=null;assert(!w.write_gpr_from_json(JSON.stringify(r)).includes('!Charges'));
 assert.throws(()=>w.parse_gpr_json('x'.repeat(1000001)),/^GPR line /);
 assert.throws(()=>w.parse_gpr_json('!Header gpr 100\n!Info 1\n!Atoms 10001\n'),/^GPR line /);
 assert.throws(()=>w.parse_gpr_json('!Header gpr 100\n!Info 1\n!Atoms 1\n0 6\n!Bonds 20001\n'),/^GPR line /);
 r.atoms=Array(10001).fill(r.atoms[0]);assert.throws(()=>w.write_gpr_from_json(JSON.stringify(r)),/^GPR line /);
});
