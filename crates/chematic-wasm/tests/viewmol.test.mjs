import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
import test from 'node:test';
const require=createRequire(import.meta.url);
const w=require(process.env.CHEMATIC_WASM_MODULE??fileURLToPath(new URL('../pkg-node/chematic_wasm.js',import.meta.url)));
const p=JSON.parse(readFileSync(new URL('../../../validation/viewmol_preservation/cases.json',import.meta.url)));
test('ViewMol raw units, title, flags and binary64 projection survive three cycles',()=>{
 for(const c of p.valid){let r=JSON.parse(w.parse_viewmol_json(c.text));assert.deepEqual(r,c.expected,c.id);
  for(let i=0;i<3;i++){assert.deepEqual(JSON.parse(w.viewmol_positions_angstrom_json(JSON.stringify(r))),c.expected_angstrom,c.id);r=JSON.parse(w.parse_viewmol_json(w.write_viewmol_from_json(JSON.stringify(r))));assert.deepEqual(r,c.expected,c.id);}
 }
});
test('ViewMol refuses unsupported input, unsafe projections, JSON violations and limits',()=>{
 for(const c of p.refused)assert.throws(()=>w.parse_viewmol_json(c.text),/^ViewMol line /,c.id);
 for(const c of p.projection_refused){const r=w.parse_viewmol_json(c.text);assert.deepEqual(JSON.parse(r),c.expected);assert.deepEqual(JSON.parse(w.parse_viewmol_json(w.write_viewmol_from_json(r))),c.expected);assert.throws(()=>w.viewmol_positions_angstrom_json(r),/^ViewMol line /);}
 const r=structuredClone(p.valid[0].expected);r.title='edited';r.factor_to_angstrom='.5';r.atoms[0].position[0]='-.12345678901234567';r.atoms[0].fixed_axes='zy';const s=JSON.stringify(r);
 assert.deepEqual(JSON.parse(w.parse_viewmol_json(w.write_viewmol_from_json(s))),r);
 for(const bad of [s.replace('"title":','"title":null,"title":'),s.replace('"title":','"t\\u0069tle":null,"title":'),s.replace('"element":','"charge":1,"element":'),s.replace('"element":','"element":"N","element":'),s.replace('".5"','0.5'),s.replace('"zy"','""')])assert.throws(()=>w.write_viewmol_from_json(bad),/^ViewMol line /);
 assert.throws(()=>w.parse_viewmol_json('x'.repeat(1000001)),/^ViewMol line /);
 assert.throws(()=>w.parse_viewmol_json('$coord 1\n'+'.1 .2 .3 C\n'.repeat(10001)+'$end\n'),/^ViewMol line /);
 r.atoms=Array(10001).fill(r.atoms[0]);assert.throws(()=>w.write_viewmol_from_json(JSON.stringify(r)),/^ViewMol line /);
});
