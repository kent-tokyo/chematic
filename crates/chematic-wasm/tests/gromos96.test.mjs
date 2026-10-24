import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import {createRequire} from 'node:module';
const require=createRequire(import.meta.url);
const w=require(process.env.CHEMATIC_WASM_MODULE||'../pkg-node/chematic_wasm.js');
const p=JSON.parse(fs.readFileSync(new URL('../../../validation/gromos96_preservation/cases.json',import.meta.url),'utf8'));
test('GROMOS96 retains every block, source field and caller-owned unit through three cycles',()=>{
 for(const c of p.valid){let r=w.parse_gromos96_json(c.text,c.record.length_unit);
  for(let i=0;i<3;i++){assert.deepEqual(JSON.parse(r),c.record,c.id);assert(Array.isArray(JSON.parse(w.gromos96_geometry_angstrom_json(r))));r=w.parse_gromos96_json(w.write_gromos96_from_json(r),c.record.length_unit);}
  assert.deepEqual(JSON.parse(r),c.record);
 }
 for(const c of p.refused)assert.throws(()=>w.parse_gromos96_json(c.text,'unspecified'),e=>String(e).startsWith('GROMOS96 line '),c.id);
 for(const c of p.projection_refused){const r=w.parse_gromos96_json(c.text,c.record.length_unit);assert.deepEqual(JSON.parse(w.parse_gromos96_json(w.write_gromos96_from_json(r),c.record.length_unit)),c.record);assert.throws(()=>w.gromos96_geometry_angstrom_json(r));}
});
test('GROMOS96 edits, strict JSON and resource limits',()=>{
 const r=structuredClone(p.valid[0].record);const l=r.blocks[1].lines[0];r.blocks[1].lines[0]=l.slice(0,24)+'.5'.padStart(15)+l.slice(39);
 assert.deepEqual(JSON.parse(w.gromos96_geometry_angstrom_json(JSON.stringify(r)))[0].positions_angstrom[0],[5,2.5,5]);
 const s=JSON.stringify(r);
 for(const bad of [s.replace('"blocks":','"blocks":[],"blocks":'),s.replace('"lines":','"l\\u0069nes":[],"lines":'),JSON.stringify({...r,unknown:1})])assert.throws(()=>w.write_gromos96_from_json(bad));
 assert.throws(()=>w.parse_gromos96_json(' '.repeat(1000001),'nm'));
 assert.throws(()=>w.parse_gromos96_json(p.valid[0].text,'bohr'));
 const huge=structuredClone(r);huge.blocks=[r.blocks[0],...Array.from({length:1025},()=>r.blocks[1])];assert.throws(()=>w.write_gromos96_from_json(JSON.stringify(huge)));
 const red=p.valid.find(c=>c.id==='reduced').record;const many=structuredClone(red);many.blocks[1].lines=Array.from({length:10001},()=>red.blocks[1].lines[0]);assert.throws(()=>w.write_gromos96_from_json(JSON.stringify(many)));
});
