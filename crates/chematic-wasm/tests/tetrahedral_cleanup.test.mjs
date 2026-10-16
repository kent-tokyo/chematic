import assert from 'node:assert/strict';
import fs from 'node:fs';
import {createRequire} from 'node:module';
const require=createRequire(import.meta.url);
const w=require(process.argv[2]);
const fixtures=JSON.parse(fs.readFileSync(new URL('../../../validation/rdkit_issues/tetrahedral-cleanup/cases.json',import.meta.url),'utf8')).cases;
const mirror=s=>s.replaceAll('@@','TEMP').replaceAll('@','@@').replaceAll('TEMP','@');
let decisions=0;
for(const c of fixtures) for(const s of [c.smiles,mirror(c.smiles)]) {
 const m=w.parse_smiles(s);const before=m.canonical_smiles();
 const r=w.remove_redundant_tetrahedral_stereo(m,c.atom_index,1000000);
 const report=JSON.parse(r.report_json());assert.equal(r.removed,c.removed,c.id);assert.equal(report.removed,c.removed);
 assert.equal(r.atom_index,c.atom_index);assert.equal(r.search_visits,report.search_visits);
 const out=r.molecule;assert.equal(out.atom_count(),m.atom_count());assert.equal(out.bond_count(),m.bond_count());
 assert.equal(m.canonical_smiles(),before);
 if(report.witness) {assert.equal(report.witness[c.atom_index],c.atom_index);assert.deepEqual([...report.witness].sort((a,b)=>a-b),Array.from({length:m.atom_count()},(_,i)=>i));}
 else assert.equal(out.canonical_smiles(),before);
 out.free();r.free();m.free();decisions++;
}
const m=w.parse_smiles('F[C@H](F)Cl');
for(const value of [true,false,'1',null,undefined,NaN,Infinity,-1,1.5,2**32]) assert.throws(()=>w.remove_redundant_tetrahedral_stereo(m,value,1000000));
for(const value of [true,false,'1000',null,undefined,NaN,Infinity,-1,1.5,0,10000001,2**32]) assert.throws(()=>w.remove_redundant_tetrahedral_stereo(m,1,value));
assert.throws(()=>w.remove_redundant_tetrahedral_stereo(m,1,1),/budget|visits|exhaust/);
assert.throws(()=>w.remove_redundant_tetrahedral_stereo(m,99,1000000),/outside/);
assert.throws(()=>w.remove_redundant_tetrahedral_stereo(w.parse_smiles('C'.repeat(65)),1,1000000),/64 atoms/);
for(const [s,idx] of [['F[C@H](F)C(Cl)(Br)I',1],['F[C@H](F)*',1],['FC=[C@AL1]=C(Cl)Br',2],['F[C@H](F)N->[Cu]',1],['C/[C@H](C)F',1]]) assert.throws(()=>w.remove_redundant_tetrahedral_stereo(w.parse_smiles(s),idx,1000000));
m.free();console.log(`Tetrahedral cleanup Node/WASM: ${decisions} decisions and boundary contracts passed`);
