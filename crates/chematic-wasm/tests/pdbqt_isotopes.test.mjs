import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
const require = createRequire(import.meta.url);
const wasm = require(process.argv[2] ?? process.env.CHEMATIC_WASM_MODULE
  ?? fileURLToPath(new URL('../pkg-node/chematic_wasm.js', import.meta.url)));
const cases = JSON.parse(readFileSync(new URL('../../../validation/openbabel_issues/pdbqt-isotopes/cases.json', import.meta.url)));
const project = (s, c, q, name='LIG') => JSON.parse(wasm.smiles_to_pdbqt_isotope_projection(s, JSON.stringify(c), JSON.stringify(q), name));
const checked = (s, c, q) => wasm.smiles_to_pdbqt_checked(s, JSON.stringify(c), JSON.stringify(q), 'LIG');
for (const c of cases) {
  const coords = c.atomic_numbers.map((_, i) => [i, -2.345, 6.789]);
  const charges = coords.map(() => -.123);
  if (c.isotopes_omitted.length) assert.throws(() => checked(c.smiles, coords, charges), /isotope/);
  const report = project(c.smiles, coords, charges);
  assert.deepEqual(report.isotopes_omitted, c.isotopes_omitted, c.name);
  assert.equal(report.pdbqt, wasm.smiles_to_pdbqt(c.smiles, JSON.stringify(coords), JSON.stringify(charges), 'LIG'));
  if (!c.isotopes_omitted.length) assert.equal(report.pdbqt, checked(c.smiles, coords, charges));
  const lines = report.pdbqt.split('\n').filter(s => s.startsWith('ATOM'));
  assert.equal(lines.length, coords.length);
  for (let i=0; i<lines.length; i++) {
    assert.equal(lines[i].length, 79);
    assert.equal(Number(lines[i].slice(6,11)), i+1);
    assert.deepEqual([[30,38],[38,46],[46,54]].map(([a,b]) => Number(lines[i].slice(a,b))), coords[i]);
    assert.equal(Number(lines[i].slice(70,76)), charges[i]);
  }
}
const coords = [[0,0,0],[1,0,0],[2,0,0],[3,0,0]], charges = [-.2,.3,-.4,.1];
assert.equal(project('[2H]O[13CH2][3H]',coords,charges).pdbqt, checked('[H]O[CH2][H]',coords,charges));
const d8 = cases.find(c => c.name === 'acetophenone-d8').smiles;
const d8Coords = Array(17).fill([1,2,3]), d8Charges = Array.from({length:17},(_,i)=>i*.01);
const d8Report = project(d8,d8Coords,d8Charges);
assert.equal(d8Report.isotopes_omitted.length,8);
assert.equal(d8Report.pdbqt,checked(d8.replaceAll('[2H]','[H]'),d8Coords,d8Charges));
const d8Types = d8Report.pdbqt.split('\n').filter(s=>s.startsWith('ATOM')).map(s=>s.slice(77,79).trim());
assert.equal(d8Types.length,17);
assert.equal(d8Types.filter(s=>s==='H').length,8);
assert.ok(!d8Types.includes('HD'));
assert.ok(!d8Report.pdbqt.includes('TORSDOF'));
for (const bad of [10000, -1000]) assert.throws(() => project('[2H]',[[bad,0,0]],[0]));
for (const [c,q,n] of [[[],[0],'LIG'], [[[0,0,0]],[],'LIG'], [[[0,0,0]],[100],'LIG'], [[[0,0,0]],[0],'\nATOM']]) {
  assert.throws(() => project('[2H]',c,q,n));
}
assert.throws(() => project('[10B]',[[0,0,0]],[0]), /AutoDock type/);
for (const [c,q] of [['[[NaN,0,0]]','[0]'],['[[0,0]]','[0]'],['[[0,0,0]]','[null]'],[' '.repeat(100001),'[0]'],['[[1e400,0,0]]','[0]']]) {
  assert.throws(() => wasm.smiles_to_pdbqt_isotope_projection('[2H]',c,q,'LIG'));
}
console.log(`PDBQT isotope Node/WASM contract passed: ${cases.length} cases`);
