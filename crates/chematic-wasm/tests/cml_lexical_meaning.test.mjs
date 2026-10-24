import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const w = await import(process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL('../pkg-node/chematic_wasm.js', import.meta.url).href);
const packet = JSON.parse(readFileSync(new URL('../../../validation/cml_lexical_meaning/cases.json', import.meta.url), 'utf8'));
function fields(text) {
  const v = JSON.parse(text), ids = new Map(v.atoms.map((a,i) => [a.id,i]));
  return {
    numbers:v.atoms.map(a => ({C:6,N:7,O:8})[a.element]),
    charges:v.atoms.map(a => a.charge ?? 0), isotopes:v.atoms.map(a => a.isotope ?? 0),
    hydrogens:v.atoms.map(a => a.hydrogens),
    bonds:v.bonds.map(b => [Math.min(ids.get(b.source_id),ids.get(b.target_id)),
      Math.max(ids.get(b.source_id),ids.get(b.target_id)),b.order]).sort(),
  };
}
for (const c of packet.cases) {
  const {coords2d, atom_maps, ...expected} = c.expected;
  for (const reader of [w.mol_from_cml,w.mol_from_cml_strict]) {
    const mol = reader(c.input);
    try { assert.deepEqual(fields(w.to_moljson(mol)),expected,c.id); }
    finally { mol.free(); }
  }
  assert.deepEqual(fields(w.convert_common_format(c.input,'cml','moljson')),expected,c.id);
  const report = JSON.parse(w.convert_with_report_json(c.input,'cml','cml','preserve_source'));
  assert.equal(report.output,c.input);
  assert.equal(report.candidate_graph_verified,false);
}
for (const source of ["<molecule><!--\n<atom elementType='C'/>\n--></molecule>",
  "<molecule><![CDATA[\n<atom elementType='C'/>\n]]></molecule>"]) {
  const mol=w.mol_from_cml(source);
  try { assert.equal(mol.atom_count(),0); } finally { mol.free(); }
  assert.throws(() => w.mol_from_cml_strict(source));
}
console.log('CML lexical meaning: 16 graph projections, source retention and opaque-only refusals passed');
