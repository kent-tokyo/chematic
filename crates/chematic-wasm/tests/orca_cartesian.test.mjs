import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(process.env.CHEMATIC_WASM_MODULE?pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href:new URL('../pkg-node/chematic_wasm.js',import.meta.url).href);
const cases=JSON.parse(readFileSync(new URL('../../../validation/orca-cartesian-2026-10-22/cases.json',import.meta.url),'utf8')).cases;
const bits=v=>{const b=Buffer.alloc(8);b.writeDoubleBE(v);return b.toString('hex');};
for(const c of cases){
  if(c.expected==='refuse'){
    assert.throws(()=>wasm.orca_cartesian_geometry_json(c.source),c.id);
    assert.throws(()=>wasm.orca_input_coords_json(c.source),c.id);
    if(c.core_refuse)assert.throws(()=>wasm.orca_input_to_json(c.source),c.id);
    continue;
  }
  const actual=JSON.parse(wasm.orca_cartesian_geometry_json(c.source)),expected=c.fields;
  assert.equal(actual.schema,'chematic.orca-cartesian.v1');
  for(const key of ['coordinate_units','units_explicit','elements','charge','multiplicity','frozen','atom_extras'])assert.deepEqual(actual[key],expected[key],`${c.id} ${key}`);
  assert.deepEqual(actual.coords.map(row=>row.map(bits)),expected.coordinate_bits,c.id);
  if(expected.coordinate_units==='bohr')assert.throws(()=>wasm.orca_input_coords_json(c.source),/explicit units/);
  else assert.deepEqual(JSON.parse(wasm.orca_input_coords_json(c.source)).coords.map(row=>row.map(bits)),expected.coordinate_bits);
}
const prefix='! Angs\n* xyz 0 1\n';
// The existing 1,024-line budget binds before the 10,000-atom budget.
const source=prefix+'He 0 0 0\n'.repeat(1_021)+'*\n';
assert.equal(JSON.parse(wasm.orca_cartesian_geometry_json(source)).coords.length,1_021);
assert.throws(()=>wasm.orca_cartesian_geometry_json(source.slice(0,-2)+'He 0 0 0\n*\n'),/lines/);
const base='! Bohrs\n* xyz 0 1\nHe 2 0 0\n*\n#';
const exact=base+'x'.repeat(1_000_000-Buffer.byteLength(base));
assert.equal(JSON.parse(wasm.orca_cartesian_geometry_json(exact)).coordinate_units,'bohr');
assert.throws(()=>wasm.orca_cartesian_geometry_json(exact+'x'));
assert.throws(()=>wasm.orca_cartesian_geometry_json(base+'é'.repeat(500_000)));
console.log('ORCA Cartesian units: 16 tagged geometries, 29 refusals, legacy Angstrom guard and resource boundaries passed');
