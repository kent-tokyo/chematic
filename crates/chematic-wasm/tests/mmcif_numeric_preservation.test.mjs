import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const w = await import(process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL('../pkg-node/chematic_wasm.js', import.meta.url).href);
const cases=JSON.parse(readFileSync(new URL('../../../validation/mmcif_numeric_preservation/cases.json', import.meta.url),'utf8')).cases;
function bits(values) {
  return values.map(value => {const buf=Buffer.alloc(8);buf.writeDoubleBE(value);return buf.toString('hex');});
}
for (const c of cases) {
  const before=JSON.parse(w.mmcif_to_json(c.input));
  const output=w.write_mmcif_json(JSON.stringify(before.atoms),JSON.stringify(before.cell),before.space_group ?? '',c.id);
  const after=JSON.parse(w.mmcif_to_json(output));
  assert.deepEqual(after.atoms,c.expected.atoms,c.id);
  assert.deepEqual(before.atoms,after.atoms,c.id);
  assert.deepEqual(after.cell,c.expected.cell,c.id);
  assert.equal(after.space_group,before.space_group,c.id);
  for(let i=0;i<after.atoms.length;i++) {
    const keys=['x','y','z','occupancy','b_iso'];
    assert.deepEqual(bits(keys.map(k=>after.atoms[i][k])),bits(keys.map(k=>c.expected.atoms[i][k])),c.id);
  }
  if(after.cell) {
    const keys=['a','b','c','alpha','beta','gamma'];
    assert.deepEqual(bits(keys.map(k=>after.cell[k])),bits(keys.map(k=>c.expected.cell[k])),c.id);
  }
}
console.log(`mmCIF numeric preservation: ${cases.length} record round trips passed`);
