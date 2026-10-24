import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { test } from 'node:test';
const require = createRequire(import.meta.url);
const wasm = require(process.env.CHEMATIC_WASM_MODULE || '../pkg-node/chematic_wasm.js');
const packet = JSON.parse(readFileSync(new URL('../../../validation/nwchem_geometry/cases.json', import.meta.url), 'utf8'));
test('NWChem geometry original shared valid and refusal packet', () => {
  for (const c of packet.valid) {
    const record = JSON.parse(wasm.nwchem_geometry_json(c.input));
    assert.deepEqual(record,c.expected,c.id);
    assert.deepEqual(JSON.parse(wasm.nwchem_geometry_json(wasm.nwchem_geometry_from_json(JSON.stringify(record)))),record,c.id);
  }
  for (const c of packet.invalid) assert.throws(() => wasm.nwchem_geometry_json(c.input),undefined,c.id);
});
test('NWChem geometry edits, strict fields and browser limits', () => {
  const record = structuredClone(packet.valid[0].expected);
  record.centers[0].nuclear_charge = -0.98765;
  assert.deepEqual(JSON.parse(wasm.nwchem_geometry_json(wasm.nwchem_geometry_from_json(JSON.stringify(record)))),record);
  record.centers[0].mass_amu = -1;
  assert.throws(() => wasm.nwchem_geometry_from_json(JSON.stringify(record)));
  record.centers[0].mass_amu = 1;
  record.centers[0].nuclear_charge = 0;
  record.centers[0].unknown = true;
  assert.throws(() => wasm.nwchem_geometry_from_json(JSON.stringify(record)));
  assert.throws(() => wasm.nwchem_geometry_json(' '.repeat(1_000_001)));
  assert.throws(() => wasm.nwchem_geometry_from_json(' '.repeat(16 * 1024 * 1024 + 1)));
});
test('NWChem browser center count boundary', () => {
  const record = structuredClone(packet.valid[0].expected);
  record.centers = Array.from({ length: 10_000 }, () => structuredClone(record.centers[0]));
  const text = wasm.nwchem_geometry_from_json(JSON.stringify(record));
  assert.equal(JSON.parse(wasm.nwchem_geometry_json(text)).centers.length, 10_000);
  record.centers.push(structuredClone(record.centers[0]));
  assert.throws(() => wasm.nwchem_geometry_from_json(JSON.stringify(record)));
  assert.throws(() => wasm.nwchem_geometry_json('geometry\n' + 'C 0 0 0\n'.repeat(10_001) + 'end\n'));
});
