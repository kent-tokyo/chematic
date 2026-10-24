import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { test } from 'node:test';
const require = createRequire(import.meta.url);
const wasm = require(process.env.CHEMATIC_WASM_MODULE || '../pkg-node/chematic_wasm.js');
const packet = JSON.parse(readFileSync(new URL('../../../validation/bgf/cases.json', import.meta.url), 'utf8'));
test('BGF original shared valid and refusal packet', () => {
  for (const c of packet.valid) {
    const record = JSON.parse(wasm.bgf_record_json(c.input));
    assert.deepEqual(record,c.expected,c.id);
    assert.deepEqual(JSON.parse(wasm.bgf_record_json(wasm.bgf_from_record_json(JSON.stringify(record)))),record,c.id);
  }
  for (const c of packet.invalid) assert.throws(() => wasm.bgf_record_json(c.input),undefined,c.id);
});
test('BGF edits, strict fields and browser limits', () => {
  const record = structuredClone(packet.valid[0].expected);
  record.atoms[0].partial_charge = -0.98765;
  assert.deepEqual(JSON.parse(wasm.bgf_record_json(wasm.bgf_from_record_json(JSON.stringify(record)))),record);
  record.atoms[0].partial_charge = 0.123456;
  assert.throws(() => wasm.bgf_from_record_json(JSON.stringify(record)));
  record.atoms[0].partial_charge = 0;
  record.atoms[0].unknown = true;
  assert.throws(() => wasm.bgf_from_record_json(JSON.stringify(record)));
  assert.throws(() => wasm.bgf_record_json(' '.repeat(1_000_001)));
  assert.throws(() => wasm.bgf_from_record_json(' '.repeat(16 * 1024 * 1024 + 1)));
});
