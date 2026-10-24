import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { test } from 'node:test';
const require = createRequire(import.meta.url);
const wasm = require(process.env.CHEMATIC_WASM_MODULE || '../pkg-node/chematic_wasm.js');
const packet = JSON.parse(readFileSync(new URL('../../../validation/dlpoly_history/cases.json', import.meta.url), 'utf8'));
test('DL-POLY HISTORY original shared valid and refusal packet', () => {
  for (const c of packet.valid) {
    const record = JSON.parse(wasm.dlpoly_history_json(c.input));
    assert.deepEqual(record,c.expected,c.id);
    assert.deepEqual(JSON.parse(wasm.dlpoly_history_json(wasm.dlpoly_history_from_json(JSON.stringify(record)))),record,c.id);
  }
  for (const c of packet.invalid) assert.throws(() => wasm.dlpoly_history_json(c.input),undefined,c.id);
});
test('DL-POLY HISTORY edits, strict fields and browser limits', () => {
  const record = structuredClone(packet.valid[0].expected);
  record.frames[0].atoms[0].charge = -0.98765;
  assert.deepEqual(JSON.parse(wasm.dlpoly_history_json(wasm.dlpoly_history_from_json(JSON.stringify(record)))),record);
  record.frames[0].atoms[0].charge = 0.123456789012345;
  assert.throws(() => wasm.dlpoly_history_from_json(JSON.stringify(record)));
  record.frames[0].atoms[0].charge = 0;
  record.frames[0].atoms[0].unknown = true;
  assert.throws(() => wasm.dlpoly_history_from_json(JSON.stringify(record)));
  assert.throws(() => wasm.dlpoly_history_json(' '.repeat(1_000_001)));
  assert.throws(() => wasm.dlpoly_history_from_json(' '.repeat(16 * 1024 * 1024 + 1)));
});
test('HISTORY browser frame and atom boundaries', () => {
  const record = structuredClone(packet.valid.find(c => c.id === 'level-0').expected);
  record.atom_count = 1;
  const frame = structuredClone(record.frames[0]);
  frame.atoms = frame.atoms.slice(0, 1);
  record.frames = Array.from({ length: 1000 }, () => structuredClone(frame));
  const text = wasm.dlpoly_history_from_json(JSON.stringify(record));
  assert.equal(JSON.parse(wasm.dlpoly_history_json(text)).frames.length, 1000);
  record.frames.push(structuredClone(frame));
  assert.throws(() => wasm.dlpoly_history_from_json(JSON.stringify(record)));
  assert.throws(() => wasm.dlpoly_history_json('title\n0 0 1\n' + 'timestep 1 1 0 0 0.002\nC 1 12.0 0.0\n0.0 0.0 0.0\n'.repeat(1001)));
  record.frames = [frame];
  frame.atoms = Array.from({ length: 10000 }, (_,i) => ({ ...structuredClone(frame.atoms[0]), index: i+1 }));
  record.atom_count = 10000;
  const large = wasm.dlpoly_history_from_json(JSON.stringify(record));
  assert.equal(JSON.parse(wasm.dlpoly_history_json(large)).atom_count, 10000);
  record.frames = [frame, structuredClone(frame)];
  assert.throws(() => wasm.dlpoly_history_from_json(JSON.stringify(record)));
  record.frames = [frame];
  frame.atoms.push({ ...structuredClone(frame.atoms[0]), index:10001 });
  record.atom_count = 10001;
  assert.throws(() => wasm.dlpoly_history_from_json(JSON.stringify(record)));
});
