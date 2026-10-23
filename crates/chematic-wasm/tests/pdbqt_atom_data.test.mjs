import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
import { createHash } from 'node:crypto';
const require = createRequire(import.meta.url);
const wasm = require(process.argv[2] ?? process.env.CHEMATIC_WASM_MODULE
  ?? fileURLToPath(new URL('../pkg-node/chematic_wasm.js', import.meta.url)));
const packet = JSON.parse(readFileSync(new URL('../../../validation/pdbqt_atom_data/cases.json', import.meta.url)));
for (const c of packet.cases) {
  assert.equal(createHash('sha256').update(c.input.replaceAll('\r\n', '\n')).digest('hex'), c.sha256);
  assert.deepEqual(JSON.parse(wasm.pdbqt_atom_data_json(c.input)), c.expected, c.id);
}
const atom = packet.cases[2].input.split('\n').find(s => s.startsWith('ATOM'));
for (const charge of ['NaN', 'inf', 'bad']) {
  assert.throws(() => wasm.pdbqt_atom_data_json(atom.slice(0, 66) + `    ${charge} Zn`));
}
for (const x of ['NaN', 'bad']) {
  assert.throws(() => wasm.pdbqt_atom_data_json(atom.slice(0, 30) + x.padStart(8) + atom.slice(38)));
}
assert.throws(() => wasm.pdbqt_atom_data_json(atom.slice(0,66) + '    1 Unknown'));
assert.deepEqual(JSON.parse(wasm.pdbqt_atom_data_json('REMARK no atoms')), {atomic_numbers: [], coords: [], partial_charges: []});
assert.throws(() => wasm.pdbqt_atom_data_json(' '.repeat(1000001)), /input bytes/);
assert.equal(JSON.parse(wasm.pdbqt_atom_data_json((atom+'\n').repeat(10000))).coords.length, 10000);
assert.throws(() => wasm.pdbqt_atom_data_json((atom+'\n').repeat(10001)), /atom records/);
console.log('PDBQT numeric data: six fixtures and malformed/limit contracts passed');
