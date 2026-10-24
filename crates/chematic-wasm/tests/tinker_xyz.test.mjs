import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
const require = createRequire(import.meta.url);
const __dirname = fileURLToPath(new URL('.', import.meta.url));
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const wasm = require(process.env.CHEMATIC_WASM_MODULE || '../pkg-node/chematic_wasm.js');
const packet = JSON.parse(fs.readFileSync(path.resolve(__dirname, '../../../validation/tinker_xyz/cases.json'), 'utf8'));
for (const c of packet.cases) {
    const record = JSON.parse(wasm.tinker_xyz_record_json(c.input));
    assert.deepEqual(record, c.expected, c.id);
    const output = wasm.tinker_xyz_from_record_json(JSON.stringify(record));
    assert.deepEqual(JSON.parse(wasm.tinker_xyz_record_json(output)), record, c.id);
}
for (const c of packet.refusals) {
    assert.throws(() => wasm.tinker_xyz_record_json(c.input), undefined, c.id);
}
const record = structuredClone(packet.cases[0].expected);
record.atoms[0].name = 'OW';
record.atoms[0].atom_type = 90001;
record.atoms[0].position_angstrom = [.12345678901234567, 0, 0];
record.cell = {lengths_angstrom:[10,11,12], angles_degrees:[70,80,90]};
assert.deepEqual(JSON.parse(wasm.tinker_xyz_record_json(wasm.tinker_xyz_from_record_json(JSON.stringify(record)))), record);
record.atoms[0].isotope = 18;
assert.throws(() => wasm.tinker_xyz_from_record_json(JSON.stringify(record)));
assert.throws(() => wasm.tinker_xyz_record_json(' '.repeat(1000001)));
// The editable JSON may be larger than the compact text input. The writer
// must accept records returned by the reader within its JSON output bound.
const largeText = '10000 compact\n' + Array.from({length:10000}, (_,i) => `${i+1} ABC 1.2345678901234567 2.345678901234567 3.456789012345678 90123\n`).join('');
const largeJson = wasm.tinker_xyz_record_json(largeText);
assert.ok(largeJson.length > 1000000);
const largeOutput = wasm.tinker_xyz_from_record_json(largeJson);
assert.deepEqual(JSON.parse(wasm.tinker_xyz_record_json(largeOutput)), JSON.parse(largeJson));
console.log(`Tinker XYZ: ${packet.cases.length} records, ${packet.refusals.length} refusals, edit and boundary checks passed`);
