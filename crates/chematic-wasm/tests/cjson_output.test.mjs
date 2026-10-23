import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';
const require = createRequire(import.meta.url);
const wasm = require(process.argv[2] ?? process.env.CHEMATIC_WASM_MODULE ?? fileURLToPath(new URL('../pkg-node/chematic_wasm.js', import.meta.url)));
const packet = JSON.parse(readFileSync(new URL('../../../validation/cjson_output/cases.json', import.meta.url)));
for (const c of packet.cases) for (const source of [c.input, c.input.replaceAll('><', '>\n<')]) {
  if (c.accepted) {
    const value = JSON.parse(wasm.convert_common_format(source, 'cml', 'cjson'));
    assert.deepEqual(value.atoms.elements.number, c.expected.numbers);
    assert.deepEqual(value.bonds.connections.index, [0, 1]);
    assert.deepEqual(value.bonds.order, [c.expected.bonds[0][2]]);
    assert.deepEqual(value.atoms.coords['3d'], c.expected.coordinates.flat());
  } else {
    assert.throws(() => wasm.convert_common_format(source, 'cml', 'cjson'), new RegExp('cannot preserve bond order ' + c.error_order));
  }
}
for (const source of ['F/C=C/F', 'N[C@@](F)(Cl)Br', '[C]', '[CH4]', 'c1ccccc1', '[CH4:7]', '*C', '[Pt@SP1](Cl)(Br)(I)F']) {
  assert.throws(() => wasm.convert_common_format(source, 'smiles', 'cjson'), /cannot preserve/);
}
assert.deepEqual(JSON.parse(wasm.convert_common_format('CCO', 'smiles', 'cjson')).bonds.order, [1, 1]);
console.log('CJSON checked output: 16 CML projections/refusals, 8 unsupported SMILES, supported topology passed');
