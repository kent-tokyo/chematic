import assert from 'node:assert/strict';
import { createRequire } from 'node:module';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

const require = createRequire(import.meta.url);
const wasm = require(process.env.CHEMATIC_WASM_MODULE
  ?? fileURLToPath(new URL('../pkg-node/chematic_wasm.js', import.meta.url)));
const source = readFileSync(new URL('../../../validation/openbabel-2429-2026-10-24/either-single.mol', import.meta.url), 'utf8');
let checked = 0;
for (const side of ['left', 'right']) for (const reverse of [false, true]) {
  for (const mirror of [false, true]) for (const threeD of [false, true]) {
    let text = source.replace(' CFG=2', '');
    const [id, a, b] = side === 'left' ? ['50', '10', '20'] : ['70', '30', '40'];
    text = text.replace(`${id} 1 ${a} ${b}`, `${id} 1 ${reverse ? b : a} ${reverse ? a : b} CFG=2`);
    if (mirror) text = text.replace('40 C 2 -1', '40 C 2 1');
    if (threeD) text = text.replace(' 2D', ' 3D').replace('10 C -1 1 0', '10 C -1 1 0.5');
    const mol = wasm.mol_from_v3000_block(text);
    const smiles = mol.canonical_smiles();
    assert(!/[\/\\]/.test(smiles), smiles);
    const output = wasm.to_mol_v3000_block(mol);
    const reread = wasm.mol_from_v3000_block(output);
    assert.equal(reread.canonical_smiles(), smiles);
    assert.equal(wasm.roundtrip_mol_v3000_block(text), text);
    mol.free(); reread.free();
    checked++;
  }
}
for (const cfg of ['', ' CFG=2']) {
  const text = source.replace(' CFG=2', '').replace('60 2 20 30', `60 2 20 30${cfg}`);
  const mol = wasm.mol_from_v3000_block(text);
  assert.equal(/[\/\\]/.test(mol.canonical_smiles()), cfg === '');
  mol.free(); checked++;
}
console.log(`V3000 either-single WASM contract passed: ${checked} cases`);
