import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { pathToFileURL } from 'node:url';
const wasm = await import(process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL('../pkg-node/chematic_wasm.js', import.meta.url).href);
const text = readFileSync(new URL('../../../validation/rdkit_issues/pdb-conect/sulfur.pdb', import.meta.url), 'utf8');
const expected = {
  atoms: [{atomic_number:16,charge:0}, {atomic_number:8,charge:0}, {atomic_number:8,charge:0}, {atomic_number:8,charge:-1}, {atomic_number:8,charge:-1}, {atomic_number:1,charge:0}],
  bonds: [[0,1,2],[0,2,2],[0,3,1],[0,4,1]],
  coords: [[0,0,0],[1.4,0,0],[-1.4,0,0],[0,1.4,0],[0,-1.4,0],[0,0,0.8]],
  serials: [10,20,30,40,50,90],
};
for (const input of [text, text.replaceAll('\n','\r\n')]) {
  assert.deepEqual(JSON.parse(wasm.pdb_conect_frame_json(input, false)), expected);
  const mol = wasm.mol_from_pdb_conect(input, false);
  assert.equal(mol.atom_count(), 6);
  const block = wasm.to_mol_v3000_block(mol);
  const rows = block.split('M  V30 BEGIN BOND\n')[1].split('M  V30 END BOND')[0].trim().split('\n');
  assert.deepEqual(rows.map(s => s.split(/\s+/).slice(3,6).map(Number)), [[2,1,2],[2,1,3],[1,1,4],[1,1,5]]);
  mol.free();
}
const extended = text.replace('CONECT   10   20   20   30   30\nCONECT   10   40   50', 'CONECT   10   20   20   30   30   40   50');
assert.throws(() => wasm.pdb_conect_frame_json(extended, false));
assert.deepEqual(JSON.parse(wasm.pdb_conect_frame_json(extended, true)), expected);
for (const input of [
  text.replace('HETATM   20', 'HETATM   10'),
  text.replace('CONECT   20   10   10', 'CONECT   20   10'),
  text.replace('CONECT   50   10', 'CONECT   50   99'),
  text.replace('CONECT   50   10', 'CONECT   50   50'),
  'MODEL        1\n'+text, 'END\n'+text,
  text.slice(0,16)+'A'+text.slice(17), text.slice(0,76)+'  '+text.slice(78),
  text.slice(0,78)+'+1'+text.slice(80), text.slice(0,30)+'     NaN'+text.slice(38),
  text+'CONECT   10   20\n', 'ATOM\n',
]) {
  assert.throws(() => wasm.pdb_conect_frame_json(input, false));
  assert.throws(() => wasm.mol_from_pdb_conect(input, false));
}
assert.equal(JSON.parse(wasm.pdb_conect_frame_json(text.split('\n').filter(s => s.startsWith('HETATM')).join('\n'), false)).bonds.length, 0);
assert.deepEqual(JSON.parse(wasm.pdb_conect_frame_json('', false)), {atoms:[],bonds:[],coords:[],serials:[]});
assert.throws(() => wasm.pdb_conect_frame_json(' '.repeat(1_000_001), false));
const atom = text.split('\n')[0];
const tooManyAtoms = Array.from({length:10_001}, (_, i) => atom.slice(0,6)+String(i+1).padStart(5)+atom.slice(11)).join('\n');
assert.throws(() => wasm.pdb_conect_frame_json(tooManyAtoms, false));
console.log('pdb_conect: graph, charges, coordinates, explicit H, limits profile and 12 rejection cases passed');
