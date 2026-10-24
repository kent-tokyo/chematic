import assert from 'node:assert/strict';
import { createRequire } from 'node:module';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
const require = createRequire(import.meta.url);
const wasm = require(process.env.CHEMATIC_WASM_MODULE ?? fileURLToPath(new URL('../pkg-node/chematic_wasm.js', import.meta.url)));
const packet = JSON.parse(readFileSync(new URL('../../../validation/rdkit-9732-2026-10-24/cases.json', import.meta.url)));
function membership(text) {
  const [base, extension=''] = text.split('|');
  const atoms=[];
  for (const [role,side] of base.trim().split(/(?<!-)>/).entries()) {
    for (const component of side.split('.').filter(Boolean)) {
      const mol=wasm.parse_smiles(component);
      const cip=new Map(JSON.parse(mol.assign_cip_json()).centers.map(c=>[c.atom,c.code]));
      for(let i=0;i<mol.atom_count();i++) atoms.push([role,mol.canonical_smiles(),cip.get(i)]);
      mol.free();
    }
  }
  const result=[];
  for(const group of extension.matchAll(/(a|&[0-9]+|o[0-9]+):([0-9]+(?:,[0-9]+)*)/g)) {
    for(const index of group[2].split(',').map(Number)) {
      assert(index<atoms.length);
      const [role,component,cip]=atoms[index];
      assert(cip,`${text}: group ${index} is not a stereocentre`);
      result.push(JSON.stringify([role,component,group[1],cip]));
    }
  }
  return result.sort();
}
for(const c of packet.positive) for(const canonical of [false,true]) {
  const out=wasm.normalize_reaction_cxsmiles(c.input,canonical);
  assert.deepEqual(membership(out),membership(c.input),c.id);
  assert.equal(wasm.normalize_reaction_cxsmiles(out,canonical),out,c.id);
}
for(const c of packet.negative) assert.throws(()=>wasm.normalize_reaction_cxsmiles(c.input,true),undefined,c.id);
assert.throws(()=>wasm.normalize_reaction_cxsmiles(' '.repeat(1000001),true),/input bytes/);
console.log('Reaction CXSMILES: 40 positive mode/case runs, 22 refusals and input limit passed');
