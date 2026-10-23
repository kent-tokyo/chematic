import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(process.env.CHEMATIC_WASM_MODULE?pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href:new URL('../pkg-node/chematic_wasm.js',import.meta.url).href);
const cases=JSON.parse(readFileSync(new URL('../../../validation/openbabel-1464-2026-10-24/cases.json',import.meta.url),'utf8')).cases;
for(const c of cases){
  const mol=wasm.parse_smiles(c.smiles);
  try{
    const legacy=mol.formula(),before=wasm.to_moljson(mol);
    assert.equal(mol.formula_with_isotopes(),c.expected_formula,c.id);
    assert.equal(mol.formula(),legacy,c.id);
    assert.equal(wasm.to_moljson(mol),before,c.id);
  }finally{mol.free();}
}
console.log('Isotope-aware total formula: 47 authored compositions, preserved legacy output and unchanged graphs passed');
