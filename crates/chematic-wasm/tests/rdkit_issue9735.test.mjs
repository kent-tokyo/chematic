import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const wasm=await import(process.env.CHEMATIC_WASM_MODULE?pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href:new URL('../pkg-node/chematic_wasm.js',import.meta.url).href);
const packet=JSON.parse(readFileSync(new URL('../../../validation/rdkit-9735-2026-10-24/cases.json',import.meta.url),'utf8'));
for(const c of packet.cases){
  const mol=wasm.parse_smiles(c.smiles);
  try{
    if(c.expected==='valence_error')assert.throws(()=>mol.rdkit_smiles(),/greater than permitted/,c.id);
    else{
      const first=mol.rdkit_smiles();const back=wasm.parse_smiles(first);
      try{assert.equal(back.rdkit_smiles(),first,c.id);}finally{back.free();}
    }
  }finally{mol.free();}
}
console.log('RDKit #9735: 13 strict valence refusals, 24 accepted/stable controls passed');
