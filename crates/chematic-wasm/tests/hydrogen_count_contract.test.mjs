import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import {pathToFileURL} from "node:url";
const wasm=await import(process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js",import.meta.url).href);
const rows=readFileSync(new URL("../../../validation/openbabel_issues/hydrogen-counts/cases.tsv",import.meta.url),"utf8").trimEnd().split("\n").filter(v=>!v.startsWith("#"));
for(const row of rows){
  const [smiles,count]=row.split("\t"),m=wasm.parse_smiles(smiles);
  assert.equal(wasm.num_hydrogens(m),Number(count),smiles);
  const again=wasm.parse_smiles(m.canonical_smiles());
  assert.equal(wasm.num_hydrogens(again),Number(count),smiles);
  again.free();m.free();
}
for(const smiles of ["C","[CH4]","[NH4+]","[OH2].[OH2]","CCO","c1ccccc1","[nH]1cccc1","N[C@@H](C)O","[2H]O[3H]","[H-][Al+3]([H-])[H-]"]){
  const m=wasm.parse_smiles(smiles),expanded=wasm.add_hydrogens(m),folded=wasm.remove_hydrogens(expanded);
  for(const other of [expanded,folded]){
    assert.equal(wasm.num_hydrogens(other),wasm.num_hydrogens(m),smiles);
    assert.equal(other.formula(),m.formula(),smiles);
    other.free();
  }
  m.free();
}
for(const [smiles,expected] of [["[NH4+]",[4]],["[H][O][H]",[0,2,0]],["[2H]O[3H]",[0,2,0]],["[Al][H-][Al]",[1,0,1]],["[H]",[0]],["[HH]",[1]],["[*H2]",[2]]]){
  const m=wasm.parse_smiles(smiles);
  for(let i=0;i<expected.length;i++)assert.equal(JSON.parse(wasm.get_atom_info(m,i)).totalHydrogens,expected[i],`${smiles}:${i}`);
  assert.equal(wasm.get_atom_info(m,999),"null");m.free();
}
console.log(`Hydrogen counts: ${rows.length} original cases, expand/fold, attached-H atom info and invalid indices passed`);
