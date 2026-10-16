import assert from "node:assert/strict";
import { pathToFileURL } from "node:url";
const url = process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href;
const wasm = await import(url);
for (const source of ["[CH:1]1=[CH:2][C@H:3]2[O:4][C@@H:5]1[CH2:6][CH2:7]2", "[13CH3:1][C@H:2]([OH:3])[F:4]", "[F:1]/[CH:2]=[CH:3]/[Cl:4]", "[N:1]->[Cu:2]"]) {
  const mol = wasm.parse_smiles(source);
  const baseline = mol.canonical_smiles();
  const first = JSON.parse(mol.smiles_rooted_checked_json(0));
  for (let root = 0; root < first.atom_order.length; root++) {
    const result = JSON.parse(mol.smiles_rooted_checked_json(root));
    assert.equal(result.atom_order[0], root);
    assert.deepEqual([...result.atom_order].sort((a,b)=>a-b), first.atom_order.map((_,i)=>i));
    const after = wasm.parse_smiles(result.smiles);
    assert.equal(after.canonical_smiles(), baseline);
    after.free();
  }
  assert.equal(mol.canonical_smiles(), baseline);
  for (const root of [-1, 0.5, NaN, Infinity, 2**40, first.atom_order.length]) {
    assert.throws(() => mol.smiles_rooted_checked_json(root));
  }
  mol.free();
}
for (const source of ["F[Pt@SP1](Cl)(Br)I", "FC=[C@AL1]=C(Cl)Br", "C[S@](=O)CC", "C~C", "C".repeat(257)]) {
  const mol = wasm.parse_smiles(source);
  assert.throws(() => mol.smiles_rooted_checked_json(0));
  mol.free();
}
const cap = wasm.parse_smiles("C".repeat(256));
assert.equal(JSON.parse(cap.smiles_rooted_checked_json(0)).atom_order.length, 256);
cap.free();
console.log("Checked rooted SMILES: all-root stereo/identity, correspondence, input and scope contracts passed");
