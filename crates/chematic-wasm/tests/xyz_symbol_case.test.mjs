import assert from "node:assert/strict";
import { pathToFileURL } from "node:url";
const url = process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href;
const wasm = await import(url);
for (const [symbol,canonical] of [["CL","Cl"],["BR","Br"],["CA","Ca"],["FE","Fe"],["NA","Na"],["OG","Og"]]) {
  const text = `1\ncase\n${symbol} 1 -2.5 1.25\n`;
  for (const reader of [wasm.mol_from_xyz, wasm.mol_from_extxyz]) {
    const mol = reader(text);
    assert.equal(mol.atom_count(), 1);
    assert.equal(wasm.to_xyz(mol).split("\n")[2].trim().split(/\s+/)[0],canonical);
    mol.free();
  }
  assert.deepEqual(JSON.parse(wasm.extxyz_frame_json(text)).coords,[[1,-2.5,1.25]]);
}
for (const token of ["cl","cL","ZZ","CL1","CLO","17","D","T","*","ＣＬ"]) {
  const text = `1\ninvalid\n${token} 0 0 0\n`;
  assert.throws(() => wasm.mol_from_xyz(text));
  assert.throws(() => wasm.mol_from_extxyz(text));
}
const ext = JSON.parse(wasm.extxyz_frame_json('1\nProperties=species:S:1:pos:R:3:name:S:1 label=CL\nCL 1 2 3 CA\n'));
assert.deepEqual(ext.properties.name, [["CA"]]);
assert.equal(ext.info.label,"CL");
assert.throws(() => wasm.parse_smiles("CL"));
console.log("XYZ species: uppercase aliases, strict invalid tokens, string metadata and SMILES boundary passed");
