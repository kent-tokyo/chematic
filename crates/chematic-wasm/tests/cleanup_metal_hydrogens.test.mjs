import assert from "node:assert/strict";
import { pathToFileURL } from "node:url";

const moduleUrl = process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href;
const wasm = await import(moduleUrl);

for (const input of [
  "[NH2][Fe@@H]([CH2]NCP(=O)(O)O)[C](=[O:27])O.[NaH]",
  "CCO[C](=O)[Pt@@H]1[C@H](O)CC[C@@H]2CN3CCc4cc5c(cc4[C@@H]3C[C@@H]21)OCO5",
]) {
  const source = wasm.parse_smiles(input);
  const cleaned = wasm.cleanup_mol(source);
  const repeated = wasm.cleanup_mol(cleaned);
  try {
    assert.equal(cleaned.formula(), source.formula(), input);
    assert.equal(repeated.formula(), source.formula(), input);
    const reparsed = wasm.parse_smiles(cleaned.rdkit_smiles());
    try {
      assert.equal(reparsed.formula(), source.formula(), input);
    } finally {
      reparsed.free();
    }
  } finally {
    repeated.free();
    cleaned.free();
    source.free();
  }
}

for (const input of [
  "C/N=C(\\S)N/N=C/C(C)=[N]/[SnH]=[CH]/N",
  "C/N=C(\\S)N/N=C/C(C)=[N]/[Sn]([H])=[CH]/N",
]) {
  const ambiguous = wasm.parse_smiles(input);
  try {
    assert.throws(() => wasm.cleanup_mol(ambiguous), /cleanup_metal_hydrogen/);
  } finally {
    ambiguous.free();
  }
}

console.log("Cleanup preserves metal hydrogen formula across the WASM API");
