import assert from "node:assert/strict";
import { pathToFileURL } from "node:url";

const moduleUrl = process.env.CHEMATIC_WASM_MODULE
  ? pathToFileURL(process.env.CHEMATIC_WASM_MODULE).href
  : new URL("../pkg-node/chematic_wasm.js", import.meta.url).href;
const wasm = await import(moduleUrl);

let tinOutput;
for (const input of [
  "[NH2][Fe@@H]([CH2]NCP(=O)(O)O)[C](=[O:27])O.[NaH]",
  "CCO[C](=O)[Pt@@H]1[C@H](O)CC[C@@H]2CN3CCc4cc5c(cc4[C@@H]3C[C@@H]21)OCO5",
  "C[N]([SnH])C",
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

for (const input of ["O->[Mg](<-O)(<-O)<-O", "C(=[O][Li])O"]) {
  const source = wasm.parse_smiles(input);
  const expanded = wasm.add_hydrogens(source);
  try {
    for (const variant of [source, expanded]) {
      assert.throws(() => wasm.cleanup_mol(variant), /cleanup_unstable_metal/);
    }
  } finally {
    expanded.free();
    source.free();
  }
}

{
  const compact = wasm.parse_smiles("Cc1c[c]([Pt])c(Oc2ccc(O)cc2CC(=O)O)s(Cl)c1");
  const explicit = wasm.parse_smiles("C(C1=C([H])[C]([Pt])=C(Oc2c([H])c([H])c(O[H])c([H])c2C(C(=O)O[H])([H])[H])S(Cl)=C1[H])([H])([H])[H]");
  const compactResult = wasm.cleanup_mol(compact);
  const explicitResult = wasm.cleanup_mol(explicit);
  try {
    assert.equal(compactResult.formula(), compact.formula());
    assert.equal(explicitResult.formula(), compact.formula());
    assert.equal(compactResult.rdkit_smiles(), explicitResult.rdkit_smiles());
  } finally {
    compactResult.free();
    explicitResult.free();
    compact.free();
    explicit.free();
  }
}

{
  const compact = wasm.parse_smiles("CN(C)C1=[C]([Cu][C@@H](CC2=CC=C(NC(=O)C3=C(Cl)C=NC=C3Cl)C=C2)C(=O)O)C(=O)C1=O");
  const explicit = wasm.parse_smiles("C(N(C([H])([H])[H])c1[c]([Cu][C@@](C(c2c([H])c([H])c(N(C(=O)c3c(Cl)c([H])nc([H])c3Cl)[H])c([H])c2[H])([H])[H])(C(=O)O[H])[H])c(=O)c1=O)([H])([H])[H]");
  const compactResult = wasm.cleanup_mol(compact);
  const explicitResult = wasm.cleanup_mol(explicit);
  try {
    assert.equal(compactResult.formula(), compact.formula());
    assert.equal(explicitResult.formula(), compact.formula());
    assert.equal(compactResult.rdkit_smiles(), explicitResult.rdkit_smiles());
  } finally {
    compactResult.free();
    explicitResult.free();
    compact.free();
    explicit.free();
  }
}

for (const input of [
  "O=C1C[CH]([Pt]2=[CH]C=C(Cl)C=[CH]2)CN1",
  "O=C(O)C1=[C]2C=C(F)C=[CH][Fe]2=[N]C(C2=CC=C(C3=CC=CC=C3)C=C2)=C1",
  "C=[C](C1=CC=CC=C1OCCC1=CC=CC=C1)[Fe]1[CH]=CN=[CH]1",
  "[Na]12=[C](C=C(F)C=[CH]1)N(C1=C[C+]=C(C=C1)Cl)C1C(=[N]2)C(=O)N(C)C(=O)N=1",
]) {
  const source = wasm.parse_smiles(input);
  const expanded = wasm.add_hydrogens(source);
  const compactResult = wasm.cleanup_mol(source);
  const expandedResult = wasm.cleanup_mol(expanded);
  try {
    assert.equal(expanded.formula(), source.formula(), input);
    assert.equal(expandedResult.formula(), source.formula(), input);
    assert.equal(expandedResult.rdkit_smiles(), compactResult.rdkit_smiles(), input);
  } finally {
    expandedResult.free();
    compactResult.free();
    expanded.free();
    source.free();
  }
}

for (const input of [
  "C/N=C(\\S)N/N=C/C(C)=[N]/[SnH]=[CH]/N",
  "C/N=C(\\S)N/N=C/C(C)=[N]/[Sn]([H])=[CH]/N",
]) {
  const source = wasm.parse_smiles(input);
  const result = wasm.cleanup_mol(source);
  try {
    assert.equal(result.formula(), source.formula());
    if (tinOutput !== undefined) assert.equal(result.rdkit_smiles(), tinOutput);
    tinOutput = result.rdkit_smiles();
  } finally {
    result.free();
    source.free();
  }
}

const ambiguous = wasm.parse_smiles("[Na]N=[N+]=[N-]");
try {
  assert.throws(() => wasm.cleanup_mol(ambiguous), /cleanup_metal_hydrogen/);
} finally {
  ambiguous.free();
}

for (const input of [
  "O=C(O)[C@@H]1([SnH])C[C@@H](/C=C/CP(=O)(O)O)CCN1",
  "O=C(O[H])[C@]1(->[Sn][H])([H])C([H])([H])[C@](/C(=C(/C(P(=O)(O[H])O[H])([H])[H])[H])[H])([H])C([H])([H])C([H])([H])N1[H]",
]) {
  const source = wasm.parse_smiles(input);
  try {
    assert.throws(() => wasm.cleanup_mol(source), /cleanup_metal_stereo/);
  } finally {
    source.free();
  }
}

console.log("Cleanup preserves metal hydrogen formula across the WASM API");
