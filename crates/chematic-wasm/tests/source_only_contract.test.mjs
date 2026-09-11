import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "..", "..", "..");
const fixture = JSON.parse(
  readFileSync(path.join(root, "validation/cross_binding_contract.json"), "utf8"),
);
const wasm = await import(path.join(root, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));

const cjson = fixture.cjson_contract;
const cjsonMol = wasm.mol_from_cjson(cjson.input);
assert.equal(cjsonMol.atom_count(), cjson.expected.atom_count);
assert.equal(cjsonMol.bond_count(), cjson.expected.bond_count);
const cjsonRoundtrip = wasm.mol_from_cjson(
  wasm.convert_common_format(cjson.input, "cjson", "cjson"),
);
assert.equal(cjsonRoundtrip.atom_count(), cjson.expected.atom_count);
cjsonMol.free();
cjsonRoundtrip.free();
assert.throws(() => wasm.mol_from_cjson(cjson.malformed_input));

const strictPdb = wasm.mol_from_pdb_strict(fixture.pdb_strict_contract.input);
assert.equal(strictPdb.atom_count(), fixture.pdb_strict_contract.expected.atom_count);
strictPdb.free();
assert.throws(() => wasm.mol_from_pdb_strict(fixture.pdb_strict_contract.malformed_input));

const strictPdbFields = [
  [6, 11, " nope", "serial"],
  [22, 26, "nope", "residue sequence"],
  [30, 38, "   nope", "x"],
  [38, 46, "   nope", "y"],
  [46, 54, "   nope", "z"],
];
for (const [start, end, replacement, field] of strictPdbFields) {
  const malformedField =
    fixture.pdb_strict_contract.input.slice(0, start) +
    replacement +
    fixture.pdb_strict_contract.input.slice(end);
  assert.throws(
    () => wasm.mol_from_pdb_strict(malformedField),
    new RegExp(`invalid PDB ${field} field`),
  );
}

const pdbqt = wasm.mol_from_pdbqt(fixture.pdbqt_contract.input);
assert.equal(pdbqt.atom_count(), fixture.pdbqt_contract.expected.atom_count);
pdbqt.free();
assert.throws(() => wasm.mol_from_pdbqt(fixture.pdbqt_contract.malformed_input));

const molecule = wasm.parse_smiles("CC");
for (const smarts of fixture.smarts_validity_contract.accepted) {
  assert.doesNotThrow(() => wasm.smarts_match_atoms(smarts, molecule), smarts);
}
for (const smarts of fixture.smarts_validity_contract.rejected) {
  assert.throws(() => wasm.smarts_match_atoms(smarts, molecule), smarts);
}
molecule.free();

for (const testCase of fixture.reaction_smarts_contract.cases) {
  assert.equal(
    wasm.reaction_smarts_match(testCase.smarts, testCase.reaction),
    testCase.matches,
    `${testCase.smarts} / ${testCase.reaction}`,
  );
}

console.log("source-only contract promotion: ok");
