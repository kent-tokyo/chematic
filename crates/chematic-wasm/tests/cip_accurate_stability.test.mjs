import assert from "node:assert/strict";
import test from "node:test";

import * as wasm from "../pkg-node/chematic_wasm.js";

const CYCLOPHOSPHAZENE = "N[P@]1(Cl)=NP(N2CC2)(N2CC2)=N[P@](N)(Cl)=N1";

test("accurate CIP labels cyclophosphazene phosphorus as RDKit, marked spelling-dependent", () => {
  const mol = wasm.parse_smiles(CYCLOPHOSPHAZENE);
  try {
    // RDKit 2026.03.6 rdCIPLabeler: atom 1 S, atom 12 R.
    assert.deepEqual(JSON.parse(wasm.cip_assignments_accurate_json(mol)), [
      { atomIdx: 1, cipCode: "S", kekuleDependent: true },
      { atomIdx: 12, cipCode: "R", kekuleDependent: true },
    ]);
    assert.deepEqual(JSON.parse(wasm.cip_unresolved_json(mol)), []);
  } finally {
    mol.free();
  }
});
