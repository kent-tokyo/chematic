import assert from "node:assert/strict";
import test from "node:test";

import * as wasm from "../pkg-node/chematic_wasm.js";

const UNSTABLE_PHOSPHORUS = "N[P@]1(Cl)=NP(N2CC2)(N2CC2)=N[P@](N)(Cl)=N1";

test("accurate CIP keeps representation-unstable phosphorus fail-closed", () => {
  const mol = wasm.parse_smiles(UNSTABLE_PHOSPHORUS);
  try {
    assert.deepEqual(JSON.parse(wasm.cip_assignments_accurate_json(mol)), []);
    assert.deepEqual(JSON.parse(wasm.cip_unresolved_json(mol)), [
      { atomIdx: 1, reason: "oracleUnstable" },
      { atomIdx: 12, reason: "oracleUnstable" },
    ]);
  } finally {
    mol.free();
  }
});
