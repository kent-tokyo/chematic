import assert from "node:assert/strict";
import {
  brenk_matches_json,
  brenk_report_json,
  parse_smiles,
} from "../pkg-node/chematic_wasm.js";

const clean = JSON.parse(brenk_report_json("CCO"));
assert.deepEqual(clean, { passes: true, alerts: [] });

const phenol = JSON.parse(brenk_report_json("Oc1ccccc1"));
assert.equal(phenol.passes, false);
assert.ok(phenol.alerts.includes("phenol"));

const mol = parse_smiles("Oc1ccccc1");
assert.equal(mol.brenk_passes(), false);
assert.deepEqual(JSON.parse(brenk_matches_json(mol)), phenol.alerts);
mol.free();

assert.throws(() => brenk_report_json("C("));

console.log("Brenk WASM API: OK");
