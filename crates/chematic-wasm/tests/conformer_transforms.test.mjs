import assert from "node:assert/strict";
import {
  get_bond_angle_from_coords_json,
  get_bond_length_from_coords_json,
  parse_smiles,
  set_bond_angle_from_coords_json,
  set_bond_length_from_coords_json,
  set_dihedral_from_coords_json,
} from "../pkg-node/chematic_wasm.js";

const butane = parse_smiles("CCCC");
try {
  const bondInput = JSON.stringify([
    [0, 0, 0],
    [1, 0, 0],
    [2, 1, 0],
    [3, 1, 1],
  ]);
  const bondOutput = set_bond_length_from_coords_json(butane, bondInput, 0, 1, 2.0);
  assert.deepEqual(JSON.parse(bondOutput), [
    [0, 0, 0],
    [2, 0, 0],
    [3, 1, 0],
    [4, 1, 1],
  ]);
  assert.equal(get_bond_length_from_coords_json(butane, bondOutput, 0, 1), 2.0);

  const angleInput = JSON.stringify([
    [1, 0, 0],
    [0, 0, 0],
    [0, 1, 0],
    [0, 2, 0],
  ]);
  const angleOutput = set_bond_angle_from_coords_json(butane, angleInput, 0, 1, 2, 60.0);
  assert.ok(
    Math.abs(get_bond_angle_from_coords_json(butane, angleOutput, 0, 1, 2) - 60.0) < 1e-10,
  );
} finally {
  butane.free();
}

const ring = parse_smiles("C1CCCCC1");
try {
  const ringCoords = JSON.stringify([
    [1, 0, 0],
    [0.5, 0.866, 0],
    [-0.5, 0.866, 0],
    [-1, 0, 0],
    [-0.5, -0.866, 0],
    [0.5, -0.866, 0],
  ]);
  assert.throws(
    () => set_bond_length_from_coords_json(ring, ringCoords, 0, 1, 2.0),
    /must not belong to a ring/,
  );
  assert.throws(
    () => set_bond_angle_from_coords_json(ring, ringCoords, 0, 1, 2, 90.0),
    /both bonds/,
  );
  assert.throws(
    () => set_dihedral_from_coords_json(ring, ringCoords, 0, 1, 2, 3, 60.0),
    /must not belong to a ring/,
  );
} finally {
  ring.free();
}

console.log("conformer_transforms.test.mjs: all assertions passed");
