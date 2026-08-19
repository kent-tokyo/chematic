// Runnable example: 3 of the v0.17.0 materials-science file-format
// bindings (Gaussian Cube, mmCIF, LAMMPS dump) via chematic-wasm, mirroring
// examples/materials_formats_quickstart.py's Python walkthrough for a
// JS/Node audience. Doubles as a CI-executed test (real assertions below,
// not just console output) -- see
// crates/chematic-wasm/tests/format_parity.test.mjs for the dedicated
// 3-language parity fixtures this file is deliberately NOT duplicating.
//
// Picked up automatically by CI's `for f in crates/chematic-wasm/tests/*.test.mjs`
// loop (.github/workflows/ci.yml, "Test (WASM)" job) -- no separate wiring
// needed.
//
// Dependencies: build the Node target first --
//
//   wasm-pack build crates/chematic-wasm --target nodejs --out-dir pkg-node --release
//   node crates/chematic-wasm/tests/materials_formats_example.test.mjs

import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");

const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const {
  mol_from_cube,
  cube_values_f64,
  cube_shape_u32,
  mol_from_mmcif,
  mmcif_to_json,
  lammps_dump_frame_to_json_str,
  lammps_dump_cartesian_positions_json,
} = wasm;

// --- Gaussian Cube: a tiny 2x2x2 density grid around one oxygen atom. ---

const CUBE_TEXT = `Water density
Generated for chematic examples
1    0.000000    0.000000    0.000000
2    1.000000    0.000000    0.000000
2    0.000000    1.000000    0.000000
2    0.000000    0.000000    1.000000
8    8.000000    0.500000    0.500000    0.500000
0.0 1.0 2.0 3.0
4.0 5.0 6.0 7.0
`;

console.log("=== Gaussian Cube ===");
{
  const mol = mol_from_cube(CUBE_TEXT);
  const shape = cube_shape_u32(CUBE_TEXT);
  const values = cube_values_f64(CUBE_TEXT);
  console.log(`atom_count=${mol.atom_count()} shape=[${shape.join(",")}] value_count=${values.length}`);
  assert.equal(mol.atom_count(), 1);
  assert.deepEqual(Array.from(shape), [2, 2, 2]);
  assert.equal(values.length, 8);
}

// --- mmCIF: 2 atoms, no bond table. ---

const MMCIF_TEXT = `data_EXAMPLE
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM   1  O  O1  HOH A 1.000 2.000 3.000
ATOM   2  H  H1  HOH A 1.500 2.500 3.500
`;

console.log("=== mmCIF ===");
{
  const mol = mol_from_mmcif(MMCIF_TEXT);
  const result = JSON.parse(mmcif_to_json(MMCIF_TEXT));
  console.log(`atom_count=${result.atoms.length}`);
  for (const atom of result.atoms) {
    console.log(`  ${atom.element} (${atom.x}, ${atom.y}, ${atom.z})`);
  }
  assert.equal(mol.atom_count(), 2);
  assert.equal(result.atoms.length, 2);
  assert.equal(result.atoms[0].element, "O");
}

// --- LAMMPS dump: one orthogonal frame, resolve real Cartesian positions. ---

const LAMMPS_DUMP_TEXT =
  "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n" +
  "0.0 10.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id x y z\n1 1.0 2.0 3.0\n2 4.0 5.0 6.0\n";

console.log("=== LAMMPS dump ===");
{
  const frameJson = lammps_dump_frame_to_json_str(LAMMPS_DUMP_TEXT);
  const frame = JSON.parse(frameJson);
  const positions = JSON.parse(lammps_dump_cartesian_positions_json(frameJson));
  console.log(`timestep=${frame.timestep} num_atoms=${frame.num_atoms}`);
  positions.forEach((p, i) => console.log(`  atom ${i}: (${p[0]}, ${p[1]}, ${p[2]})`));
  assert.equal(frame.num_atoms, 2);
  assert.deepEqual(positions, [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
  ]);
}

console.log("materials_formats_example.test.mjs: all assertions passed");
