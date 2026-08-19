// Cross-language parity fixtures for the v0.18.0 "Binding Quality Pack".
//
// Same 4 fixtures (verbatim text/values) and the same hardcoded expected
// numbers as crates/chematic-mol/tests/format_binding_parity.rs (Rust) and
// crates/chematic-py/tests/test_new_formats.py's `test_parity_*` tests
// (Python) -- each computed once, independently, not by trusting another
// binding's output. See the Rust file's module doc comment for why this
// independently-hardcoded approach still proves 3-way parity without a
// shared fixture file. Goal: catch a binding silently doing its own unit
// conversion, coordinate reordering, or lossy write that diverges from
// chematic-mol, the source of truth.
//
// Both the JSON-returning functions (cube_grid_json/opendx_grid_json/
// mmcif_to_json/lammps_dump_cartesian_positions_json) and their
// typed-array siblings added in this same PR (cube_values_f64/
// cube_shape_u32/opendx_values_f64/opendx_shape_u32/
// lammps_dump_cartesian_positions_f64) are checked here, cross-checked
// against each other where both exist.
//
// Picked up automatically by CI's `for f in crates/chematic-wasm/tests/*.test.mjs`
// loop (.github/workflows/ci.yml, "Test (WASM)" job) -- no separate wiring
// needed. Run manually after building the Node target:
//
//   wasm-pack build crates/chematic-wasm --target nodejs --out-dir pkg-node --release
//   node crates/chematic-wasm/tests/format_parity.test.mjs

import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");

const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const {
  cube_grid_json,
  cube_values_f64,
  cube_shape_u32,
  opendx_grid_json,
  opendx_values_f64,
  opendx_shape_u32,
  mmcif_to_json,
  lammps_dump_cartesian_positions_json,
  lammps_dump_cartesian_positions_f64,
} = wasm;

// ---------------------------------------------------------------------------
// Gaussian Cube -- shared verbatim with chematic-mol's/chematic-py's
// CUBE_FIXTURE/CUBE_2X2X2.
// ---------------------------------------------------------------------------

const CUBE_FIXTURE = `Water density
Generated for chematic tests
1    0.000000    0.000000    0.000000
2    1.000000    0.000000    0.000000
2    0.000000    1.000000    0.000000
2    0.000000    0.000000    1.000000
8    8.000000    0.500000    0.500000    0.500000
0.0 1.0 2.0 3.0
4.0 5.0 6.0 7.0
`;

{
  const grid = JSON.parse(cube_grid_json(CUBE_FIXTURE));
  assert.deepEqual(grid.shape, [2, 2, 2]);
  assert.equal(grid.units, "bohr");
  assert.equal(grid.values.length, 8);
  assert.equal(grid.atoms.length, 1);
  // First/last (reversed-flatten tripwire) plus two interior values only a
  // correctly k-fastest-ordered flatten reproduces.
  assert.equal(grid.values[0], 0.0);
  assert.equal(grid.values[grid.values.length - 1], 7.0);
  assert.equal(grid.values[2], 2.0); // get(0,1,0)
  assert.equal(grid.values[4], 4.0); // get(1,0,0)

  // Typed-array sibling must carry the exact same values/shape.
  const values = cube_values_f64(CUBE_FIXTURE);
  assert.ok(values instanceof Float64Array);
  assert.deepEqual(Array.from(values), grid.values);
  const shape = cube_shape_u32(CUBE_FIXTURE);
  assert.ok(shape instanceof Uint32Array);
  assert.deepEqual(Array.from(shape), grid.shape);
}

// ---------------------------------------------------------------------------
// OpenDX -- shared verbatim with chematic-mol's/chematic-py's
// OPENDX_FIXTURE/OPENDX_2X2X2.
// ---------------------------------------------------------------------------

const OPENDX_FIXTURE = `object 1 class gridpositions counts 2 2 2
origin -1.0 -1.0 -1.0
delta 0.5 0.0 0.0
delta 0.0 0.5 0.0
delta 0.0 0.0 0.5
object 2 class gridconnections counts 2 2 2
object 3 class array type double rank 0 items 8 data follows
0.0 1.0 2.0
3.0 4.0 5.0
6.0 7.0
attribute "dep" string "positions"
object "regular positions regular connections" class field
component "positions" value 1
component "connections" value 2
component "data" value 3
`;

{
  const grid = JSON.parse(opendx_grid_json(OPENDX_FIXTURE));
  assert.deepEqual(grid.shape, [2, 2, 2]);
  assert.equal(grid.units, "angstrom");
  assert.equal(grid.values.length, 8);
  assert.equal(grid.atoms.length, 0);
  // axes: a Bohr<->Angstrom unit-conversion bug would show up here directly.
  assert.deepEqual(grid.axes, [
    [0.5, 0.0, 0.0],
    [0.0, 0.5, 0.0],
    [0.0, 0.0, 0.5],
  ]);
  assert.equal(grid.values[2], 2.0); // get(0,1,0)
  assert.equal(grid.values[4], 4.0); // get(1,0,0)

  const values = opendx_values_f64(OPENDX_FIXTURE);
  assert.ok(values instanceof Float64Array);
  assert.deepEqual(Array.from(values), grid.values);
  const shape = opendx_shape_u32(OPENDX_FIXTURE);
  assert.ok(shape instanceof Uint32Array);
  assert.deepEqual(Array.from(shape), grid.shape);
}

// ---------------------------------------------------------------------------
// mmCIF -- shared verbatim with chematic-mol's/chematic-py's MMCIF_FIXTURE.
// ---------------------------------------------------------------------------

const MMCIF_FIXTURE = `data_TEST
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

{
  const result = JSON.parse(mmcif_to_json(MMCIF_FIXTURE));
  assert.equal(result.atoms.length, 2);
  // No occupancy column at all in this fixture -- exercises both the
  // coordinate field-mapping and the spec-mandated occupancy default (1.0).
  assert.equal(result.atoms[0].element, "O");
  assert.equal(result.atoms[0].x, 1.0);
  assert.equal(result.atoms[0].y, 2.0);
  assert.equal(result.atoms[0].z, 3.0);
  assert.equal(result.atoms[0].occupancy, 1.0);
}

// ---------------------------------------------------------------------------
// LAMMPS dump (triclinic) -- exact same box/tilt/xs values as
// chematic_mol::lammps_dump's own triclinic_frame() test fixture, passed
// as an already-resolved TRUE-box frame JSON object (never derived from
// hand-written dump-file "ITEM: BOX BOUNDS" text, which carries the BOUND
// box, not the true box -- see the Rust parity file's fixture comment).
// This is the single highest-value parity check in this pack: PR #343 had
// a real bug in exactly this triclinic-shear-term resolution, caught and
// fixed during that PR's own review.
// ---------------------------------------------------------------------------

const LAMMPS_TRICLINIC_FRAME_JSON = JSON.stringify({
  timestep: 2000,
  num_atoms: 1,
  box_bounds: { lo: [0.0, 0.0, 0.0], hi: [10.0, 10.0, 10.0], tilt: [2.0, 1.0, 0.5] },
  boundary_flags: ["pp", "ff", "ss"],
  column_names: ["id", "xs", "ys", "zs"],
  rows: [[1.0, 0.5, 0.5, 0.5]],
});

{
  // Hand-computed (see the Rust parity file for the derivation):
  //   x = 0 + 0.5*10 + 0.5*2 + 0.5*1 = 6.5
  //   y = 0 + 0.5*10 + 0.5*0.5       = 5.25
  //   z = 0 + 0.5*10                 = 5.0
  const positionsJson = JSON.parse(
    lammps_dump_cartesian_positions_json(LAMMPS_TRICLINIC_FRAME_JSON),
  );
  assert.equal(positionsJson.length, 1);
  assert.ok(Math.abs(positionsJson[0][0] - 6.5) < 1e-9);
  assert.ok(Math.abs(positionsJson[0][1] - 5.25) < 1e-9);
  assert.ok(Math.abs(positionsJson[0][2] - 5.0) < 1e-9);

  const flat = lammps_dump_cartesian_positions_f64(LAMMPS_TRICLINIC_FRAME_JSON);
  assert.ok(flat instanceof Float64Array);
  assert.equal(flat.length, 3);
  assert.ok(Math.abs(flat[0] - 6.5) < 1e-9);
  assert.ok(Math.abs(flat[1] - 5.25) < 1e-9);
  assert.ok(Math.abs(flat[2] - 5.0) < 1e-9);
}

console.log("format_parity.test.mjs: all assertions passed");
