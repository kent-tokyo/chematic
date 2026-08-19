// Tests for the 5 js_sys::Float64Array/Uint32Array-returning "typed-array
// accessor" functions added to crates/chematic-wasm/src/format_io.rs
// alongside their existing JSON-returning siblings (cube_grid_json/
// opendx_grid_json/lammps_dump_frame_to_json_str/
// lammps_dump_cartesian_positions_json -- all still present, unchanged):
// cube_values_f64, cube_shape_u32, opendx_values_f64, opendx_shape_u32,
// lammps_dump_rows_f64, lammps_dump_cartesian_positions_f64.
//
// cube_values_f64/cube_shape_u32/opendx_values_f64/opendx_shape_u32/
// lammps_dump_cartesian_positions_f64 are also cross-checked against their
// JSON siblings as part of the 4-format parity fixtures in
// format_parity.test.mjs -- this file covers what that one doesn't:
// lammps_dump_rows_f64 (no JSON-sibling parity fixture of its own) and the
// lammps_dump_cartesian_positions_f64 error-on-unresolvable-frame case
// (the JSON sibling returns `null` there; a Float64Array has no `null`, so
// this is a real, disclosed API-shape difference, not an oversight -- see
// both functions' Rust doc comments).
//
// Picked up automatically by CI's `for f in crates/chematic-wasm/tests/*.test.mjs`
// loop (.github/workflows/ci.yml, "Test (WASM)" job) -- no separate wiring
// needed. Run manually after building the Node target:
//
//   wasm-pack build crates/chematic-wasm --target nodejs --out-dir pkg-node --release
//   node crates/chematic-wasm/tests/typed_array_accessors.test.mjs

import assert from "node:assert/strict";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(__dirname, "..", "..", "..");

const wasm = await import(path.join(repoRoot, "crates/chematic-wasm/pkg-node/chematic_wasm.js"));
const {
  lammps_dump_frame_to_json_str,
  lammps_dump_rows_f64,
  lammps_dump_cartesian_positions_json,
  lammps_dump_cartesian_positions_f64,
} = wasm;

// Same LAMMPS dump text as crates/chematic-py/tests/test_new_formats.py's
// LAMMPS_DUMP_FIXTURE and chematic-mol's lammps_dump.rs orthogonal_frame
// (in spirit -- a plain orthogonal x/y/z frame, one atom).
const LAMMPS_DUMP_TEXT =
  "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n" +
  "0.0 10.0\n0.0 10.0\n0.0 10.0\nITEM: ATOMS id type x y z\n1 1 5.0 5.0 5.0\n";

// ---------------------------------------------------------------------------
// lammps_dump_rows_f64: flattens `rows` (one row per atom, column_names.length
// values per row), row-major, same data as the JSON sibling's "rows" field.
// ---------------------------------------------------------------------------

{
  const frameJson = lammps_dump_frame_to_json_str(LAMMPS_DUMP_TEXT);
  const frame = JSON.parse(frameJson);
  assert.deepEqual(frame.column_names, ["id", "type", "x", "y", "z"]);
  assert.deepEqual(frame.rows, [[1, 1, 5.0, 5.0, 5.0]]);

  const flat = lammps_dump_rows_f64(frameJson);
  assert.ok(flat instanceof Float64Array);
  const expected = frame.rows.flat();
  assert.deepEqual(Array.from(flat), expected);
  // Row length is column_names.length, which the caller already has from
  // the JSON-based function -- documented behavior, not a missing API.
  assert.equal(flat.length, frame.rows.length * frame.column_names.length);
}

// Multi-atom frame, to make sure row-major flattening order (not
// column-major) is what's actually produced.
{
  const multiAtomText =
    "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n" +
    "0.0 10.0\n0.0 20.0\n0.0 30.0\nITEM: ATOMS id x y z\n1 1.0 2.0 3.0\n2 4.0 5.0 6.0\n";
  const frameJson = lammps_dump_frame_to_json_str(multiAtomText);
  const flat = lammps_dump_rows_f64(frameJson);
  // Row-major: atom 0's columns (id,x,y,z)=(1,1,2,3), then atom 1's.
  assert.deepEqual(Array.from(flat), [1, 1.0, 2.0, 3.0, 2, 4.0, 5.0, 6.0]);
}

// ---------------------------------------------------------------------------
// lammps_dump_cartesian_positions_f64: Err (not null) when unresolvable.
// ---------------------------------------------------------------------------

{
  const unwrappedOnlyFrameJson = JSON.stringify({
    timestep: 1000,
    num_atoms: 2,
    box_bounds: { lo: [0.0, 0.0, 0.0], hi: [10.0, 20.0, 30.0], tilt: null },
    boundary_flags: ["pp", "pp", "pp"],
    column_names: ["id", "xu", "yu", "zu"],
    rows: [
      [1.0, 100.0, 200.0, 300.0],
      [2.0, 1.0, 2.0, 3.0],
    ],
  });

  // JSON sibling: null, not an error, not [].
  const jsonResult = JSON.parse(lammps_dump_cartesian_positions_json(unwrappedOnlyFrameJson));
  assert.equal(jsonResult, null);

  // Typed-array sibling: throws, since a Float64Array can't represent null.
  assert.throws(
    () => lammps_dump_cartesian_positions_f64(unwrappedOnlyFrameJson),
    (err) => {
      const message = typeof err === "string" ? err : String(err);
      assert.match(message, /no recognized coordinate columns/);
      return true;
    },
  );
}

console.log("typed_array_accessors.test.mjs: all assertions passed");
