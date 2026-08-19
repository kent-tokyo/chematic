//! Cross-language parity fixtures for the v0.18.0 "Binding Quality Pack".
//!
//! Goal: prove the Rust source of truth, the Python bindings
//! (`crates/chematic-py/tests/test_new_formats.py`), and the WASM bindings
//! (`crates/chematic-wasm/tests/format_parity.test.mjs`) all agree on the
//! SAME small fixture's meaning for 4 formats (Gaussian Cube, OpenDX,
//! mmCIF, LAMMPS dump) -- not full-format-spec coverage, just enough to
//! catch a binding silently doing its own unit conversion, coordinate
//! reordering, or lossy write that diverges from this crate.
//!
//! This file is independent ground truth: every expected value here is
//! either computed directly from the fixture text (Cube/OpenDX/mmCIF) or
//! hand-derived in a doc comment (the LAMMPS triclinic transform, shared
//! with `chematic_mol::lammps_dump`'s own
//! `scaled_coordinates_triclinic_hand_computed` test) -- never copied from
//! what another language binding happened to output. The Python/WASM test
//! files hardcode these same numbers independently rather than importing
//! this file, so a bug shared by all 3 entry points (unlikely, since each
//! calls into this same `chematic_mol` crate, but the whole point of this
//! file existing at all) is still visible as "3 tests independently
//! assert the same known-correct number", not just "3 tests call the same
//! Rust function".

use chematic_mol::{GridUnits, LammpsBox, LammpsDumpFrame, parse_cube, parse_mmcif, parse_opendx};

// ---------------------------------------------------------------------------
// Gaussian Cube
// ---------------------------------------------------------------------------

/// Shared verbatim with `crates/chematic-wasm/src/format_io.rs`'s
/// `CUBE_2X2X2` and `crates/chematic-py/tests/test_new_formats.py`'s
/// `CUBE_FIXTURE`.
const CUBE_FIXTURE: &str = "Water density\n\
Generated for chematic tests\n\
1    0.000000    0.000000    0.000000\n\
2    1.000000    0.000000    0.000000\n\
2    0.000000    1.000000    0.000000\n\
2    0.000000    0.000000    1.000000\n\
8    8.000000    0.500000    0.500000    0.500000\n\
0.0 1.0 2.0 3.0\n\
4.0 5.0 6.0 7.0\n";

#[test]
fn cube_fixture_summary() {
    let grid = parse_cube(CUBE_FIXTURE).unwrap();
    assert_eq!(grid.shape, [2, 2, 2]);
    assert_eq!(grid.units, GridUnits::Bohr);
    assert_eq!(grid.point_count().unwrap(), 8);
    assert_eq!(grid.atoms.len(), 1);
    // Format-specific check: first/last `values` entries (catches a
    // reversed/transposed flatten order immediately), plus two interior
    // values only a correctly k-fastest-ordered flatten reproduces (a
    // first/last-only check alone survives some transposes on this
    // fixture, since 0.0 and 7.0 sit at the same flat position under a
    // few wrong orderings too).
    assert_eq!(grid.values.first().copied(), Some(0.0));
    assert_eq!(grid.values.last().copied(), Some(7.0));
    assert_eq!(grid.get(0, 1, 0), Some(2.0));
    assert_eq!(grid.get(1, 0, 0), Some(4.0));
}

// ---------------------------------------------------------------------------
// OpenDX
// ---------------------------------------------------------------------------

/// Shared verbatim with `crates/chematic-wasm/src/format_io.rs`'s
/// `OPENDX_2X2X2` and `crates/chematic-py/tests/test_new_formats.py`'s
/// `OPENDX_FIXTURE`.
const OPENDX_FIXTURE: &str = "object 1 class gridpositions counts 2 2 2\n\
origin -1.0 -1.0 -1.0\n\
delta 0.5 0.0 0.0\n\
delta 0.0 0.5 0.0\n\
delta 0.0 0.0 0.5\n\
object 2 class gridconnections counts 2 2 2\n\
object 3 class array type double rank 0 items 8 data follows\n\
0.0 1.0 2.0\n\
3.0 4.0 5.0\n\
6.0 7.0\n\
attribute \"dep\" string \"positions\"\n\
object \"regular positions regular connections\" class field\n\
component \"positions\" value 1\n\
component \"connections\" value 2\n\
component \"data\" value 3\n";

#[test]
fn opendx_fixture_summary() {
    let grid = parse_opendx(OPENDX_FIXTURE).unwrap();
    assert_eq!(grid.shape, [2, 2, 2]);
    assert_eq!(grid.units, GridUnits::Angstrom);
    assert_eq!(grid.point_count().unwrap(), 8);
    assert_eq!(grid.atoms.len(), 0);
    // Format-specific check: `axes` values -- a units-conversion bug
    // (Bohr vs Angstrom, ~1.89x) would show up directly here as a wrong
    // diagonal value.
    assert_eq!(
        grid.axes,
        [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.5]]
    );
    assert_eq!(grid.get(0, 1, 0), Some(2.0));
    assert_eq!(grid.get(1, 0, 0), Some(4.0));
}

// ---------------------------------------------------------------------------
// mmCIF
// ---------------------------------------------------------------------------

/// Shared verbatim with `crates/chematic-py/tests/test_new_formats.py`'s
/// `MMCIF_FIXTURE`. Deliberately has NO `_atom_site.occupancy` column, so
/// this fixture's occupancy check exercises the *default* path (1.0), not
/// just a column-value passthrough -- catching both a field-mapping bug
/// (wrong column read for `x`) and a defaulting bug in one assertion pair.
const MMCIF_FIXTURE: &str = "data_TEST\n\
loop_\n\
_atom_site.group_PDB\n\
_atom_site.id\n\
_atom_site.type_symbol\n\
_atom_site.label_atom_id\n\
_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n\
_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n\
ATOM   1  O  O1  HOH A 1.000 2.000 3.000\n\
ATOM   2  H  H1  HOH A 1.500 2.500 3.500\n";

#[test]
fn mmcif_fixture_summary() {
    let result = parse_mmcif(MMCIF_FIXTURE).unwrap();
    assert_eq!(result.atoms.len(), 2);
    // Format-specific check: one specific `_atom_site` field per atom --
    // the first atom's coordinate (catches a field-mapping/column-shift
    // bug) and its occupancy (catches a wrong default; this fixture has no
    // occupancy column at all, so the spec-mandated default is 1.0).
    assert_eq!(result.atoms[0].element.symbol(), "O");
    assert!((result.atoms[0].x - 1.0).abs() < 1e-9);
    assert!((result.atoms[0].y - 2.0).abs() < 1e-9);
    assert!((result.atoms[0].z - 3.0).abs() < 1e-9);
    assert!((result.atoms[0].occupancy - 1.0).abs() < 1e-9);
}

// ---------------------------------------------------------------------------
// LAMMPS dump (triclinic)
// ---------------------------------------------------------------------------

/// The exact triclinic fixture from `chematic_mol::lammps_dump`'s own
/// `triclinic_frame()` test helper and
/// `scaled_coordinates_triclinic_hand_computed` test -- reused verbatim
/// (not re-derived) here and in the Python/WASM parity tests, since this
/// is the single highest-value parity check in this pack: PR #343 (the
/// prior WASM-bindings PR) had a real bug in exactly this
/// bound-box-vs-true-box / triclinic-shear-term resolution, independently
/// caught and fixed during that PR's review. Constructed directly (not
/// parsed from dump-file text) so this fixture can never accidentally pick
/// up a "BOX BOUNDS" bound-box value where a true-box value was intended
/// (`parse_lammps_dump_frame` applies `box_bounds_to_true` internally
/// before a frame is ever built; a hand-written dump-text fixture would
/// have to get that conversion right by hand too, which is exactly the
/// kind of step this fixture is designed to make impossible to get wrong).
fn triclinic_frame() -> LammpsDumpFrame {
    LammpsDumpFrame {
        timestep: 2000,
        num_atoms: 1,
        box_bounds: LammpsBox {
            lo: [0.0, 0.0, 0.0],
            hi: [10.0, 10.0, 10.0],
            tilt: Some([2.0, 1.0, 0.5]),
        },
        boundary_flags: ["pp".to_string(), "ff".to_string(), "ss".to_string()],
        column_names: vec![
            "id".to_string(),
            "xs".to_string(),
            "ys".to_string(),
            "zs".to_string(),
        ],
        rows: vec![vec![1.0, 0.5, 0.5, 0.5]],
    }
}

#[test]
fn lammps_dump_triclinic_fixture_summary() {
    let frame = triclinic_frame();
    assert_eq!(frame.num_atoms, 1);
    // Format-specific check: cartesian_positions() for the one atom, in a
    // fixture with a genuine triclinic tilt (all 3 shear terms nonzero).
    // Hand-computed:
    //   x = xlo + xs*(xhi-xlo) + ys*xy + zs*xz = 0 + 0.5*10 + 0.5*2 + 0.5*1 = 6.5
    //   y = ylo + ys*(yhi-ylo) + zs*yz         = 0 + 0.5*10 + 0.5*0.5       = 5.25
    //   z = zlo + zs*(zhi-zlo)                 = 0 + 0.5*10                = 5.0
    let pos = frame.cartesian_positions().unwrap();
    assert_eq!(pos.len(), 1);
    assert!((pos[0][0] - 6.5).abs() < 1e-9);
    assert!((pos[0][1] - 5.25).abs() < 1e-9);
    assert!((pos[0][2] - 5.0).abs() < 1e-9);
}
