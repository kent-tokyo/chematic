//! Open Babel #2995: coordinate labels must describe the values emitted.
use chematic_mol::{CifError, UnitCell, parse_cif, parse_xyz, write_cif};

const XYZ: &str = include_str!("../../../validation/openbabel_issues/cif/cartesian-no-cell.xyz");

fn close(a: (f64, f64, f64), b: (f64, f64, f64), tolerance: f64) {
    assert!(
        (a.0 - b.0).abs() < tolerance
            && (a.1 - b.1).abs() < tolerance
            && (a.2 - b.2).abs() < tolerance,
        "{a:?} != {b:?}"
    );
}

#[test]
fn openbabel_2995_no_cell_writer_emits_cartesian_and_round_trips() {
    let frame = parse_xyz(XYZ).unwrap();
    let mol = frame.to_molecule();
    let coords = frame.coords();
    let text = write_cif(&mol, &coords, None);
    assert!(text.contains("_atom_site_Cartn_x"));
    assert!(!text.contains("_atom_site_fract_"));
    assert!(!text.contains("_cell_length_"));
    let parsed = parse_cif(&text).unwrap();
    assert!(parsed.cell.is_none());
    assert_eq!(parsed.mol.atom_count(), mol.atom_count());
    for (a, b) in coords.iter().zip(parsed.coords) {
        close(*a, b, 1e-10);
    }
}

#[test]
fn cell_writer_keeps_fractional_values_for_orthogonal_and_triclinic_cells() {
    let frame = parse_xyz(XYZ).unwrap();
    let mol = frame.to_molecule();
    let coords = frame.coords();
    for cell in [
        UnitCell {
            a: 4.,
            b: 6.,
            c: 8.,
            alpha: 90.,
            beta: 90.,
            gamma: 90.,
        },
        UnitCell {
            a: 8.,
            b: 9.,
            c: 10.,
            alpha: 75.,
            beta: 83.,
            gamma: 64.,
        },
    ] {
        let text = write_cif(&mol, &coords, Some(&cell));
        assert!(text.contains("_atom_site_fract_x"));
        assert!(!text.contains("_atom_site_Cartn_"));
        let parsed = parse_cif(&text).unwrap();
        assert!(parsed.cell.is_some());
        for (a, b) in coords.iter().zip(parsed.coords) {
            close(*a, b, 2e-4);
        }
    }
}

#[test]
fn mixed_incomplete_coordinate_families_are_rejected() {
    for columns in [
        "_atom_site_fract_x\n_atom_site_Cartn_y\n_atom_site_Cartn_z",
        "_atom_site_Cartn_x\n_atom_site_fract_y\n_atom_site_Cartn_z",
        "_atom_site_Cartn_x\n_atom_site_Cartn_y\n_atom_site_fract_z",
    ] {
        let text = format!(
            "data_test\n_cell_length_a 5\nloop_\n_atom_site_type_symbol\n{columns}\nC 1 2 3\n"
        );
        assert!(
            matches!(parse_cif(&text), Err(CifError::MissingCoordinateColumns)),
            "{columns}"
        );
    }
}

#[test]
fn complete_cartesian_family_wins_over_stray_fractional_column() {
    let text = "data_test\nloop_\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_Cartn_z\n_atom_site_Cartn_x\n_atom_site_Cartn_y\nC 0.1 4.75 2.5 -3.25\n";
    let parsed = parse_cif(text).unwrap();
    assert!(parsed.cell.is_none());
    assert_eq!(parsed.coords, [(2.5, -3.25, 4.75)]);
}

#[test]
fn nonfinite_coordinates_are_rejected_in_either_family() {
    for family in ["Cartn", "fract"] {
        for value in ["NaN", "inf", "-inf", "1e309"] {
            let text = format!(
                "data_test\n_cell_length_a 5\n_cell_length_b 6\n_cell_length_c 7\nloop_\n_atom_site_type_symbol\n_atom_site_{family}_x\n_atom_site_{family}_y\n_atom_site_{family}_z\nC {value} 0 0\n"
            );
            assert!(
                matches!(parse_cif(&text), Err(CifError::InvalidCoordinate(_))),
                "{family}: {value}"
            );
        }
    }
    // Individually finite inputs can still overflow during the cell transform.
    let overflow = "data_test\n_cell_length_a 5\n_cell_length_b 6\n_cell_length_c 7\nloop_\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\nC 1e308 0 0\n";
    assert!(matches!(
        parse_cif(overflow),
        Err(CifError::InvalidCoordinate(_))
    ));
}
