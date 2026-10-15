//! Original minimal MIT OR Apache-2.0 EXYZ fixtures. Independent oracle recipe
//! and normalized payload hashes: scripts/exyz_openbabel_gate.py.
use chematic_core::Element;
use chematic_mol::{
    ExyzAtom, ExyzGeometry, ExyzParseLimits, parse_exyz, parse_exyz_with_limits, write_exyz,
};

const TEXT: &str = "2\n title  %VIRTUAL %PBC \nSi 0.12345678901234567 -2 3\nCl -1 0 2 VIRTUAL\n\nVector1 4 1 0\nVector2 -1 5 2\nVector3 1 0 -6\nOffset 0.5 -1 2\n";

#[test]
fn complete_geometry_preserves_basis_offset_flags_and_text() {
    let g = diagnostic(parse_exyz(TEXT)).unwrap();
    assert_eq!(g.header, " title  %VIRTUAL %PBC ");
    assert_eq!(g.source_text, TEXT);
    assert_eq!(g.atoms[0].element, Element::SI);
    assert!(!g.atoms[0].virtual_atom);
    assert!(g.atoms[1].virtual_atom);
    assert_eq!(g.cell.as_ref().unwrap().offset_angstrom, [0.5, -1.0, 2.0]);
    let out = diagnostic(write_exyz(&g)).unwrap();
    let back = diagnostic(parse_exyz(&out)).unwrap();
    assert_eq!(back.header, g.header);
    assert_eq!(back.atoms, g.atoms);
    assert_eq!(back.cell, g.cell);
}

#[test]
fn writer_uses_structured_data_after_edits() {
    let mut g = diagnostic(parse_exyz(TEXT)).unwrap();
    g.atoms[0].position_angstrom[0] = -0.125;
    g.atoms[1].virtual_atom = false;
    g.cell.as_mut().unwrap().offset_angstrom = [1.0, 2.0, 3.0];
    let back = diagnostic(parse_exyz(&diagnostic(write_exyz(&g)).unwrap())).unwrap();
    assert_eq!(g.atoms, back.atoms);
    assert_eq!(g.cell, back.cell);
}

#[test]
fn all_elements_and_finite_extremes_roundtrip_exactly() {
    let mut atoms = Vec::new();
    for number in 1..=118 {
        atoms.push(ExyzAtom {
            element: Element::from_atomic_number(number).unwrap(),
            position_angstrom: [f64::MIN_POSITIVE, f64::MAX, -f64::from(number)],
            virtual_atom: number % 2 == 0,
        });
    }
    atoms[0].position_angstrom = [-0.0, f64::from_bits(1), -f64::from_bits(1)];
    let g = ExyzGeometry {
        header: "%VIRTUAL".into(),
        atoms,
        cell: None,
        source_text: String::new(),
    };
    let back = diagnostic(parse_exyz(&diagnostic(write_exyz(&g)).unwrap())).unwrap();
    assert_eq!(g.atoms, back.atoms);
    assert_eq!(
        g.atoms[0].position_angstrom[0].to_bits(),
        back.atoms[0].position_angstrom[0].to_bits()
    );
}

#[test]
fn plain_xyz_subset_crlf_and_final_newline_are_supported() {
    for text in ["1\n\nC 0 0 0", "1\r\nplain title\r\nC 0 0 0\r\n\r\n"] {
        let g = diagnostic(parse_exyz(text)).unwrap();
        assert!(g.cell.is_none());
        assert!(!g.atoms[0].virtual_atom);
    }
}

#[test]
fn count_preflight_and_caller_limits_fail_before_allocation() {
    for text in [
        "18446744073709551615\n\n",
        "1000000\n\nC 0 0 0\n",
        "-1\n\n",
        "0\n\n",
        "2\n\nC 0 0 0\n",
    ] {
        assert!(diagnostic(parse_exyz(text)).is_err(), "{text:?}");
    }
    for limits in [
        ExyzParseLimits {
            max_input_bytes: TEXT.len() - 1,
            ..Default::default()
        },
        ExyzParseLimits {
            max_atoms: 1,
            ..Default::default()
        },
        ExyzParseLimits {
            max_line_len: 5,
            ..Default::default()
        },
    ] {
        assert!(diagnostic(parse_exyz_with_limits(TEXT, &limits)).is_err());
    }
}

#[test]
fn unsupported_headers_and_atom_fields_fail() {
    for text in [
        "1\n%UNKNOWN\nC 0 0 0\n",
        "1\n%VIRTUAL %VIRTUAL\nC 0 0 0\n",
        "1\nProperties=species:S:1:pos:R:3\nC 0 0 0\n",
        "1\nLattice=\"1 0 0 0 1 0 0 0 1\"\nC 0 0 0\n",
        "1\n\nC 0 0 0 VIRTUAL\n",
        "1\n%VIRTUAL\nC 0 0 0 EXTRA\n",
        "1\n\nc 0 0 0\n",
        "1\n\n13C 0 0 0\n",
        "1\n\nX 0 0 0\n",
        "1\n\nC NaN 0 0\n",
        "1\n\nC 0 inf 0\n",
        "1\n\nC 0 0 1e999\n",
        "1\n\nC 0 0\n",
        "1\n\nC 0 0 0 1\n",
        "1\n\nC 0 0 0\nextra\n",
    ] {
        assert!(diagnostic(parse_exyz(text)).is_err(), "{text:?}");
    }
}

#[test]
fn strict_cell_block_rejects_incomplete_reordered_or_unflagged_data() {
    for text in [
        TEXT.replace("%PBC", ""),
        TEXT.replace("%PBC", "%PBC %PBC"),
        TEXT.replace("\n\nVector1", "\nVector1"),
        TEXT.replace("Vector2", "Vector1"),
        TEXT.replace("Offset 0.5 -1 2\n", ""),
        TEXT.replace("Offset 0.5 -1 2", "Offset 0.5 -1"),
        TEXT.replace("Vector1 4 1 0", "Vector1 0 0 0"),
        TEXT.replace("Vector3 1 0 -6", "Vector3 4 1 0"),
        TEXT.replace("Vector3 1 0 -6", "Vector3 4 1 1e-12"),
        TEXT.replace("Offset 0.5 -1 2", "Offset NaN -1 2"),
    ] {
        assert!(diagnostic(parse_exyz(&text)).is_err(), "{text:?}");
    }
}

#[test]
fn errors_identify_rejected_source_line() {
    assert_eq!(
        diagnostic(parse_exyz("1\n%NO\nC 0 0 0")).unwrap_err().line,
        2
    );
    assert_eq!(
        diagnostic(parse_exyz("1\n\nC 0 NaN 0")).unwrap_err().line,
        3
    );
    assert_eq!(
        diagnostic(parse_exyz(&TEXT.replace("Vector2", "bad")))
            .unwrap_err()
            .line,
        7
    );
}

#[test]
fn mutable_invalid_geometry_cannot_be_written() {
    for change in 0..7 {
        let mut g = diagnostic(parse_exyz(TEXT)).unwrap();
        match change {
            0 => g.atoms.clear(),
            1 => g.atoms[0].position_angstrom[0] = f64::INFINITY,
            2 => g.header = "%PBC".into(),
            3 => g.cell = None,
            4 => g.header.push('\n'),
            5 => g.cell.as_mut().unwrap().offset_angstrom[0] = f64::NAN,
            _ => g.cell.as_mut().unwrap().lattice_angstrom[0] = [0.0; 3],
        }
        assert!(diagnostic(write_exyz(&g)).is_err());
    }
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
