//! Original synthetic MIT OR Apache-2.0 fixtures; no upstream structure copied.
use chematic_mol::{PdbSelectionLimits, select_pdb_atoms, select_pdb_atoms_with_limits};
const TEXT: &str = include_str!("../../../validation/openbabel_issues/pdb-selection/rich.pdb");
fn atoms(text: &str) -> Vec<&str> {
    text.lines()
        .filter(|s| s.starts_with("ATOM  ") || s.starts_with("HETATM"))
        .collect()
}

#[test]
fn removing_ligand_keeps_protein_identity_and_exact_records() {
    let r = diagnostic(select_pdb_atoms(TEXT, Some(&[0, 1, 2, 3]), false)).unwrap();
    assert_eq!(atoms(&r.pdb), atoms(TEXT)[..4]);
    assert_eq!(r.source_text, TEXT);
    assert_eq!(r.source_atom_count, 6);
    assert_eq!(
        r.atoms
            .iter()
            .map(|a| a.residue_name.as_str())
            .collect::<Vec<_>>(),
        ["VAL", "VAL", "LEU", "LEU"]
    );
    assert_eq!(r.atoms[0].atom_name, " CA ");
    assert_eq!(r.atoms[0].alt_loc, 'A');
    assert_eq!(r.atoms[2].residue_number, -2);
    assert_eq!(r.atoms[2].insertion_code, 'I');
    assert_eq!(
        r.old_to_new,
        [Some(0), Some(1), Some(2), Some(3), None, None]
    );
    assert!(r.retargeted_ter_lines.is_empty());
    assert_eq!(r.pdb.lines().filter(|s| s.starts_with("ANISOU")).count(), 4);
    assert_eq!(r.pdb.lines().filter(|s| s.starts_with("TER")).count(), 2);
    assert!(
        !r.pdb
            .lines()
            .filter(|s| s.starts_with("CONECT"))
            .any(|s| s.contains("    9") || s.contains("   10"))
    );
}

#[test]
fn unordered_selection_retains_source_order_and_retargets_ter() {
    let r = diagnostic(select_pdb_atoms(TEXT, Some(&[4, 2, 0]), false)).unwrap();
    assert_eq!(r.selected_source_indices, [0, 2, 4]);
    assert_eq!(r.old_to_new, [Some(0), None, Some(1), None, Some(2), None]);
    assert_eq!(r.retargeted_ter_lines, [5, 10, 15]);
    assert_eq!(
        r.pdb
            .lines()
            .filter(|s| s.starts_with("TER"))
            .map(|s| s[6..11].trim())
            .collect::<Vec<_>>(),
        ["2", "6", "10"]
    );
    assert_eq!(r.pdb.lines().filter(|s| s.starts_with("CONECT")).count(), 0);
    assert_eq!(r.pdb.lines().filter(|s| s.starts_with("ANISOU")).count(), 3);
    let again = diagnostic(select_pdb_atoms(&r.pdb, None, false)).unwrap();
    assert_eq!(again.pdb, r.pdb);
}

#[test]
fn all_records_and_duplicate_connectivity_multiplicity_are_preserved() {
    let r = diagnostic(select_pdb_atoms(TEXT, None, false)).unwrap();
    assert_eq!(r.pdb, TEXT);
    assert!(r.filtered_conect_lines.is_empty());
    assert!(r.pdb.contains("CONECT    1    2    2"));
    assert_eq!(r.atoms[5].record_type, "HETATM");
    assert_eq!(atoms(&r.pdb)[5][78..80], *"1-");
}

#[test]
fn partial_conect_target_filter_retains_order_and_multiplicity() {
    let source = TEXT.replace("END\n", "CONECT    1    2    5    5    9\nEND\n");
    let r = diagnostic(select_pdb_atoms(&source, Some(&[0, 1, 2]), false)).unwrap();
    assert!(r.pdb.contains("CONECT    1    2    5    5\n"));
    assert!(!r.pdb.contains("    5    5    9"));
}

#[test]
fn explicit_single_model_crlf_and_empty_selection() {
    let source = format!(
        "MODEL        7\n{}ENDMDL\nEND\n",
        TEXT.lines()
            .filter(|s| !s.starts_with("CONECT") && *s != "END")
            .collect::<Vec<_>>()
            .join("\n")
            + "\n"
    )
    .replace('\n', "\r\n");
    let r = diagnostic(select_pdb_atoms(&source, Some(&[0]), false)).unwrap();
    assert_eq!(r.model_num, Some(7));
    assert!(r.pdb.starts_with("MODEL        7\n"));
    assert!(r.pdb.ends_with("ENDMDL\nEND\n"));
    assert_eq!(r.source_text, source);
    let r = diagnostic(select_pdb_atoms(&source, Some(&[]), false)).unwrap();
    assert_eq!(r.pdb, "END\n");
    assert!(r.old_to_new.iter().all(Option::is_none));
    assert_eq!(
        diagnostic(select_pdb_atoms("", None, false)).unwrap().pdb,
        "END\n"
    );
}

#[test]
fn metadata_loss_requires_consent_and_keeps_source_sidecar() {
    let source = format!("HEADER    original synthetic\n{TEXT}");
    assert!(diagnostic(select_pdb_atoms(&source, None, false)).is_err());
    let r = diagnostic(select_pdb_atoms(&source, None, true)).unwrap();
    assert_eq!(r.pdb, TEXT);
    assert_eq!(r.source_text, source);
    assert_eq!(r.losses.len(), 1);
    assert_eq!(r.losses[0].line, 1);
    assert_eq!(r.losses[0].record_name, "HEADER");
}

#[test]
fn invalid_indices_and_reference_fields_fail_before_output() {
    for selected in [&[0, 0][..], &[6][..]] {
        assert!(diagnostic(select_pdb_atoms(TEXT, Some(selected), false)).is_err());
    }
    for source in [
        TEXT.replace("CONECT    1    2    2", "CONECT    1   99"),
        TEXT.replace("CONECT    1    2    2", "CONECT    1    1"),
        TEXT.replace("ANISOU    1  CA AVAL", "ANISOU   99  CA AVAL"),
        format!("{TEXT}ATOM"),
        format!("ENDMDL\n{TEXT}"),
        format!("MODEL        1\n{TEXT}"),
        TEXT.replace("TER       3", "TER       4"),
    ] {
        assert!(
            diagnostic(select_pdb_atoms(&source, None, false)).is_err(),
            "{source}"
        );
    }
    let row = atoms(TEXT)[0];
    assert!(
        diagnostic(select_pdb_atoms(
            &format!("{row}\n{row}\nEND\n"),
            None,
            false
        ))
        .is_err()
    );
}

#[test]
fn malformed_atom_and_anisou_fields_fail_with_physical_line_numbers() {
    let mut row = atoms(TEXT)[0].to_string();
    row.replace_range(30..38, "     NaN");
    assert_eq!(
        diagnostic(select_pdb_atoms(&format!("\n{row}\nEND\n"), None, false))
            .unwrap_err()
            .line,
        2
    );
    let source = TEXT.replace("    100    200", "    bad    200");
    assert!(diagnostic(select_pdb_atoms(&source, None, false)).is_err());
    assert!(
        diagnostic(select_pdb_atoms(
            &TEXT.replace(" 23.45", "   inf"),
            None,
            false
        ))
        .is_err()
    );
    assert!(diagnostic(select_pdb_atoms(&format!("REMARK λ\n{TEXT}"), None, true)).is_err());
}

#[test]
fn resource_limits_are_checked_and_cannot_raise_hard_caps() {
    for limits in [
        PdbSelectionLimits {
            max_atoms: 5,
            ..Default::default()
        },
        PdbSelectionLimits {
            max_input_bytes: 2,
            ..Default::default()
        },
        PdbSelectionLimits {
            max_line_len: 20,
            ..Default::default()
        },
        PdbSelectionLimits {
            max_records: 2,
            ..Default::default()
        },
    ] {
        assert!(diagnostic(select_pdb_atoms_with_limits(TEXT, None, false, &limits)).is_err());
    }
    assert!(
        diagnostic(select_pdb_atoms_with_limits(
            &format!("REMARK{}\nEND\n", "a".repeat(1025)),
            None,
            true,
            &PdbSelectionLimits {
                max_line_len: usize::MAX,
                ..Default::default()
            }
        ))
        .is_err()
    );
}

#[test]
fn ter_collision_is_an_error_instead_of_renumbering_atoms() {
    let source = atoms(TEXT);
    let mut a = source[0].to_string();
    a.replace_range(6..11, "   10");
    let mut b = source[1].to_string();
    b.replace_range(6..11, "   20");
    let mut c = source[2].to_string();
    c.replace_range(6..11, "   11");
    let text = format!(
        "{a}\n{b}\nTER      21      {}\n{c}\nTER      12      {}\nEND\n",
        &b[17..27],
        &c[17..27]
    );
    assert!(diagnostic(select_pdb_atoms(&text, None, false)).is_ok());
    assert!(
        diagnostic(select_pdb_atoms(&text, Some(&[0, 2]), false))
            .unwrap_err()
            .detail
            .contains("collides")
    );
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
