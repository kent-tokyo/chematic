use chematic_mol::{
    MmcifError, MmcifParseLimits, MmcifResidueError, parse_mmcif_residues,
    parse_mmcif_residues_with_limits,
};

const FIXTURE: &str =
    include_str!("../../../validation/openbabel_issues/mmcif/branched-residues.cif");
fn rows() -> Vec<Vec<String>> {
    FIXTURE
        .lines()
        .filter(|s| s.starts_with("HETATM"))
        .map(|s| s.split_whitespace().map(str::to_owned).collect())
        .collect()
}
fn render(rows: &[Vec<String>]) -> String {
    let header = FIXTURE.split("HETATM").next().unwrap();
    format!(
        "{}{}\n#\n",
        header,
        rows.iter()
            .map(|r| r.join(" "))
            .collect::<Vec<_>>()
            .join("\n")
    )
}

#[test]
fn nonpolymer_sequence_placeholders_do_not_merge_three_sugars() {
    let doc = diagnostic(parse_mmcif_residues(FIXTURE)).unwrap();
    assert_eq!(doc.residues.len(), 3);
    assert_eq!(doc.atom_to_residue, vec![0, 0, 1, 1, 2, 2]);
    let names: Vec<_> = doc
        .residues
        .iter()
        .map(|r| r.identity.label_comp_id.as_deref().unwrap())
        .collect();
    assert_eq!(names, vec!["NAG", "NAG", "BMA"]);
    for (i, r) in doc.residues.iter().enumerate() {
        assert_eq!(
            r.identity.auth_seq_id.as_deref(),
            Some((i + 1).to_string().as_str())
        );
        assert_eq!(r.identity.label_seq_id, None);
        assert_eq!(r.atom_indices, vec![i * 2, i * 2 + 1]);
    }
    let (mol, coords) = doc.parsed.to_molecule();
    assert_eq!(mol.atom_count(), 6);
    assert_eq!(mol.bond_count(), 0);
    assert_eq!(coords[4], (60.0, 0.0, 0.0));
}

fn permute(order: &mut [usize], start: usize, visit: &mut impl FnMut(&[usize])) {
    if start == order.len() {
        visit(order);
        return;
    }
    for i in start..order.len() {
        order.swap(start, i);
        permute(order, start + 1, visit);
        order.swap(start, i);
    }
}
#[test]
fn all_720_row_permutations_keep_partition_and_source_mapping() {
    let original = rows();
    let mut count = 0;
    permute(&mut [0, 1, 2, 3, 4, 5], 0, &mut |order| {
        let reordered: Vec<_> = order.iter().map(|&i| original[i].clone()).collect();
        let doc = diagnostic(parse_mmcif_residues(&render(&reordered))).unwrap();
        assert_eq!(doc.residues.len(), 3);
        for (i, &old) in order.iter().enumerate() {
            let r = &doc.residues[doc.atom_to_residue[i]];
            assert_eq!(
                r.identity.auth_seq_id.as_deref(),
                Some((old / 2 + 1).to_string().as_str())
            );
            assert_eq!(r.atom_indices.len(), 2);
            assert!(r.atom_indices.contains(&i));
            assert_eq!(doc.parsed.atoms[i].serial, (old + 1) as i64);
        }
        count += 1;
    });
    assert_eq!(count, 720);
}

#[test]
fn keeps_both_chain_namespaces_when_author_chain_collides() {
    let mut r = rows();
    r.truncate(2);
    r[1][3] = "C1".into();
    r[1][18] = "C1".into();
    r[1][6] = "F".into();
    let doc = diagnostic(parse_mmcif_residues(&render(&r))).unwrap();
    assert_eq!(doc.residues.len(), 2);
    assert_eq!(
        doc.residues[0].identity.auth_asym_id,
        doc.residues[1].identity.auth_asym_id
    );
    assert_ne!(
        doc.residues[0].identity.label_asym_id,
        doc.residues[1].identity.label_asym_id
    );
}

#[test]
fn model_entity_insertion_and_components_are_distinct_identity_fields() {
    for (column, value) in [
        (19, "2"),
        (7, "2"),
        (9, "AB"),
        (5, "MAN"),
        (16, "LIG"),
        (8, "99"),
    ] {
        let mut r = rows();
        r.truncate(2);
        r[1][3] = "C1".into();
        r[1][18] = "C1".into();
        r[1][column] = value.into();
        assert_eq!(
            diagnostic(parse_mmcif_residues(&render(&r)))
                .unwrap()
                .residues
                .len(),
            2,
            "column {column}"
        );
    }
}

#[test]
fn alternate_locations_stay_atoms_of_the_same_residue_and_keep_full_codes() {
    let mut r = rows();
    r.truncate(2);
    r[1][3] = "C1".into();
    r[1][18] = "C1".into();
    r[0][4] = "AA".into();
    r[1][4] = "BB".into();
    r[0][9] = "AB".into();
    r[1][9] = "AB".into();
    let d = diagnostic(parse_mmcif_residues(&render(&r))).unwrap();
    assert_eq!(d.residues.len(), 1);
    assert_eq!(d.atom_ids[0].alt_id.as_deref(), Some("AA"));
    assert_eq!(d.atom_ids[1].alt_id.as_deref(), Some("BB"));
    assert_eq!(d.residues[0].identity.insertion_code.as_deref(), Some("AB"));
}

#[test]
fn duplicate_atom_names_are_not_silently_merged() {
    let mut r = rows();
    r.truncate(2);
    r[1][3] = "C1".into();
    r[1][18] = "AUTH_DIFFERENT".into();
    assert!(matches!(
        diagnostic(parse_mmcif_residues(&render(&r))),
        Err(MmcifResidueError::DuplicateAtom {
            atom: 1,
            residue: 0
        })
    ));
}

#[test]
fn missing_identity_errors_and_explicit_author_zero_remains_valid() {
    for (columns, field) in [
        (&[6, 17][..], "chain identifier"),
        (&[8, 15][..], "residue sequence identifier"),
        (&[5, 16][..], "component identifier"),
        (&[3, 18][..], "atom identifier"),
    ] {
        let mut r = rows();
        r.truncate(1);
        for &i in columns {
            r[0][i] = "?".into();
        }
        assert_eq!(
            diagnostic(parse_mmcif_residues(&render(&r))).err(),
            Some(MmcifResidueError::MissingIdentity { atom: 0, field })
        );
    }
    let mut r = rows();
    r.truncate(1);
    r[0][15] = "0".into();
    assert_eq!(
        diagnostic(parse_mmcif_residues(&render(&r)))
            .unwrap()
            .residues[0]
            .identity
            .auth_seq_id
            .as_deref(),
        Some("0")
    );
}

#[test]
fn author_and_label_fallbacks_are_explicit_in_identity() {
    for cols in [&[3, 5, 6, 8][..], &[15, 16, 17, 18][..]] {
        let mut r = rows();
        r.truncate(1);
        r[0][8] = "17".into();
        for &i in cols {
            r[0][i] = ".".into();
        }
        let d = diagnostic(parse_mmcif_residues(&render(&r))).unwrap();
        assert_eq!(d.residues.len(), 1);
    }
}

#[test]
fn leading_zero_author_identifiers_remain_distinct() {
    let mut r = rows();
    r.truncate(2);
    r[1][3] = "C1".into();
    r[1][18] = "C1".into();
    r[1][15] = "01".into();
    assert_eq!(
        diagnostic(parse_mmcif_residues(&render(&r)))
            .unwrap()
            .residues
            .len(),
        2
    );
}

#[test]
fn long_component_names_and_atom_aliases_are_preserved() {
    let mut r = rows();
    r.truncate(1);
    r[0][5] = "A1AJ1".into();
    r[0][16] = "AUTHOR_LONG".into();
    r[0][18] = "NAME_LONG".into();
    let d = diagnostic(parse_mmcif_residues(&render(&r))).unwrap();
    assert_eq!(
        d.residues[0].identity.label_comp_id.as_deref(),
        Some("A1AJ1")
    );
    assert_eq!(
        d.residues[0].identity.auth_comp_id.as_deref(),
        Some("AUTHOR_LONG")
    );
    assert_eq!(d.atom_ids[0].label_atom_id.as_deref(), Some("C1"));
    assert_eq!(d.atom_ids[0].auth_atom_id.as_deref(), Some("NAME_LONG"));
}

#[test]
fn incomplete_rows_and_empty_loops_are_errors() {
    let mut r = rows();
    r.truncate(1);
    r[0].pop();
    assert!(matches!(
        diagnostic(parse_mmcif_residues(&render(&r))),
        Err(MmcifResidueError::IncompleteRow { atom: 0 })
    ));
    assert!(matches!(
        diagnostic(parse_mmcif_residues("data_empty\n")),
        Err(MmcifResidueError::Parse(MmcifError::NoAtomSiteLoop))
    ));
}

#[test]
fn limits_apply_before_atom_or_token_growth() {
    for l in [
        MmcifParseLimits {
            max_input_bytes: 8,
            ..Default::default()
        },
        MmcifParseLimits {
            max_atoms: 2,
            ..Default::default()
        },
        MmcifParseLimits {
            max_line_len: 8,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            diagnostic(parse_mmcif_residues_with_limits(FIXTURE, &l)),
            Err(MmcifResidueError::Parse(_))
        ));
    }
    let exact = MmcifParseLimits {
        max_input_bytes: FIXTURE.len(),
        max_atoms: 6,
        max_line_len: 8192,
    };
    assert!(diagnostic(parse_mmcif_residues_with_limits(FIXTURE, &exact)).is_ok());
}

#[test]
fn model_range_is_checked_instead_of_wrapping() {
    for v in ["0", "-1", "2147483648", "4294967297"] {
        let mut r = rows();
        r.truncate(1);
        r[0][19] = v.into();
        assert!(matches!(
            diagnostic(parse_mmcif_residues(&render(&r))),
            Err(MmcifResidueError::InvalidField {
                field: "model number",
                ..
            })
        ));
    }
}

#[test]
fn unsupported_author_strings_and_invalid_numbers_fail_explicitly() {
    for (column, v) in [(15, "12A"), (10, "NaN"), (2, "NotAnElement")] {
        let mut r = rows();
        r.truncate(1);
        r[0][column] = v.into();
        assert!(matches!(
            diagnostic(parse_mmcif_residues(&render(&r))),
            Err(MmcifResidueError::Parse(_))
        ));
    }
}

#[test]
fn preserves_cell_space_group_and_reports_unhandled_columns() {
    let input=FIXTURE.replacen("loop_\n","_cell.length_a 10\n_cell.length_b 11\n_cell.length_c 12\n_cell.angle_alpha 90\n_cell.angle_beta 90\n_cell.angle_gamma 90\n_symmetry.space_group_name_H-M 'P 1'\nloop_\n",1);
    let d = diagnostic(parse_mmcif_residues(&input)).unwrap();
    assert_eq!(d.parsed.cell.unwrap().a, 10.0);
    assert_eq!(d.parsed.space_group.as_deref(), Some("P 1"));
    let input = FIXTURE
        .replacen(
            "_atom_site.group_PDB",
            "_atom_site.custom_tag\n_atom_site.group_PDB",
            1,
        )
        .lines()
        .map(|s| {
            if s.starts_with("HETATM") {
                format!("opaque {s}")
            } else {
                s.into()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert_eq!(
        diagnostic(parse_mmcif_residues(&input))
            .unwrap()
            .parsed
            .unhandled_columns,
        vec!["_atom_site.custom_tag"]
    );
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
