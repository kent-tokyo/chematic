//! Reused token storage versus the frozen per-line allocation path.
use super::*;

fn frozen_per_line_tokens(input: &str) -> Result<MolReadReport, MolParseError> {
    let lines = validate_v3000_physical_lines(input)?;
    let (metadata, line2) = v3000_metadata(&lines)?;
    let logical_lines = collect_v30_lines(&lines);
    let mut parser = V3000Parser::new();
    for line in &logical_lines {
        let tokens = line.payload.split_whitespace().collect::<Vec<_>>();
        parser.consume(line, &tokens)?;
    }
    finish_v3000_report(parser, metadata, line2)
}

fn compare(input: &str) -> bool {
    match (
        frozen_per_line_tokens(input),
        read_mol_v3000_with_diagnostics(input),
    ) {
        (Err(old), Err(new)) => {
            assert_eq!(old.to_string(), new.to_string());
            false
        }
        (Ok(old), Ok(new)) => {
            assert_eq!(
                old.mol.atoms().collect::<Vec<_>>(),
                new.mol.atoms().collect::<Vec<_>>()
            );
            assert_eq!(
                old.mol.bonds().collect::<Vec<_>>(),
                new.mol.bonds().collect::<Vec<_>>()
            );
            assert_eq!(old.mol.stereo_groups(), new.mol.stereo_groups());
            assert_eq!(
                old.coords
                    .iter()
                    .map(|(x, y)| (x.to_bits(), y.to_bits()))
                    .collect::<Vec<_>>(),
                new.coords
                    .iter()
                    .map(|(x, y)| (x.to_bits(), y.to_bits()))
                    .collect::<Vec<_>>()
            );
            for (i, _) in old.mol.atoms() {
                assert_eq!(old.mol.neighbor_slice(i), new.mol.neighbor_slice(i));
                assert_eq!(
                    old.mol.stereo_neighbor_order(i),
                    new.mol.stereo_neighbor_order(i)
                );
                assert_eq!(old.mol.r_group_label(i), new.mol.r_group_label(i));
                assert_eq!(old.mol.atom_tag(i), new.mol.atom_tag(i));
            }
            for i in 0..old.mol.stereo_groups().len() {
                assert_eq!(
                    old.mol.stereo_group_bond_indices(i),
                    new.mol.stereo_group_bond_indices(i)
                );
            }
            assert_eq!(
                format!(
                    "{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}",
                    old.metadata,
                    old.stereo_diagnostics,
                    old.ez_diagnostics,
                    old.coordinate_dimension,
                    old.geometry_rank,
                    old.stereo3d_diagnostics,
                    old.square_planar_diagnostics,
                    old.conformer
                ),
                format!(
                    "{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}",
                    new.metadata,
                    new.stereo_diagnostics,
                    new.ez_diagnostics,
                    new.coordinate_dimension,
                    new.geometry_rank,
                    new.stereo3d_diagnostics,
                    new.square_planar_diagnostics,
                    new.conformer
                )
            );
            assert_eq!(
                write_mol_v3000(&old.mol, &old.metadata, &old.coords),
                write_mol_v3000(&new.mol, &new.metadata, &new.coords)
            );
            true
        }
        _ => panic!("token storage changed parse success/error boundary"),
    }
}

fn chain(n: usize, space: &str, continuation: bool) -> String {
    let mut text = format!(
        "token reuse\n  chematic          2D\nhand-generated regression\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS {n} {} 0 0 0\nM  V30 BEGIN ATOM\n",
        n - 1
    );
    for i in 0..n {
        let mut row = format!("{} C {} -0.0 0 0", i + 1, i).replace(' ', space);
        if i == 0 {
            let props = (0..80)
                .map(|p| format!("VENDOR{p}=value{p}"))
                .collect::<Vec<_>>()
                .join(space);
            if continuation {
                row.push_str(" -\nM  V30 ");
            } else {
                row.push_str(space);
            }
            row.push_str(&props);
        }
        text.push_str(&format!("M  V30 {row}\nM  V30 \t\n"));
    }
    text.push_str("M  V30 END ATOM\nM  V30 BEGIN BOND\n");
    for i in 1..n {
        text.push_str(&format!("M  V30 {} 1 {} {}\n", i, i, i + 1));
    }
    text.push_str("M  V30 END BOND\nM  V30 END CTAB\nM  END\n");
    text
}

#[test]
fn long_short_empty_and_continued_rows_keep_full_reports() {
    for n in [1, 2, 8, 64] {
        for space in [" ", "\t", "\u{2003}"] {
            for continuation in [false, true] {
                assert!(compare(&chain(n, space, continuation)));
            }
        }
    }
}

#[test]
fn malformed_physical_and_logical_rows_keep_error_details() {
    let input = chain(3, " ", false);
    let lines = input.lines().collect::<Vec<_>>();
    let mut errors = 0;
    for i in 0..lines.len() {
        let mut changed = lines.clone();
        changed.remove(i);
        errors += usize::from(!compare(&(changed.join("\n") + "\n")));
        for replacement in [
            "",
            "M V30 1 C 0 0 0 0",
            "M  V30 ",
            "M  V30 \u{2003}",
            "M  V30 1",
            "M  V30 END CTAB",
            "M  END",
        ] {
            let mut changed = lines.clone();
            changed[i] = replacement;
            errors += usize::from(!compare(&(changed.join("\n") + "\n")));
        }
    }
    assert!(errors > 100);
    for suffix in [
        "M  END ",
        "M  V30 trailing -",
        "M  V30 -\nnot V30",
        "\n\n",
        "M  V30 COUNTS 99999999999999999999999999999 0 0",
    ] {
        compare(&(input.clone() + suffix));
    }
}

#[test]
fn attributes_and_bad_numbers_keep_first_match_and_line_errors() {
    let input = chain(3, " ", false);
    for tail in [
        "CHG=1 CHG=-1 MASS=13 HCOUNT=2",
        "RGROUPS=(1 7) SUBST=2 RBCNT=3 UNSAT=1",
        "RGROUPS=(2 1 2)",
        "CFG=2 ENDPTS=(2 1 2) ATTACH=ANY",
        "CHG=invalid MASS=99999 HCOUNT=-1",
        "",
    ] {
        let changed = input.replace("1 C 0 -0.0 0 0", &format!("1 C 0 -0.0 0 0 {tail}"));
        compare(&changed);
    }
    for value in ["nan", "inf", "-inf", "bad", "1e999", "-0", "1e-300"] {
        compare(&input.replace("1 C 0 -0.0 0 0", &format!("1 C {value} -0.0 0 0")));
    }
}
