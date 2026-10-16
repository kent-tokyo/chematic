//! Frozen main reader: full reports and exact errors must be preserved.
use super::*;

fn frozen_collected_rows(input: &str) -> Result<MolReadReport, MolParseError> {
    let lines = validate_v3000_physical_lines(input)?;
    let (metadata, line2) = v3000_metadata(&lines)?;
    let logical_lines = collect_v30_lines(&lines);
    let mut parser = V3000Parser::new();
    parser.allocation_limit = logical_lines.len();
    // Logical payloads live until the entire parse ends, including owned
    // continuation text. Reuse only the token vector's storage; clear every
    // row before borrowing that row's exact Unicode-whitespace tokens.
    let mut tokens = Vec::new();
    for line in &logical_lines {
        tokens.clear();
        tokens.extend(line.payload.split_whitespace());
        parser.consume(line, &tokens)?;
    }
    finish_v3000_report(parser, metadata, line2)
}

fn compare(input: &str) -> bool {
    match (
        frozen_collected_rows(input),
        read_mol_v3000_with_diagnostics(input),
    ) {
        (Err(old), Err(new)) => {
            assert_eq!(old, new);
            false
        }
        (Ok(old), Ok(new)) => {
            let points = |report: &MolReadReport| {
                report.conformer.as_ref().map(|coords| {
                    coords
                        .points
                        .iter()
                        .map(|p| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()])
                        .collect::<Vec<_>>()
                })
            };
            assert_eq!(points(&old), points(&new));
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
        _ => panic!("borrowed rows changed parse success/error boundary"),
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
fn plain_and_continued_rows_keep_full_reports() {
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

#[test]
fn all_pinned_v3000_rows_preserve_full_reports() {
    use flate2::read::GzDecoder;
    use std::io::Read;
    let corpora: [(&[u8], usize); 3] = [
        (
            include_bytes!(
                "../../../../validation/benchmark_corpora/openbabel-v3000-speed-slice0-1000-v1.molblocks.gz"
            ),
            1000,
        ),
        (
            include_bytes!(
                "../../../../validation/benchmark_corpora/openbabel-v3000-speed-slice5000-1000-v1.molblocks.gz"
            ),
            1000,
        ),
        (
            include_bytes!(
                "../../../../validation/benchmark_corpora/openbabel-v3000-speed-chains64-v1.molblocks.gz"
            ),
            64,
        ),
    ];
    let mut rows = 0;
    let mut plain = 0;
    for (corpus, expected) in corpora {
        let mut input = String::new();
        GzDecoder::new(corpus).read_to_string(&mut input).unwrap();
        let mut count = 0;
        for block in input
            .split("$$$$\n")
            .filter(|block| !block.trim().is_empty())
        {
            let lines = validate_v3000_physical_lines(block).unwrap();
            if let Some(n) = plain_v3000_line_count(&lines) {
                assert_eq!(n, collect_v30_lines(&lines).len());
                plain += 1;
            }
            for _ in 0..2 {
                assert!(compare(block));
            }
            count += 1;
            rows += 1;
        }
        assert_eq!(count, expected);
    }
    assert_eq!(rows, 2064);
    assert!(plain > 2000);
}

#[test]
fn allocation_hints_and_continuation_selection_match_original_rows() {
    for continuation in [false, true] {
        let input = chain(3, " ", continuation);
        let lines = validate_v3000_physical_lines(&input).unwrap();
        assert_eq!(plain_v3000_line_count(&lines).is_some(), !continuation);
        let mut old = V3000Parser::new();
        let logical = collect_v30_lines(&lines);
        old.allocation_limit = logical.len();
        for line in &logical {
            let tokens = line.payload.split_whitespace().collect::<Vec<_>>();
            old.consume(line, &tokens).unwrap();
        }
        let mut new = V3000Parser::new();
        consume_v3000_rows(&mut new, &lines).unwrap();
        assert_eq!(new.allocation_limit, old.allocation_limit);
        assert_eq!(new.coords.capacity(), old.coords.capacity());
        assert_eq!(new.raw_z.capacity(), old.raw_z.capacity());
        assert_eq!(
            new.atom_idx_map.entries.capacity(),
            old.atom_idx_map.entries.capacity()
        );
        assert_eq!(
            new.bond_idx_map.entries.capacity(),
            old.bond_idx_map.entries.capacity()
        );
    }
    let input = chain(3, " ", false);
    for header in ["M  V30 title", "M  V30 title-", "M  V30 ", "plain title-"] {
        let (_, tail) = input.split_once('\n').unwrap();
        compare(&format!("{header}\n{tail}"));
    }
    for counts in ["0 0 0 0 0", "999999999 999999999 0 0 0", "3 2 0 0 0"] {
        compare(&input.replace("COUNTS 3 2 0 0 0", &format!("COUNTS {counts}")));
    }
}

#[test]
fn three_dimensional_stereo_and_collections_keep_full_reports() {
    for smiles in ["CCO", "C[C@H](F)Cl", "F/C=C/F", "N->[Pt@SP1](<-N)(Cl)Cl"] {
        let mut mol = chematic_smiles::parse(smiles).unwrap();
        mol.set_stereo_groups(vec![
            StereoGroup::new(StereoGroupKind::Absolute, vec![AtomIdx(0)]),
            StereoGroup::new(StereoGroupKind::And(7), vec![AtomIdx(1)]),
        ]);
        for three_d in [false, true] {
            let coords = Coords3D {
                points: (0..mol.atom_count())
                    .map(|i| {
                        Point3::new(
                            i as f64,
                            (i % 2) as f64,
                            if three_d { (i % 3) as f64 } else { -0.0 },
                        )
                    })
                    .collect(),
            };
            let text = write_mol_v3000_with_conformer(&mol, &MolMetadata::default(), &coords);
            assert!(compare(&text));
            assert!(compare(&text.replace("\n", "\r\n")));
        }
    }
}
