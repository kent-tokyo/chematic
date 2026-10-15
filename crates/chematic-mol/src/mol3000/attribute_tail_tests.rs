//! Complete old bond reader and endpoint normalizer frozen before editing.
use super::*;

fn reference_parse_v3000_bond_line(
    tokens: &[&str],
    line: usize,
    builder: &mut MoleculeBuilder,
    atom_idx_map: &V3000IndexMap<AtomIdx>,
    bond_idx_map: &mut V3000IndexMap<BondIdx>,
    explicitly_unspecified_ez: &mut std::collections::HashSet<BondIdx>,
    bond_properties: &mut Vec<(u32, String)>,
) -> Result<(), MolParseError> {
    // Bond line: <idx> <type> <atom1> <atom2> [KEY=VAL ...]
    if tokens.len() < 4 {
        return Err(MolParseError::InvalidBondLine {
            line,
            detail: format!(
                "V3000 bond line needs at least 4 fields, got {}",
                tokens.len()
            ),
        });
    }

    let v3k_idx = tokens[0]
        .parse::<u32>()
        .map_err(|_| MolParseError::InvalidBondLine {
            line,
            detail: format!("cannot parse bond index from '{}'", tokens[0]),
        })?;
    if bond_idx_map.resolve(v3k_idx).is_some() {
        return Err(MolParseError::InvalidBondLine {
            line,
            detail: format!("duplicate V3000 bond index {v3k_idx}"),
        });
    }
    let btype_raw = tokens[1]
        .parse::<u8>()
        .map_err(|_| MolParseError::InvalidBondLine {
            line,
            detail: format!("cannot parse bond type from '{}'", tokens[1]),
        })?;
    let a1_v3k = tokens[2]
        .parse::<u32>()
        .map_err(|_| MolParseError::InvalidBondLine {
            line,
            detail: format!("cannot parse atom1 index from '{}'", tokens[2]),
        })?;
    let a2_v3k = tokens[3]
        .parse::<u32>()
        .map_err(|_| MolParseError::InvalidBondLine {
            line,
            detail: format!("cannot parse atom2 index from '{}'", tokens[3]),
        })?;
    let a1 = atom_idx_map
        .resolve(a1_v3k)
        .ok_or_else(|| MolParseError::InvalidBondLine {
            line,
            detail: format!("atom index {a1_v3k} not found in atom block"),
        })?;
    let a2 = atom_idx_map
        .resolve(a2_v3k)
        .ok_or_else(|| MolParseError::InvalidBondLine {
            line,
            detail: format!("atom index {a2_v3k} not found in atom block"),
        })?;

    let kv_tokens = tokens.get(4..).unwrap_or(&[]);
    let order = match btype_raw {
        0 => BondOrder::Zero,
        1 => match parse_kv(kv_tokens, "CFG") {
            Some("1") => BondOrder::Up,
            Some("3") => BondOrder::Down,
            _ => BondOrder::Single,
        },
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        5 => BondOrder::QuerySingleOrDouble,
        6 => BondOrder::QuerySingleOrAromatic,
        7 => BondOrder::QueryDoubleOrAromatic,
        8 => BondOrder::QueryAny,
        9 => BondOrder::Dative,
        _ => BondOrder::Single,
    };

    let bidx = builder
        .add_bond(a1, a2, order)
        .map_err(|e| MolParseError::InvalidBondLine {
            line,
            detail: e.to_string(),
        })?;
    bond_idx_map.push((v3k_idx, bidx));
    if btype_raw == 2 && parse_kv(kv_tokens, "CFG") == Some("2") {
        explicitly_unspecified_ez.insert(bidx);
    }
    // Preserve richer V3000 bond attributes without pretending to interpret
    // them. Wedge CFG is regenerated from BondOrder; double-bond CFG=2
    // has no core tag, so retain it explicitly for lossless rewriting. ENDPTS is the one supported reference-bearing attribute: its
    // atom ids must follow the normalized output order when source ids are
    // non-contiguous. ATTACH and future keys remain opaque.
    let opaque_properties = kv_tokens
        .iter()
        .filter(|token| !token.starts_with("CFG=") || (btype_raw == 2 && token.trim() == "CFG=2"))
        .copied()
        .collect::<Vec<_>>()
        .join(" ");
    let opaque_properties = normalize_v3000_endpoints(&opaque_properties, atom_idx_map);
    if !opaque_properties.is_empty() {
        // As with atom metadata, normalize a source bond id to the writer's
        // sequential id. V3000 permits non-contiguous ids, whereas the
        // writer always emits `BondIdx + 1`.
        bond_properties.push((bidx.0 + 1, opaque_properties));
    }
    Ok(())
}

/// Convert an `ENDPTS=(n atom-id...)` property from source V3000 ids to the
/// sequential ids emitted by this writer. Leave malformed or unknown-shaped
/// values unchanged: they remain opaque metadata rather than becoming a
/// guessed semantic edit.
fn normalize_v3000_endpoints(properties: &str, atom_idx_map: &[(u32, AtomIdx)]) -> String {
    let Some(start) = properties.find("ENDPTS=(") else {
        return properties.to_owned();
    };
    let value_start = start + "ENDPTS=(".len();
    let Some(close_offset) = properties[value_start..].find(')') else {
        return properties.to_owned();
    };
    let value_end = value_start + close_offset;
    let values: Vec<&str> = properties[value_start..value_end]
        .split_whitespace()
        .collect();
    let Some((declared, atom_ids)) = values.split_first() else {
        return properties.to_owned();
    };
    let Ok(declared) = declared.parse::<usize>() else {
        return properties.to_owned();
    };
    if atom_ids.len() != declared {
        return properties.to_owned();
    }
    let Some(writer_ids) = atom_ids
        .iter()
        .map(|id| {
            id.parse::<u32>()
                .ok()
                .and_then(|source_id| resolve_atom_idx(source_id, atom_idx_map))
                .map(|idx| (idx.0 + 1).to_string())
        })
        .collect::<Option<Vec<_>>>()
    else {
        return properties.to_owned();
    };
    let normalized = format!("ENDPTS=({declared} {})", writer_ids.join(" "));
    format!(
        "{}{}{}",
        &properties[..start],
        normalized,
        &properties[value_end + 1..]
    )
}

fn reference_read(input: &str) -> Result<MolReadReport, MolParseError> {
    let lines = validate_v3000_physical_lines(input)?;
    let (metadata, line2) = v3000_metadata(&lines)?;
    let logical = collect_v30_lines(&lines);
    let mut parser = V3000Parser::new();
    for line in &logical {
        let tokens: Vec<_> = line.payload.split_whitespace().collect();
        if !tokens.is_empty()
            && matches!(parser.state, V3000State::InBondBlock)
            && !is_marker(&tokens, "END", "BOND")
        {
            reference_parse_v3000_bond_line(
                &tokens,
                line.line_num,
                &mut parser.builder,
                &parser.atom_idx_map,
                &mut parser.bond_idx_map,
                &mut parser.unspecified_ez,
                &mut parser.bond_properties,
            )?;
            parser.parsed_bonds += 1;
        } else {
            parser.consume(line, &tokens)?;
        }
    }
    finish_v3000_report(parser, metadata, line2)
}

fn compare_read(input: &str) -> bool {
    match (
        reference_read(input),
        read_mol_v3000_with_diagnostics(input),
    ) {
        (Ok(a), Ok(b)) => {
            assert_eq!(
                a.mol.atoms().collect::<Vec<_>>(),
                b.mol.atoms().collect::<Vec<_>>()
            );
            assert_eq!(
                a.mol.bonds().collect::<Vec<_>>(),
                b.mol.bonds().collect::<Vec<_>>()
            );
            assert_eq!(a.mol.stereo_groups(), b.mol.stereo_groups());
            assert_eq!(
                a.coords
                    .iter()
                    .map(|(x, y)| (x.to_bits(), y.to_bits()))
                    .collect::<Vec<_>>(),
                b.coords
                    .iter()
                    .map(|(x, y)| (x.to_bits(), y.to_bits()))
                    .collect::<Vec<_>>()
            );
            for (i, _) in a.mol.atoms() {
                assert_eq!(a.mol.neighbor_slice(i), b.mol.neighbor_slice(i));
                assert_eq!(
                    a.mol.stereo_neighbor_order(i),
                    b.mol.stereo_neighbor_order(i)
                );
                assert_eq!(a.mol.r_group_label(i), b.mol.r_group_label(i));
                assert_eq!(a.mol.atom_tag(i), b.mol.atom_tag(i));
            }
            for (i, _) in a.mol.bonds() {
                assert_eq!(a.mol.bond_direction(i), b.mol.bond_direction(i));
                assert_eq!(
                    a.mol.bond_direction_anchor(i),
                    b.mol.bond_direction_anchor(i)
                );
            }
            for i in 0..a.mol.stereo_groups().len() {
                assert_eq!(
                    a.mol.stereo_group_bond_indices(i),
                    b.mol.stereo_group_bond_indices(i)
                );
            }
            assert_eq!(
                format!(
                    "{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}",
                    a.metadata,
                    a.stereo_diagnostics,
                    a.ez_diagnostics,
                    a.coordinate_dimension,
                    a.geometry_rank,
                    a.stereo3d_diagnostics,
                    a.square_planar_diagnostics,
                    a.conformer
                ),
                format!(
                    "{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}:{:?}",
                    b.metadata,
                    b.stereo_diagnostics,
                    b.ez_diagnostics,
                    b.coordinate_dimension,
                    b.geometry_rank,
                    b.stereo3d_diagnostics,
                    b.square_planar_diagnostics,
                    b.conformer
                )
            );
            assert_eq!(
                write_mol_v3000(&a.mol, &a.metadata, &a.coords),
                write_mol_v3000(&b.mol, &b.metadata, &b.coords)
            );
            true
        }
        (Err(a), Err(b)) => {
            assert_eq!(a.to_string(), b.to_string());
            false
        }
        _ => panic!("bond attributes changed parse success or error boundary"),
    }
}

#[test]
fn pinned_v3000_reports_and_rewrite_bytes_match_the_frozen_bond_reader() {
    use std::io::Read;
    let corpora: [&[u8]; 3] = [
        include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/benchmark_corpora/openbabel-v3000-speed-slice0-1000-v1.molblocks.gz"
        )),
        include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/benchmark_corpora/openbabel-v3000-speed-slice5000-1000-v1.molblocks.gz"
        )),
        include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/benchmark_corpora/openbabel-v3000-speed-chains64-v1.molblocks.gz"
        )),
    ];
    let mut rows = 0;
    for bytes in corpora {
        let mut text = String::new();
        flate2::read::GzDecoder::new(bytes)
            .read_to_string(&mut text)
            .unwrap();
        for block in text
            .split("$$$$\n")
            .filter(|block| !block.trim().is_empty())
        {
            assert!(compare_read(block));
            rows += 1;
        }
    }
    assert_eq!(rows, 2064);
    eprintln!("frozen V3000 bond reader: {rows} complete reports and rewrites matched");
}

fn fixture(bond_rows: &str, ids: [u32; 3], declared_bonds: usize) -> String {
    format!(
        "attributes\n  chematic          2D\nmetadata oracle\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 3 {declared_bonds} 0 0 0\nM  V30 BEGIN ATOM\nM  V30 {} C -0.0000 0.0 0.0 0\nM  V30 {} C 1.5 0.0 0.0 0\nM  V30 {} O 3.0 1.0 0.0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\n{bond_rows}M  V30 END BOND\nM  V30 END CTAB\nM  END\n",
        ids[0], ids[1], ids[2]
    )
}

#[test]
fn absent_cfg_opaque_sparse_and_malformed_attributes_keep_reports_and_errors() {
    let tails = [
        "",
        " CFG=1",
        " CFG=2",
        " CFG=3",
        " CFG=oops",
        " CFG=2 CFG=1",
        " ATTACH=ANY ENDPTS=(3 1 2 3)",
        " VENDOR=\"left right\" UNKNOWN=α",
        " ENDPTS=(2 55 7) ATTACH=ANY",
        " ENDPTS=(2 1)",
        " ENDPTS=(bad 1)",
        " ENDPTS=(2 1 99)",
        " ENDPTS=(2 1 2",
        " ENDPTS=(0 )",
        " CFG=2 ENDPTS=(2 1 2)",
    ];
    let mut successes = 0;
    let mut failures = 0;
    for ids in [[1u32, 2, 3], [55, 7, 101]] {
        for kind in 0..=12 {
            for tail in tails {
                let rows = format!(
                    "M  V30 1 {kind} {} {}{tail}\nM  V30 2 1 {} {}\n",
                    ids[0], ids[1], ids[1], ids[2]
                );
                let input = fixture(&rows, ids, 2);
                if compare_read(&input) {
                    successes += 1;
                } else {
                    failures += 1;
                }
                for bad in [
                    input.replace("M  V30 2 1", "M  V30 1 1"),
                    input.replace("COUNTS 3 2", "COUNTS 3 3"),
                    input.replace("M  V30 2 1", "M  V30 2 nope"),
                    input.replace("M  V30 2 1", "M  V30 2 1 99999999"),
                ] {
                    assert!(!compare_read(&bad));
                    failures += 1;
                }
            }
        }
    }
    assert!(successes > 0 && failures > 0);
    eprintln!(
        "bond-attribute report/error oracle: {successes} successes, {failures} errors matched"
    );
}

#[test]
fn moved_endpoint_properties_match_the_complete_original_string_transform() {
    let maps = [
        vec![],
        vec![(1, AtomIdx(0)), (2, AtomIdx(1)), (3, AtomIdx(2))],
        vec![(55, AtomIdx(0)), (7, AtomIdx(1)), (101, AtomIdx(2))],
        vec![
            (0, AtomIdx(2)),
            (7, AtomIdx(0)),
            (7, AtomIdx(1)),
            (u32::MAX, AtomIdx(3)),
        ],
    ];
    let pieces = [
        "",
        "CFG=2",
        " ENDPTS=(",
        "ENDPTS=(2 1 2)",
        "ENDPTS=(2 55 7)",
        "ENDPTS=(0 )",
        "ENDPTS=(2 0 7)",
        "ENDPTS=(2 7 4294967295)",
        "ENDPTS=(2 7 7)",
        "ENDPTS=(x 1)",
        "ENDPTS=(2 1)",
        "ENDPTS=(1 x)",
        "ENDPTS=(1 1)",
        "ENDPTS=(2\u{2003}1\t2)",
        ")",
        " ATTACH=ANY",
        " VENDOR=α",
        "ENDPTS=(2 1 2) ENDPTS=(1 3)",
    ];
    let mut checked = 0;
    for map in &maps {
        for a in pieces {
            for b in pieces {
                for c in pieces {
                    let value = format!("{a}{b}{c}");
                    assert_eq!(
                        super::normalize_v3000_endpoints(value.clone(), map),
                        normalize_v3000_endpoints(&value, map)
                    );
                    checked += 1;
                }
            }
        }
    }
    assert_eq!(checked, 23328);
    eprintln!("endpoint normalization oracle: {checked} strings matched");
}
