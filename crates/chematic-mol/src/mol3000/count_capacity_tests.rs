//! Frozen COUNTS transition oracle; see the count-capacity fixture manifest.
use super::*;

fn reference_in_ctab(
    parser: &mut V3000Parser,
    tokens: &[&str],
    line: usize,
) -> Result<(), MolParseError> {
    if tokens[0] == "COUNTS" {
        if tokens.len() < 4 {
            return Err(v3k_err(line, "COUNTS line has fewer than 3 values"));
        }
        parser.expected_atoms = parse_v3000_count(tokens[1], line, "atom")?;
        parser.expected_bonds = parse_v3000_count(tokens[2], line, "bond")?;
        parser.expected_sgroups = Some(parse_v3000_count(tokens[3], line, "SGROUP")?);
    } else if is_marker(tokens, "BEGIN", "ATOM") {
        parser.state = V3000State::InAtomBlock;
    } else if is_marker(tokens, "END", "CTAB") {
        parser.state = V3000State::Done;
    }
    Ok(())
}

fn recipe() -> serde_json::Value {
    serde_json::from_str(include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/benchmark_corpora/v3000-count-capacity-v2.recipe.json"
    )))
    .unwrap()
}

fn compare_counts(tokens: &[&str], limit: usize, line: usize) {
    let mut expected = V3000Parser::new();
    let mut actual = V3000Parser::new();
    expected.allocation_limit = limit;
    actual.allocation_limit = limit;
    let prime = ["COUNTS", "7", "9", "11"];
    reference_in_ctab(&mut expected, &prime, 1).unwrap();
    actual.in_ctab(&prime, 1).unwrap();
    let left = reference_in_ctab(&mut expected, tokens, line);
    let right = actual.in_ctab(tokens, line);
    assert_eq!(format!("{left:?}"), format!("{right:?}"));
    assert_eq!(actual.expected_atoms, expected.expected_atoms);
    assert_eq!(actual.expected_bonds, expected.expected_bonds);
    assert_eq!(actual.expected_sgroups, expected.expected_sgroups);
    assert_eq!(
        std::mem::discriminant(&actual.state),
        std::mem::discriminant(&expected.state)
    );
    assert!(actual.coords.capacity() <= limit);
    assert!(actual.raw_z.capacity() <= limit);
    assert!(actual.atom_idx_map.entries.capacity() <= limit);
    assert!(actual.bond_idx_map.entries.capacity() <= limit);
    assert_eq!(actual.builder.atom_count(), 0);
}

#[test]
fn count_transitions_errors_and_allocation_bounds_match_original() {
    let recipe = recipe();
    let values: Vec<_> = recipe["count_lexemes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap())
        .collect();
    let mut cases = 0;
    for limit in recipe["allocation_limits"].as_array().unwrap() {
        for line in recipe["source_lines"].as_array().unwrap() {
            for &a in &values {
                for &b in &values {
                    for &sg in &values {
                        compare_counts(
                            &["COUNTS", a, b, sg],
                            limit.as_u64().unwrap() as usize,
                            line.as_u64().unwrap() as usize,
                        );
                        cases += 1;
                    }
                }
            }
        }
    }
    assert_eq!(cases, 27000);
}

#[test]
fn incomplete_counts_and_block_markers_match_original() {
    let recipe = recipe();
    for row in recipe["markers"].as_array().unwrap() {
        let tokens: Vec<_> = row
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();
        // consume() skips empty rows before dispatching to in_ctab().
        if !tokens.is_empty() {
            compare_counts(&tokens, 5, 42);
        }
    }
}

#[test]
fn allocation_waits_for_atom_block_and_is_bounded_by_real_rows() {
    let recipe = recipe();
    for limit in recipe["allocation_limits"].as_array().unwrap() {
        let limit = limit.as_u64().unwrap() as usize;
        let mut parser = V3000Parser::new();
        parser.allocation_limit = limit;
        let large = usize::MAX.to_string();
        for _ in 0..16 {
            parser.in_ctab(&["COUNTS", &large, &large, "0"], 1).unwrap();
            assert_eq!(parser.coords.capacity(), 0);
            assert_eq!(parser.raw_z.capacity(), 0);
            assert_eq!(parser.atom_idx_map.entries.capacity(), 0);
            assert_eq!(parser.bond_idx_map.entries.capacity(), 0);
        }
        parser.in_ctab(&["BEGIN", "ATOM"], 2).unwrap();
        assert_eq!(parser.coords.capacity(), limit);
        assert_eq!(parser.raw_z.capacity(), limit);
        assert_eq!(parser.atom_idx_map.entries.capacity(), limit);
        assert_eq!(parser.bond_idx_map.entries.capacity(), limit);
    }
}

#[test]
fn rounded_capacity_matches_original_growth_with_sufficient_row_budget() {
    for count in 0..=1024 {
        let mut original = Vec::<(f64, f64)>::new();
        for _ in 0..count {
            original.push((0.0, 0.0));
        }
        assert_eq!(
            v3000_storage_capacity(count, count * 2 + 7),
            original.capacity()
        );
    }
    assert_eq!(v3000_storage_capacity(usize::MAX, 0), 0);
    assert_eq!(v3000_storage_capacity(usize::MAX, 5), 5);
    assert_eq!(v3000_storage_capacity(usize::MAX, usize::MAX), usize::MAX);
}
