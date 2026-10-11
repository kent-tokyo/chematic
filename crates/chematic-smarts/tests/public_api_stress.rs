//! Exercise matching modes over the checked-in special-chemistry corpus.

use chematic_smarts::*;

const STRESS: &str = include_str!("../../../validation/cosmolkit_comparison/corpora/stress-v1.smi");

#[test]
fn matching_surfaces_handle_special_chemistry() {
    let queries: Vec<_> = [
        "*",
        "[#6]",
        "[!#6]",
        "[C,c]",
        "[N;H0,H1]",
        "[+,-]",
        "[13C]",
        "[R]",
        "[R0]",
        "[R1]",
        "[r5]",
        "[r6]",
        "[x2]",
        "[D3]",
        "[v4]",
        "[^1]",
        "[^2]",
        "[$([N;H1](C)C)]",
        "C@C",
        "C!@C",
        "C~[O,N,S]",
        "c1ccccc1",
        "[C:1](=[O:2])[O,N:3]",
    ]
    .into_iter()
    .map(|source| parse_smarts(source).unwrap())
    .collect();
    let mut exercised = 0usize;
    let mut mcs_pool = Vec::new();

    for line in STRESS.lines() {
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = chematic_smiles::parse(smiles) else {
            continue;
        };
        if mol.atom_count() == 0 || mol.atom_count() > 120 {
            continue;
        }
        if mcs_pool.len() < 128 {
            mcs_pool.push(mol.clone());
        }
        let rings = chematic_perception::find_sssr(&mol);
        for query in &queries {
            let config = MatchConfig {
                max_matches: Some(32),
                max_visit_budget: Some(100_000),
                ..Default::default()
            };
            let _ = find_matches(query, &mol);
            let _ = find_matches_with_config(query, &mol, &config);
            let _ = has_match_with_config(query, &mol, &config);
            let _ = has_match_bounded(query, &mol, &rings, &config);
            let _ = find_matches_perceived(query, &mol, &config);
            let _ = find_match_atom_sets_perceived(query, &mol, &config);
            let _ = has_match_perceived(query, &mol, &config);

            let parity = RdkitParityConfig {
                base: config,
                ..Default::default()
            };
            let _ = find_matches_rdkit_parity(query, &mol, &parity);
            let _ = find_match_atom_sets_rdkit_parity(query, &mol, &parity);
            let _ = has_match_rdkit_parity_bounded(query, &mol, &parity);
        }
        exercised += 1;
    }

    for pair in mcs_pool.windows(2) {
        let refs = [&pair[0], &pair[1]];
        let _ = find_mcs(&refs);
        let config = McsConfig {
            timeout_ms: Some(50),
            max_search_nodes: Some(20_000),
            ..Default::default()
        };
        let _ = find_mcs_with_config_checked(&refs, &config);
    }
    assert!(
        exercised >= 4_000,
        "stress corpus unexpectedly shrank: {exercised}"
    );
}
