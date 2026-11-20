use chematic_core::AtomIdx;
use chematic_fp::{RdkitRdkConfig as Config, RdkitRdkConfigError as Error, rdkit_rdk_config};
use chematic_smiles::parse;

#[test]
fn linear_ring_closures_depend_on_the_maximum_path_length() {
    let mol = parse("C1CC1C").unwrap();
    let mut config = Config {
        min_path: 3,
        max_path: 3,
        nbits: 31,
        bits_per_feature: 1,
        branched_paths: false,
        ..Default::default()
    };
    // RDKit 2026.09.2 counts the triangle once and both open three-edge
    // paths. Rooting at the terminal carbon drops the triangle.
    let all = rdkit_rdk_config(&mol, &config).unwrap();
    assert_eq!(all.sparse_counts, [(2172292514, 1), (2911990635, 2)]);
    config.from_atoms = vec![AtomIdx(3), AtomIdx(3)];
    let rooted = rdkit_rdk_config(&mol, &config).unwrap();
    assert_eq!(rooted.sparse_counts, [(2911990635, 4)]);
    config.max_path = 4;
    let longer = rdkit_rdk_config(&mol, &config).unwrap();
    assert_eq!(longer.sparse_counts, [(2911990635, 4), (3011716336, 2)]);
}

#[test]
fn configured_rng_and_count_simulation_match_rdkit() {
    let mol = parse("C=C").unwrap();
    let config = Config {
        nbits: 13,
        bits_per_feature: 7,
        use_bond_order: false,
        atom_invariants: Some(vec![17, 31]),
        count_bounds: Some(vec![1, 2]),
        ..Default::default()
    };
    let result = rdkit_rdk_config(&mol, &config).unwrap();
    assert_eq!(result.fingerprint.to_le_bytes(), [0x45, 0x0f]);
    assert_eq!(
        result.folded_counts,
        [(0, 1), (5, 2), (6, 1), (9, 1), (10, 2)]
    );
    assert_eq!(
        result.sparse_counts,
        [
            (25469856, 1),
            (345776735, 1),
            (654385672, 1),
            (727540759, 1),
            (1113457810, 1),
            (1483898877, 1),
            (2306184725, 1),
        ]
    );
}

#[test]
fn skips_query_bonds_and_applies_ignore_atoms_to_both_path_modes() {
    let mol = parse("C~CC").unwrap();
    let mut config = Config {
        nbits: 17,
        bits_per_feature: 1,
        ..Default::default()
    };
    for branched_paths in [false, true] {
        config.branched_paths = branched_paths;
        assert_eq!(
            rdkit_rdk_config(&mol, &config).unwrap().sparse_counts,
            [(4275705116, 1)]
        );
        config.ignore_atoms = vec![AtomIdx(2)];
        assert!(
            rdkit_rdk_config(&mol, &config)
                .unwrap()
                .sparse_counts
                .is_empty()
        );
        config.ignore_atoms.clear();
    }
}

#[test]
fn malformed_indices_and_path_budget_return_errors() {
    let mol = parse("CCO").unwrap();
    let mut config = Config {
        from_atoms: vec![AtomIdx(3)],
        ..Default::default()
    };
    assert_eq!(
        rdkit_rdk_config(&mol, &config).unwrap_err(),
        Error::AtomIndices
    );
    config.from_atoms.clear();
    config.max_visits = 1;
    assert_eq!(
        rdkit_rdk_config(&mol, &config).unwrap_err(),
        Error::VisitLimit
    );
    config.max_visits = 100;
    config.max_paths = 1;
    assert_eq!(
        rdkit_rdk_config(&mol, &config).unwrap_err(),
        Error::PathLimit
    );
}
