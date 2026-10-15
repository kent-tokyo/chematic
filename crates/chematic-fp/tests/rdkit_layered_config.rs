//! Original minimal tests; manifest: validation/rdkit_issues/fingerprints/layered-config.json.
use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};
use chematic_fp::{RdkitLayeredConfig as Config, rdkit_layered_config, rdkit_layered_fp};
fn mol(s: &str) -> chematic_core::Molecule {
    chematic_smiles::parse(s).unwrap()
}
#[test]
fn defaults_preserve_existing_packed_fingerprints() {
    for s in [
        "CCO",
        "CC(C)O",
        "C1CC1",
        "c1ccccc1O",
        "c1ccc2ccccc2c1",
        "[2H]C([2H])([2H])NC=O",
        "C1C2CC3CC1CC(C2)C3",
    ] {
        let m = mol(s);
        let r = diagnostic(rdkit_layered_config(&m, &Config::default())).unwrap();
        assert_eq!(r.bits, rdkit_layered_fp(&m).to_le_bytes());
    }
}
#[test]
fn root_multiplicity_and_seeded_counter_contract() {
    let m = mol("CC(C)O");
    let mut c = Config {
        from_atoms: vec![AtomIdx(0)],
        ..Default::default()
    };
    let a = diagnostic(rdkit_layered_config(&m, &c)).unwrap();
    assert_eq!(a.atom_counts, [4, 4, 2, 2]);
    c.from_atoms.push(AtomIdx(0));
    c.atom_counts = Some(vec![10, 20, 30, 40, 99]);
    let b = diagnostic(rdkit_layered_config(&m, &c)).unwrap();
    assert_eq!(b.atom_counts, [18, 28, 34, 44, 99]);
    assert_eq!(b.bits, a.bits);
    assert_eq!(c.atom_counts.unwrap(), [10, 20, 30, 40, 99]);
}
#[test]
fn bit_mask_controls_counters_even_after_hash_collisions() {
    let m = mol("CC(C)O");
    let mut c = Config {
        fp_size: 1,
        ..Default::default()
    };
    let r = diagnostic(rdkit_layered_config(&m, &c)).unwrap();
    assert_eq!(r.bits, [1]);
    assert_eq!(r.atom_counts, [4, 7, 4, 4]);
    c.set_only_bits = Some(vec![0]);
    let r = diagnostic(rdkit_layered_config(&m, &c)).unwrap();
    assert_eq!(r.bits, [0]);
    assert_eq!(r.atom_counts, [0; 4]);
}
#[test]
fn selected_layers_union_and_zero_layers() {
    let m = mol("c1ccccc1O");
    let all = diagnostic(rdkit_layered_config(&m, &Config::default())).unwrap();
    let mut union = vec![0; 256];
    for layer in 0..6 {
        let c = Config {
            layer_flags: 1 << layer,
            ..Default::default()
        };
        let r = diagnostic(rdkit_layered_config(&m, &c)).unwrap();
        for (u, b) in union.iter_mut().zip(r.bits) {
            *u |= b;
        }
    }
    assert_eq!(union, all.bits);
    let c = Config {
        layer_flags: 0,
        atom_counts: Some(vec![9; 7]),
        ..Default::default()
    };
    let r = diagnostic(rdkit_layered_config(&m, &c)).unwrap();
    assert!(r.bits.iter().all(|b| *b == 0));
    assert_eq!(r.atom_counts, [9; 7]);
}
#[test]
fn rooted_linear_uses_endpoint_semantics_and_final_range_closures() {
    let m = mol("CC(C)O");
    let c = Config {
        branched_paths: false,
        from_atoms: vec![AtomIdx(1)],
        ..Default::default()
    };
    let r = diagnostic(rdkit_layered_config(&m, &c)).unwrap();
    assert_eq!(r.atom_counts, [1, 3, 1, 1]);
    let m = mol("C1CC1");
    let a = Config {
        branched_paths: false,
        from_atoms: vec![AtomIdx(0)],
        max_path: 3,
        ..Default::default()
    };
    let mut b = a.clone();
    b.max_path = 7;
    assert!(
        diagnostic(rdkit_layered_config(&m, &a))
            .unwrap()
            .selected_paths
            > diagnostic(rdkit_layered_config(&m, &b))
                .unwrap()
                .selected_paths
    );
}
#[test]
fn counters_wrap_like_unsigned_rdkit_counters() {
    let c = Config {
        atom_counts: Some(vec![u32::MAX; 3]),
        ..Default::default()
    };
    assert_eq!(
        diagnostic(rdkit_layered_config(&mol("CCO"), &c))
            .unwrap()
            .atom_counts,
        [1, 2, 1]
    );
}
#[test]
fn invalid_options_and_known_invalid_upstream_branch_are_checked() {
    let m = mol("CCO");
    for c in [
        Config {
            fp_size: 0,
            ..Default::default()
        },
        Config {
            min_path: 0,
            ..Default::default()
        },
        Config {
            min_path: 3,
            max_path: 2,
            ..Default::default()
        },
        Config {
            max_path: 13,
            ..Default::default()
        },
        Config {
            from_atoms: vec![AtomIdx(3)],
            ..Default::default()
        },
        Config {
            atom_counts: Some(vec![0; 2]),
            ..Default::default()
        },
        Config {
            fp_size: 9,
            set_only_bits: Some(vec![0, 2]),
            ..Default::default()
        },
        Config {
            set_only_bits: Some(vec![]),
            ..Default::default()
        },
    ] {
        assert!(diagnostic(rdkit_layered_config(&m, &c)).is_err());
    }
    let c = Config {
        branched_paths: false,
        ..Default::default()
    };
    assert_eq!(
        diagnostic(rdkit_layered_config(&m, &c)).unwrap_err().code,
        "unrooted_linear_unsupported"
    );
}
#[test]
fn exhausted_budgets_return_no_partial_mutation() {
    let m = mol("CC(C)O");
    let mut c = Config {
        atom_counts: Some(vec![7; 4]),
        ..Default::default()
    };
    let result = diagnostic(rdkit_layered_config(&m, &c)).unwrap();
    c.max_visits = result.search_visits - 1;
    assert_eq!(
        diagnostic(rdkit_layered_config(&m, &c)).unwrap_err().code,
        "visit_limit"
    );
    assert_eq!(c.atom_counts.as_ref().unwrap(), &[7; 4]);
    c.max_visits = 1_000_000;
    c.max_paths = 1;
    assert_eq!(
        diagnostic(rdkit_layered_config(&m, &c)).unwrap_err().code,
        "path_limit"
    );
}
#[test]
fn query_zero_and_dative_bonds_are_rejected() {
    for order in [BondOrder::QueryAny, BondOrder::Zero, BondOrder::Dative] {
        let mut b = MoleculeBuilder::new();
        let a = b.add_atom(Atom::new(Element::C));
        let z = b.add_atom(Atom::new(Element::C));
        b.add_bond(a, z, order).unwrap();
        assert_eq!(
            diagnostic(rdkit_layered_config(&b.build(), &Config::default()))
                .unwrap_err()
                .code,
            "unsupported_bond"
        );
    }
}
#[test]
fn explicit_hydrogen_roots_keep_original_indices_but_set_no_bits() {
    let m = mol("[2H]OC");
    let c = Config {
        from_atoms: vec![AtomIdx(0)],
        atom_counts: Some(vec![7, 8, 9]),
        ..Default::default()
    };
    let r = diagnostic(rdkit_layered_config(&m, &c)).unwrap();
    assert!(r.bits.iter().all(|b| *b == 0));
    assert_eq!(r.atom_counts, [7, 8, 9]);
}

#[test]
fn rooted_linear_frontier_and_visit_caps_fail_atomically() {
    let m = mol("CC(C)O");
    let mut c = Config {
        branched_paths: false,
        from_atoms: vec![AtomIdx(1)],
        atom_counts: Some(vec![7; 4]),
        max_paths: 1,
        ..Default::default()
    };
    assert_eq!(
        diagnostic(rdkit_layered_config(&m, &c)).unwrap_err().code,
        "frontier_limit"
    );
    c.max_paths = 100000;
    c.max_visits = 1;
    assert_eq!(
        diagnostic(rdkit_layered_config(&m, &c)).unwrap_err().code,
        "visit_limit"
    );
    assert_eq!(c.atom_counts.unwrap(), [7; 4]);
}
#[test]
fn graph_limits_are_checked_before_perception() {
    let mut b = MoleculeBuilder::new();
    for _ in 0..513 {
        b.add_atom(Atom::new(Element::C));
    }
    assert_eq!(
        diagnostic(rdkit_layered_config(&b.build(), &Config::default()))
            .unwrap_err()
            .code,
        "graph_limit"
    );
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
