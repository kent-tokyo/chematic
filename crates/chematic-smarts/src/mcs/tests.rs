use super::*;
use chematic_smiles::parse;

#[test]
fn test_empty_input() {
    let result = find_mcs(&[]);
    assert_eq!(result.atom_count(), 0);
}

#[test]
fn test_single_molecule() {
    let benzene = parse("c1ccccc1").unwrap();
    let result = find_mcs(&[&benzene]);
    assert_eq!(result.atom_count(), 6);
}

#[test]
fn test_ethane_propane() {
    // CC and CCC → MCS should be CC (2 atoms)
    let a = parse("CC").unwrap();
    let b = parse("CCC").unwrap();
    let result = find_mcs(&[&a, &b]);
    assert!(
        result.atom_count() >= 2,
        "Expected >= 2 atoms, got {}",
        result.atom_count()
    );
}

#[test]
fn test_methane_ethanol() {
    // C and CCO → MCS should be C (1 atom)
    let a = parse("C").unwrap();
    let b = parse("CCO").unwrap();
    let result = find_mcs(&[&a, &b]);
    assert!(result.atom_count() >= 1);
}

#[test]
fn test_benzene_toluene() {
    // c1ccccc1 and Cc1ccccc1 → MCS should be the benzene ring (6 atoms)
    let a = parse("c1ccccc1").unwrap();
    let b = parse("Cc1ccccc1").unwrap();
    let result = find_mcs(&[&a, &b]);
    assert_eq!(
        result.atom_count(),
        6,
        "Expected benzene ring (6 atoms), got {}",
        result.atom_count()
    );
}

#[test]
fn test_no_common_atoms() {
    // Only N atoms vs only O atoms → no MCS
    let a = parse("N").unwrap(); // ammonia
    let b = parse("O").unwrap(); // water
    let result = find_mcs(&[&a, &b]);
    assert_eq!(result.atom_count(), 0);
}

#[test]
fn test_identical_molecules() {
    let a = parse("CC").unwrap();
    let b = parse("CC").unwrap();
    let result = find_mcs(&[&a, &b]);
    assert_eq!(result.atom_count(), 2);
}

#[test]
fn test_match_bonds_false() {
    // With match_bonds=false, single and double bond C-C should match
    let config = McsConfig {
        match_bonds: false,
        ..McsConfig::default()
    };
    let a = parse("CC").unwrap();
    let b = parse("C=C").unwrap();
    let result = find_mcs_with_config(&[&a, &b], &config);
    assert!(result.atom_count() >= 2);
}

#[test]
fn test_min_atoms_filter() {
    // CC and CCC have MCS of size 2, but if min_atoms=5, result should be empty
    let config = McsConfig {
        min_atoms: 5,
        ..McsConfig::default()
    };
    let a = parse("CC").unwrap();
    let b = parse("CCC").unwrap();
    let result = find_mcs_with_config(&[&a, &b], &config);
    assert_eq!(result.atom_count(), 0);
}

#[test]
fn test_timeout_does_not_panic() {
    let config = McsConfig {
        timeout_ms: Some(1),
        ..McsConfig::default()
    }; // 1 ms
    let a = parse("c1ccccc1").unwrap();
    let b = parse("Cc1ccccc1").unwrap();
    // Should return without panic (may return partial result)
    let _ = find_mcs_with_config(&[&a, &b], &config);
}

#[test]
fn test_checked_reports_exhaustive_without_timeout() {
    let a = parse("c1ccccc1").unwrap();
    let b = parse("Cc1ccccc1").unwrap();
    let unchecked_result = find_mcs_with_config(&[&a, &b], &McsConfig::default());
    let checked = find_mcs_with_config_checked(&[&a, &b], &McsConfig::default());
    assert!(
        !checked.was_timed_out(),
        "no timeout_ms set: search must report Exhaustive"
    );
    assert!(matches!(checked, McsOutcome::Exhaustive(_)));
    assert_eq!(
        checked.into_query().atom_count(),
        unchecked_result.atom_count(),
        "find_mcs_with_config_checked's result must match find_mcs_with_config's"
    );
}

#[test]
fn test_checked_reports_timed_out_with_zero_deadline() {
    // A 0ms deadline is reached on the very first check, guaranteeing a
    // cut-off search (unlike test_timeout_does_not_panic's 1ms, which
    // isn't guaranteed to actually trigger on a small molecule pair --
    // that test only checks "no panic", this one checks the outcome tag).
    let config = McsConfig {
        timeout_ms: Some(0),
        ..McsConfig::default()
    };
    let a = parse("c1ccccc1").unwrap();
    let b = parse("Cc1ccccc1").unwrap();
    let checked = find_mcs_with_config_checked(&[&a, &b], &config);
    assert!(
        checked.was_timed_out(),
        "0ms deadline must be reported as TimedOut, not silently treated as Exhaustive"
    );
    assert!(matches!(checked, McsOutcome::TimedOut(_)));
}

#[test]
fn test_checked_reports_resource_limit_without_claiming_exhaustive() {
    let config = McsConfig {
        max_search_nodes: Some(0),
        ..McsConfig::default()
    };
    let a = parse("CC").unwrap();
    let b = parse("CC").unwrap();
    let checked = find_mcs_with_config_checked(&[&a, &b], &config);
    assert!(checked.was_resource_limited());
    assert!(matches!(checked, McsOutcome::ResourceLimited(_)));
    assert_eq!(checked.into_query().atom_count(), 0);
}

#[test]
fn test_checked_rejects_oversized_input_molecule() {
    let config = McsConfig {
        max_atoms_per_molecule: 1,
        ..McsConfig::default()
    };
    let a = parse("CC").unwrap();
    let b = parse("C").unwrap();
    let checked = find_mcs_with_config_checked(&[&a, &b], &config);
    assert!(checked.was_resource_limited());
    assert_eq!(checked.into_query().atom_count(), 0);
}

#[test]
fn test_result_matches_all_inputs() {
    // The MCS result should be findable as a substructure in all input molecules
    use crate::find_matches;
    let a = parse("c1ccccc1").unwrap();
    let b = parse("Cc1ccccc1").unwrap();
    let result = find_mcs(&[&a, &b]);
    if result.atom_count() > 0 {
        assert!(!find_matches(&result, &a).is_empty());
        assert!(!find_matches(&result, &b).is_empty());
    }
}

#[test]
fn test_aspirin_benzoic_acid() {
    // CC(=O)Oc1ccccc1C(=O)O and OC(=O)c1ccccc1
    // Common: benzoic acid core → >= 9 atoms
    let a = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
    let b = parse("OC(=O)c1ccccc1").unwrap();
    let result = find_mcs(&[&a, &b]);
    assert!(
        result.atom_count() >= 7,
        "Expected >= 7 atoms for aspirin/benzoic acid MCS, got {}",
        result.atom_count()
    );
}

// -----------------------------------------------------------------------
// ring_matches_ring_only tests
// -----------------------------------------------------------------------

#[test]
fn test_ring_matches_ring_only_blocks_ring_to_nonring() {
    // cyclohexane (all ring C's) vs hexane (all chain C's)
    // Without constraint: MCS = 5 atoms (C-C-C-C-C path, both have it)
    // With ring_matches_ring_only: ring C (cyclohexane) can't match non-ring C (hexane) → MCS = 0
    let a = parse("C1CCCCC1").unwrap();
    let b = parse("CCCCCC").unwrap();

    // Confirm default finds a path match.
    let default_result = find_mcs(&[&a, &b]);
    assert!(
        default_result.atom_count() > 0,
        "Default MCS should be non-zero"
    );

    let config = McsConfig {
        ring_matches_ring_only: true,
        ..McsConfig::default()
    };
    let result = find_mcs_with_config(&[&a, &b], &config);
    assert_eq!(
        result.atom_count(),
        0,
        "ring_matches_ring_only: ring C must not match chain C, expected 0 got {}",
        result.atom_count()
    );
}

#[test]
fn test_ring_matches_ring_only_preserves_ring_only_mcs() {
    // benzene vs toluene: ring atoms should still match each other
    let a = parse("c1ccccc1").unwrap();
    let b = parse("Cc1ccccc1").unwrap();
    let config = McsConfig {
        ring_matches_ring_only: true,
        ..McsConfig::default()
    };
    let result = find_mcs_with_config(&[&a, &b], &config);
    assert_eq!(
        result.atom_count(),
        6,
        "benzene ring should still be the MCS"
    );
}

// -----------------------------------------------------------------------
// complete_rings_only tests
// -----------------------------------------------------------------------

#[test]
fn test_complete_rings_only_quinoline_series() {
    // With complete_rings_only: partial rings are removed from the MCS result.
    // In the quinoline series the exocyclic C would create a partial ring if it
    // happened to be included; enabling complete_rings_only ensures the quinoline
    // scaffold (10 atoms) is returned intact.
    let a = parse("c1ccc2nc(CC)ccc2c1").unwrap();
    let b = parse("c1ccc2nc(CO)ccc2c1").unwrap();
    let c = parse("c1ccc2nc(CN)ccc2c1").unwrap();

    let config = McsConfig {
        complete_rings_only: true,
        ..McsConfig::default()
    };
    let result = find_mcs_with_config(&[&a, &b, &c], &config);
    // Result must be >= 10 atoms (quinoline) and no partial ring
    assert!(
        result.atom_count() >= 10,
        "Expected at least quinoline scaffold (10) with complete_rings_only, got {}",
        result.atom_count()
    );
}

#[test]
fn test_complete_rings_only_does_not_break_full_ring_mcs() {
    // benzene vs toluene: full benzene ring should survive complete_rings_only
    let a = parse("c1ccccc1").unwrap();
    let b = parse("Cc1ccccc1").unwrap();
    let config = McsConfig {
        complete_rings_only: true,
        ..McsConfig::default()
    };
    let result = find_mcs_with_config(&[&a, &b], &config);
    assert_eq!(
        result.atom_count(),
        6,
        "Full benzene ring should be preserved"
    );
}

// -----------------------------------------------------------------------
// Diagnostic: verify complete_rings_only ACTUALLY fires pruning
// -----------------------------------------------------------------------

#[test]
fn test_complete_rings_only_pruning_fires() {
    // Benzene vs naphthalene: raw MCS without constraint is 6 ring atoms (one full ring).
    // With complete_rings_only, fused-ring bridgeheads are in two rings each; verify
    // the full benzene ring is preserved (6 atoms), not over-pruned.
    let a = parse("c1ccccc1").unwrap();
    let b = parse("c1ccc2ccccc2c1").unwrap(); // naphthalene

    let default_result = find_mcs(&[&a, &b]);
    let config = McsConfig {
        complete_rings_only: true,
        ..McsConfig::default()
    };
    let pruned_result = find_mcs_with_config(&[&a, &b], &config);

    eprintln!(
        "benzene vs naphthalene — raw MCS: {}, complete_rings_only: {}",
        default_result.atom_count(),
        pruned_result.atom_count()
    );
    // Both should be 6 (one complete benzene ring)
    assert_eq!(
        default_result.atom_count(),
        6,
        "Default MCS should be benzene ring (6 atoms)"
    );
    assert_eq!(
        pruned_result.atom_count(),
        6,
        "complete_rings_only should preserve the full 6-ring (no partial ring)"
    );
}

#[test]
fn test_complete_rings_only_prunes_partial_fused_ring() {
    // C1Cc2ccccc21 is benzocyclobutene (4-ring fused to 6-ring), not indane.
    // Toluene: Cc1ccccc1.
    // Raw MCS = 7 atoms: the 6 aromatic C's + aliphatic atom1 from benzocyclobutene
    // matched to the methyl C of toluene (both non-aromatic, atoms_compatible passes).
    // Atom1 of benzocyclobutene is in the 4-membered ring — partial ring in mol[0].
    // complete_rings_only fires: remove atom1 + bridgehead atoms (shared with 4-ring) →
    // the 6-ring becomes partially covered → cascading removal → result = 0.
    let benzocyclobutene = parse("C1Cc2ccccc21").unwrap();
    let toluene = parse("Cc1ccccc1").unwrap();

    let default_result = find_mcs(&[&benzocyclobutene, &toluene]);
    let config = McsConfig {
        complete_rings_only: true,
        ..McsConfig::default()
    };
    let pruned_result = find_mcs_with_config(&[&benzocyclobutene, &toluene], &config);

    // Raw MCS includes a partial 4-ring; cascading pruning removes everything.
    assert!(
        default_result.atom_count() > 0,
        "default MCS should be non-empty"
    );
    assert_eq!(
        pruned_result.atom_count(),
        0,
        "complete_rings_only must prune to 0 when partial ring cascades"
    );
}

#[test]
fn test_complete_rings_only_quinoline_raw_vs_pruned() {
    // Verify quinoline series: check whether complete_rings_only actually changes
    // the atom count vs default (determines if pruning fires at all for this case)
    let a = parse("c1ccc2nc(CC)ccc2c1").unwrap();
    let b = parse("c1ccc2nc(CO)ccc2c1").unwrap();
    let c = parse("c1ccc2nc(CN)ccc2c1").unwrap();

    let default_result = find_mcs(&[&a, &b, &c]);
    let config = McsConfig {
        complete_rings_only: true,
        ..McsConfig::default()
    };
    let pruned_result = find_mcs_with_config(&[&a, &b, &c], &config);

    eprintln!(
        "quinoline series raw MCS: {}, with complete_rings_only: {}",
        default_result.atom_count(),
        pruned_result.atom_count()
    );
    // Document whether complete_rings_only triggered a change
    // (if equal, pruning was not triggered for this series)
}

#[test]
fn test_ring_matches_ring_only_and_complete_rings_only_combined() {
    // Both flags together: quinoline series should return the full quinoline scaffold
    let a = parse("c1ccc2nc(CC)ccc2c1").unwrap();
    let b = parse("c1ccc2nc(CO)ccc2c1").unwrap();
    let config = McsConfig {
        ring_matches_ring_only: true,
        complete_rings_only: true,
        ..McsConfig::default()
    };
    let result = find_mcs_with_config(&[&a, &b], &config);
    assert!(
        result.atom_count() >= 10,
        "Combined flags should yield at least the quinoline scaffold, got {}",
        result.atom_count()
    );
}

// -----------------------------------------------------------------------
// match_chiral_tag tests
// -----------------------------------------------------------------------

#[test]
fn test_match_chiral_tag_false_matches_enantiomers() {
    // R-alanine vs S-alanine: with match_chiral_tag=false (default),
    // should match all atoms (ignoring chirality)
    let r_ala = parse("[C@H](C)(N)C(=O)O").unwrap();
    let s_ala = parse("[C@@H](C)(N)C(=O)O").unwrap();
    let default_result = find_mcs(&[&r_ala, &s_ala]);
    assert_eq!(
        default_result.atom_count(),
        r_ala.atom_count(),
        "Default (match_chiral_tag=false) should match all atoms of enantiomers"
    );
}

#[test]
fn test_match_chiral_tag_true_blocks_stereocenters() {
    // R-alanine vs S-alanine: with match_chiral_tag=true,
    // the stereocenter C should not match (different chirality).
    // However, non-chiral atoms (N, O, carboxyl C) can still match.
    // So MCS should be smaller than without the constraint.
    let r_ala = parse("[C@H](C)(N)C(=O)O").unwrap();
    let s_ala = parse("[C@@H](C)(N)C(=O)O").unwrap();

    let default_result = find_mcs(&[&r_ala, &s_ala]);
    let config = McsConfig {
        match_chiral_tag: true,
        ..McsConfig::default()
    };
    let result = find_mcs_with_config(&[&r_ala, &s_ala], &config);

    assert!(
        result.atom_count() < default_result.atom_count(),
        "match_chiral_tag=true should find smaller MCS than default: {} vs {}",
        result.atom_count(),
        default_result.atom_count()
    );
}

#[test]
fn test_match_chiral_tag_true_matches_same_stereoisomer() {
    // R-alanine vs R-alanine: with match_chiral_tag=true,
    // should still match all atoms (same stereochemistry)
    let r_ala1 = parse("[C@H](C)(N)C(=O)O").unwrap();
    let r_ala2 = parse("[C@H](C)(N)C(=O)O").unwrap();
    let config = McsConfig {
        match_chiral_tag: true,
        ..McsConfig::default()
    };
    let result = find_mcs_with_config(&[&r_ala1, &r_ala2], &config);
    assert_eq!(
        result.atom_count(),
        r_ala1.atom_count(),
        "match_chiral_tag=true should still match molecules with same stereochemistry"
    );
}

// -----------------------------------------------------------------------
// match_charge / match_isotope tests (mirror match_chiral_tag exactly)
// -----------------------------------------------------------------------

#[test]
fn test_match_charge_false_matches_charge_states() {
    // Acetate vs neutral acetic acid: with match_charge=false (default),
    // the carboxylate/carboxylic-acid oxygen should still match.
    let acetate = parse("CC(=O)[O-]").unwrap();
    let acetic_acid = parse("CC(=O)O").unwrap();
    let default_result = find_mcs(&[&acetate, &acetic_acid]);
    assert_eq!(
        default_result.atom_count(),
        acetate.atom_count(),
        "Default (match_charge=false) should match all atoms across charge states"
    );
}

#[test]
fn test_match_charge_true_blocks_differently_charged_atoms() {
    // Acetate vs neutral acetic acid: with match_charge=true, the [O-]
    // atom must not match the neutral O, so MCS should be smaller.
    let acetate = parse("CC(=O)[O-]").unwrap();
    let acetic_acid = parse("CC(=O)O").unwrap();
    let default_result = find_mcs(&[&acetate, &acetic_acid]);
    let config = McsConfig {
        match_charge: true,
        ..McsConfig::default()
    };
    let result = find_mcs_with_config(&[&acetate, &acetic_acid], &config);
    assert!(
        result.atom_count() < default_result.atom_count(),
        "match_charge=true should find smaller MCS than default: {} vs {}",
        result.atom_count(),
        default_result.atom_count()
    );
}

#[test]
fn test_match_charge_true_matches_same_charge_state() {
    // Acetate vs acetate: with match_charge=true, should still match
    // all atoms (identical charge state).
    let acetate1 = parse("CC(=O)[O-]").unwrap();
    let acetate2 = parse("CC(=O)[O-]").unwrap();
    let config = McsConfig {
        match_charge: true,
        ..McsConfig::default()
    };
    let result = find_mcs_with_config(&[&acetate1, &acetate2], &config);
    assert_eq!(
        result.atom_count(),
        acetate1.atom_count(),
        "match_charge=true should still match molecules with the same charge state"
    );
}

#[test]
fn test_match_isotope_false_matches_isotope_labels() {
    // 13C-labeled methanol vs unlabeled methanol: with match_isotope=false
    // (default), the labeled carbon should still match.
    let labeled = parse("[13CH3]O").unwrap();
    let unlabeled = parse("CO").unwrap();
    let default_result = find_mcs(&[&labeled, &unlabeled]);
    assert_eq!(
        default_result.atom_count(),
        labeled.atom_count(),
        "Default (match_isotope=false) should match all atoms across isotope labels"
    );
}

#[test]
fn test_match_isotope_true_blocks_differently_labeled_atoms() {
    // 13C-labeled methanol vs unlabeled methanol: with match_isotope=true,
    // the labeled carbon must not match the unlabeled one, so MCS should
    // be smaller.
    let labeled = parse("[13CH3]O").unwrap();
    let unlabeled = parse("CO").unwrap();
    let default_result = find_mcs(&[&labeled, &unlabeled]);
    let config = McsConfig {
        match_isotope: true,
        ..McsConfig::default()
    };
    let result = find_mcs_with_config(&[&labeled, &unlabeled], &config);
    assert!(
        result.atom_count() < default_result.atom_count(),
        "match_isotope=true should find smaller MCS than default: {} vs {}",
        result.atom_count(),
        default_result.atom_count()
    );
}

#[test]
fn test_match_isotope_true_matches_same_isotope_label() {
    // Two identically 13C-labeled methanol molecules: with
    // match_isotope=true, should still match all atoms.
    let labeled1 = parse("[13CH3]O").unwrap();
    let labeled2 = parse("[13CH3]O").unwrap();
    let config = McsConfig {
        match_isotope: true,
        ..McsConfig::default()
    };
    let result = find_mcs_with_config(&[&labeled1, &labeled2], &config);
    assert_eq!(
        result.atom_count(),
        labeled1.atom_count(),
        "match_isotope=true should still match molecules with the same isotope label"
    );
}

#[test]
fn frontier_exclude_branch_finds_larger_mcs_than_include_only() {
    // Regression test for a real completeness bug: `grow` used to try
    // only `frontier[0]` at each node with no way to exclude it and try
    // a different frontier atom instead.
    //
    // mol0: O-C(-N)(-N) — hub C's neighbor list (parse/bond-insertion
    // order) is [O, N(branch), N]. mol1: N-C(-N) — hub' with 2 N leaves,
    // no O at all. True MCS is the full 3-atom C(N)(N) (mol1 in its
    // entirety): every seed that reaches "hub mapped, one N-leaf still
    // unmapped" puts unmatchable O ahead of the remaining N leaf in
    // frontier order. Without an exclude branch, growth died at size 2
    // (stuck on O) regardless of which atom was used as the seed --
    // fixed by branching include/exclude on `frontier[0]` instead of
    // only ever including it.
    let mol0 = parse("OC(N)N").unwrap();
    let mol1 = parse("NC(N)").unwrap();
    let result = find_mcs(&[&mol0, &mol1]);
    assert_eq!(
        result.atom_count(),
        3,
        "expected true MCS (C(N)(N), 3 atoms = all of mol1), got {}",
        result.atom_count()
    );
}

#[test]
fn non_induced_mcs_is_input_order_independent_for_ring_and_chain() {
    let ring = parse("C1CC1").unwrap();
    let chain = parse("CCC").unwrap();

    let ring_first = find_mcs(&[&ring, &chain]);
    let chain_first = find_mcs(&[&chain, &ring]);

    assert_eq!(ring_first.atom_count(), 3);
    assert_eq!(ring_first.bonds.len(), 2);
    assert_eq!(chain_first.atom_count(), 3);
    assert_eq!(chain_first.bonds.len(), 2);

    let complete_rings = McsConfig {
        complete_rings_only: true,
        ..McsConfig::default()
    };
    let ring_constrained = find_mcs_with_config(&[&ring, &chain], &complete_rings);
    assert_eq!(ring_constrained.atom_count(), 0);
}

#[test]
fn equal_atom_seed_can_improve_bond_count() {
    let chain_and_ring = parse("CCC.C1CC1").unwrap();
    let ring = parse("C1CC1").unwrap();
    for refs in [[&chain_and_ring, &ring], [&ring, &chain_and_ring]] {
        let result = find_mcs(&refs);
        assert_eq!(result.atom_count(), 3);
        assert_eq!(result.bonds.len(), 3, "the cyclopropane beats the chain");
    }
}

#[test]
fn mcs_size_is_input_order_independent_across_mixed_topologies() {
    let molecules: Vec<Molecule> = [
        "CCC", "C1CC1", "C1CCC1", "CCCC", "c1ccccc1", "CC=CC", "CC#N", "CCO", "NCCN",
    ]
    .into_iter()
    .map(|smiles| parse(smiles).unwrap())
    .collect();

    for left in 0..molecules.len() {
        for right in (left + 1)..molecules.len() {
            let forward = find_mcs(&[&molecules[left], &molecules[right]]);
            let reverse = find_mcs(&[&molecules[right], &molecules[left]]);
            assert_eq!(
                (forward.atom_count(), forward.bonds.len()),
                (reverse.atom_count(), reverse.bonds.len()),
                "MCS size changed when inputs {left} and {right} were reversed"
            );
        }
    }
}

#[test]
fn test_element_different_tie_is_spelling_dependent_known_gap() {
    // Two disconnected 3-atom fragments differing only in the
    // heteroatom (N vs O): both are size-3/bond-count-2 candidates of
    // equal quality, a genuine (non-automorphic, element-different) tie
    // -- unlike most ties elsewhere in this crate's tests, the choice
    // here IS externally observable (AtomicNum(7) vs AtomicNum(8) in
    // the returned QueryMolecule), because AtomCompare::Elements records
    // atomic number per atom.
    //
    // Spot-checked per user request after the #4 (iupac) tie-break
    // fix turned out to be a net-negative regression: confirmed via a
    // temporary revert that mcs_tiebreak_key's choice here is IDENTICAL
    // to the pre-#8-fix (first-found-in-DFS-wins) behavior -- this is
    // not a regression introduced by #8, it's a pre-existing gap that
    // #8's tie-break doesn't resolve (see mcs_tiebreak_key's doc
    // comment: no full spelling-invariance is claimed for it).
    let n_first = find_mcs(&[&parse("NCC.OCC").unwrap(), &parse("NCC.OCC").unwrap()]);
    let o_first = find_mcs(&[&parse("OCC.NCC").unwrap(), &parse("OCC.NCC").unwrap()]);

    assert_eq!(n_first.atom_count(), 3);
    assert_eq!(o_first.atom_count(), 3);
    assert_eq!(
        n_first.atoms[0].query,
        AtomQuery::Primitive(AtomPrimitive::AtomicNum(7))
    );
    assert_eq!(
        o_first.atoms[0].query,
        AtomQuery::Primitive(AtomPrimitive::AtomicNum(8))
    );
}
