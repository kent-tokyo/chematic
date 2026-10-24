use super::*;
use crate::parser::parse_smarts;
use chematic_smiles::parse;

#[test]
fn organometallic_ring_closure_is_single_in_raw_graph() {
    let mol = parse("CN(C)C[C-]12C3=C4C5=C1[Fe++]23456789[C-]%10C6=C7C8=C9%10").unwrap();
    let (bond_idx, bond) = mol.bond_between(AtomIdx(4), AtomIdx(9)).unwrap();
    assert_eq!(bond.order, BondOrder::Single);
    assert_eq!(mol.atom(AtomIdx(4)).charge, -1);
    assert_eq!(mol.atom(AtomIdx(9)).element.atomic_number(), 26);
    assert_eq!(
        mol.neighbors(AtomIdx(4))
            .map(|(_, idx)| mol.bond(idx).order.order_int() as usize)
            .sum::<usize>(),
        4
    );
    assert_eq!(mol.bond(bond_idx).order, BondOrder::Single);

    let view = rdkit_parity_iron_carbanion_view(&mol).unwrap();
    assert_eq!(view.bond(bond_idx).order, BondOrder::Dative);
    assert_eq!(mol.bond(bond_idx).order, BondOrder::Single);
    let (ordinary_fe_bond, _) = mol.bond_between(AtomIdx(9), AtomIdx(10)).unwrap();
    assert_eq!(view.bond(ordinary_fe_bond).order, BondOrder::Single);
    let rings = chematic_perception::find_sssr(&view);
    assert_eq!(
        rings
            .rings()
            .iter()
            .filter(|ring| ring.contains(&AtomIdx(4)))
            .count(),
        1
    );

    for (pattern, should_contain) in [
        ("[R1]", true),
        ("[R2]", false),
        ("[x2]", true),
        ("[x3]", false),
    ] {
        let query = parse_smarts(pattern).unwrap();
        let (matches, exhausted) =
            find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
        assert!(!exhausted);
        assert_eq!(
            matches
                .iter()
                .any(|mapping| mapping.values().any(|idx| *idx == AtomIdx(4))),
            should_contain,
            "{pattern}"
        );
    }

    for (pattern, should_contain) in [("*@*", false), ("*!@*", true)] {
        let query = parse_smarts(pattern).unwrap();
        let (matches, exhausted) =
            find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
        assert!(!exhausted);
        assert_eq!(
            atom_sets(&matches).contains(&vec![4, 9]),
            should_contain,
            "{pattern}"
        );
    }
}

fn atom_sets(matches: &[FxHashMap<usize, AtomIdx>]) -> Vec<Vec<u32>> {
    let mut sets: Vec<Vec<u32>> = matches
        .iter()
        .map(|m| {
            let mut v: Vec<u32> = m.values().map(|a| a.0).collect();
            v.sort_unstable();
            v
        })
        .collect();
    sets.sort();
    sets
}

#[test]
fn matches_default_on_ordinary_query() {
    // A query with no [RN] at all: opt-in mode must produce the exact
    // same match set as the default matcher (no ring model even built).
    let mol = parse("c1ccc2ccccc2c1").unwrap(); // naphthalene
    let query = parse_smarts("c1ccccc1").unwrap();
    let default = crate::find_matches(&query, &mol);
    let (parity, exhausted) =
        find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
    assert!(!exhausted);
    assert_eq!(atom_sets(&default), atom_sets(&parity));
}

#[test]
fn r2_matches_extra_ring_bridgeheads_on_adamantane() {
    // Adamantane: RDKit ground truth (rdkit==2026.03.3, live oracle,
    // atom indices 0-9 in the same left-to-right SMILES parse order as
    // chematic's) gives exactly 4 atoms (indices 1,3,5,7) a ring count
    // of 3, the rest (0,2,4,6,8,9) a ring count of 2.
    //
    // chematic's raw SSSR (cycle_rank 3, no augmentation) happens to
    // pick a *different but equally valid* basis where atom 1 alone
    // already sits in all 3 basis rings -- a concrete illustration of
    // the "genuine SSSR-basis-cardinality disagreement" `docs/
    // rdkit_compat.md`'s SMARTS-R0/R2 sections describe: different
    // (both valid) basis choices produce different naive [R3] answers
    // even before any symmetrization. See `dbg_print_sssr` below for
    // the raw per-atom counts this depends on.
    let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap();
    let query = parse_smarts("[R3]").unwrap();
    let (parity, _) =
        find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
    let mut parity_atoms: Vec<u32> = parity.iter().map(|m| m[&0].0).collect();
    parity_atoms.sort_unstable();
    assert_eq!(parity_atoms, vec![1, 3, 5, 7], "RDKit-parity [R3] atom set");

    // Default matcher (plain SSSR) only finds atom 1 -- the one atom
    // chematic's own basis choice happens to already give count 3 to,
    // without needing the 4th "extra" ring at all. This is exactly the
    // partial/coincidental overlap the mode difference is meant to fix
    // (the other 3 -- atoms 3, 5, 7 -- are missed on plain SSSR).
    let default = crate::find_matches(&query, &mol);
    let mut default_atoms: Vec<u32> = default.iter().map(|m| m[&0].0).collect();
    default_atoms.sort_unstable();
    assert_eq!(default_atoms, vec![1], "default (plain-SSSR) [R3] atom set");
}

#[test]
fn k6_uses_symmetry_equivalent_ring_in_opt_in_mode() {
    // In this bridged fragment, atom 46 belongs to the six-membered
    // ring selected by RDKit 2026.03.6, but not to chematic's plain
    // SSSR basis. The default matcher must retain its existing answer.
    let mol = parse("C=CC[C@H](NC(=O)[C@@H]1C[C@@H](CCCc2cccc3ccccc23)c2c(Cl)nc(NCc3cccc(OC)c3)c(=O)n21)B1OC2CC3CC(C3(C)C)[C@@]2(C)O1").unwrap();
    let query = parse_smarts("[k6]").unwrap();
    let default = crate::find_matches(&query, &mol);
    let (parity, exhausted) =
        find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
    assert!(!exhausted);
    assert!(!default.iter().any(|m| m[&0] == AtomIdx(46)));
    assert!(parity.iter().any(|m| m[&0] == AtomIdx(46)));
}

#[test]
fn k6_corpus_3498_uses_mixed_size_replacement() {
    // RDKit 2026.03.6 includes atom 22 in the replacement six-ring
    // 1-2-20-21-22-23. Perception's symmetrized set keeps only the
    // four-ring through atom 22, while the SMARTS-specific bounded
    // selector independently recovers the six-ring.
    let mol = parse("C=C1[C@H]2Oc3cc(C(C)(C)CCCCCC)cc(O)c3[C@H]2[C@H]2C[C@@H]1C2(C)C").unwrap();
    let query = parse_smarts("[k6]").unwrap();
    let (parity, _) =
        find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
    let mut atoms: Vec<u32> = parity.iter().map(|m| m[&0].0).collect();
    atoms.sort_unstable();
    assert_eq!(atoms, vec![1, 2, 4, 5, 6, 16, 17, 19, 20, 21, 22, 23, 24]);
}

#[test]
fn k6_ring_model_budget_exceeded_is_typed() {
    let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap();
    let query = parse_smarts("[k6]").unwrap();
    let config = RdkitParityConfig {
        ring_model_budget: RdkitRingModelBudget { max_candidates: 0 },
        ..RdkitParityConfig::default()
    };
    assert!(matches!(
        find_matches_rdkit_parity(&query, &mol, &config),
        Err(RdkitParityError::RingModelBudgetExceeded { cap: 0, .. })
    ));
}

#[test]
fn ring_model_budget_exceeded_is_typed_not_silent() {
    let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap(); // adamantane
    let query = parse_smarts("[R3]").unwrap();
    let config = RdkitParityConfig {
        ring_model_budget: RdkitRingModelBudget { max_candidates: 0 },
        ..RdkitParityConfig::default()
    };
    let result = find_matches_rdkit_parity(&query, &mol, &config);
    assert!(matches!(
        result,
        Err(RdkitParityError::RingModelBudgetExceeded {
            candidates_examined: 1,
            cap: 0,
        })
    ));
}

#[test]
fn shared_ring_model_propagates_measured_budget_exhaustion() {
    let mol = parse("C12C3C4C1C5C4C3C25").unwrap();
    let query = parse_smarts("[R3]").unwrap();
    let config = RdkitParityConfig {
        use_shared_symmetrized_sssr: true,
        ring_model_budget: RdkitRingModelBudget { max_candidates: 0 },
        ..RdkitParityConfig::default()
    };
    let result = find_matches_rdkit_parity(&query, &mol, &config);
    assert!(matches!(
        result,
        Err(RdkitParityError::RingModelBudgetExceeded {
            candidates_examined: 1,
            cap: 0,
        })
    ));
}

#[test]
fn order_dependent_ring_model_uses_exact_rdkit_ring_order() {
    // This highly charged fused polyaromatic scaffold was the sole
    // SMARTS parity regression in the 5,021-molecule A4 census: the
    // generated replacement rings produced many false-positive [R3]
    // matches. RDKit's symmetrized SSSR depends on atom order here; the
    // exact port of its ring order gives RDKit 2026.03's answer (atoms 3
    // and 4) for the molecule in its given order.
    let mol = parse("c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)-c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1")
        .unwrap();
    let query = parse_smarts("[R3]").unwrap();
    let (matches, _) =
        find_matches_rdkit_parity(&query, &mol, &RdkitParityConfig::default()).unwrap();
    let mut atoms: Vec<u32> = matches.iter().map(|m| m[&0].0).collect();
    atoms.sort_unstable();
    assert_eq!(atoms, vec![3, 4]);
}

#[test]
fn non_ring_count_query_unaffected_by_zero_ring_budget() {
    // A query that never touches [RN] must not even attempt to build
    // the ring model, so a starved budget must not affect it.
    let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap(); // adamantane
    let query = parse_smarts("[R]").unwrap();
    let config = RdkitParityConfig {
        ring_model_budget: RdkitRingModelBudget { max_candidates: 0 },
        ..RdkitParityConfig::default()
    };
    let (parity, _) = find_matches_rdkit_parity(&query, &mol, &config).unwrap();
    assert_eq!(parity.len(), 10); // every atom in adamantane is in a ring
}

#[test]
fn has_match_rdkit_parity_bounded_three_way_outcome() {
    let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap();
    let query = parse_smarts("[R3]").unwrap();
    let found =
        has_match_rdkit_parity_bounded(&query, &mol, &RdkitParityConfig::default()).unwrap();
    assert_eq!(found, MatchOutcome::Found);

    let query_none = parse_smarts("[R5]").unwrap();
    let not_found =
        has_match_rdkit_parity_bounded(&query_none, &mol, &RdkitParityConfig::default()).unwrap();
    assert_eq!(not_found, MatchOutcome::NotFound);
}

// -- SMARTS-A0 bridgehead-N bucket: this pipeline vs. RdkitLike re-perception --
//
// This reproducer guards the bridgehead-N ring-fusion aromaticity case.
// RDKit ground truth is that only atoms 2-7, the benzo ring, are
// aromatic. Both the direct parser and explicit RdkitLike re-perception
// must preserve that result.

#[test]
fn smarts_a0_does_not_fire_on_direct_parse_no_reperception() {
    // Neither `find_matches` nor `find_matches_rdkit_parity` ever calls
    // any aromaticity re-perception on their own -- they match whatever
    // flags the input molecule already carries. The reproducer's mixed
    // aromatic/Kekule SMILES already encodes the *correct* (RDKit-
    // matching) flags directly from parsing, so SMARTS-A0's precondition
    // (a `RdkitLike` re-perception pass) is simply never reached on this
    // pipeline as used by this crate's default entry points.
    let mol = parse("C1=Cc2ccccc2C2=NCCCN12").unwrap();
    let c_query = parse_smarts("c").unwrap();
    let default = crate::find_matches(&c_query, &mol);
    let mut atoms: Vec<u32> = default.iter().map(|m| m[&0].0).collect();
    atoms.sort_unstable();
    // RDKit ground truth (rdkit==2026.03.3, live oracle): exactly atoms
    // 2-7 (the benzo ring), never atoms 0, 1, or 13.
    assert_eq!(atoms, vec![2, 3, 4, 5, 6, 7]);
}

#[test]
fn smarts_a0_remains_rdkit_compatible_after_explicit_reperception() {
    // Callers may explicitly re-perceive with
    // `AromaticityAlgorithm::RdkitLike` before handing the molecule to
    // this crate. The fixed aromaticity path must retain the RDKit
    // ground-truth benzo ring rather than extend it into the fused ring.
    let mol = parse("C1=Cc2ccccc2C2=NCCCN12").unwrap();
    let kekulized = chematic_core::kekulize(&mol)
        .map(|k| chematic_core::apply_kekule(&mol, &k))
        .unwrap_or(mol);
    let reperceived = chematic_perception::apply_aromaticity_ex(
        &kekulized,
        chematic_perception::AromaticityAlgorithm::RdkitLike,
    );
    let c_query = parse_smarts("c").unwrap();
    let default = crate::find_matches(&c_query, &reperceived);
    let mut atoms: Vec<u32> = default.iter().map(|m| m[&0].0).collect();
    atoms.sort_unstable();
    assert_eq!(atoms, vec![2, 3, 4, 5, 6, 7]);
}

#[test]
fn use_rdkit_parity_aromaticity_flag_avoids_the_rdkit_like_bug_on_this_reproducer() {
    // This crate's OTHER opt-in flag, `use_rdkit_parity_aromaticity`
    // (a from-scratch port of RDKit's actual aromaticity algorithm, see
    // `chematic_perception::apply_aromaticity_rdkit_parity_experimental`'s
    // own doc comment) is a materially different engine from the
    // heuristic `AromaticityAlgorithm::RdkitLike` used above. On this
    // specific bare-core reproducer it does NOT reproduce SMARTS-A0's
    // over-extension -- checked directly, not assumed, since the two
    // engines are unrelated code paths and this crate must not imply
    // one fixes the other without evidence.
    let mol = parse("C1=Cc2ccccc2C2=NCCCN12").unwrap();
    let c_query = parse_smarts("c").unwrap();
    let config = RdkitParityConfig {
        use_rdkit_parity_aromaticity: true,
        ..RdkitParityConfig::default()
    };
    let (parity, _) = find_matches_rdkit_parity(&c_query, &mol, &config).unwrap();
    let mut atoms: Vec<u32> = parity.iter().map(|m| m[&0].0).collect();
    atoms.sort_unstable();
    assert_eq!(
        atoms,
        vec![2, 3, 4, 5, 6, 7],
        "use_rdkit_parity_aromaticity=true should match RDKit's real answer \
             (benzo ring only) on this reproducer, not RdkitLike's over-extension"
    );
}
