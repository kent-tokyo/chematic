use super::*;
#[test]
fn group_candidates_partitions_by_canonical_smiles() {
    use chematic_smiles::parse;
    let mols = vec![
        parse("CCO").unwrap(),
        parse("OCC").unwrap(), // same molecule, different spelling
        parse("CCN").unwrap(),
    ];
    let groups = group_candidates(&mols);
    assert_eq!(groups.len(), 2, "{groups:?}");
    // CCO and OCC (indices 0,1) share a canonical SMILES; CCN (index 2) doesn't.
    let group_with_zero = groups.iter().find(|g| g.contains(&0)).unwrap();
    assert_eq!(group_with_zero, &vec![0, 1]);
}

#[cfg(not(feature = "native-inchi"))]
#[test]
fn compare_without_native_inchi_is_always_unavailable() {
    use chematic_smiles::parse;
    let a = parse("c1ccccc1").unwrap();
    let b = parse("c1ccccc1").unwrap();
    assert_eq!(
        compare_molecules(&a, &b, IdentityPolicy::StandardInchiString),
        DedupRelation::VerificationUnavailable
    );
}

#[cfg(not(feature = "native-inchi"))]
#[test]
fn deduplicate_verified_without_native_inchi_reports_all_unavailable() {
    use chematic_smiles::parse;
    let mols = vec![parse("CCO").unwrap(), parse("OCC").unwrap()];
    let report = deduplicate_verified(&mols, IdentityPolicy::StandardInchiString);
    assert!(report.groups.is_empty());
    assert_eq!(report.verification_unavailable, vec![0, 1]);
}
// ── Accurate-CIP preflight (issue #161) ────────────────────────────────

/// The original 4663/4664 diastereomer pair PR #156 fixed must NEVER
/// become `VerifiedDuplicate` under the accurate-CIP preflight either --
/// per `validation/results/dedup_stereo_guard_corpus_audit.jsonl`, the
/// accurate engine resolves both previously-unresolved centres (atoms 13
/// and 32) for both molecules, but to genuinely DIFFERENT codes (4663:
/// S,R; 4664: R,S) -- so the stereo supplement must differ and the pair
/// must resolve to something other than `VerifiedDuplicate`. This is the
/// single most important regression guard in this file: a design that
/// merely relaxed the guard (without folding the accurate codes into the
/// comparison key) would silently reopen this exact false merge.
#[cfg(feature = "native-inchi")]
#[test]
fn accurate_preflight_never_reopens_the_4663_4664_false_duplicate() {
    use chematic_smiles::parse;
    let a = parse(
        "O=C(Oc1c(O)cc(C(=O)O[C@@H]2C[C@](O)(C(=O)O)C[C@@H](OC(=O)c3cc(O)c(O)c(O)c3)[C@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
    )
    .unwrap();
    let b = parse(
        "O=C(Oc1c(O)cc(C(=O)O[C@@H]2C[C@@](O)(C(=O)O)C[C@@H](OC(=O)c3cc(O)c(O)c(O)c3)[C@@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
    )
    .unwrap();

    // The plain path must never merge this known diastereomer pair. It
    // may either fail closed (the original PR #156 behavior) or now
    // distinguish the pair directly if the legacy ranking path resolves
    // the centres; both outcomes preserve the safety invariant.
    assert_ne!(
        compare_molecules(&a, &b, IdentityPolicy::StandardInchiString),
        DedupRelation::VerifiedDuplicate,
        "the plain path must not reopen PR #156's false duplicate"
    );

    let with_preflight =
        compare_molecules_with_accurate_cip_preflight(&a, &b, IdentityPolicy::StandardInchiString);
    assert_ne!(
        with_preflight,
        DedupRelation::VerifiedDuplicate,
        "accurate-CIP preflight must never claim these two real diastereomers \
         are the same molecule: {with_preflight:?}"
    );
}

/// Case B from the corpus audit (tricyclic-cage bridgehead, corpus
/// idx=443/590): unresolved by *both* chematic CIP engines (and by
/// RDKit's own modern CIP labeler) -- must stay `VerificationUnavailable`
/// even with the preflight engaged, never silently promoted.
#[cfg(feature = "native-inchi")]
#[test]
fn accurate_preflight_still_fails_closed_on_genuine_ties() {
    use chematic_smiles::parse;
    let idx_443 =
        parse("CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@]23C[C@H]4C[C@H](C[C@H](C4)C2)C3)co1")
            .unwrap();
    assert_eq!(
        compare_molecules_with_accurate_cip_preflight(
            &idx_443,
            &idx_443,
            IdentityPolicy::StandardInchiString
        ),
        DedupRelation::VerificationUnavailable,
        "genuine tricyclic-cage tie (case B) must still fail closed under the preflight"
    );
}

/// One of the 3 "case A blocked by tied rank" molecules (corpus
/// idx=1609): a 3-fold-symmetric adamantane-cage substituent where all 3
/// legacy-unresolved centres share one `morgan_ranks` value. This directly
/// exercises `accurate_stereo_supplement` (not just the outcome via
/// `compare_molecules_with_accurate_cip_preflight`) to confirm the
/// *mechanism* is really the rank-collision guard firing --
/// `AmbiguousStereoRankCorrespondence` -- and not some other failure path
/// (e.g. the accurate engine tying or erroring on these atoms, which
/// would make the module doc comment's causal claim wrong).
#[cfg(feature = "native-inchi")]
#[test]
fn accurate_stereo_supplement_blocks_via_rank_collision_not_some_other_path() {
    use chematic_smiles::parse;
    let mol =
        parse("O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3").unwrap();
    assert_eq!(
        accurate_stereo_supplement(&mol),
        Err(IdentityDiagnostic::AmbiguousStereoRankCorrespondence),
        "idx=1609 must be blocked specifically by the rank-collision guard, \
         not by AccurateCipTied/AccurateCipBudgetExceeded/AccurateCipEngineError"
    );
}

/// Case A from the corpus audit (corpus idx=1567): legacy fails on the
/// ring-fusion atom, accurate resolves it to a plain 'R' agreeing with
/// RDKit -- the preflight must recover verified-comparison capability
/// here (comparing the molecule against itself must no longer be
/// `VerificationUnavailable`).
#[cfg(feature = "native-inchi")]
#[test]
fn accurate_preflight_recovers_case_a_self_comparison() {
    use chematic_smiles::parse;
    let mol = parse("NS(=O)(=O)OC[C@@]12CCCC[C@@H]1CCC2").unwrap();
    assert_eq!(
        compare_molecules(&mol, &mol, IdentityPolicy::StandardInchiString),
        DedupRelation::VerificationUnavailable,
        "sanity: plain path is unavailable for this legacy-unresolved centre"
    );
    assert_eq!(
        compare_molecules_with_accurate_cip_preflight(
            &mol,
            &mol,
            IdentityPolicy::StandardInchiString
        ),
        DedupRelation::VerifiedDuplicate,
        "accurate-CIP preflight should recover a molecule compared against itself"
    );
}

/// Every EXISTING public entry point (`compare`, `compare_molecules`,
/// `deduplicate_verified`) must be completely unaffected by the new
/// preflight machinery existing in this file -- verified here by
/// checking a molecule with no stereo ambiguity at all behaves
/// identically either way, and (above) that the guarded molecules'
/// PLAIN-path behavior is unchanged from PR #156.
#[test]
fn preflight_is_opt_in_existing_paths_unchanged() {
    use chematic_smiles::parse;
    let a = parse("CCO").unwrap();
    let b = parse("OCC").unwrap();
    assert_eq!(
        compare_molecules(&a, &b, IdentityPolicy::StandardInchiString),
        compare_molecules_with_accurate_cip_preflight(&a, &b, IdentityPolicy::StandardInchiString),
        "an ordinary molecule with no legacy-CIP-unresolved centre must behave \
         identically under both entry points"
    );
}
// ── Indexed graph relation (requires matching atom-index correspondence) ──

#[test]
fn indexed_graph_relation_same_molecule_is_identical() {
    use chematic_smiles::parse;
    let a = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
    let b = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
    let mode = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::RawGraphExact,
        atom_maps: AtomMapPolicy::Ignore,
    };
    let evidence = compare_indexed_graph_relation(&a, &b, mode);
    assert_eq!(evidence.graph_relation, Some(GraphRelation::Identical));
}

#[test]
fn indexed_graph_relation_different_formula_is_distinct() {
    use chematic_smiles::parse;
    let a = parse("CCO").unwrap();
    let b = parse("CCN").unwrap();
    let mode = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::RawGraphExact,
        atom_maps: AtomMapPolicy::Ignore,
    };
    let evidence = compare_indexed_graph_relation(&a, &b, mode);
    assert_eq!(evidence.graph_relation, Some(GraphRelation::Distinct));
}

/// Enantiomers: same atom order, same connectivity, opposite tetrahedral
/// parity -- must be `Distinct` under exact graph identity (unlike
/// `IdentityPolicy::StereoIgnored`, which deliberately collapses them).
#[test]
fn indexed_graph_relation_enantiomers_are_distinct() {
    use chematic_smiles::parse;
    let l = parse("N[C@@H](C)C(=O)O").unwrap();
    let d = parse("N[C@H](C)C(=O)O").unwrap();
    let mode = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::RawGraphExact,
        atom_maps: AtomMapPolicy::Ignore,
    };
    let evidence = compare_indexed_graph_relation(&l, &d, mode);
    assert_eq!(evidence.graph_relation, Some(GraphRelation::Distinct));
}

/// E/Z isomers: same atom order/connectivity, opposite double-bond
/// direction markers -- must be `Distinct`.
#[test]
fn indexed_graph_relation_ez_isomers_are_distinct() {
    use chematic_smiles::parse;
    let e = parse("C/C=C/C").unwrap();
    let z = parse("C/C=C\\C").unwrap();
    let mode = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::RawGraphExact,
        atom_maps: AtomMapPolicy::Ignore,
    };
    let evidence = compare_indexed_graph_relation(&e, &z, mode);
    assert_eq!(evidence.graph_relation, Some(GraphRelation::Distinct));
}

/// `RawGraphExact` treats aromatic-order-spelled benzene and an explicit
/// Kekulé respelling as literally different bond orders.
#[test]
fn indexed_graph_relation_raw_graph_exact_distinguishes_aromatic_spelling() {
    use chematic_smiles::parse;
    let aromatic = parse("c1ccccc1").unwrap();
    let kekule = parse("C1=CC=CC=C1").unwrap();
    let mode = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::RawGraphExact,
        atom_maps: AtomMapPolicy::Ignore,
    };
    let evidence = compare_indexed_graph_relation(&aromatic, &kekule, mode);
    assert_ne!(
        evidence.graph_relation,
        Some(GraphRelation::Identical),
        "RawGraphExact must not equate aromatic-order and Kekule-order spelling"
    );
}

/// `AtomMapPolicy::Include` must distinguish molecules that otherwise
/// have identical structure but different atom-mapping numbers;
/// `AtomMapPolicy::Ignore` must not.
#[test]
fn indexed_graph_relation_atom_map_policy_is_respected() {
    use chematic_smiles::parse;
    let mapped_a = parse("[CH4:1]").unwrap();
    let mapped_b = parse("[CH4:2]").unwrap();

    let include = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::RawGraphExact,
        atom_maps: AtomMapPolicy::Include,
    };
    let ignore = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::RawGraphExact,
        atom_maps: AtomMapPolicy::Ignore,
    };

    assert_eq!(
        compare_indexed_graph_relation(&mapped_a, &mapped_b, include).graph_relation,
        Some(GraphRelation::Distinct)
    );
    assert_eq!(
        compare_indexed_graph_relation(&mapped_a, &mapped_b, ignore).graph_relation,
        Some(GraphRelation::Identical)
    );
}

/// When atom ordering genuinely doesn't line up between the two inputs
/// (a real respelling that starts the traversal from a different atom),
/// the bounded index-correspondence comparator must report `None`
/// (inconclusive) rather than a guess in either direction.
#[test]
fn indexed_graph_relation_reordered_atoms_is_inconclusive_not_a_guess() {
    use chematic_smiles::{parse, random_smiles};
    let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
    let mode = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::RawGraphExact,
        atom_maps: AtomMapPolicy::Ignore,
    };
    // Find a random respelling seed that actually reorders atoms (not
    // every seed necessarily does, for a small/symmetric molecule) --
    // try a handful of seeds and require at least one exercises the
    // inconclusive path.
    let mut saw_inconclusive = false;
    for seed in 0..8u64 {
        let respelled = parse(&random_smiles(&mol, seed)).unwrap();
        let evidence = compare_indexed_graph_relation(&mol, &respelled, mode);
        if evidence.graph_relation.is_none() {
            assert!(
                evidence.diagnostics.contains(
                    &IndexedGraphRelationDiagnostic::AtomOrderCorrespondenceNotEstablished
                ),
                "{:?}",
                evidence.diagnostics
            );
            saw_inconclusive = true;
        } else {
            // A seed that happens to preserve index order entirely is
            // fine too (still not a WRONG answer) -- just not the case
            // this test is trying to exercise.
            assert_eq!(evidence.graph_relation, Some(GraphRelation::Identical));
        }
    }
    assert!(
        saw_inconclusive,
        "expected at least one random-respelling seed to reorder atoms for this molecule"
    );
}

/// Never a false `Identical` and never a false `Distinct`: for a batch of
/// known-distinct small molecules compared pairwise under
/// `ChemicalGraphExact`, nothing spuriously collapses.
#[test]
fn indexed_graph_relation_chemical_graph_exact_no_false_positives() {
    use chematic_smiles::parse;
    let mols = ["CCO", "CCN", "CCC", "c1ccccc1", "CC(=O)O"]
        .iter()
        .map(|s| parse(s).unwrap())
        .collect::<Vec<_>>();
    let mode = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::ChemicalGraphExact,
        atom_maps: AtomMapPolicy::Ignore,
    };
    for i in 0..mols.len() {
        for j in (i + 1)..mols.len() {
            let evidence = compare_indexed_graph_relation(&mols[i], &mols[j], mode);
            assert_ne!(
                evidence.graph_relation,
                Some(GraphRelation::Identical),
                "{i} vs {j} must not be Identical"
            );
        }
    }
}
