//! End-to-end fixtures for `chematic_inchi::dedup`'s indexed graph relation API
//! (`IndexedGraphRelationMode`/`IndexedGraphRelationEvidence`/`compare_indexed_graph_relation`).
//! Named "indexed graph relation," not "graph identity" -- see the module
//! doc comment in `src/dedup.rs` for why: the comparator requires the two
//! input molecules to already share the same atom indexing, and is not a
//! general graph-isomorphism test.
//!
//! Unlike `IdentityPolicy`'s four native-InChI-backed policies, this API's
//! structural comparator (`same_structure_under_identity` in `src/dedup.rs`)
//! does not require the `native-inchi` feature at all -- it never uses
//! canonical SMILES or native InChI as the actual proof of identity (only as
//! optional corroborating evidence in `IndexedGraphRelationEvidence::standard_inchi`/
//! `inchikey`, which are simply `None` without the feature). So, unlike
//! `tests/dedup_fixtures.rs`, this file runs in both feature configurations.

use chematic_inchi::dedup::{
    AtomMapPolicy, GraphRelation, GraphStrictness, IndexedGraphRelationDiagnostic,
    IndexedGraphRelationMode, compare_indexed_graph_relation,
};
use chematic_smiles::{parse, random_smiles};

fn mol(smiles: &str) -> chematic_core::Molecule {
    parse(smiles).unwrap_or_else(|e| panic!("parse {smiles:?}: {e}"))
}

const RAW_IGNORE: IndexedGraphRelationMode = IndexedGraphRelationMode {
    graph_strictness: GraphStrictness::RawGraphExact,
    atom_maps: AtomMapPolicy::Ignore,
};
const CHEMICAL_IGNORE: IndexedGraphRelationMode = IndexedGraphRelationMode {
    graph_strictness: GraphStrictness::ChemicalGraphExact,
    atom_maps: AtomMapPolicy::Ignore,
};

#[test]
fn same_molecule_twice_is_identical_under_both_strictness_modes() {
    for mode in [RAW_IGNORE, CHEMICAL_IGNORE] {
        let a = mol("CC(=O)Oc1ccccc1C(=O)O");
        let b = mol("CC(=O)Oc1ccccc1C(=O)O");
        let evidence = compare_indexed_graph_relation(&a, &b, mode);
        assert_eq!(
            evidence.graph_relation,
            Some(GraphRelation::Identical),
            "{mode:?}"
        );
        // `diagnostics` may carry a `NativeInchiUnavailable` note when this
        // crate is built without the `native-inchi` feature (the InChI
        // corroboration fields are simply unavailable, not an error in the
        // structural comparator itself) -- not asserted empty here.
    }
}

#[test]
fn different_element_is_distinct() {
    let a = mol("CCO");
    let b = mol("CCN");
    let evidence = compare_indexed_graph_relation(&a, &b, RAW_IGNORE);
    assert_eq!(evidence.graph_relation, Some(GraphRelation::Distinct));
}

#[test]
fn different_charge_is_distinct() {
    let a = mol("[NH4+]");
    let b = mol("[NH3]");
    let evidence = compare_indexed_graph_relation(&a, &b, RAW_IGNORE);
    assert_eq!(evidence.graph_relation, Some(GraphRelation::Distinct));
}

#[test]
fn different_isotope_is_distinct() {
    let a = mol("[13CH4]");
    let b = mol("C");
    let evidence = compare_indexed_graph_relation(&a, &b, RAW_IGNORE);
    assert_eq!(evidence.graph_relation, Some(GraphRelation::Distinct));
}

#[test]
fn different_fragment_composition_is_distinct() {
    // Same total heavy-atom multiset, different fragmentation (one molecule
    // vs. two salt fragments) -- must not collide.
    let a = mol("CCCC"); // butane, 1 fragment
    let b = mol("CC.CC"); // two ethane fragments
    let evidence = compare_indexed_graph_relation(&a, &b, RAW_IGNORE);
    assert_eq!(evidence.graph_relation, Some(GraphRelation::Distinct));
}

#[test]
fn enantiomers_are_distinct_not_collapsed() {
    let l = mol("N[C@@H](C)C(=O)O");
    let d = mol("N[C@H](C)C(=O)O");
    for mode in [RAW_IGNORE, CHEMICAL_IGNORE] {
        let evidence = compare_indexed_graph_relation(&l, &d, mode);
        assert_eq!(
            evidence.graph_relation,
            Some(GraphRelation::Distinct),
            "{mode:?}"
        );
    }
}

#[test]
fn diastereomers_are_distinct() {
    // Tartaric acid: (2R,3R) vs meso (2R,3S).
    let rr = mol("OC(=O)[C@H](O)[C@H](O)C(=O)O");
    let meso = mol("OC(=O)[C@H](O)[C@@H](O)C(=O)O");
    let evidence = compare_indexed_graph_relation(&rr, &meso, RAW_IGNORE);
    assert_eq!(evidence.graph_relation, Some(GraphRelation::Distinct));
}

#[test]
fn ez_isomers_are_distinct() {
    let e = mol("C/C=C/C");
    let z = mol("C/C=C\\C");
    let evidence = compare_indexed_graph_relation(&e, &z, RAW_IGNORE);
    assert_eq!(evidence.graph_relation, Some(GraphRelation::Distinct));
}

#[test]
fn ez_same_configuration_written_differently_is_identical() {
    // (E)-but-2-ene, two equivalent ways to write the same E configuration
    // (marker on the other side of the double bond).
    let a = mol("C/C=C/C");
    let b = mol("C/C=C/C");
    let evidence = compare_indexed_graph_relation(&a, &b, RAW_IGNORE);
    assert_eq!(evidence.graph_relation, Some(GraphRelation::Identical));
}

#[test]
fn raw_graph_exact_distinguishes_aromatic_vs_kekule_spelling() {
    let aromatic = mol("c1ccccc1");
    let kekule = mol("C1=CC=CC=C1");
    let evidence = compare_indexed_graph_relation(&aromatic, &kekule, RAW_IGNORE);
    assert_ne!(
        evidence.graph_relation,
        Some(GraphRelation::Identical),
        "RawGraphExact must not equate aromatic-order and Kekule-order spelling: {evidence:?}"
    );
}

#[test]
fn chemical_graph_exact_equates_aromatic_vs_same_kekule_spelling() {
    // `kekulize()` is deterministic given the same underlying graph -- so an
    // aromatic-flagged benzene and an explicit Kekule respelling that
    // happens to already match `kekulize()`'s own chosen alternation must
    // compare Identical under ChemicalGraphExact.
    let aromatic = mol("c1ccccc1");
    let kekule = mol("C1=CC=CC=C1");
    let evidence = compare_indexed_graph_relation(&aromatic, &kekule, CHEMICAL_IGNORE);
    assert_eq!(
        evidence.graph_relation,
        Some(GraphRelation::Identical),
        "{evidence:?}"
    );
}

#[test]
fn chemical_graph_exact_never_false_positive_across_distinct_small_molecules() {
    let smiles = ["CCO", "CCN", "CCC", "c1ccccc1", "CC(=O)O", "c1ccncc1"];
    let mols: Vec<_> = smiles.iter().map(|s| mol(s)).collect();
    for i in 0..mols.len() {
        for j in (i + 1)..mols.len() {
            let evidence = compare_indexed_graph_relation(&mols[i], &mols[j], CHEMICAL_IGNORE);
            assert_ne!(
                evidence.graph_relation,
                Some(GraphRelation::Identical),
                "{} vs {}",
                smiles[i],
                smiles[j],
            );
        }
    }
}

#[test]
fn atom_map_include_distinguishes_include_ignores_dont() {
    let a = mol("[CH4:1]");
    let b = mol("[CH4:2]");
    let include = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::RawGraphExact,
        atom_maps: AtomMapPolicy::Include,
    };
    assert_eq!(
        compare_indexed_graph_relation(&a, &b, include).graph_relation,
        Some(GraphRelation::Distinct)
    );
    assert_eq!(
        compare_indexed_graph_relation(&a, &b, RAW_IGNORE).graph_relation,
        Some(GraphRelation::Identical)
    );
}

#[test]
fn atom_map_axis_is_orthogonal_to_strictness_axis() {
    // Same molecule, different atom maps, aromatic vs Kekule spelling: all
    // 4 combinations of the two axes must be independently selectable and
    // give the expected, axis-specific answer.
    let a = mol("[cH:1]1ccccc1");
    let b_same_map_kekule = mol("[CH:1]1=CC=CC=C1");

    let raw_include = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::RawGraphExact,
        atom_maps: AtomMapPolicy::Include,
    };
    let chemical_include = IndexedGraphRelationMode {
        graph_strictness: GraphStrictness::ChemicalGraphExact,
        atom_maps: AtomMapPolicy::Include,
    };

    // RawGraphExact still distinguishes the aromatic/Kekule spelling
    // regardless of atom-map policy.
    assert_ne!(
        compare_indexed_graph_relation(&a, &b_same_map_kekule, raw_include).graph_relation,
        Some(GraphRelation::Identical)
    );
    // ChemicalGraphExact + Include atom maps: same map number, same
    // chemistry once Kekulized -- Identical.
    assert_eq!(
        compare_indexed_graph_relation(&a, &b_same_map_kekule, chemical_include).graph_relation,
        Some(GraphRelation::Identical)
    );
}

#[test]
fn reordered_atoms_are_inconclusive_never_a_wrong_guess() {
    let base = mol("CC(=O)Oc1ccccc1C(=O)O");
    let mut saw_inconclusive = false;
    for seed in 0..8u64 {
        let respelled = mol(&random_smiles(&base, seed));
        let evidence = compare_indexed_graph_relation(&base, &respelled, RAW_IGNORE);
        match evidence.graph_relation {
            None => {
                assert!(
                    evidence.diagnostics.contains(
                        &IndexedGraphRelationDiagnostic::AtomOrderCorrespondenceNotEstablished
                    ),
                    "{:?}",
                    evidence.diagnostics
                );
                saw_inconclusive = true;
            }
            Some(GraphRelation::Identical) => {} // a seed that happened to preserve index order -- still not wrong
            Some(GraphRelation::Distinct) => {
                panic!(
                    "a respelling of the SAME molecule must never be reported Distinct: seed={seed}"
                )
            }
        }
    }
    assert!(
        saw_inconclusive,
        "expected at least one seed to actually reorder atoms for this molecule"
    );
}

#[test]
fn evidence_never_claims_identical_and_distinct_diagnostics_simultaneously() {
    // Sanity/shape check: whenever `graph_relation` is `Some`, it is
    // decisive -- no lingering "why inconclusive" diagnostic should be
    // attached to a conclusive result (those diagnostics are reserved for
    // `graph_relation: None`).
    let a = mol("CCO");
    let b = mol("CCN");
    let evidence = compare_indexed_graph_relation(&a, &b, RAW_IGNORE);
    assert_eq!(evidence.graph_relation, Some(GraphRelation::Distinct));
    assert!(
        !evidence
            .diagnostics
            .contains(&IndexedGraphRelationDiagnostic::AtomOrderCorrespondenceNotEstablished)
    );
}
