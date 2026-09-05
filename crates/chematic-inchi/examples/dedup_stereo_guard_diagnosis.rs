//! One-off diagnostic (not part of the dedup public API) for the 12 corpus
//! molecules newly flagged `VerificationUnavailable` by
//! `has_unresolved_specified_tetrahedral_stereo`, beyond the originally
//! disclosed 4663/4664 pair. Answers, per specified-chiral atom: does the
//! legacy CIP engine (`chematic_chem::assign_cip`) resolve it, does the
//! accurate engine (`CipMode::Accurate`) resolve it, and does
//! `tetrahedral_stereo_neighbors` (what native InChI conversion actually
//! depends on) resolve it. Feeds the A/B/C/D/E classification in the PR
//! discussion; this file is diagnosis-only, it does not change any
//! production code path.
//!
//! Usage: `cargo run -p chematic-inchi --features native-inchi --example
//! dedup_stereo_guard_diagnosis`

#![cfg(feature = "native-inchi")]

use chematic_chem::{
    CipMode, CipUnresolvedReason, assign_cip, assign_cip_with_mode, tetrahedral_stereo_neighbors,
};
use chematic_core::Chirality;
use chematic_smiles::{canonical_smiles, parse, random_smiles};

const CASES: &[(usize, &str)] = &[
    (
        196,
        "CCCCc1cn([C@H]2[C@H](C)CCC[C@@H]2C)c(=O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)nc1",
    ),
    (
        443,
        "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@]23C[C@H]4C[C@H](C[C@H](C4)C2)C3)co1",
    ),
    (
        590,
        "O=c1cc(COC2CCOCC2)occ1OC(=O)[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2",
    ),
    (1567, "NS(=O)(=O)OC[C@@]12CCCC[C@@H]1CCC2"),
    (
        1609,
        "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
    ),
    (
        1643,
        "COc1ccccc1N1CCN(CCCCNC(=O)C23C[C@H]4C[C@@H](C2)C[C@@H](C3)C4)CC1",
    ),
    (
        4047,
        "O=C(Oc1c(O)cc(C(=O)O[C@H]2[C@H](OC(=O)c3cc(O)c(O)c(O)c3)C[C@](O)(C(=O)O)C[C@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
    ),
    (
        4178,
        "O=C(NCCCCN1CCCC(/C=C\\c2ccccc2)C1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
    ),
    (
        4412,
        "O=C(O[C@H]1[C@H](O)C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
    ),
    (
        4413,
        "O=C(O[C@H]1[C@H](O)C[C@](O)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
    ),
    (
        4509,
        "O=C(O[C@@H]1C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1OC(=O)c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1",
    ),
    (4745, "N[C@]1(C(=O)O)C[C@H](C(=O)O)[C@H](C(=O)O)C1"),
    // Reference: the originally disclosed pair.
    (
        4663,
        "O=C(Oc1c(O)cc(C(=O)O[C@@H]2C[C@](O)(C(=O)O)C[C@@H](OC(=O)c3cc(O)c(O)c(O)c3)[C@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
    ),
    (
        4664,
        "O=C(Oc1c(O)cc(C(=O)O[C@@H]2C[C@@](O)(C(=O)O)C[C@@H](OC(=O)c3cc(O)c(O)c(O)c3)[C@@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
    ),
];

fn main() {
    for &(idx, smiles) in CASES {
        let mol = parse(smiles).unwrap_or_else(|e| panic!("idx={idx} parse failed: {e}"));
        let legacy = assign_cip(&mol);
        let accurate = assign_cip_with_mode(&mol, CipMode::Accurate);
        let canon = canonical_smiles(&mol);

        println!("MOLECULE idx={idx} smiles={smiles:?} canonical={canon:?}");

        for (aidx, atom) in mol.atoms() {
            if !matches!(
                atom.chirality,
                Chirality::CounterClockwise | Chirality::Clockwise
            ) {
                continue;
            }
            let tsn = tetrahedral_stereo_neighbors(&mol, aidx);
            let legacy_code = legacy.get(aidx);
            let (accurate_code, accurate_unresolved) = match &accurate {
                Ok(a) => (
                    a.get(aidx),
                    a.unresolved
                        .iter()
                        .find(|(i, _)| *i == aidx)
                        .map(|(_, r)| *r),
                ),
                Err(e) => {
                    println!("  ACCURATE_ENGINE_ERROR idx={idx} atom={aidx:?} err={e:?}");
                    (None, None)
                }
            };
            let accurate_reason = accurate_unresolved.map(|r| match r {
                CipUnresolvedReason::Tied => "Tied",
                CipUnresolvedReason::BudgetExceeded => "BudgetExceeded",
                CipUnresolvedReason::OracleUnstable => "OracleUnstable",
            });

            println!(
                "  ATOM idx={idx} atom={aidx:?} chirality={:?} tetrahedral_stereo_neighbors_resolved={} legacy_cip={:?} accurate_cip={:?} accurate_unresolved_reason={:?}",
                atom.chirality,
                tsn.is_some(),
                legacy_code,
                accurate_code,
                accurate_reason,
            );
        }

        // Renumbering stability (8 respellings): confirm the SET of
        // (unresolved-by-tetrahedral_stereo_neighbors) PHYSICAL atoms is
        // renumbering-invariant, tracked via each respelling's own
        // canonical SMILES agreeing with the reference (same molecule) and
        // the *count* of chirality-tagged/unresolved atoms staying fixed.
        let mut unresolved_counts = Vec::new();
        let mut tagged_counts = Vec::new();
        let mut canonical_mismatches = Vec::new();
        for seed in 0..8u64 {
            let s = random_smiles(&mol, seed);
            let m2 = parse(&s).unwrap_or_else(|e| panic!("idx={idx} reparse failed: {e}"));
            let c2 = canonical_smiles(&m2);
            if c2 != canon {
                canonical_mismatches.push((seed, c2.clone()));
            }
            let tagged = m2
                .atoms()
                .filter(|(_, a)| {
                    matches!(
                        a.chirality,
                        Chirality::CounterClockwise | Chirality::Clockwise
                    )
                })
                .count();
            let unresolved = m2
                .atoms()
                .filter(|(aidx, a)| {
                    matches!(
                        a.chirality,
                        Chirality::CounterClockwise | Chirality::Clockwise
                    ) && tetrahedral_stereo_neighbors(&m2, *aidx).is_none()
                })
                .count();
            tagged_counts.push(tagged);
            unresolved_counts.push(unresolved);
        }
        println!(
            "  RENUMBERING_STABILITY idx={idx} tagged_counts={tagged_counts:?} unresolved_counts={unresolved_counts:?} canonical_mismatches={:?}",
            canonical_mismatches
                .iter()
                .map(|(s, _)| *s)
                .collect::<Vec<_>>()
        );
        for (seed, c2) in &canonical_mismatches {
            println!("    CANONICAL_MISMATCH idx={idx} seed={seed} got={c2:?}");
        }
    }
}
