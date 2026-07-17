//! Kekule-S0: `chematic_core::apply_kekule` must not change a molecule's
//! represented stereochemistry.
//!
//! Root cause (found while validating Aromaticity-A1-1b-1's canonical
//! round-trip numbers): `apply_kekule` rebuilt its output via
//! `MoleculeBuilder` without copying `stereo_neighbor_order` (the
//! SMILES-text-order neighbor sequence `@`/`@@` is defined relative to).
//! `atom.chirality` alone survived the rebuild, but without the neighbor
//! order to interpret it against, downstream canonical-SMILES writing could
//! serialize the same tetrahedral center with the wrong `@`/`@@` — silently,
//! since `chirality` itself never became `None`.
//!
//! These tests check *chemical* invariance (same molecule, before vs after
//! an explicit `kekulize`+`apply_kekule` round trip produces the *same*
//! canonical SMILES), not merely that `@`/`@@` characters appear somewhere.
//! The direct-vs-kekulize-then-rearomatize comparisons need
//! `chematic_perception::apply_aromaticity` (Huckel re-perception) purely so
//! both sides land back in aromatic-form notation before being compared as
//! strings -- kekulize/apply_kekule themselves only ever change bond-order
//! *representation*, never aromaticity perception, so this is a fair,
//! notation-neutral comparison, not a different claim about a different
//! function.

use chematic_core::{apply_kekule, kekulize};
use chematic_perception::apply_aromaticity;
use chematic_smiles::{canonical_smiles, parse};

/// Canonicalize `smi` directly, and again after routing it through an
/// explicit `kekulize` -> `apply_kekule` -> re-aromatize round trip. Both
/// must produce the identical canonical string -- kekulization only changes
/// bond-order representation, never the molecule's stereochemistry, and
/// re-aromatizing returns both sides to the same (aromatic-form) notation
/// for a fair string comparison.
fn check_kekulize_invariant(smi: &str) -> Result<(), String> {
    let mol = parse(smi).map_err(|e| format!("PARSE_FAIL '{smi}': {e}"))?;
    let direct = canonical_smiles(&mol);

    let k = kekulize(&mol).map_err(|e| format!("KEKULIZE_FAIL '{smi}': {e}"))?;
    let kekulized = apply_kekule(&mol, &k);
    let reperceived = apply_aromaticity(&kekulized);
    let via_kekulize = canonical_smiles(&reperceived);

    if direct != via_kekulize {
        Err(format!(
            "KEKULE_STEREO_DRIFT '{smi}': direct='{direct}' via_kekulize='{via_kekulize}'"
        ))
    } else {
        Ok(())
    }
}

#[track_caller]
fn assert_all_kekulize_invariant(cases: &[&str]) {
    let failures: Vec<_> = cases
        .iter()
        .filter_map(|s| check_kekulize_invariant(s).err())
        .collect();
    assert!(
        failures.is_empty(),
        "kekulize stereo invariance failures:\n{}",
        failures.join("\n")
    );
}

#[test]
fn simple_chiral_centers_survive_kekulization() {
    assert_all_kekulize_invariant(&[
        "N[C@@H](C)C(=O)O",   // L-alanine
        "N[C@H](C)C(=O)O",    // D-alanine
        "[C@@H]1(N)CCCC1",    // aminocyclopentane
        "F[C@H]1CCCCC1",      // fluorocyclohexane
        "C[C@H](O)c1ccccc1",  // 1-phenylethanol
        "C[C@@H](O)c1ccccc1", // its enantiomer
    ]);
}

#[test]
fn chiral_substituent_on_aromatic_ring_survives_kekulization() {
    // A stereocenter directly exocyclic to an aromatic ring -- exactly the
    // shape that exposed the bug: the ring needs kekulization, the
    // substituent's stereo descriptor is defined relative to
    // stereo_neighbor_order recorded when the substituent was parsed.
    assert_all_kekulize_invariant(&[
        "c1ccc(cc1)[C@H](F)Cl",
        "c1ccc(cc1)[C@@H](F)Cl",
        "Cc1cn([C@H]2CCCC[C@@H]2C)c(=O)n1C",
    ]);
}

#[test]
fn real_corpus_case_survives_kekulization() {
    // One of the 25 "experimental_only" A1-1b-1 canonical-round-trip
    // instability cases -- three stereocenters, an aromatic ring directly
    // bonded to two of them. Confirmed via direct inspection (before this
    // fix) that stereo_neighbor_order was `None` for all three centers
    // immediately after `apply_kekule`, and canonical SMILES flipped
    // `@`/`@@` between the direct and kekulize-first paths.
    assert_all_kekulize_invariant(&[
        "CCCCc1cn([C@H]2[C@H](C)CCC[C@@H]2C)c(=O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)nc1",
    ]);
}

#[test]
fn enantiomer_pairs_remain_mirror_images_after_kekulization() {
    // For a true enantiomer pair (every `@`/`@@` swapped, topology
    // identical), routing both through kekulize+apply_kekule and
    // canonicalizing must still produce a string pair that is the same
    // mirror-image relationship: identical modulo `@`<->`@@` swaps, not
    // collapsed to the same string or diverging into an unrelated one. Both
    // sides stay in Kekule-form notation here (no re-aromatization needed)
    // since the comparison is symmetric -- kekule-form-a vs kekule-form-b.
    //
    // NOTE: the real_corpus_case_survives_kekulization molecule (3
    // stereocenters on a ring whose substitution pattern is constitutionally
    // symmetric -- A-Me-CH2-CH2-CH2-Me-A reads the same both directions) is
    // deliberately NOT used as a pair here: flipping all 3 centers coincides
    // with that ring's own topological automorphism, so naively-constructed
    // "flip every @" does not produce a distinguishable enantiomer for it
    // (confirmed: even the *direct*, no-kekulize-at-all canonicalization
    // collapses both spellings to the same string -- a canonicalizer/graph-
    // symmetry property unrelated to apply_kekule). Use asymmetric
    // stereocenters here instead, where "flip every @" is unambiguous.
    let pairs: &[(&str, &str)] = &[
        ("C[C@H](O)c1ccccc1", "C[C@@H](O)c1ccccc1"),
        ("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"),
        ("c1ccc(cc1)[C@H](F)Cl", "c1ccc(cc1)[C@@H](F)Cl"),
    ];

    fn via_kekulize(smi: &str) -> String {
        let mol = parse(smi).expect("valid SMILES");
        let k = kekulize(&mol).expect("kekulizable");
        canonical_smiles(&apply_kekule(&mol, &k))
    }

    /// Swap every `@`/`@@` in a canonical SMILES string. `@@` must be
    /// handled before lone `@` or the second `@` of `@@` gets treated as a
    /// separate lone `@`.
    fn swap_at_symbols(s: &str) -> String {
        s.replace("@@", "\u{0}")
            .replace('@', "@@")
            .replace('\u{0}', "@")
    }

    for (a, b) in pairs {
        let ca = via_kekulize(a);
        let cb = via_kekulize(b);
        assert_eq!(
            swap_at_symbols(&ca),
            cb,
            "enantiomer pair diverged from a clean mirror relationship: \
             {a} -> {ca}, {b} -> {cb}"
        );
    }
}
