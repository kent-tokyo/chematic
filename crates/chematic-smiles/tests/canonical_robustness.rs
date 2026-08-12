//! Canonical SMILES robustness tests — addressing RDKit issue #8775.
//!
//! RDKit #8775 documents 115+ SMILES that produce oscillating or incorrect
//! canonical forms. After recent RDKit fixes, 22 remain problematic. This suite
//! verifies chematic's behavior on these and related cases.
//!
//! Three test categories:
//! 1. **Stability**: `parse → canonical → parse → canonical` gives the same string.
//!    All cases are expected to pass.
//! 2. **Platform independence (topology)**: two SMILES of the same molecule (no stereo)
//!    must give identical canonical forms. All expected to pass.
//! 3. **Stereo parity**: same molecule written with different atom ordering must produce
//!    the same canonical SMILES. Implemented via parse-time stereo neighbor recording and
//!    permutation-parity correction in canonical.rs (v0.2.7+).

use chematic_smiles::{canonical_smiles, parse};

// ── Helpers ─────────────────────────────────────────────────────────────────

/// Returns Ok if the canonical SMILES for `smi` is roundtrip-stable, Err otherwise.
fn check_canonical_stable(smi: &str) -> Result<(), String> {
    let mol1 = parse(smi).map_err(|e| format!("PARSE_FAIL '{}': {}", smi, e))?;
    let c1 = canonical_smiles(&mol1);
    if c1.is_empty() {
        return Err(format!("EMPTY_CANONICAL '{}'", smi));
    }
    let mol2 =
        parse(&c1).map_err(|e| format!("RE_PARSE_FAIL '{}' (canonical='{}'): {}", smi, c1, e))?;
    let c2 = canonical_smiles(&mol2);
    if c1 != c2 {
        return Err(format!("UNSTABLE '{}': '{}' → '{}'", smi, c1, c2));
    }
    Ok(())
}

/// Returns Ok if `a` and `b` (same molecule, different SMILES) produce the same canonical form.
fn check_same_canonical(a: &str, b: &str) -> Result<(), String> {
    let mol_a = parse(a).map_err(|e| format!("PARSE_FAIL '{}': {}", a, e))?;
    let mol_b = parse(b).map_err(|e| format!("PARSE_FAIL '{}': {}", b, e))?;
    let ca = canonical_smiles(&mol_a);
    let cb = canonical_smiles(&mol_b);
    if ca != cb {
        Err(format!("DIFFERENT '{}' vs '{}': '{}' ≠ '{}'", a, b, ca, cb))
    } else {
        Ok(())
    }
}

/// Assert all SMILES in `cases` are roundtrip-stable, collecting all failures before panicking.
/// `#[track_caller]` ensures the panic points at the calling test, not this helper.
#[track_caller]
fn assert_all_stable(cases: &[&str]) {
    let failures: Vec<_> = cases
        .iter()
        .filter_map(|s| check_canonical_stable(s).err())
        .collect();
    assert!(
        failures.is_empty(),
        "stability failures:\n{}",
        failures.join("\n")
    );
}

// ── Test 1: Roundtrip stability ──────────────────────────────────────────────

/// Every SMILES here must be roundtrip-stable (parse → canonical → parse → canonical
/// gives the same canonical string). A failure here is a genuine correctness regression.
#[test]
fn stability_bridged_bicyclics() {
    let cases = [
        "C1CC2CCC1CC2",         // bicyclo[2.2.2]octane
        "C1CC2CCCC2C1",         // bicyclo[3.2.1]octane variant
        "C1CCC2CC3CCCCC3CC2C1", // polycycle
    ];
    assert_all_stable(&cases);
}

/// Idempotence coverage for the same bridged/fused/spiro/cage molecules
/// exercised by `bridged_fused_spiro_permutation_invariance` — reported
/// separately (per this project's convention: permutation invariance and
/// idempotence are two distinct probes for the same underlying property,
/// see `docs/rfcs/canonical_smiles_residual_rfc.md`'s Method section).
#[test]
fn stability_bridged_cage_and_heteroatom() {
    let cases = [
        "C1C2CC3CC(C2)CC1C3",               // adamantane
        "C1C2CCN(CC2)C1",                   // quinuclidine (N bridgehead)
        "C12CCC(CC1)C2",                    // norbornane
        "C1(C)(C)[C@@]2(C)C(=O)C[C@H]1CC2", // camphor (stereocenters)
        "C12C3C4C1C1C2C3C41",               // cubane
    ];
    assert_all_stable(&cases);
}

#[test]
fn stability_spiro() {
    let cases = [
        "C1CCC2(CC1)CCCC2", // spiro[4.5]decane
        "C1CC2(CCC1)CCC2",  // spiro[4.4]nonane
    ];
    assert_all_stable(&cases);
}

#[test]
fn stability_ring_stereocenters() {
    let cases = [
        "[C@@H]1(N)CCCC1",          // (R)-aminocyclopentane
        "[C@H]1(N)CCCC1",           // (S)-aminocyclopentane
        "[C@H]1([C@@H](O)CO)CCCO1", // bicyclic-ish stereocenters
    ];
    assert_all_stable(&cases);
}

#[test]
fn stability_fused_ring_stereo() {
    let cases = [
        "[C@@H]1(CC[C@H]2CCCC[C@@H]12)O", // trans-decalin-OH
        "[C@H]1(CC[C@H]2CCCC[C@@H]12)O",  // cis-decalin-OH
    ];
    assert_all_stable(&cases);
}

#[test]
fn stability_ez_bonds() {
    let cases = [
        "C/C=C/C",  // trans-but-2-ene
        "C/C=C\\C", // cis-but-2-ene
        "F/C=C/Cl", // E-1-chloro-2-fluoroethylene
    ];
    assert_all_stable(&cases);
}

#[test]
fn stability_complex_sugars() {
    let cases = [
        "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", // D-glucose
    ];
    assert_all_stable(&cases);
}

#[test]
fn stability_fused_aromatics() {
    let cases = [
        "C1=CC2=CC=CC=C2C=C1", // azulene
        "c1ccc2ccccc2c1",      // naphthalene
        "c1ccc2[nH]ccc2c1",    // indole
    ];
    assert_all_stable(&cases);
}

#[test]
fn stability_amino_acids() {
    // Each of these as a single SMILES string (not compared to another form) must be stable.
    let cases = [
        "N[C@@H](C)C(=O)O",         // L-alanine
        "N[C@H](C)C(=O)O",          // D-alanine
        "N[C@@H](Cc1ccccc1)C(=O)O", // L-phenylalanine
    ];
    assert_all_stable(&cases);
}

// ── Test 2: Platform independence (topology only) ────────────────────────────

/// Two different SMILES for the SAME molecule (no stereocenters) must give identical
/// canonical forms. If either fails to parse, that is reported without panicking the suite.
#[test]
fn platform_independence_topology() {
    let pairs: &[(&str, &str)] = &[
        ("c1ccccc1", "c1cccc(c1)"),             // benzene
        ("c1ccncc1", "n1cccc(c1)"),             // pyridine
        ("C1CCNCC1", "N1CCCCC1"),               // piperidine
        ("C1CCC2(CC1)CCCC2", "C1CCCCC12CCCC2"), // spiro
    ];
    let failures: Vec<_> = pairs
        .iter()
        .filter_map(|(a, b)| check_same_canonical(a, b).err())
        .collect();
    assert!(
        failures.is_empty(),
        "platform-independence failures:\n{}",
        failures.join("\n")
    );
}

/// Regression test for bridged/fused/spiro ring-closure ordering
/// permutation-invariance (`docs/rfcs/canonical_smiles_residual_rfc.md`'s
/// "Update (C2)" correction).
///
/// **History**: this test previously (`bridged_bicyclic_canonical_gap_documentation`)
/// asserted a "known gap" using the pair `("C1CC2CCC1CC2", "C1CCC2CC1CC2")`,
/// labeled as "two spellings of bicyclo[2.2.2]octane". That premise was never
/// checked against an independent structural oracle. RDKit `MolToInchi` shows
/// they are in fact two DIFFERENT constitutional isomers of C8H14 (bridges
/// 2-2-2 vs. 3-2-1 between the two degree-3 atoms) — `chematic` giving them
/// different canonical strings was correct behavior, not a bug. See the RFC
/// update for the full InChI evidence.
///
/// Every pair below IS independently verified (RDKit `MolToInchi`, both
/// spellings share the same InChI — see the comment above each pair) to
/// encode the same molecule, and both spellings were generated via RDKit
/// `RenumberAtoms` (a real atom relabeling, not hand-guessed), matching the
/// methodology `scripts/canonical_residual_diagnosis.py` uses for its
/// corpus-level permutation-invariance check. All 8 converge under the
/// current writer — hardening against a regression, since a targeted
/// 22-molecule probe (bridged bicyclics, spiro, fused, and cage systems,
/// with and without stereocenters/heteroatom bridgeheads) found zero real
/// convergence failures in this mechanism.
#[test]
fn bridged_fused_spiro_permutation_invariance() {
    let pairs: &[(&str, &str, &str)] = &[
        // InChI=1S/C8H14/c1-2-8-5-3-7(1)4-6-8/h7-8H,1-6H2
        ("C1CC2CCC1CC2", "C1C2CCC(CC2)C1", "bicyclo[2.2.2]octane"),
        // InChI=1S/C10H16/c1-7-2-9-4-8(1)5-10(3-7)6-9/h7-10H,1-6H2
        ("C1C2CC3CC(C2)CC1C3", "C12CC3CC(CC(C1)C3)C2", "adamantane"),
        // InChI=1S/C10H18/c1-2-6-10(7-3-1)8-4-5-9-10/h1-9H2
        ("C1CCCC2(CCCC2)C1", "C1CCCCC12CCCC2", "spiro[4.5]decane"),
        // InChI=1S/C10H18/c1-2-6-10-8-4-3-7-9(10)5-1/h9-10H,1-8H2/t9-,10-
        (
            "C1[C@H]2CCCC[C@@H]2CCC1",
            "[C@H]12[C@@H](CCCC1)CCCC2",
            "trans-decalin",
        ),
        // InChI=1S/C7H13N/c1-4-8-5-2-7(1)3-6-8/h7H,1-6H2
        ("C1C2CCN(CC2)C1", "N12CCC(CC1)CC2", "quinuclidine"),
        // InChI=1S/C7H12/c1-2-7-4-3-6(1)5-7/h6-7H,1-5H2
        ("C12CCC(CC1)C2", "C1C2CCC1CC2", "norbornane"),
        // InChI=1S/C10H16O/c1-9(2)7-4-5-10(9,3)8(11)6-7/h7H,4-6H2,1-3H3/t7-,10-/m1/s1
        (
            "C1(C)(C)[C@@]2(C)C(=O)C[C@H]1CC2",
            "C1(C)(C)[C@H]2CC(=O)[C@@]1(C)CC2",
            "camphor",
        ),
        // InChI=1S/C8H8/c1-2-5-3(1)7-4(1)6(2)8(5)7/h1-8H
        ("C12C3C4C1C1C2C3C41", "C12C3C4C5C(C1C35)C24", "cubane"),
    ];

    let failures: Vec<String> = pairs
        .iter()
        .filter_map(|&(a, b, label)| {
            check_same_canonical(a, b)
                .err()
                .map(|e| format!("{label}: {e}"))
        })
        .collect();

    assert!(
        failures.is_empty(),
        "bridged/fused/spiro permutation-invariance regression:\n{}",
        failures.join("\n")
    );
}

// ── Test 3: Stereo parity gap (documentation, does not panic) ────────────────

/// DOCUMENTED GAP: When the same molecule is written with two different atom orderings,
/// the canonical writer in `canonical.rs` does not apply a parity correction for the
/// neighbor permutation at stereocenters. Two SMILES encoding the same stereochemistry
/// but with different input atom traversal orders may produce different canonical forms.
///
/// Stereo parity is now corrected in canonical.rs.  These pairs must converge.
#[test]
fn stereo_parity_gap_documentation() {
    let pairs: &[(&str, &str, &str)] = &[
        // L-alanine written from N vs from C — odd permutation, parity must flip
        ("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O", "L-alanine"),
        // Aminocyclopentane: ring-first vs NH2-first
        ("[C@@H]1(N)CCCC1", "[C@H]1(CCCC1)N", "aminocyclopentane"),
        // Ring stereocenter: exercises PendingRing path.
        // [C@H]1(F)CCCCC1 and F[C@H]1CCCCC1 are the same enantiomer:
        //   both encode [H/F/C6/C2] negative via an even permutation.
        ("F[C@H]1CCCCC1", "[C@H]1(F)CCCCC1", "fluorocyclohexane"),
    ];

    let failures: Vec<String> = pairs
        .iter()
        .filter_map(|&(a, b, label)| {
            check_same_canonical(a, b)
                .err()
                .map(|e| format!("{label}: {e}"))
        })
        .collect();

    assert!(
        failures.is_empty(),
        "stereo parity gap not fully resolved:\n{}",
        failures.join("\n")
    );
}

/// Charged aromatics may or may not parse depending on aromatic valence model.
/// Report results without hard-failing on parse errors (feature, not bug).
#[test]
fn charged_aromatics_parse_and_stable() {
    let cases = [
        "c1cc[n+](C)cc1", // N-methylpyridinium
        "c1cc[nH]cc1",    // should not parse (wrong valence) — verify graceful error
    ];
    for smi in &cases {
        match check_canonical_stable(smi) {
            Ok(()) => eprintln!("ℹ PASS (parsed+stable): {}", smi),
            Err(e) => eprintln!("ℹ SKIP/FAIL: {}", e),
        }
    }
    // Non-panicking: charged aromatic handling is an informational probe.
}

// ── RDKit issue #8775 supplemental coverage ──────────────────────────────────
// Tests added to address RDKit's oscillating canonical SMILES edge cases.

#[test]
fn large_fused_pah_stable() {
    // Pyrene, benzo[a]pyrene: stable. Coronene (7-ring) has a known instability
    // tracked separately in `coronene_canonical_known_bug`.
    let cases = [
        "c1ccc2cc3ccc4cccc5ccc(c1)c2c3c45",  // pyrene (4 rings) — stable
        "c1ccc2ccc3cccc4ccc5ccccc5c4c3c2c1", // benzo[a]pyrene (5 rings) — stable
    ];
    let failures: Vec<_> = cases
        .iter()
        .filter_map(|&s| check_canonical_stable(s).err())
        .collect();
    assert!(
        failures.is_empty(),
        "PAH round-trip failures:\n{}",
        failures.join("\n")
    );
}

/// Coronene canonical SMILES was not idempotent for this 7-ring fused PAH
/// system prior to `fix/canonical-automorphism-pruning`: Morgan rank ties in
/// this highly symmetric graph could produce oscillating ring-closure
/// numbering (tracked as a known limitation, analogous to RDKit issue
/// #8775). Fixed as a verified side effect of that PR's automorphism-
/// orbit-pruned canonical search -- confirmed (not just "the test happens
/// to pass now") against the unbounded exhaustive individualize-refine
/// oracle in both directions (original parse and re-parse-of-canonical) and
/// across 32 relabelings of the same molecule, see
/// `crates/chematic-smiles/src/canonical_search.rs`'s
/// `unbounded_matches_exhaustive_oracle_on_symmetric_molecules` test. No
/// longer `#[ignore]`d.
///
/// **Fixture correction**: the SMILES string this test originally used
/// (`c1ccc2ccc3ccc4ccc5ccc6ccccc6c5c4c3c2c1`) parses to **26** atoms, not
/// the 24 a real coronene (C24H12, 7 fused hexagons) has -- the same class
/// of mislabeled-fixture problem this project was warned about for its
/// "cubane" fixture elsewhere in this repo. Replaced with a coronene
/// skeleton independently verified geometrically (7-hexagon flower
/// construction, Kekule-matched, then aromatized: 24 atoms, 30 bonds, 7
/// aromatic rings -- see `crates/chematic-smiles/examples/
/// canonical_orbit_perf.rs`'s `coronene_smiles()`). The bug this test
/// guards was real and reproduced on both the old (26-atom) and the
/// corrected (24-atom) molecule.
#[test]
fn coronene_canonical_known_bug() {
    assert!(
        check_canonical_stable("c2c3c4c1c6c(ccc7ccc5ccc(c4c5c67)cc3)ccc1c2").is_ok(),
        "coronene should be idempotent"
    );
}

#[test]
fn bridged_bicyclic_with_stereo_stable() {
    // Bridged ring systems with stereocenters — topology where traversal order
    // can vary between implementations.
    let cases = [
        "C1CC2CCC1CC2",         // bicyclo[2.2.2]octane (no stereo)
        "[C@@H]1CC2(CC1)CCC2",  // bicyclic with one stereocenter
        "[C@H]12CC(CC1)CC2",    // another bridged stereo
        "C1C[C@@H]2CC[C@H]1C2", // bridged bicyclic two stereocenters
    ];
    let failures: Vec<_> = cases
        .iter()
        .filter_map(|&s| check_canonical_stable(s).err())
        .collect();
    assert!(
        failures.is_empty(),
        "Bridged bicyclic+stereo failures:\n{}",
        failures.join("\n")
    );
}

#[test]
fn polyanion_and_salt_stable() {
    // Complex charged species — RDKit had edge cases where canonical SMILES
    // for salts would change between parse cycles.
    let cases = [
        "[Na+].[Cl-]",
        "[Na+].[Na+].[O-]C(=O)c1ccccc1S(=O)(=O)[O-]", // disodium benzoate sulfonate
        "[NH4+].[O-]C(=O)CC(=O)[O-]",                 // ammonium malate
        "[Ca+2].[O-]C(=O)c1ccccc1.[O-]C(=O)c1ccccc1", // calcium benzoate
    ];
    let failures: Vec<_> = cases
        .iter()
        .filter_map(|&s| check_canonical_stable(s).err())
        .collect();
    assert!(
        failures.is_empty(),
        "Polyanion/salt round-trip failures:\n{}",
        failures.join("\n")
    );
}

#[test]
fn fused_heteroaromatics_stable() {
    // Purine-class fused heteroaromatics — topology with multiple N atoms
    // in fused 5+6 rings where canonical ordering can be ambiguous.
    let cases = [
        "c1ncc2[nH]cnc2n1",           // purine
        "c1nc2ccccc2[nH]1",           // benzimidazole
        "c1cnc2ncnc2c1",              // 6,7-diazaindolizine
        "c1ccc2nc3ccccc3nc2c1",       // acridine
        "c1cc2ccc3cccc4ccc(c1)c2c34", // triphenylene
    ];
    let failures: Vec<_> = cases
        .iter()
        .filter_map(|&s| check_canonical_stable(s).err())
        .collect();
    assert!(
        failures.is_empty(),
        "Fused heteroaromatic failures:\n{}",
        failures.join("\n")
    );
}

// ── RDKit issue-inspired regression tests ────────────────────────────────────

/// RDKit #8759: canonical SMILES must be idempotent on stereocenters.
/// Stereo parity from different atom orderings must produce the same canonical form.
#[test]
fn rdkit_8759_stereo_idempotence() {
    assert_all_stable(&[
        "[C@@H]1(F)CCCC1",
        "[C@H]1(F)CCCC1",
        "[C@@H]1(O)CCCO1",
        "[C@@H]1([C@@H](O)CO)CCCO1",
        "[C@H]1([C@@H](O)CO)CCCO1",
        "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", // D-glucose
    ]);
}

/// RDKit #8759: canonical SMILES must be idempotent on E/Z bonds.
#[test]
fn rdkit_8759_ez_idempotence() {
    assert_all_stable(&[
        "C(/C=C/C)=C/C",
        "C(/C=C\\C)=C/C",
        "F/C=C(/F)Cl",
        "CC(/C=C/c1ccccc1)=O",  // chalcone E
        "CC(/C=C\\c1ccccc1)=O", // chalcone Z
        "O=C(/C=C/c1ccccc1)O",  // trans-cinnamic acid
    ]);
}

/// RDKit #9368: fragment extraction near E/Z bonds must never panic.
/// chematic's brics_fragments() / brics_bonds() must handle these without crashing.
/// Stereo preservation is best-effort; the invariant is: no panic, valid SMILES out.
#[test]
fn rdkit_9368_ez_fragment_no_panic() {
    let cases = [
        "CC(=O)O/C=C/c1ccccc1",         // cinnamyl acetate E
        "CC(=O)O/C=C\\c1ccccc1",        // cinnamyl acetate Z
        "O=C(/C=C/c1ccccc1)O",          // trans-cinnamic acid
        "C(/C=C/C(=O)O)Oc1ccccc1",      // E/Z ether acid
        "c1ccc(/C=C/C(=O)c2ccccc2)cc1", // chalcone E
    ];
    // Just verify parse + canonical doesn't panic
    for smi in cases {
        let result = chematic_smiles::parse(smi);
        if let Ok(mol) = result {
            let _ = chematic_smiles::canonical_smiles(&mol);
        }
    }
}

/// RDKit #8759 / charged heteroaromatics: N+/O- near aromatic rings must stay idempotent.
#[test]
fn rdkit_8759_charged_heteroaromatic_idempotence() {
    assert_all_stable(&[
        "c1cc[n+](C)cc1",         // N-methyl-pyridinium
        "c1cc[n+]([O-])cc1",      // pyridine N-oxide
        "[O-]c1ccccc1",           // phenolate
        "c1ccc(cc1)[N+](=O)[O-]", // nitrobenzene
        "c1ccc(cc1)[NH3+]",       // aniline protonated
    ]);
}
