//! Canonical-Stereo-D0: `assign_cip`'s (legacy) CIP element multiset must be
//! stable across a `canonical_smiles` round trip, for each of the three
//! parser bond-creation paths that can carry a stashed aromatic direction
//! (see `chematic_smiles::parser`'s `resolve_aromatic_direction_stash`).
//!
//! Before this fix, the ring-closure and branch-attachment paths stored a
//! `/`/`\` between two aromatic atoms as the bond's own literal `Up`/`Down`
//! order instead of stashing it -- inconsistent with the (already-correct)
//! plain chain-edge path. That inconsistency meant the SAME physical bond's
//! visibility to `assign_ez` (which reads `bond.order` directly, never the
//! stash) could change between the direct parse and a re-parse of
//! chematic's own canonical output, purely depending on which of the three
//! paths a given canonical DFS happened to route the bond through.
//!
//! This is a stability fix, not a correctness fix for E/Z detection itself:
//! `assign_ez` still doesn't read `bond_direction` at all, so these
//! molecules correctly show an EMPTY CIP set both before and after the
//! round trip here -- consistently, not as "no stereo" (an independent
//! RDKit check confirms real E/Z exists on these bonds; see EZ-A0/EZ-S1,
//! `docs/verification_coverage.md`). These tests intentionally do NOT
//! assert what the "right" E/Z value should be -- only that the round trip
//! doesn't change what's visible.

use chematic_chem::assign_cip;
use chematic_core::Molecule;
use chematic_smiles::{canonical_smiles, parse};

fn sorted_codes(mol: &Molecule) -> Vec<String> {
    let mut codes: Vec<String> = assign_cip(mol)
        .assignments
        .iter()
        .map(|(_, c)| format!("{c:?}"))
        .collect();
    codes.sort();
    codes
}

#[track_caller]
fn assert_cip_multiset_stable_across_round_trip(smi: &str, path_label: &str) {
    let mol = parse(smi).unwrap_or_else(|e| panic!("{path_label}: parse '{smi}': {e}"));
    let before = sorted_codes(&mol);

    let c1 = canonical_smiles(&mol);
    let mol2 = parse(&c1).unwrap_or_else(|e| panic!("{path_label}: re-parse '{c1}': {e}"));
    let after = sorted_codes(&mol2);

    assert_eq!(
        before, after,
        "{path_label}: CIP element multiset must be unchanged across a canonical \
         round trip ('{smi}' -> '{c1}')"
    );
}

#[test]
fn chain_edge_path_cip_multiset_stable() {
    assert_cip_multiset_stable_across_round_trip(r"N=c1\c(O)c(O)\c1=N", "path1(chain-edge)");
}

#[test]
fn ring_closure_path_cip_multiset_stable() {
    assert_cip_multiset_stable_across_round_trip(r"C/N=c1ccccc/1", "path2(ring-closure)");
}

#[test]
fn branch_attachment_path_cip_multiset_stable() {
    assert_cip_multiset_stable_across_round_trip(r"Cc1ccc(/c1)N", "path3(branch-attachment)");
}

#[test]
fn real_corpus_case_cip_multiset_stable_and_existing_stereocenter_unaffected() {
    // The actual molecule that surfaced this bug: a real, unrelated
    // tetrahedral stereocenter (`[C@H]`, CIP S) plus an aromatic mancude
    // ring flanked by an exocyclic imine. The S center must survive
    // unchanged; no extra element may appear or disappear alongside it.
    let smi = "O=C(Nc1ccc(C[C@H](/N=c2/c(N3CCSCC3)c(O)c2=O)C(=O)O)cc1)c1c(Cl)cncc1Cl";
    let mol = parse(smi).unwrap();
    let before = sorted_codes(&mol);
    assert_eq!(
        before,
        vec!["S".to_string()],
        "test setup sanity: exactly one real stereocenter (S), no E/Z visible yet"
    );

    let c1 = canonical_smiles(&mol);
    let mol2 = parse(&c1).unwrap();
    let after = sorted_codes(&mol2);
    assert_eq!(
        before, after,
        "the real stereocenter must survive, and no CIP element may appear or \
         disappear as a side effect of the round trip"
    );
}
