//! Canonical-Stereo-D0: `assign_cip`'s (legacy) CIP element multiset must be
//! stable across a `canonical_smiles` round trip, for each of the three
//! parser bond-creation paths that can carry a stashed aromatic direction
//! (see `chematic_smiles::parser`'s `resolve_aromatic_direction_stash`).
//!
//! Before the D0 fix, the ring-closure and branch-attachment paths stored a
//! `/`/`\` between two aromatic atoms as the bond's own literal `Up`/`Down`
//! order instead of stashing it -- inconsistent with the (already-correct)
//! plain chain-edge path. That inconsistency meant the SAME physical bond's
//! visibility to `assign_ez` (which reads `bond.order` directly, never the
//! stash) could change between the direct parse and a re-parse of
//! chematic's own canonical output, purely depending on which of the three
//! paths a given canonical DFS happened to route the bond through.
//!
//! At the time D0 landed, `assign_ez` didn't read `bond_direction` at all,
//! so every case here showed an EMPTY CIP set both before and after the
//! round trip -- these tests only ever checked that the round trip didn't
//! change what's (not) visible. EZ-S1 (`crates/chematic-chem/src/cip.rs`'s
//! `substituent_is_up`) closed that gap, so `chain_edge_path...` and
//! `branch_attachment_path...` below are unaffected (their molecules have
//! no substituent, or no double bond at all, for `assign_ez` to act on --
//! not a case this file exists to cover), while `ring_closure_path...` and
//! `real_corpus_case...` now assert the real, RDKit-confirmed E/Z value
//! survives the round trip, not just "still empty."
//!
//! `ring_closure_path_cip_multiset_stable`'s original D0-era SMILES
//! (`C/N=c1ccccc/1`, an exocyclic imine on an *unsubstituted* benzo ring)
//! is not RDKit-parseable and, worse, turned out to be a genuine CIP
//! priority tie between its two ring branches (the unsubstituted ring is
//! symmetric across the ipso/para axis) -- confirmed via
//! `compare_branches(mol, alkene_end, subs[0], subs[1]) ==
//! Ordering::Equal` in both argument orders. That's a missing tie guard in
//! `highest_stereo_sub` (fixed in EZ-S1 alongside the stash read -- see
//! `crates/chematic-chem/src/cip.rs`'s
//! `test_highest_stereo_sub_symmetric_ring_is_not_stereogenic`), not the
//! deeper shell-pooling-comparator instability it initially looked like.
//! Swapped here for an asymmetric ring (`[nH]` breaks the tie, and has a
//! real RDKit-confirmed E) so this file keeps testing what it's for --
//! round-trip stability of an actually-stereogenic bond through the
//! ring-closure parser path -- rather than duplicating the tie-guard
//! regression test.

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
    // RDKit rdCIPLabeler confirms E for this exact SMILES; `assert_eq!` on
    // the exact value below (not just before==after) also pins that value.
    let smi = r"C/N=c1cccc[nH]\1";
    let mol = chematic_smiles::parse(smi).unwrap();
    assert_eq!(sorted_codes(&mol), vec!["E".to_string()]);
    assert_cip_multiset_stable_across_round_trip(smi, "path2(ring-closure)");
}

#[test]
fn branch_attachment_path_cip_multiset_stable() {
    assert_cip_multiset_stable_across_round_trip(r"Cc1ccc(/c1)N", "path3(branch-attachment)");
}

#[test]
fn real_corpus_case_cip_multiset_stable_and_existing_stereocenter_unaffected() {
    // The actual molecule that surfaced this bug (see EZ-A0/EZ-S1): a real,
    // unrelated tetrahedral stereocenter (`[C@H]`, CIP S) plus an aromatic
    // mancude ring flanked by an exocyclic imine, RDKit-confirmed Z. Both
    // must survive the round trip; no extra element may appear or
    // disappear alongside them.
    let smi = "O=C(Nc1ccc(C[C@H](/N=c2/c(N3CCSCC3)c(O)c2=O)C(=O)O)cc1)c1c(Cl)cncc1Cl";
    let mol = parse(smi).unwrap();
    let before = sorted_codes(&mol);
    assert_eq!(
        before,
        vec!["S".to_string(), "Z".to_string()],
        "test setup sanity: the real stereocenter (S) plus the now-correctly-read \
         exocyclic-imine E/Z (Z)"
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
