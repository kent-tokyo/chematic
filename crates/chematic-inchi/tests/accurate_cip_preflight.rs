//! Corpus regression fixture for issue #161 (accurate-CIP dedup preflight).
//!
//! The 12 molecules below (SMILES quoted verbatim) are taken from
//! `validation/results/dedup_stereo_guard_corpus_audit.jsonl` /
//! `_summary.json` -- PR #156's own full 5,000-molecule corpus audit of
//! every molecule `has_unresolved_specified_tetrahedral_stereo` newly fails
//! closed. That audit's own classification (case A = legacy fails, accurate
//! resolves; case B = both chematic engines tie, matching RDKit's own
//! `rdCIPLabeler`'s inability to assign a label either) is **not** taken on
//! faith here: re-running these exact molecules against the current
//! `chematic-inchi`/`chematic-chem`/`chematic-cip` HEAD (this file does not
//! re-parse the mutable `~/Downloads/SMILES.csv` at test time -- no
//! JSON-parsing dev-dependency is available to this crate under its
//! file-ownership scope for this PR, so these SMILES are pinned literals,
//! same approach as `examples/dedup_stereo_guard_diagnosis.rs`) confirmed the
//! classification, with one correction: the audit JSONL's own `corpus_index`
//! is uniformly off by one from the molecule's 0-based position in the
//! corpus file (audit idx=N sits at 0-based line N+1 -- confirmed for 11/12
//! rows by canonicalizing both the audit's `input_smiles` and the corpus
//! line at `N+1` and finding they match), and the row labeled `idx=4412`
//! additionally has a bad `input_smiles` transcription that matches no
//! molecule anywhere in the current, SHA-256-pinned 5,000-molecule corpus
//! (`validation/manifests/dataset_provenance.json`). The corpus file itself
//! has **not** drifted since the audit -- an earlier draft of this file
//! claimed otherwise and was wrong. The real molecule at 0-based corpus line
//! 4413 (sourced by position, not by trusting that audit row's string) DOES
//! still trigger the guard and IS recovered by this preflight, same as the
//! other 6 originally-recovered molecules -- see `CASE_A_FULLY_RECOVERED`
//! below, now 7 entries, not 6. See the PR body for the full re-derivation,
//! the independent live RDKit 2026.03.3 cross-check (own isolated venv), and
//! the fresh, current 5,000-molecule corpus rerun this file's groupings are
//! consistent with.
//!
//! Requires the `native-inchi` feature.

#![cfg(feature = "native-inchi")]

use chematic_inchi::dedup::{
    DedupRelation, IdentityPolicy, compare_molecules, compare_molecules_with_accurate_cip_preflight,
};
use chematic_smiles::parse;

fn mol(smiles: &str) -> chematic_core::Molecule {
    parse(smiles).unwrap_or_else(|e| panic!("parse {smiles:?}: {e}"))
}

/// 7 of the audit's original 10 "case A" molecules (corrected from a
/// previous, mistaken count of 6 -- see the module doc comment): legacy CIP
/// fails to rank at least one specified tetrahedral centre, `CipMode::Accurate`
/// resolves every such centre, and the accurate-CIP preflight fully
/// recovers verified-comparison capability. Each recovered atom's code was
/// independently re-checked against a live RDKit 2026.03.3 oracle
/// (`rdCIPLabeler` + `FindMolChiralCenters(includeUnassigned=True,
/// useLegacyImplementation=False)`, own isolated venv) -- 15/15 atoms agree
/// exactly for the first 6 entries (see the PR body for the full table);
/// the 7th entry below (corpus line 4413) is the corrected fixture for what
/// the audit intended as `idx=4412` (see the module doc comment for why its
/// original `input_smiles` was wrong and how the replacement was sourced).
const CASE_A_FULLY_RECOVERED: &[(usize, &str)] = &[
    (
        196,
        "CCCCc1cn([C@H]2[C@H](C)CCC[C@@H]2C)c(=O)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)nc1",
    ),
    (1567, "NS(=O)(=O)OC[C@@]12CCCC[C@@H]1CCC2"),
    (
        4047,
        "O=C(Oc1c(O)cc(C(=O)O[C@H]2[C@H](OC(=O)c3cc(O)c(O)c(O)c3)C[C@](O)(C(=O)O)C[C@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
    ),
    (
        // Corrected fixture for the audit's mistranscribed "idx=4412" row.
        // Sourced directly from 0-based corpus line 4413 of
        // `~/Downloads/SMILES.csv` (SHA-256
        // 1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e,
        // matching `validation/manifests/dataset_provenance.json`'s
        // `sha256_at_baseline`) -- NOT the audit JSONL's own `input_smiles`
        // for that row, which does not canonically match this or any other
        // molecule in the corpus. A di-galloylquinic-acid family member,
        // same family as the 4047/4413/4509 entries here.
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
];

/// The 2 "case B" molecules from the audit: tricyclic-cage bridgeheads,
/// unresolved by *both* chematic CIP engines. Re-confirmed directly against
/// a live RDKit 2026.03.3 oracle in this PR's own venv: RDKit's own modern
/// CIP labeler leaves the exact same bridgehead atom unlabeled
/// (`Tet_CCW`/`None` -- defined parity, no CIP descriptor) -- a genuine,
/// shared CIP-ranking limit, not fixable by this preflight.
const CASE_B_GENUINE_TIE: &[(usize, &str)] = &[
    (
        443,
        "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@]23C[C@H]4C[C@H](C[C@H](C4)C2)C3)co1",
    ),
    (
        590,
        "O=c1cc(COC2CCOCC2)occ1OC(=O)[C@]12C[C@H]3C[C@H](C[C@H](C3)C1)C2",
    ),
];

/// 3 of the audit's original 10 "case A" molecules that this preflight's
/// implementation does **not** recover, for a reason the audit itself did
/// not anticipate: each has a 3-fold-symmetric substituent (an
/// adamantane-cage-like scaffold) whose 3 legacy-unresolved stereocentres
/// all land on the SAME `chematic_smiles::morgan_ranks` value. The accurate
/// CIP engine resolves all 3 individually (each independently re-confirmed
/// against RDKit's `rdCIPLabeler` here: 9/9 atoms agree, all labeled `s`),
/// but `accurate_stereo_supplement`'s cross-molecule correspondence
/// mechanism deliberately fails closed on ANY rank collision among flagged
/// atoms in the same molecule (`IdentityDiagnostic::AmbiguousStereoRankCorrespondence`)
/// -- see that function's doc comment in `src/dedup.rs` for why a bijection
/// can't be safely assumed here without a real automorphism-consistency
/// check (out of scope for this PR, flagged as a follow-up). This is a
/// deliberate, conservative recall/safety trade-off, not a bug: it never
/// produces a wrong answer, it simply declines to recover these 3
/// molecules rather than risk mispairing a symmetric stereocentre.
const CASE_A_BLOCKED_BY_TIED_RANK: &[(usize, &str)] = &[
    (
        1609,
        "O=C(NCCCCN1CCN(c2cccc3ccccc23)CC1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
    ),
    (
        1643,
        "COc1ccccc1N1CCN(CCCCNC(=O)C23C[C@H]4C[C@@H](C2)C[C@@H](C3)C4)CC1",
    ),
    (
        4178,
        "O=C(NCCCCN1CCCC(/C=C\\c2ccccc2)C1)C12C[C@H]3C[C@@H](C1)C[C@@H](C2)C3",
    ),
];

/// For every audited molecule (case A fully-recovered, case A
/// blocked-by-tied-rank, and case B alike -- all 12, with no "no longer
/// applicable" exception; see the module doc comment for why an earlier
/// version of this file wrongly carved one molecule out here), the PLAIN
/// (legacy-only) path is unchanged from PR #156's shipped behavior:
/// comparing the molecule against itself is `VerificationUnavailable`. This
/// is the "before" half of issue #161's before/after claim.
#[test]
fn plain_path_unchanged_unavailable_for_every_audited_molecule() {
    for &(idx, smiles) in CASE_A_FULLY_RECOVERED
        .iter()
        .chain(CASE_B_GENUINE_TIE)
        .chain(CASE_A_BLOCKED_BY_TIED_RANK)
    {
        let m = mol(smiles);
        assert_eq!(
            compare_molecules(&m, &m, IdentityPolicy::StandardInchiString),
            DedupRelation::VerificationUnavailable,
            "idx={idx}: plain path must still fail closed (PR #156 unchanged)"
        );
    }
}

/// Case A (fully recovered): the accurate-CIP preflight recovers
/// verified-comparison capability -- a molecule compared against itself
/// becomes `VerifiedDuplicate` instead of `VerificationUnavailable`.
#[test]
fn case_a_fully_recovered_molecules_become_verified_duplicate_under_preflight() {
    for &(idx, smiles) in CASE_A_FULLY_RECOVERED {
        let m = mol(smiles);
        assert_eq!(
            compare_molecules_with_accurate_cip_preflight(
                &m,
                &m,
                IdentityPolicy::StandardInchiString
            ),
            DedupRelation::VerifiedDuplicate,
            "idx={idx}: case A molecule should recover under the accurate-CIP preflight"
        );
    }
}

/// Case B: even with the preflight engaged, a genuine tie in *both* chematic
/// engines (matching RDKit's own inability to label the same centre) must
/// stay `VerificationUnavailable` -- never silently promoted.
#[test]
fn case_b_molecules_still_fail_closed_under_preflight() {
    for &(idx, smiles) in CASE_B_GENUINE_TIE {
        let m = mol(smiles);
        assert_eq!(
            compare_molecules_with_accurate_cip_preflight(
                &m,
                &m,
                IdentityPolicy::StandardInchiString
            ),
            DedupRelation::VerificationUnavailable,
            "idx={idx}: genuine tie (case B) must still fail closed under the preflight"
        );
    }
}

/// Case A blocked by tied morgan-rank correspondence: must ALSO still fail
/// closed under the preflight (this is the conservative, safe direction --
/// see `CASE_A_BLOCKED_BY_TIED_RANK`'s doc comment). A future, more capable
/// correspondence mechanism might recover these; this preflight must never
/// promote them via a mispaired guess in the meantime.
#[test]
fn case_a_blocked_by_tied_rank_still_fails_closed_under_preflight() {
    for &(idx, smiles) in CASE_A_BLOCKED_BY_TIED_RANK {
        let m = mol(smiles);
        assert_eq!(
            compare_molecules_with_accurate_cip_preflight(
                &m,
                &m,
                IdentityPolicy::StandardInchiString
            ),
            DedupRelation::VerificationUnavailable,
            "idx={idx}: tied-morgan-rank case must stay VerificationUnavailable, never a guess"
        );
    }
}

/// `IdentityPolicy::StandardInchiKey` and `IdentityPolicy::IsotopeIgnored`
/// are the other two policies the corpus audit found this guard fires
/// under (`StereoIgnored` deliberately does not check the guard at all --
/// see the module docs) -- spot-check one fully-recovered case-A molecule
/// under each to confirm the preflight isn't accidentally
/// `StandardInchiString`-only.
#[test]
fn case_a_recovery_also_applies_under_standard_inchi_key_and_isotope_ignored() {
    let m = mol("NS(=O)(=O)OC[C@@]12CCCC[C@@H]1CCC2");
    for policy in [
        IdentityPolicy::StandardInchiKey,
        IdentityPolicy::IsotopeIgnored,
    ] {
        assert_eq!(
            compare_molecules(&m, &m, policy),
            DedupRelation::VerificationUnavailable,
            "sanity ({policy:?}): plain path unavailable"
        );
        assert_eq!(
            compare_molecules_with_accurate_cip_preflight(&m, &m, policy),
            DedupRelation::VerifiedDuplicate,
            "{policy:?}: should recover under the accurate-CIP preflight"
        );
    }
}
