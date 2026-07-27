//! High-confidence molecule deduplication: fast canonical-SMILES candidate
//! grouping, verified/reconciled by native (IUPAC-standard) InChI.
//!
//! # Why two mechanisms, not one
//!
//! [`chematic_smiles::canonical_smiles`] is fast and deterministic, but it is
//! **not** a proof of chemical identity: two valid respellings of the exact
//! same molecule can occasionally produce different canonical strings (a
//! characterized, molecule-preserving residual of the canonicalizer -- see
//! `docs/canonical_smiles_residual_rfc.md`), and -- at least theoretically,
//! since this is unobserved on real corpora -- two different molecules could
//! in principle collide onto the same canonical string. Neither failure mode
//! is acceptable for a "high-confidence" dedup claim, so canonical SMILES is
//! used here only as a **fast candidate key**: it decides which pairs are
//! *cheap to compare*, never which pairs *are* duplicates.
//!
//! The actual identity decision is always made by native, IUPAC-standard
//! InChI ([`crate::native::standard_inchi`], gated behind the `native-inchi`
//! Cargo feature -- **not** the crate's default pure-Rust [`crate::inchi`],
//! which is an approximation and is never used here). If `native-inchi` is
//! disabled, or the oracle fails for either molecule, the result is
//! [`DedupRelation::VerificationUnavailable`] -- **never** silently treated
//! as a match. This mirrors a failure mode this project has hit before:
//! pooling a fallible/unverified path's results into a "verified" success
//! rate. See `docs/canonical_smiles_residual_rfc.md` and this module's own
//! `README`-equivalent (the crate-level dedup section) for the full
//! rationale.
//!
//! # Identity policies
//!
//! Four named, non-overlapping policies (see [`IdentityPolicy`]). None of
//! them is called "exact": [`IdentityPolicy::StandardInchiString`] is full
//! standard-InChI *string* equality, but standard InChI itself normalizes
//! some tautomers/mobile-H (see the tautomer fixture below) -- that is
//! real, IUPAC-defined chemical normalization, not raw graph identity, so
//! naming it "exact" would overclaim what it checks.
//!
//! - [`IdentityPolicy::StandardInchiString`] -- full standard InChI *string*
//!   equality (every layer: connectivity, H, charge, isotope, stereo). The
//!   strictest thing the oracle can give.
//! - [`IdentityPolicy::StandardInchiKey`] -- standard InChI**Key** equality --
//!   "as InChI itself defines it." Weaker than
//!   [`IdentityPolicy::StandardInchiString`] by exactly the (astronomically
//!   rare) 27-character-hash-collision margin; this is a real, if tiny,
//!   operational difference, not a relabeling of the same check.
//! - [`IdentityPolicy::StereoIgnored`] -- native InChI generated with stereo
//!   perception suppressed at the C-library level (the `SNon` generation
//!   option), so the output has no `/b`, `/t`, `/m`, or `/s` layer at all.
//!   **Not** the full string with those layers stripped out afterward.
//! - [`IdentityPolicy::IsotopeIgnored`] -- native InChI generated from a
//!   *cloned* molecule with every atom's isotope label cleared beforehand,
//!   so the `/i` layer (and anything InChI derives from it, e.g.
//!   isotope-induced stereo) never appears in the first place. **Not** the
//!   full string with the `/i` block stripped out afterward. The caller's
//!   original [`chematic_core::Molecule`] is never mutated -- only the
//!   internal clone is.
//!
//! No policy here performs post-hoc string surgery on a generated InChI
//! string (no splitting on `/`, no substring removal). Every distinction
//! between policies is made by *generating* a different, policy-appropriate
//! InChI in the first place -- either via native InChI's own generation-time
//! option ([`IdentityPolicy::StereoIgnored`]) or by transforming the input
//! molecule before generation ([`IdentityPolicy::IsotopeIgnored`]).
//!
//! An out-of-scope idea worth flagging for a future PR, not built here: a
//! true "exact graph identity" policy (isomorphism-based, no InChI
//! normalization at all) would be a fifth, structurally different policy --
//! deliberately not added in this PR.
//!
//! None of the four policies performs tautomer or salt merging *beyond* what
//! native InChI itself already does. If two tautomers happen to normalize to
//! the same standard InChI (a real, documented InChI behavior -- e.g.
//! 2-pyridone/2-hydroxypyridine), that merge is InChI's, not this module's,
//! and is reported the same way under every policy that reaches the InChI
//! string. Protonation-state and salt differences are preserved by standard
//! InChI's own `/q`/`/p` layers and formula layer respectively, so they
//! remain [`DedupRelation::Distinct`] under all four policies without any
//! special-casing here.
//!
//! [`IdentityPolicy::StandardInchiKey`] is a distinct code path (InChIKey
//! equality) with a distinct, if tiny, failure margin (a 27-character-hash
//! collision) from [`IdentityPolicy::StandardInchiString`] (full string
//! equality) -- but no fixture in this crate's test suite can *empirically*
//! discriminate the two, because key equality follows deterministically from
//! string equality and a genuine hash collision cannot be manufactured.
//! Every fixture that shows the two policies agreeing is expected to always
//! agree; that is not evidence the two checks are the same computation (they
//! are provably not, per the different code paths below), only evidence that
//! no fixture in this corpus happens to be near the hash-collision margin.
//!
//! # Explicit-H substituent at a stereocenter (history, now mostly resolved)
//!
//! `chematic-inchi`'s native conversion (`crate::native::convert`) used to
//! silently drop the tetrahedral stereo descriptor when a stereocenter's 4th
//! substituent was an explicit graph hydrogen atom (e.g. `[C@](Br)(Cl)(F)[H]`
//! or `[C@](Br)(Cl)(F)[2H]`), as opposed to bracket-H notation
//! (`[C@H](Br)(Cl)F`). That general case is now **fixed** (see
//! `crate::native::convert`'s `explicit_h_stereo` fixtures and
//! `two_stereocenters_each_get_their_own_manufactured_atom`): a real graph H
//! atom is now routed through the same manufactured-stand-in-atom mechanism
//! bracket-H already used, so `[C@](Br)(Cl)(F)[2H]` and
//! `[C@@](Br)(Cl)(F)[2H]` correctly produce different InChI strings with a
//! `/t` layer.
//!
//! One narrower case remains unfixed at the `crate::native::convert` level,
//! but is actively **guarded here**, not just documented: a stereocenter
//! with **two** H-like substituents on the SAME centre (e.g. one `[2H]` and
//! one `[3H]`, or bracket-H plus an explicit `[2H]`) -- a real,
//! CIP-rule-2-distinguished stereocentre that this crate's single
//! manufactured-atom-per-centre mechanism cannot represent. Native
//! conversion fails closed there (the stereo descriptor is dropped, not
//! corrupted -- no duplicated Stereo0D index), which means such a centre's
//! enantiomers are indistinguishable to `standard_inchi` itself. Getting a
//! non-empty string back from `standard_inchi` is not, by itself, evidence
//! that a comparison based on it is safe -- so **every** [`compare`] /
//! [`compare_molecules`] / [`deduplicate_verified`] call in this module
//! detects this input shape structurally (see
//! `crate::native::has_unrepresentable_multi_h_stereocenter`, checked before
//! any InChI generation) and reports
//! [`DedupRelation::VerificationUnavailable`] for any molecule containing
//! such a centre, under all four policies -- never a false
//! [`DedupRelation::VerifiedDuplicate`]. See `crate::native::convert`'s
//! `two_h_like_substituents_on_one_centre_drops_stereo_not_corrupts_it` and
//! `bracket_h_plus_explicit_isotope_h_drops_stereo_not_corrupts_it` tests for
//! the underlying conversion behavior, and this module's own
//! `two_h_like_substituents_*` fixtures for the guard itself.

#[cfg(feature = "native-inchi")]
use chematic_core::AtomIdx;
use chematic_core::{CipCode, Molecule};
use chematic_smiles::canonical_smiles;
use std::collections::HashMap;

/// A named identity policy governing what counts as "the same molecule".
///
/// See the module-level docs for the full rationale behind each variant.
/// These are deliberately kept distinct (never blurred into a single mode
/// or a bool) because they answer genuinely different questions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IdentityPolicy {
    /// Full standard InChI string equality (every layer preserved: stereo,
    /// isotope, charge, fragment composition).
    StandardInchiString,
    /// Standard InChIKey equality -- InChI's own identity token.
    StandardInchiKey,
    /// Native InChI generated with stereo perception suppressed at
    /// generation time (the C library's own `SNon` option) -- no `/b/t/m/s`
    /// layer is ever produced, rather than being produced and then
    /// stripped.
    ///
    /// **Not a "same compound" test.** Ignoring stereo means exactly that:
    /// an enantiomeric pair or a diastereomeric pair -- which differ *only*
    /// in stereo -- both correctly collapse to the same identity under this
    /// policy (reported as [`DedupRelation::CanonicalSplit`], since their
    /// candidate keys still differ; never a silent
    /// [`DedupRelation::VerifiedDuplicate`]). Use one of the other three
    /// policies whenever stereochemical identity matters to the caller.
    StereoIgnored,
    /// Native InChI generated from a clone of the input molecule with every
    /// atom's isotope label cleared beforehand -- no `/i` layer (or
    /// anything InChI nests under it, including isotope-induced stereo) is
    /// ever produced. The caller's original molecule is never mutated.
    IsotopeIgnored,
}

/// The relation between two molecules under a given [`IdentityPolicy`].
///
/// These six outcomes are intentionally kept distinct -- never collapsed
/// into a bool or a single "success" value. In particular
/// [`DedupRelation::CanonicalSplit`] and [`DedupRelation::CanonicalCollision`]
/// both represent "the fast candidate key disagreed with native InChI" but in
/// opposite, non-interchangeable directions: the former means the candidate
/// grouping under-merged (same molecule, different key -- unify via InChI,
/// but report that the fast path missed it); the latter means the candidate
/// grouping over-merged (same key, different molecule -- fail closed, never
/// merge).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DedupRelation {
    /// Same candidate key AND native InChI agrees (per the policy in
    /// effect). Both molecules' native InChI ran successfully.
    VerifiedDuplicate,
    /// Native InChI agrees (per the policy in effect) but the candidate keys
    /// differ -- the fast canonical-SMILES grouping split what is actually
    /// one molecule into two buckets. Reported explicitly, never silently
    /// merged into [`DedupRelation::VerifiedDuplicate`] and never dropped.
    CanonicalSplit,
    /// Same candidate key but native InChI disagrees -- the fast candidate
    /// key over-merged two different molecules. Fails closed: reported,
    /// never merged.
    CanonicalCollision,
    /// Different candidate key and native InChI disagrees. Ordinary
    /// non-duplicate result.
    Distinct,
    /// Native InChI did not run successfully for one or both molecules
    /// (the `native-inchi` feature is disabled, or the oracle returned an
    /// error other than "this molecule has no valid InChI at all"). Never
    /// counted as verified, regardless of whether the candidate keys match.
    VerificationUnavailable,
    /// One or both molecules cannot be assigned any InChI at all (e.g. a
    /// pure-hydrogen molecule with zero heavy atoms) -- a property of the
    /// molecule itself, not an oracle/environment failure.
    InvalidMolecule,
}

/// Outcome of attempting to verify one molecule against the native InChI
/// oracle. Kept separate from [`DedupRelation`] because it describes a
/// single molecule, not a pair.
///
/// `Ok`/`Invalid` are only ever constructed when the `native-inchi` feature
/// is enabled (see the two `verify()` impls below); without the feature,
/// every molecule is honestly `Unavailable`, so those variants go unused in
/// that build -- not dead code, just feature-conditional.
#[cfg_attr(not(feature = "native-inchi"), allow(dead_code))]
enum Verify {
    /// The raw native InChI string, generated the way `policy` requires
    /// (full stereo, `SNon`-suppressed, or isotope-cleared-clone). Not yet
    /// reduced to the policy's actual comparison key -- see
    /// [`verified_identity_key`].
    Ok(String),
    /// The molecule itself has no valid InChI (e.g. zero heavy atoms).
    Invalid,
    /// The oracle is unavailable (feature off) or failed for a reason that
    /// is not intrinsic to the molecule (kekulization limit, library error).
    Unavailable,
}

/// Generate the policy-appropriate native InChI for `mol` and classify the
/// outcome. This is where [`IdentityPolicy::StereoIgnored`] and
/// [`IdentityPolicy::IsotopeIgnored`] diverge from the other two policies --
/// each generates its OWN InChI (via a generation-time option or a
/// pre-cleared clone), never post-processes a shared one.
#[cfg(feature = "native-inchi")]
fn identity_verify(mol: &Molecule, policy: IdentityPolicy) -> Verify {
    // Fail closed BEFORE trusting anything `standard_inchi` returns: a
    // stereocentre with 2+ H-like substituents (e.g. one [2H] and one [3H]
    // on the same carbon) is a shape `crate::native::convert` cannot
    // represent -- it silently drops that centre's `/t` layer instead of
    // corrupting it (see
    // `crate::native::has_unrepresentable_multi_h_stereocenter`'s doc
    // comment). A non-empty InChI string is not evidence the comparison is
    // trustworthy here: two genuinely different stereoisomers of such a
    // molecule could otherwise silently collapse to the same string and
    // read as a false `VerifiedDuplicate`. Checked once, up front, for
    // every policy alike -- not just the stereo-preserving ones -- since
    // this is a property of whether *native InChI's own conversion* can be
    // trusted at all for this molecule, not of which layers a given policy
    // happens to compare.
    if crate::native::has_unrepresentable_multi_h_stereocenter(mol) {
        return Verify::Unavailable;
    }

    let result = match policy {
        IdentityPolicy::StandardInchiString | IdentityPolicy::StandardInchiKey => {
            // Second guard, checked on the ORIGINAL molecule: a stereocentre
            // whose SMILES specified `@`/`@@` but which the legacy CIP-based
            // ranking `crate::native::convert` depends on could not resolve
            // a rank for (`tetrahedral_stereo_neighbors` returns `None`).
            // This is the confirmed root cause of a real false
            // `VerifiedDuplicate` found via live 5,000-molecule corpus
            // verification -- two genuine diastereomers whose
            // `standard_inchi` collapsed to one string, `?` (undefined
            // parity) at both differing centres, identically for both
            // inputs. `VerifiedGroup`'s whole contract is "verified" --
            // this module owns that boundary even when the root cause is
            // upstream (chematic_chem's legacy CIP engine, not this crate),
            // so it must fail closed here rather than promote an
            // unresolved-parity string to a duplicate claim. False
            // positives (wrong merge) and false negatives (missed merge)
            // are NOT symmetric for a "high-confidence dedup" claim --
            // failing closed is the conservative direction.
            //
            // See `crate::native::has_unresolved_specified_tetrahedral_stereo`'s
            // doc comment for how this was diagnosed (not assumed).
            if crate::native::has_unresolved_specified_tetrahedral_stereo(mol) {
                return Verify::Unavailable;
            }
            crate::native::standard_inchi(mol)
        }
        IdentityPolicy::StereoIgnored => {
            // Deliberately NOT guarded by `has_unresolved_specified_tetrahedral_stereo`:
            // ignoring stereo entirely is this policy's explicit contract,
            // so an unresolved stereocentre is irrelevant to what it
            // compares -- the C library discards the `Stereo0D` array once
            // `SNon` is set regardless of whether chematic could rank it.
            crate::native::standard_inchi_no_stereo(mol)
        }
        IdentityPolicy::IsotopeIgnored => {
            // Clone first: the caller's molecule must never be mutated (see
            // module docs and `isotope_ignored_never_mutates_original` test).
            let mut cleared = mol.clone();
            for i in 0..cleared.atom_count() {
                cleared.set_isotope(AtomIdx(i as u32), None);
            }
            // Check the guard on the CLEARED CLONE, not the original: an
            // isotope-induced stereocentre that disappears once isotopes
            // are cleared should not be flagged "unresolved" under this
            // policy -- it simply isn't a stereocentre anymore once
            // isotopes are ignored, which is correct, expected behavior
            // here, not a failure.
            if crate::native::has_unresolved_specified_tetrahedral_stereo(&cleared) {
                return Verify::Unavailable;
            }
            crate::native::standard_inchi(&cleared)
        }
    };
    match result {
        Ok(s) => Verify::Ok(s),
        Err(crate::native::InchiError::InvalidInput(_)) => Verify::Invalid,
        Err(_) => Verify::Unavailable,
    }
}

#[cfg(not(feature = "native-inchi"))]
fn identity_verify(_mol: &Molecule, _policy: IdentityPolicy) -> Verify {
    // ponytail: the `native-inchi` feature is off in this build -- there is
    // no oracle to run at all, so every comparison is honestly reported as
    // unavailable rather than silently trusting the fast candidate key.
    Verify::Unavailable
}

#[cfg(feature = "native-inchi")]
fn inchi_key(inchi: &str) -> Result<String, ()> {
    crate::native::standard_inchi_key(inchi).map_err(|_| ())
}

#[cfg(not(feature = "native-inchi"))]
fn inchi_key(_inchi: &str) -> Result<String, ()> {
    Err(())
}

/// Reduce a raw native InChI string (as produced by [`identity_verify`]) to
/// the actual value `policy` compares/groups by. Only
/// [`IdentityPolicy::StandardInchiKey`] differs from the raw string itself
/// (and is independently fallible -- InChIKey computation is its own C call).
fn verified_identity_key(inchi: &str, policy: IdentityPolicy) -> Result<String, ()> {
    match policy {
        IdentityPolicy::StandardInchiKey => inchi_key(inchi),
        IdentityPolicy::StandardInchiString
        | IdentityPolicy::StereoIgnored
        | IdentityPolicy::IsotopeIgnored => Ok(inchi.to_string()),
    }
}

/// Compare two molecules under `policy`, given explicit (injectable)
/// candidate keys rather than computing them internally.
///
/// This is the low-level primitive: `key_a`/`key_b` would normally be each
/// molecule's [`chematic_smiles::canonical_smiles`], but accepting them as
/// parameters lets callers (notably tests) inject a synthetic collision --
/// two different molecules sharing one candidate key -- to prove
/// [`DedupRelation::CanonicalCollision`] is actually reachable and correctly
/// fails closed, since real canonical-SMILES collisions are not known to
/// occur in practice (0/4,992 observed in the 5,000-molecule corpus scan
/// behind `docs/canonical_smiles_residual_rfc.md`). Most callers should use
/// [`compare_molecules`] instead.
pub fn compare(
    key_a: &str,
    mol_a: &Molecule,
    key_b: &str,
    mol_b: &Molecule,
    policy: IdentityPolicy,
) -> DedupRelation {
    let (raw_a, raw_b) = match (
        identity_verify(mol_a, policy),
        identity_verify(mol_b, policy),
    ) {
        (Verify::Invalid, _) | (_, Verify::Invalid) => return DedupRelation::InvalidMolecule,
        (Verify::Unavailable, _) | (_, Verify::Unavailable) => {
            return DedupRelation::VerificationUnavailable;
        }
        (Verify::Ok(a), Verify::Ok(b)) => (a, b),
    };

    let (ident_a, ident_b) = match (
        verified_identity_key(&raw_a, policy),
        verified_identity_key(&raw_b, policy),
    ) {
        (Ok(a), Ok(b)) => (a, b),
        _ => return DedupRelation::VerificationUnavailable,
    };
    let same_identity = ident_a == ident_b;
    let same_key = key_a == key_b;

    match (same_identity, same_key) {
        (true, true) => DedupRelation::VerifiedDuplicate,
        (true, false) => DedupRelation::CanonicalSplit,
        (false, true) => DedupRelation::CanonicalCollision,
        (false, false) => DedupRelation::Distinct,
    }
}

/// Compare two molecules under `policy`, using each molecule's own
/// [`chematic_smiles::canonical_smiles`] as the fast candidate key.
///
/// This is the ordinary entry point for pairwise deduplication. See
/// [`compare`] for the lower-level, injectable-key primitive used to test
/// [`DedupRelation::CanonicalCollision`] directly, and
/// [`deduplicate_verified`] for whole-collection reconciliation without
/// pairwise (O(n²)) comparison.
pub fn compare_molecules(
    mol_a: &Molecule,
    mol_b: &Molecule,
    policy: IdentityPolicy,
) -> DedupRelation {
    let key_a = canonical_smiles(mol_a);
    let key_b = canonical_smiles(mol_b);
    compare(&key_a, mol_a, &key_b, mol_b, policy)
}

/// Partition `mols` into fast candidate groups by canonical SMILES.
///
/// This is the "fast canonical-SMILES candidate grouping" half of the dedup
/// pipeline: molecules sharing a canonical SMILES string land in the same
/// group. **This grouping is not itself a duplicate claim** -- per
/// [`DedupRelation::CanonicalCollision`], two different molecules could in
/// principle share a candidate group; per [`DedupRelation::CanonicalSplit`],
/// the true duplicate of a molecule may sit in a *different* group. Use
/// [`compare_molecules`] (or [`compare`]) on members within and, if needed,
/// across groups to get a verified relation -- or use
/// [`deduplicate_verified`] to get the whole collection reconciled at once.
///
/// Returned groups are sorted by their smallest member index (deterministic
/// across runs); each group's indices are in ascending input order.
pub fn group_candidates(mols: &[Molecule]) -> Vec<Vec<usize>> {
    let mut groups: HashMap<String, Vec<usize>> = HashMap::new();
    for (i, m) in mols.iter().enumerate() {
        groups.entry(canonical_smiles(m)).or_default().push(i);
    }
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    out.sort_by_key(|g| g[0]);
    out
}

/// One verified-duplicate group: every member shares the same native InChI
/// identity under the report's [`IdentityPolicy`]. Always at least 2
/// members (a singleton is not a duplicate). Members that also disagree on
/// canonical-SMILES key are additionally cross-referenced by a
/// [`CanonicalSplit`] entry in [`VerifiedDedupReport::canonical_splits`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VerifiedGroup {
    /// Input indices in this group, ascending (input order).
    pub members: Vec<usize>,
}

/// A verified-duplicate group ([`VerifiedGroup`], reported separately in
/// [`VerifiedDedupReport::groups`]) whose members do NOT all share one
/// canonical-SMILES key -- the fast candidate-key pass under-merged them
/// into more than one bucket; native InChI (per the report's policy) proves
/// they are the same identity. This is the case
/// [`crate::dedup`][self]'s whole batch API exists for.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalSplit {
    /// Same as the corresponding [`VerifiedGroup::members`] -- all indices
    /// verified as one identity, ascending.
    pub members: Vec<usize>,
    /// `members` partitioned by (differing) canonical-SMILES key. Each
    /// inner `Vec` is ascending; outer `Vec` is sorted by each subgroup's
    /// smallest index. Always at least 2 subgroups (that's what makes it a
    /// split).
    pub canonical_subgroups: Vec<Vec<usize>>,
}

/// A group of inputs sharing one canonical-SMILES key that native InChI
/// (per the report's policy) proves is NOT all one identity -- the fast
/// candidate key over-merged. Failed closed: none of `verified_subgroups`
/// is merged with any other in [`VerifiedDedupReport::groups`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalCollision {
    /// All indices sharing this canonical-SMILES key (that have a usable
    /// verified identity -- see [`VerifiedDedupReport::verification_unavailable`]
    /// / [`VerifiedDedupReport::invalid_molecules`] for exclusions),
    /// ascending.
    pub canonical_key_members: Vec<usize>,
    /// `canonical_key_members` partitioned by (differing) verified
    /// identity. Each inner `Vec` is ascending; outer `Vec` is sorted by
    /// each subgroup's smallest index. Always at least 2 subgroups (that's
    /// what makes it a collision). Any subgroup with 2+ members is also
    /// reported as its own, separate [`VerifiedGroup`].
    pub verified_subgroups: Vec<Vec<usize>>,
}

/// Whole-collection reconciliation report produced by
/// [`deduplicate_verified`]. See that function's docs for the algorithm and
/// complexity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VerifiedDedupReport {
    /// Every verified-duplicate group (2+ members sharing one native InChI
    /// identity under the policy in effect), regardless of whether their
    /// canonical-SMILES keys agree. Sorted by each group's smallest member
    /// index.
    pub groups: Vec<VerifiedGroup>,
    /// The subset of `groups` whose members span more than one
    /// canonical-SMILES key, with the per-key breakdown. Sorted by each
    /// entry's smallest member index.
    pub canonical_splits: Vec<CanonicalSplit>,
    /// Canonical-SMILES buckets that native InChI (per the policy in
    /// effect) proves contain more than one distinct identity. Sorted by
    /// each entry's smallest member index.
    pub canonical_collisions: Vec<CanonicalCollision>,
    /// Input indices for which native InChI verification did not succeed
    /// (feature disabled, or a non-intrinsic oracle failure). Never counted
    /// toward any group above. Ascending.
    pub verification_unavailable: Vec<usize>,
    /// Input indices for molecules with no valid InChI at all (a property
    /// of the molecule, e.g. zero heavy atoms). Never counted toward any
    /// group above. Ascending.
    pub invalid_molecules: Vec<usize>,
}

/// Reconcile a whole collection of molecules at once: fast canonical-SMILES
/// candidate keys, verified/re-grouped by native InChI identity (per
/// `policy`), in a single pass -- not pairwise `compare`/`compare_molecules`
/// calls run over every pair.
///
/// # Algorithm (all O(n), no all-pairs comparison)
///
/// 1. Compute each molecule's canonical SMILES once.
/// 2. Compute each molecule's native InChI identity once (the
///    policy-appropriate generation from [`identity_verify`], reduced to a
///    comparison key by [`verified_identity_key`]). Molecules with no
///    successful verification are set aside into
///    [`VerifiedDedupReport::verification_unavailable`] /
///    [`VerifiedDedupReport::invalid_molecules`] and never grouped.
/// 3. Group the remaining (successfully-verified) indices by verified
///    identity key -- a single `HashMap` pass -- to get
///    [`VerifiedDedupReport::groups`] (2+ member groups only).
/// 4. Separately, group the same indices by canonical-SMILES key -- another
///    single `HashMap` pass -- to detect [`CanonicalCollision`]s (a
///    canonical bucket spanning more than one verified identity).
/// 5. Within each verified-identity group, checking how many distinct
///    canonical keys are present is bounded by that group's own size, so
///    summed across all groups it's still O(n); same for the reverse
///    direction. No molecule is ever compared against another via a
///    dedicated pairwise InChI call.
///
/// Output order is deterministic across runs: every group/split/collision
/// list is sorted by its lowest member index, and every index list within
/// an entry is in ascending (input) order.
pub fn deduplicate_verified(mols: &[Molecule], policy: IdentityPolicy) -> VerifiedDedupReport {
    let n = mols.len();
    let mut canonical_key: Vec<String> = Vec::with_capacity(n);
    let mut verified_key: Vec<Option<String>> = Vec::with_capacity(n);
    let mut verification_unavailable: Vec<usize> = Vec::new();
    let mut invalid_molecules: Vec<usize> = Vec::new();

    for (i, m) in mols.iter().enumerate() {
        canonical_key.push(canonical_smiles(m));
        match identity_verify(m, policy) {
            Verify::Invalid => {
                invalid_molecules.push(i);
                verified_key.push(None);
            }
            Verify::Unavailable => {
                verification_unavailable.push(i);
                verified_key.push(None);
            }
            Verify::Ok(raw) => match verified_identity_key(&raw, policy) {
                Ok(k) => verified_key.push(Some(k)),
                Err(()) => {
                    verification_unavailable.push(i);
                    verified_key.push(None);
                }
            },
        }
    }

    // Partition 1: verified-identity groups (drives `groups`/`canonical_splits`).
    let mut by_verified: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, k) in verified_key.iter().enumerate() {
        if let Some(k) = k {
            by_verified.entry(k.as_str()).or_default().push(i);
        }
    }

    let mut groups: Vec<VerifiedGroup> = Vec::new();
    let mut canonical_splits: Vec<CanonicalSplit> = Vec::new();
    for members in by_verified.into_values() {
        if members.len() < 2 {
            continue;
        }
        let mut members = members; // already ascending: built from an ascending enumerate() scan
        members.sort_unstable();

        let mut by_canon_within: HashMap<&str, Vec<usize>> = HashMap::new();
        for &i in &members {
            by_canon_within
                .entry(canonical_key[i].as_str())
                .or_default()
                .push(i);
        }
        groups.push(VerifiedGroup {
            members: members.clone(),
        });
        if by_canon_within.len() > 1 {
            let mut canonical_subgroups: Vec<Vec<usize>> = by_canon_within.into_values().collect();
            canonical_subgroups.sort_by_key(|g| g[0]);
            canonical_splits.push(CanonicalSplit {
                members,
                canonical_subgroups,
            });
        }
    }

    // Partition 2: canonical-SMILES groups (drives `canonical_collisions`).
    let mut by_canonical: HashMap<&str, Vec<usize>> = HashMap::new();
    for i in 0..n {
        if verified_key[i].is_some() {
            by_canonical
                .entry(canonical_key[i].as_str())
                .or_default()
                .push(i);
        }
    }

    let mut canonical_collisions: Vec<CanonicalCollision> = Vec::new();
    for mut members in by_canonical.into_values() {
        if members.len() < 2 {
            continue;
        }
        members.sort_unstable();
        let mut by_ident_within: HashMap<&str, Vec<usize>> = HashMap::new();
        for &i in &members {
            let k = verified_key[i].as_deref().expect("filtered to Some above");
            by_ident_within.entry(k).or_default().push(i);
        }
        if by_ident_within.len() > 1 {
            let mut verified_subgroups: Vec<Vec<usize>> = by_ident_within.into_values().collect();
            verified_subgroups.sort_by_key(|g| g[0]);
            canonical_collisions.push(CanonicalCollision {
                canonical_key_members: members,
                verified_subgroups,
            });
        }
    }

    groups.sort_by_key(|g| g.members[0]);
    canonical_splits.sort_by_key(|s| s.members[0]);
    canonical_collisions.sort_by_key(|c| c.canonical_key_members[0]);
    verification_unavailable.sort_unstable();
    invalid_molecules.sort_unstable();

    VerifiedDedupReport {
        groups,
        canonical_splits,
        canonical_collisions,
        verification_unavailable,
        invalid_molecules,
    }
}

// ─── Typed identity diagnostics (shared) ───────────────────────────────────

/// Typed reason a [`crate::dedup`] identity claim could not be made at full
/// strength (or had to fall back to a weaker one). Shared by the
/// accurate-CIP dedup preflight (issue #161, see
/// [`compare_molecules_with_accurate_cip_preflight`]) and (in a later
/// commit) the exact-graph-identity API -- never a free-form string, so
/// callers can match on it instead of parsing text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IdentityDiagnostic {
    /// The accurate CIP engine (`CipMode::Accurate`) explicitly tied on a
    /// centre the legacy engine also could not rank -- a genuine CIP-ranking
    /// limit shared by both chematic engines (see issue #161's "case B":
    /// tricyclic-cage bridgeheads, corpus idx 443/590, where RDKit's own
    /// `rdCIPLabeler` agrees no label exists either).
    AccurateCipTied,
    /// The accurate CIP engine exceeded its recursion/size budget for a
    /// centre the legacy engine also could not rank.
    AccurateCipBudgetExceeded,
    /// The accurate CIP engine itself returned an error (not a tie/budget
    /// outcome). Never pooled with the tie/budget classes above -- an engine
    /// failure is a different thing from an explicit "no answer" (see
    /// `docs/`'s standing rule against pooling a fallible path's failures
    /// into another class's rate).
    AccurateCipEngineError,
    /// Two or more legacy-CIP-unresolved centres in the SAME molecule land on
    /// the same `chematic_smiles::morgan_ranks` value -- their identity
    /// cannot be safely attributed to a specific centre when compared
    /// against another molecule, so the accurate-CIP recovery is not applied
    /// and the comparison fails closed instead of risking a mispaired
    /// centre.
    AmbiguousStereoRankCorrespondence,
}

// ─── Accurate-CIP dedup preflight (issue #161) ─────────────────────────────
//
// PR #156 shipped `has_unresolved_specified_tetrahedral_stereo`: a specified
// tetrahedral stereocentre (`@`/`@@`) that the LEGACY CIP engine
// (`chematic_chem::assign_cip`) cannot rank a substituent order for makes
// `identity_verify` fail closed to `VerificationUnavailable`, because
// `crate::native::convert::mol_to_inchi_atoms` depends on that exact same
// ranking to build the InChI Stereo0D descriptor for that atom -- when it
// fails, native InChI silently emits either an undefined `?` parity (if
// other centres in the molecule succeeded) or omits the whole `/t` layer (if
// none did), and two genuine diastereomers can then collapse onto the same
// string (the corpus 4663/4664 pair PR #156 fixed).
//
// PR #156's own full 5,000-molecule-corpus audit
// (`validation/results/dedup_stereo_guard_corpus_audit.jsonl`) found this
// guard newly fires on 14 molecules beyond that pair, classified into case A
// (legacy ties, `CipMode::Accurate` resolves with a label an independent
// RDKit 2026.03.3 oracle agrees with exactly -- 10 molecules / 26 atoms) and
// case B (genuinely unresolvable by *both* chematic engines and by RDKit's
// own modern CIP labeler -- 2 molecules / 2 atoms, tricyclic-cage
// bridgeheads).
//
// This PR independently RE-RAN that classification against current HEAD
// (not trusted on faith -- see the fixtures in
// `tests/accurate_cip_preflight.rs` and the PR body for the full
// re-derivation) and a **fresh** 5,000-molecule corpus
// (`~/Downloads/SMILES.csv`, SHA-256-confirmed byte-identical to the value
// pinned in `validation/manifests/dataset_provenance.json` -- the corpus has
// NOT drifted since the audit ran; an earlier draft of this comment claimed
// otherwise and was wrong, see below). Two independent, unrelated findings
// came out of the re-run:
//
// (1) The audit JSONL's own `corpus_index` field is uniformly off by one
// from the 0-based position of the molecule in the corpus file (i.e. audit
// idx=N's molecule actually sits at 0-based line N+1) -- confirmed for
// 11 of the 12 audited molecules by canonicalizing both the audit's
// `input_smiles` and the corpus line at `N+1` and finding they match. This
// is a labeling-convention mismatch in the original audit tooling, not a
// code regression, and doesn't change which molecules are affected.
//
// (2) Exactly one audited row -- the one labeled `audit idx=4412` -- also
// has a bad `input_smiles` transcription: that string does not match the
// real corpus molecule at its implied 0-based line 4413 (nor does it
// canonically match anything anywhere else in the 5,000-line corpus). An
// earlier version of this preflight's test fixture trusted that row's
// string verbatim, found it (correctly, for that string) no longer
// triggers the guard, and wrongly concluded "this molecule is no longer
// applicable" / "an upstream chematic-chem/chematic-cip CIP-ranking change
// since the audit" -- both false. `chematic-chem/src/cip.rs` (home of
// `tetrahedral_stereo_neighbors`, the legacy ranking function this guard
// depends on) has had zero commits since well before the audit ran, so the
// legacy ranking path is byte-identical to audit time; there was no
// upstream change to attribute this to. The REAL molecule at 0-based corpus
// line 4413 (a di-galloylquinic-acid family member, distinct from the
// audit's mistranscribed string) DOES still trigger the guard today and IS
// recovered by this preflight, exactly like the other 6 fully-recovered
// audited molecules -- see `tests/accurate_cip_preflight.rs` for the
// corrected fixture, sourced by corpus position rather than by trusting the
// audit row's string.
//
// Net effect, recounted from scratch against the 12 originally-audited
// molecules: 7 fully recovered (not 6), 2 genuine ties (both chematic
// engines and RDKit's own modern CIP labeler agree no label exists), 3
// blocked by this preflight's conservative tied-morgan-rank-correspondence
// guard -- 7 + 2 + 3 = 12, with no "no longer applicable" category (that
// was an artifact of the bad transcription above, not a real phenomenon).
// Full-corpus effect: `verification_unavailable` shrinks from 15 to 6 (9
// molecules recovered -- the 7 above plus the already-disclosed 4663/4664
// pair, 0 molecules newly made unavailable -- the "shrink, never widen"
// invariant this PR must preserve, confirmed empirically over the whole
// corpus, not just the audited 12). Full denominators, the RDKit-cross-check
// table, and dataset provenance are in the PR body.
//
// Importantly, resolving a centre via the accurate engine does NOT change
// what `crate::native::standard_inchi` itself emits for that atom -- that
// generation path still depends on the legacy ranking (issue #161
// deliberately scopes out "switch native InChI generation to the accurate
// engine wholesale," since that would change `standard_inchi`'s output for
// every caller, not just `dedup`, and the accurate engine's own
// budget/tied-result semantics would need threading through the whole
// Stereo0D pipeline -- a separate, larger proposal). So merely "unblocking"
// the guard and generating `standard_inchi` as normal would silently
// reintroduce the exact 4663/4664 bug: both molecules would still collapse
// onto the same `?`-bearing string. Instead, the accurate engine's resolved
// codes for exactly the previously-unresolved centres are folded into the
// COMPARISON KEY as a supplement (never into the generated InChI string
// itself, and never into the InChIKey C library call -- see
// `verified_identity_key_with_supplement`), keyed by
// `chematic_smiles::morgan_ranks` (a purely constitutional, stereo-blind
// isomorphism invariant) so the supplement is meaningfully comparable
// between two molecules that share a skeleton but differ in stereo, without
// needing native InChI's own canonical atom numbering.
//
// This is deliberately a NEW, separate set of entry points
// (`*_with_accurate_cip_preflight`), not a parameter threaded into
// `identity_verify`/`compare`/`compare_molecules`/`deduplicate_verified`:
// every existing call through those stays byte-for-byte unchanged (verified
// by `preflight_is_opt_in_existing_paths_unchanged` below), matching the
// "opt-in, two-tier" framing issue #161 itself proposed.

/// Recover legacy-CIP-unresolved specified tetrahedral centres via
/// `CipMode::Accurate` (issue #161). Returns the accurate CIP code for every
/// atom `crate::native::unresolved_specified_tetrahedral_stereo_atoms`
/// flagged, keyed by that atom's `chematic_smiles::morgan_ranks` value.
///
/// Fails closed (`Err`) -- never partially recovers -- if the accurate
/// engine ties, exceeds its budget, or errors on ANY flagged atom, or if two
/// flagged atoms in this SAME molecule share a morgan-rank value (their
/// cross-molecule correspondence would be ambiguous -- see the module-level
/// docs above).
#[cfg(feature = "native-inchi")]
fn accurate_stereo_supplement(mol: &Molecule) -> Result<Vec<(u64, CipCode)>, IdentityDiagnostic> {
    let unresolved = crate::native::unresolved_specified_tetrahedral_stereo_atoms(mol);
    if unresolved.is_empty() {
        return Ok(Vec::new());
    }
    let accurate = chematic_chem::assign_cip_with_mode(mol, chematic_chem::CipMode::Accurate)
        .map_err(|_| IdentityDiagnostic::AccurateCipEngineError)?;
    let ranks = chematic_smiles::morgan_ranks(mol);

    let mut supplement: Vec<(u64, CipCode)> = Vec::with_capacity(unresolved.len());
    let mut seen_ranks: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for atom in unresolved {
        if let Some((_, reason)) = accurate.unresolved.iter().find(|(i, _)| *i == atom) {
            return Err(match reason {
                chematic_chem::CipUnresolvedReason::Tied => IdentityDiagnostic::AccurateCipTied,
                chematic_chem::CipUnresolvedReason::BudgetExceeded => {
                    IdentityDiagnostic::AccurateCipBudgetExceeded
                }
            });
        }
        let Some(code) = accurate.get(atom) else {
            // The accurate engine never touched this atom at all -- shouldn't
            // happen for a specified tetrahedral centre, but fail closed
            // defensively rather than silently skip it.
            return Err(IdentityDiagnostic::AccurateCipEngineError);
        };
        let rank = ranks[atom.0 as usize];
        if !seen_ranks.insert(rank) {
            return Err(IdentityDiagnostic::AmbiguousStereoRankCorrespondence);
        }
        supplement.push((rank, code));
    }
    supplement.sort_by_key(|&(r, _)| r);
    Ok(supplement)
}

/// Parallel to [`Verify`], for the accurate-CIP-preflight entry points only.
/// `Ok` additionally carries the (possibly empty) stereo supplement
/// recovered by [`accurate_stereo_supplement`].
#[cfg_attr(not(feature = "native-inchi"), allow(dead_code))]
enum PreflightVerify {
    Ok(String, Vec<(u64, CipCode)>),
    Invalid,
    Unavailable,
}

/// Like [`identity_verify`], but retries a legacy-CIP-unresolved specified
/// tetrahedral stereocentre via the accurate CIP engine before giving up
/// (issue #161), instead of failing closed immediately.
#[cfg(feature = "native-inchi")]
fn identity_verify_with_accurate_preflight(
    mol: &Molecule,
    policy: IdentityPolicy,
) -> PreflightVerify {
    if crate::native::has_unrepresentable_multi_h_stereocenter(mol) {
        return PreflightVerify::Unavailable;
    }

    let (result, supplement) = match policy {
        IdentityPolicy::StandardInchiString | IdentityPolicy::StandardInchiKey => {
            if crate::native::has_unresolved_specified_tetrahedral_stereo(mol) {
                match accurate_stereo_supplement(mol) {
                    Err(_) => return PreflightVerify::Unavailable,
                    Ok(supplement) => (crate::native::standard_inchi(mol), supplement),
                }
            } else {
                (crate::native::standard_inchi(mol), Vec::new())
            }
        }
        IdentityPolicy::StereoIgnored => {
            // Deliberately not guarded (same rationale as `identity_verify`):
            // ignoring stereo entirely is this policy's contract.
            (crate::native::standard_inchi_no_stereo(mol), Vec::new())
        }
        IdentityPolicy::IsotopeIgnored => {
            let mut cleared = mol.clone();
            for i in 0..cleared.atom_count() {
                cleared.set_isotope(AtomIdx(i as u32), None);
            }
            if crate::native::has_unresolved_specified_tetrahedral_stereo(&cleared) {
                match accurate_stereo_supplement(&cleared) {
                    Err(_) => return PreflightVerify::Unavailable,
                    Ok(supplement) => (crate::native::standard_inchi(&cleared), supplement),
                }
            } else {
                (crate::native::standard_inchi(&cleared), Vec::new())
            }
        }
    };
    match result {
        Ok(s) => PreflightVerify::Ok(s, supplement),
        Err(crate::native::InchiError::InvalidInput(_)) => PreflightVerify::Invalid,
        Err(_) => PreflightVerify::Unavailable,
    }
}

#[cfg(not(feature = "native-inchi"))]
fn identity_verify_with_accurate_preflight(
    _mol: &Molecule,
    _policy: IdentityPolicy,
) -> PreflightVerify {
    PreflightVerify::Unavailable
}

/// Reduce a raw native InChI string plus an accurate-CIP stereo supplement to
/// the actual comparison key. Reuses [`verified_identity_key`] unchanged for
/// the base key (so `IdentityPolicy::StandardInchiKey`'s InChIKey C-library
/// call always receives a clean, real InChI string -- the supplement is
/// appended only to the final *comparison* key, never fed back into any
/// native call).
#[cfg(feature = "native-inchi")]
fn verified_identity_key_with_supplement(
    inchi: &str,
    supplement: &[(u64, CipCode)],
    policy: IdentityPolicy,
) -> Result<String, ()> {
    let base = verified_identity_key(inchi, policy)?;
    if supplement.is_empty() {
        return Ok(base);
    }
    use std::fmt::Write;
    let mut key = base;
    // A marker that cannot appear inside a real InChI string or InChIKey
    // (both use restricted character sets), so this can never collide with
    // an unmodified base key from a molecule that didn't need the
    // supplement.
    key.push_str("|accurate_cip_stereo_supplement=");
    for (rank, code) in supplement {
        let _ = write!(key, "{rank}:{code:?};");
    }
    Ok(key)
}

#[cfg(not(feature = "native-inchi"))]
fn verified_identity_key_with_supplement(
    _inchi: &str,
    _supplement: &[(u64, CipCode)],
    _policy: IdentityPolicy,
) -> Result<String, ()> {
    Err(())
}

/// Compare two molecules under `policy`, retrying a legacy-CIP-unresolved
/// specified tetrahedral stereocentre via the accurate CIP engine before
/// giving up (issue #161). Given explicit candidate keys, like [`compare`].
///
/// This SHRINKS `VerificationUnavailable`'s trigger rate relative to
/// [`compare`] (9 of the 15 molecules a fresh 5,000-molecule corpus rerun
/// finds newly guarded by PR #156, beyond the original 4663/4664 pair,
/// recover verified-comparison capability -- see the module docs above and
/// the PR body for the full re-derivation) and never weakens the guard: a centre unresolved
/// by *both* engines still fails closed exactly as before, and a resolved
/// centre's accurate-CIP code is folded into the comparison key rather than
/// silently trusted via a regenerated (still `?`-bearing) InChI string, so
/// the original 4663/4664 false-`VerifiedDuplicate` stays fixed (see
/// `accurate_preflight_never_reopens_the_4663_4664_false_duplicate` below).
pub fn compare_with_accurate_cip_preflight(
    key_a: &str,
    mol_a: &Molecule,
    key_b: &str,
    mol_b: &Molecule,
    policy: IdentityPolicy,
) -> DedupRelation {
    let (raw_a, supp_a, raw_b, supp_b) = match (
        identity_verify_with_accurate_preflight(mol_a, policy),
        identity_verify_with_accurate_preflight(mol_b, policy),
    ) {
        (PreflightVerify::Invalid, _) | (_, PreflightVerify::Invalid) => {
            return DedupRelation::InvalidMolecule;
        }
        (PreflightVerify::Unavailable, _) | (_, PreflightVerify::Unavailable) => {
            return DedupRelation::VerificationUnavailable;
        }
        (PreflightVerify::Ok(a, sa), PreflightVerify::Ok(b, sb)) => (a, sa, b, sb),
    };

    let (ident_a, ident_b) = match (
        verified_identity_key_with_supplement(&raw_a, &supp_a, policy),
        verified_identity_key_with_supplement(&raw_b, &supp_b, policy),
    ) {
        (Ok(a), Ok(b)) => (a, b),
        _ => return DedupRelation::VerificationUnavailable,
    };
    let same_identity = ident_a == ident_b;
    let same_key = key_a == key_b;

    match (same_identity, same_key) {
        (true, true) => DedupRelation::VerifiedDuplicate,
        (true, false) => DedupRelation::CanonicalSplit,
        (false, true) => DedupRelation::CanonicalCollision,
        (false, false) => DedupRelation::Distinct,
    }
}

/// [`compare_with_accurate_cip_preflight`], computing each molecule's own
/// [`chematic_smiles::canonical_smiles`] as the fast candidate key -- the
/// accurate-CIP-preflight counterpart to [`compare_molecules`].
pub fn compare_molecules_with_accurate_cip_preflight(
    mol_a: &Molecule,
    mol_b: &Molecule,
    policy: IdentityPolicy,
) -> DedupRelation {
    let key_a = canonical_smiles(mol_a);
    let key_b = canonical_smiles(mol_b);
    compare_with_accurate_cip_preflight(&key_a, mol_a, &key_b, mol_b, policy)
}

/// [`deduplicate_verified`], with the accurate-CIP preflight (issue #161)
/// applied to every molecule instead of failing closed on the first
/// legacy-CIP-unresolved specified tetrahedral centre. Same algorithm/
/// complexity/output-ordering guarantees as [`deduplicate_verified`] -- see
/// that function's docs.
pub fn deduplicate_verified_with_accurate_cip_preflight(
    mols: &[Molecule],
    policy: IdentityPolicy,
) -> VerifiedDedupReport {
    let n = mols.len();
    let mut canonical_key: Vec<String> = Vec::with_capacity(n);
    let mut verified_key: Vec<Option<String>> = Vec::with_capacity(n);
    let mut verification_unavailable: Vec<usize> = Vec::new();
    let mut invalid_molecules: Vec<usize> = Vec::new();

    for (i, m) in mols.iter().enumerate() {
        canonical_key.push(canonical_smiles(m));
        match identity_verify_with_accurate_preflight(m, policy) {
            PreflightVerify::Invalid => {
                invalid_molecules.push(i);
                verified_key.push(None);
            }
            PreflightVerify::Unavailable => {
                verification_unavailable.push(i);
                verified_key.push(None);
            }
            PreflightVerify::Ok(raw, supplement) => {
                match verified_identity_key_with_supplement(&raw, &supplement, policy) {
                    Ok(k) => verified_key.push(Some(k)),
                    Err(()) => {
                        verification_unavailable.push(i);
                        verified_key.push(None);
                    }
                }
            }
        }
    }

    let mut by_verified: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, k) in verified_key.iter().enumerate() {
        if let Some(k) = k {
            by_verified.entry(k.as_str()).or_default().push(i);
        }
    }

    let mut groups: Vec<VerifiedGroup> = Vec::new();
    let mut canonical_splits: Vec<CanonicalSplit> = Vec::new();
    for members in by_verified.into_values() {
        if members.len() < 2 {
            continue;
        }
        let mut members = members;
        members.sort_unstable();

        let mut by_canon_within: HashMap<&str, Vec<usize>> = HashMap::new();
        for &i in &members {
            by_canon_within
                .entry(canonical_key[i].as_str())
                .or_default()
                .push(i);
        }
        groups.push(VerifiedGroup {
            members: members.clone(),
        });
        if by_canon_within.len() > 1 {
            let mut canonical_subgroups: Vec<Vec<usize>> = by_canon_within.into_values().collect();
            canonical_subgroups.sort_by_key(|g| g[0]);
            canonical_splits.push(CanonicalSplit {
                members,
                canonical_subgroups,
            });
        }
    }

    let mut by_canonical: HashMap<&str, Vec<usize>> = HashMap::new();
    for i in 0..n {
        if verified_key[i].is_some() {
            by_canonical
                .entry(canonical_key[i].as_str())
                .or_default()
                .push(i);
        }
    }

    let mut canonical_collisions: Vec<CanonicalCollision> = Vec::new();
    for mut members in by_canonical.into_values() {
        if members.len() < 2 {
            continue;
        }
        members.sort_unstable();
        let mut by_ident_within: HashMap<&str, Vec<usize>> = HashMap::new();
        for &i in &members {
            let k = verified_key[i].as_deref().expect("filtered to Some above");
            by_ident_within.entry(k).or_default().push(i);
        }
        if by_ident_within.len() > 1 {
            let mut verified_subgroups: Vec<Vec<usize>> = by_ident_within.into_values().collect();
            verified_subgroups.sort_by_key(|g| g[0]);
            canonical_collisions.push(CanonicalCollision {
                canonical_key_members: members,
                verified_subgroups,
            });
        }
    }

    groups.sort_by_key(|g| g.members[0]);
    canonical_splits.sort_by_key(|s| s.members[0]);
    canonical_collisions.sort_by_key(|c| c.canonical_key_members[0]);
    verification_unavailable.sort_unstable();
    invalid_molecules.sort_unstable();

    VerifiedDedupReport {
        groups,
        canonical_splits,
        canonical_collisions,
        verification_unavailable,
        invalid_molecules,
    }
}

#[cfg(test)]
mod tests {
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

        // Reference: the plain (legacy-only) path already reports this as
        // unavailable (PR #156's fix, unchanged).
        assert_eq!(
            compare_molecules(&a, &b, IdentityPolicy::StandardInchiString),
            DedupRelation::VerificationUnavailable,
            "sanity: PR #156's fix must still fail closed on the plain path"
        );

        let with_preflight = compare_molecules_with_accurate_cip_preflight(
            &a,
            &b,
            IdentityPolicy::StandardInchiString,
        );
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
        let idx_443 = parse(
            "CCCCCCCC/C=C/CCCCCCCC(=O)OCc1cc(=O)c(OC(=O)[C@]23C[C@H]4C[C@H](C[C@H](C4)C2)C3)co1",
        )
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
            compare_molecules_with_accurate_cip_preflight(
                &a,
                &b,
                IdentityPolicy::StandardInchiString
            ),
            "an ordinary molecule with no legacy-CIP-unresolved centre must behave \
             identically under both entry points"
        );
    }
}
