//! RDKit-equivalent Morgan environment suppression — additive, experimental.
//!
//! chematic's original Morgan expansion ([`crate::ecfp::ecfp`] and friends)
//! emits every atom at every radius unconditionally. RDKit's default
//! (`includeRedundantEnvironments=False`) instead suppresses an atom once
//! its cumulative bond-index-set duplicates one already emitted by another
//! atom at the same or an earlier radius. This module reproduces that
//! algorithm, verified directly against RDKit's real source at a pinned
//! commit (not a mutable `master` reference):
//! `Code/GraphMol/Fingerprints/MorganGenerator.cpp`,
//! `MorganEnvGenerator<OutputType>::getEnvironments`, commit
//! [`0062b670640352ab63d6256be608615e87e1af53`](https://github.com/rdkit/rdkit/blob/0062b670640352ab63d6256be608615e87e1af53/Code/GraphMol/Fingerprints/MorganGenerator.cpp).
//! That commit is **not** an ancestor of release tag `Release_2026_03_4`
//! (which resolves to `8afba32ec539dcb2369bc84549d802aca3f7eb39`,
//! independently verified via the GitHub tags API during Morgan M4-A0) —
//! diverged history, not a simple predecessor. Independently diffed during
//! M4-A0's provenance audit: the file as a whole differs by one unrelated
//! line in a different function (`updateAdditionalOutput`'s `bitId`
//! parameter type), but `getEnvironments` itself — the function this module
//! actually ports — is byte-identical between the two commits, so the
//! algorithm above is unaffected. See [`crate::rdkit_morgan_hash`]'s module
//! docs and `THIRD_PARTY_NOTICES.md` for the full three-citation picture
//! (this module, `chematic-perception`'s aromaticity port, and
//! `rdkit_morgan_hash.rs` each cite a different commit under the same
//! `Release_2026_03_4` label).
//!
//! Per molecule: radius 0 is emitted unconditionally for every atom (no
//! suppression concept at round 0). For each subsequent layer (0-indexed,
//! emitting radius `layer+1`):
//! 1. Every atom not yet dead computes a cumulative `BondSet` (its own
//!    incident bonds unioned with each neighbor's *previous*-round
//!    `BondSet`) and a fresh invariant via [`expand_atom_id`]. Atoms of
//!    degree 0 die immediately (before computing anything) — this happens
//!    unconditionally, in both emission modes.
//! 2. Atoms are grouped by `BondSet`. Under
//!    [`EnvironmentEmissionMode::IncludeRdkitRedundant`] every group member
//!    emits, no death from this step. Under
//!    [`EnvironmentEmissionMode::SuppressRdkitRedundant`]: if the `BondSet`
//!    was already seen in an earlier round, the *entire* group dies
//!    (nobody emits, including what would otherwise be the winner); if not,
//!    the member with the smallest `(invariant, atom_idx)` emits, its
//!    `BondSet` is recorded as seen, and every other group member dies.
//! 3. Freshly computed invariants are written for every atom that reached
//!    step 1 this layer, *before* the group-death decision — so an atom
//!    that dies this round still contributes a real, non-zero invariant to
//!    its live neighbors one round later, and only reads as zero from two
//!    rounds after its death onward. This one-round grace period falls out
//!    automatically from writing the fresh-invariant buffer unconditionally
//!    and swapping it in after the death walk, matching RDKit exactly.

use chematic_core::{AtomIdx, Molecule};
use rustc_hash::{FxHashMap, FxHashSet};

use crate::bitvec::BitVec2048;
use crate::ecfp::{
    EcfpConfig, EcfpInvariantMode, MAX_ECFP_RADIUS, expand_atom_id, initial_atom_id, record_bit,
};

/// Whether repeated (redundant) atom environments are suppressed. See the
/// module docs for the exact algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EnvironmentEmissionMode {
    /// RDKit's `includeRedundantEnvironments=true` lifecycle: every atom's
    /// bond-environment is emitted at every radius it computes one for, with
    /// no death from bondset collision. Degree-zero atoms still stop after
    /// radius 0 — that rule is unconditional on this flag, matching RDKit —
    /// so this is **not** the same as chematic's existing legacy
    /// `ecfp()`/`ecfp4()`/`ecfp6()`/`morgan_fp_counts()` behavior, which
    /// re-hashes and emits a degree-zero atom at every radius too.
    #[allow(dead_code)] // exercised by tests; no production caller needs it yet
    IncludeRdkitRedundant,
    /// RDKit's default (`includeRedundantEnvironments=false`): once a
    /// bond-index-set has been emitted by any atom at any radius, no atom
    /// emits that same bond-index-set again.
    SuppressRdkitRedundant,
}

/// One bit per bond index in the molecule. Equality/hash is plain
/// word-vector equality/hash — sufficient because every `BondSet` compared
/// against another within one call is sized to the same molecule's
/// `bond_count()`. `pub(crate)`: pure bitset infrastructure with no
/// hash-specific behavior, reused as-is by [`crate::rdkit_morgan_hash`]'s
/// independent (RDKit-exact-hash) suppression pass — only the invariant/hash
/// computation differs between the two modules, not this plumbing.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct BondSet(Vec<u64>);

impl BondSet {
    pub(crate) fn empty(bond_count: usize) -> Self {
        BondSet(vec![0u64; bond_count.div_ceil(64)])
    }

    pub(crate) fn set(&mut self, bond: u32) {
        let word = bond as usize / 64;
        let bit = bond as usize % 64;
        self.0[word] |= 1u64 << bit;
    }

    pub(crate) fn union_with(&mut self, other: &BondSet) {
        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            *a |= b;
        }
    }
}

/// Core Morgan expansion, returning every *emitted* environment as
/// `(atom_idx, radius, raw_id)`. Radius 0 always emits every atom exactly
/// once; suppression (if enabled) only ever applies to radius >= 1.
fn ecfp_environments_emitted(
    mol: &Molecule,
    config: &EcfpConfig,
    invariant_mode: EcfpInvariantMode,
    emission_mode: EnvironmentEmissionMode,
) -> Vec<(u32, u32, u64)> {
    let n = mol.atom_count();
    let mut emitted = Vec::new();
    if n == 0 {
        return emitted;
    }

    let radius = config.radius.min(MAX_ECFP_RADIUS);
    let ring_set = chematic_perception::find_sssr(mol);
    let bond_count = mol.bond_count();

    let mut current_invariants: Vec<u64> = Vec::with_capacity(n);
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let id = initial_atom_id(mol, idx, &ring_set, config.use_chirality, invariant_mode);
        emitted.push((i as u32, 0, id));
        current_invariants.push(id);
    }

    let mut dead = vec![false; n];
    let mut atom_neighborhoods: Vec<BondSet> = (0..n).map(|_| BondSet::empty(bond_count)).collect();
    let mut seen: FxHashSet<BondSet> = FxHashSet::default();

    for layer in 0..radius {
        let mut next_invariants = vec![0u64; n];
        let mut round_atom_neighborhoods = atom_neighborhoods.clone();
        let mut groups: FxHashMap<BondSet, Vec<(u64, u32)>> = FxHashMap::default();

        for i in 0..n {
            if dead[i] {
                continue;
            }
            let idx = AtomIdx(i as u32);
            let neighbors: Vec<(AtomIdx, chematic_core::BondIdx)> = mol.neighbors(idx).collect();
            if neighbors.is_empty() {
                dead[i] = true;
                continue;
            }
            let mut bond_env = BondSet::empty(bond_count);
            for (nb_idx, bond_idx) in &neighbors {
                bond_env.set(bond_idx.0);
                bond_env.union_with(&atom_neighborhoods[nb_idx.0 as usize]);
            }
            let invariant = expand_atom_id(mol, i, layer + 1, &current_invariants);
            next_invariants[i] = invariant;
            round_atom_neighborhoods[i] = bond_env.clone();
            groups
                .entry(bond_env)
                .or_default()
                .push((invariant, i as u32));
        }

        for (bond_env, mut members) in groups {
            match emission_mode {
                EnvironmentEmissionMode::IncludeRdkitRedundant => {
                    for &(invariant, atom_idx) in &members {
                        emitted.push((atom_idx, layer + 1, invariant));
                    }
                }
                EnvironmentEmissionMode::SuppressRdkitRedundant => {
                    if seen.contains(&bond_env) {
                        for &(_, atom_idx) in &members {
                            dead[atom_idx as usize] = true;
                        }
                    } else {
                        members.sort_unstable();
                        let (winner_invariant, winner_idx) = members[0];
                        emitted.push((winner_idx, layer + 1, winner_invariant));
                        seen.insert(bond_env);
                        for &(_, atom_idx) in &members[1..] {
                            dead[atom_idx as usize] = true;
                        }
                    }
                }
            }
        }

        current_invariants = next_invariants;
        atom_neighborhoods = round_atom_neighborhoods;
    }

    emitted
}

/// Like [`crate::ecfp::ecfp_with_bitinfo_and_mode`], but applying
/// `emission_mode` on top — the suppression-aware counterpart used by the
/// `_rdkit_environment_experimental` public functions in
/// [`crate::ecfp`].
pub(crate) fn ecfp_environments(
    mol: &Molecule,
    config: &EcfpConfig,
    invariant_mode: EcfpInvariantMode,
    emission_mode: EnvironmentEmissionMode,
) -> (BitVec2048, FxHashMap<usize, Vec<(u32, u32)>>) {
    let mut fp = BitVec2048::new();
    let mut info: FxHashMap<usize, Vec<(u32, u32)>> = FxHashMap::default();
    for (atom_idx, radius, raw_id) in
        ecfp_environments_emitted(mol, config, invariant_mode, emission_mode)
    {
        record_bit(
            &mut fp,
            &mut info,
            raw_id,
            atom_idx,
            radius,
            config.nbits,
            config.use_double_fold,
        );
    }
    (fp, info)
}

/// Diagnostic-only: raw (unfolded) *suppressed* emitted environments as
/// `(atom_idx, radius, raw_environment_id)`, always using
/// `EcfpInvariantMode::RdkitMorgan` + `SuppressRdkitRedundant` — the real
/// production suppression path, not a general-purpose entry point. Exists
/// so validation tooling can measure raw-identifier *count multiplicity*
/// (sparse-count shape) against an external oracle without going through
/// 2048-bit folding, which can collide distinct raw ids into the same bit
/// and muddy a count comparison. Re-exported under the `diagnostics`
/// feature — deliberately not part of the public API (mirrors
/// [`ecfp_environments`]'s own narrow scope; `EnvironmentEmissionMode`
/// itself stays crate-private).
#[allow(dead_code)] // only reachable via the `diagnostics` feature (matches ecfp_diagnostics.rs's own convention)
pub fn suppressed_environments_diagnostic(
    mol: &Molecule,
    config: &EcfpConfig,
) -> Vec<(u32, u32, u64)> {
    ecfp_environments_emitted(
        mol,
        config,
        EcfpInvariantMode::RdkitMorgan,
        EnvironmentEmissionMode::SuppressRdkitRedundant,
    )
}

/// Raw (unfolded) emitted environments as `(atom_idx, radius, raw_id)` —
/// test-only. Production/validation code gets the same `(atom, radius)`
/// pairs by flattening [`ecfp_environments`]'s `info` map; this exists
/// because tests want the *value* too, without going through 2048-bit
/// folding (which could — for larger radii or nbits — collide two distinct
/// raw ids into the same bit, muddying a hand-verified expectation). See
/// [`ecfp_environments_emitted`].
#[cfg(test)]
pub(crate) fn ecfp_environments_sparse(
    mol: &Molecule,
    config: &EcfpConfig,
    invariant_mode: EcfpInvariantMode,
    emission_mode: EnvironmentEmissionMode,
) -> Vec<(u32, u32, u64)> {
    ecfp_environments_emitted(mol, config, invariant_mode, emission_mode)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn emitted_pairs(
        mol: &Molecule,
        radius: u32,
        mode: EnvironmentEmissionMode,
    ) -> Vec<(u32, u32)> {
        let config = EcfpConfig {
            radius,
            ..EcfpConfig::default()
        };
        let mut pairs: Vec<(u32, u32)> =
            ecfp_environments_sparse(mol, &config, EcfpInvariantMode::RdkitMorgan, mode)
                .into_iter()
                .map(|(atom, radius, _)| (atom, radius))
                .collect();
        pairs.sort_unstable();
        pairs
    }

    #[test]
    fn bondset_union_order_independent() {
        let mut a = BondSet::empty(130);
        a.set(3);
        a.set(70);
        let mut b = BondSet::empty(130);
        b.set(70);
        b.set(3);
        assert_eq!(a, b);
        let mut h = std::collections::hash_map::DefaultHasher::new();
        use std::hash::{Hash, Hasher};
        a.hash(&mut h);
        let ha = h.finish();
        let mut h2 = std::collections::hash_map::DefaultHasher::new();
        b.hash(&mut h2);
        assert_eq!(ha, h2.finish());
    }

    #[test]
    fn ethane_suppresses_duplicate_atom_at_radius_1_then_saturates_at_radius_2() {
        // CC: atom 0 and atom 1 are symmetric, share the single (only) bond,
        // so their radius-1 bond-environment is identical — exactly one
        // emits under suppression. Ethane has just 1 bond total, so the
        // survivor's radius-2 bond-environment ({bond0} again, nothing left
        // to discover) exactly repeats the radius-1 `seen` entry: radius 2
        // must have ZERO emissions under suppression, not 1 — the whole
        // molecule is already fully discovered by radius 1.
        let mol = parse("CC").unwrap();
        let suppressed = emitted_pairs(&mol, 2, EnvironmentEmissionMode::SuppressRdkitRedundant);
        let radius1_count = suppressed.iter().filter(|&&(_, r)| r == 1).count();
        assert_eq!(
            radius1_count, 1,
            "radius 1: expected exactly 1 emission, got {radius1_count}"
        );
        let radius2_count = suppressed.iter().filter(|&&(_, r)| r == 2).count();
        assert_eq!(
            radius2_count, 0,
            "radius 2: expected zero emissions (ethane fully discovered by radius 1), \
             got {radius2_count}"
        );
        let all = emitted_pairs(&mol, 2, EnvironmentEmissionMode::IncludeRdkitRedundant);
        for r in [1u32, 2] {
            let count = all.iter().filter(|&&(_, radius)| radius == r).count();
            assert_eq!(
                count, 2,
                "IncludeRdkitRedundant radius {r}: expected 2 emissions, got {count}"
            );
        }
    }

    #[test]
    fn neopentane_center_saturates_before_radius_2() {
        // CC(C)(C)C: the central atom's bond-environment stops growing after
        // radius 1 (all neighbors are leaf methyls) — under suppression it
        // must not re-emit at radius 2.
        let mol = parse("CC(C)(C)C").unwrap();
        let suppressed = emitted_pairs(&mol, 2, EnvironmentEmissionMode::SuppressRdkitRedundant);
        let center = 1u32;
        assert!(
            !suppressed.contains(&(center, 2)),
            "central atom must not re-emit at radius 2 once saturated: {suppressed:?}"
        );
        assert!(
            suppressed.contains(&(center, 1)),
            "central atom must emit at radius 1: {suppressed:?}"
        );
    }

    #[test]
    fn degree_zero_atom_never_emits_past_radius_0_in_either_mode() {
        let mol = parse("[Cl-]").unwrap();
        for mode in [
            EnvironmentEmissionMode::IncludeRdkitRedundant,
            EnvironmentEmissionMode::SuppressRdkitRedundant,
        ] {
            let pairs = emitted_pairs(&mol, 2, mode);
            assert_eq!(pairs, vec![(0, 0)], "mode {mode:?}: {pairs:?}");
        }
    }

    #[test]
    fn neopentane_all_atoms_suppressed_at_radius_2_via_cross_round_seen() {
        // Deeper than "center doesn't re-emit": at radius 2 every methyl's
        // bond-environment becomes IDENTICAL to the center's radius-1
        // bond-environment (each methyl's own single bond is already inside
        // the center's 4-bond set, and unioning with the center's *previous*
        // round bondset contributes nothing new) — that bondset was already
        // recorded in `seen` when the center emitted it at radius 1. So all
        // 5 atoms collide with a *cross-round* `seen` entry simultaneously:
        // zero emissions at radius 2, not just the center's own.
        let mol = parse("CC(C)(C)C").unwrap();
        let suppressed = emitted_pairs(&mol, 2, EnvironmentEmissionMode::SuppressRdkitRedundant);
        let radius2_count = suppressed.iter().filter(|&&(_, r)| r == 2).count();
        assert_eq!(
            radius2_count, 0,
            "neopentane: expected zero radius-2 emissions (all suppressed via cross-round \
             `seen`), got {suppressed:?}"
        );
    }

    #[test]
    fn dying_atom_invariant_has_one_round_grace_period() {
        // A methyl (atom 0) on a wide branch point (atom 1, also bonded to
        // two long 8-carbon chains) saturates fast: it emits radius 0 and 1,
        // then its radius-2 bond-environment collides with an already-`seen`
        // set and it dies *at* the round computing radius 2 — having already
        // written a real (non-zero) radius-2 invariant into `next_invariants`
        // before that group-death decision. Atom 1 keeps growing via its
        // other two branches for many more rounds, so its radius-3
        // computation (one round later) is directly observable, and must
        // read atom 0's real radius-2 value (grace period), not 0.
        //
        // This white-box test computes the expected radius-3 invariant using
        // the same production `expand_atom_id` the implementation itself
        // calls, fed with every neighbor's real radius-2 invariant
        // (collected from `IncludeRdkitRedundant` mode, where invariant computation is
        // identical to `SuppressRdkitRedundant` — modes only differ in who
        // survives to *emit*, not in the invariant formula, so atom 0's real
        // radius-2 value is directly readable there even though it never
        // emits it under suppression) — proving the implementation actually
        // used the live value, not a zeroed one.
        let mol = parse("CC(CCCCCCCC)CCCCCCCC").unwrap();
        let center = 1u32;
        let dying_neighbor = 0u32;
        let config = EcfpConfig {
            radius: 3,
            ..EcfpConfig::default()
        };

        let all_env = ecfp_environments_sparse(
            &mol,
            &config,
            EcfpInvariantMode::RdkitMorgan,
            EnvironmentEmissionMode::IncludeRdkitRedundant,
        );
        let mut ids_at_radius2 = vec![0u64; mol.atom_count()];
        for &(atom, radius, raw_id) in &all_env {
            if radius == 2 {
                ids_at_radius2[atom as usize] = raw_id;
            }
        }
        assert!(
            ids_at_radius2.iter().all(|&v| v != 0),
            "every atom must have a real (nonzero) radius-2 invariant in IncludeRdkitRedundant mode"
        );

        let suppressed_env = ecfp_environments_sparse(
            &mol,
            &config,
            EcfpInvariantMode::RdkitMorgan,
            EnvironmentEmissionMode::SuppressRdkitRedundant,
        );
        assert!(
            suppressed_env
                .iter()
                .any(|&(atom, r, _)| atom == dying_neighbor && r == 1),
            "atom 0 must emit at radius 1 (last real emission before dying at radius 2): \
             {suppressed_env:?}"
        );
        assert!(
            !suppressed_env
                .iter()
                .any(|&(atom, r, _)| atom == dying_neighbor && r == 2),
            "atom 0 must NOT emit at radius 2 (dies from a cross-round `seen` collision): \
             {suppressed_env:?}"
        );

        let actual_radius3 = suppressed_env
            .iter()
            .find(|&&(atom, r, _)| atom == center && r == 3)
            .map(|&(_, _, raw_id)| raw_id)
            .expect("the branch center must survive to emit at radius 3");

        let expected_radius3 = expand_atom_id(&mol, center as usize, 3, &ids_at_radius2);
        assert_eq!(
            actual_radius3, expected_radius3,
            "radius-3 invariant must be computed from the real (grace-period) radius-2 \
             invariant of the neighbor that died at radius 2, not a zeroed one"
        );
    }

    #[test]
    fn suppression_never_over_emits_relative_to_include_redundant() {
        // The suppressed emission set must always be a subset (by
        // (atom,radius) pair) of the IncludeRdkitRedundant set — suppression
        // only ever removes emissions, never adds new ones.
        for smi in ["CC(C)(C)C", "c1ccccc1", "CCO.c1ccccc1", "N[C@@H](C)C(=O)O"] {
            let mol = parse(smi).unwrap();
            let all: std::collections::HashSet<_> =
                emitted_pairs(&mol, 2, EnvironmentEmissionMode::IncludeRdkitRedundant)
                    .into_iter()
                    .collect();
            let suppressed: std::collections::HashSet<_> =
                emitted_pairs(&mol, 2, EnvironmentEmissionMode::SuppressRdkitRedundant)
                    .into_iter()
                    .collect();
            assert!(
                suppressed.is_subset(&all),
                "{smi}: suppressed set must be a subset of the IncludeRdkitRedundant set"
            );
        }
    }
}
