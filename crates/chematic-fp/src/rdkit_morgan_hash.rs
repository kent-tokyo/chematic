//! RDKit-exact Morgan raw-identifier reference path — Milestone M4-A0,
//! diagnostic only.
//!
//! [`crate::ecfp::EcfpInvariantMode::RdkitMorgan`] and
//! [`crate::morgan_environment`] (Phase B, PR #123) already reproduce
//! RDKit's atom-invariant *partition* and emission *lifecycle* — which atoms
//! are chemically equivalent, and which one wins/dies at each radius — but
//! neither reproduces RDKit's actual 32-bit hash *values*: both still hash
//! through chematic's own FNV-1a over a byte layout of chematic's own
//! design. This module is a from-scratch, source-verified port of RDKit's
//! real hashing machinery, so raw identifiers can be compared *numerically*
//! against a real RDKit oracle, not just partition-for-partition.
//!
//! Ported from RDKit release `Release_2026_03_4`, commit
//! [`8afba32ec539dcb2369bc84549d802aca3f7eb39`](https://github.com/rdkit/rdkit/blob/8afba32ec539dcb2369bc84549d802aca3f7eb39/Code/GraphMol/Fingerprints/MorganGenerator.cpp)
//! — independently resolved this session via the GitHub tags API
//! (`GET /repos/rdkit/rdkit/git/refs/tags/Release_2026_03_4`), not reused
//! from a prior port's citation. See `THIRD_PARTY_NOTICES.md` at the repo
//! root for the required BSD 3-Clause attribution and license text. Note:
//! two *other* SHAs already exist in this project's history under the same
//! "`Release_2026_03_4`" label —
//! `crates/chematic-perception/src/rdkit_parity.rs`/`THIRD_PARTY_NOTICES.md`
//! cite `e89c9f656a694fab4105139844cba88d2e013354`, and
//! [`crate::morgan_environment`]'s own doc comment cites
//! `0062b670640352ab63d6256be608615e87e1af53`. Both are ancestors of the
//! true tag resolution, not the tag itself; not reconciled as part of this
//! milestone (flagged, not fixed — see the M4-A0 report).
//!
//! Ported functions/formulas, each cited at its own definition below:
//! - `gboost::hash_combine` (`Code/RDGeneral/hash/hash.hpp`)
//! - `hash_value(std::pair<A,B>)`, `hash_range`/`hash_value(vector<T>)` (same file)
//! - `MorganAtomInvGenerator::getConnectivityInvariants` (`FingerprintUtil.cpp`)
//! - `MorganBondInvGenerator::getBondInvariants` (`FingerprintUtil.cpp`) — only
//!   the `useBondTypes=true, useChirality=false` branch (this diagnostic's
//!   pinned oracle config; the chirality-perturbed branch is out of scope)
//! - `MorganEnvGenerator<OutputType>::getEnvironments`'s per-round hash loop
//!   (`MorganGenerator.cpp`) — the atom-CIP chirality re-fold is likewise out
//!   of scope (`includeChirality=false` pinned throughout this workstream)
//!
//! Deliberately kept structurally separate from both the FNV-1a production
//! path (`ecfp.rs`) and the Phase B suppression path
//! (`morgan_environment.rs`), even though it duplicates some control-flow
//! shape (round loop, degree-0 death, bondset-keyed suppression) — reusing
//! only genuinely hash-independent infrastructure
//! ([`crate::morgan_environment::BondSet`], and `ecfp.rs`'s
//! `rdkit_total_degree`/`rdkit_total_h_count`/`rdkit_isotope_delta` value
//! computations) so this module was never pre-emptively entangled with
//! either existing path.
//!
//! **Numeric exactness proven (M4-A0), now promoted to a real production
//! API:** [`crate::rdkit_morgan_ecfp4`]'s `rdkit_morgan_ecfp4_experimental`
//! reuses this module's [`expand_one_pass`] and [`checked_bond_invariant`]
//! directly (`pub(crate)`) rather than re-deriving the same hash logic a
//! second time. The `rdkit_morgan_raw_trace`/`RdkitMorganRawTraceEntry`
//! diagnostic surface stays diagnostics-feature-gated and computes *both*
//! RDKit lifecycles at once for comparison purposes — the production path
//! only ever runs the single `suppress=true` ("default") lifecycle RDKit
//! itself uses.
//!
//! **Caller contract:** [`rdkit_morgan_raw_trace`] expects `mol` to already
//! have aromaticity perception applied
//! ([`chematic_perception::apply_aromaticity`]) when Kekule-spelled aromatic
//! input is possible. Hand-verified this session: a Kekule-spelled pyridine
//! (`C1=CC=NC=C1`) fed in *without* that step diverges from RDKit's real
//! identifiers at every radius >= 1 (aromatic ring bonds hash as literal
//! Single/Double instead of [`bond_invariant`]'s Aromatic=12), while the
//! *same* molecule with `apply_aromaticity` applied first matches both
//! RDKit's real oracle and chematic's own aromatic-spelled `c1ccncc1` input
//! bit-for-bit at every radius — RDKit's own `Chem.MolFromSmiles` always
//! sanitizes (perceives aromaticity) before hashing, so this is matching
//! RDKit's real precondition, not working around a chematic quirk.

#![allow(dead_code)] // only reachable via the `diagnostics` feature + this module's own tests

use chematic_core::{AtomIdx, BondIdx, BondOrder, CipCode, Molecule};
use chematic_perception::RingSet;
use rustc_hash::{FxHashMap, FxHashSet};
use smallvec::SmallVec;

use crate::ecfp::{MAX_ECFP_RADIUS, rdkit_isotope_delta, rdkit_total_degree, rdkit_total_h_count};
use crate::morgan_environment::BondSet;

/// `gboost::hash_combine`: `seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) +
/// (seed>>2)`. Source: `Code/RDGeneral/hash/hash.hpp` (`hash_combine`
/// template, golden-ratio Boost-style mixing). `hash_result_t` is a fixed
/// `std::uint32_t` (`Code/RDGeneral/hash/hash_fwd.hpp`) regardless of the
/// fingerprint's `OutputType` template parameter — all mixing happens in
/// 32-bit arithmetic; only the final identifier is later widened (never
/// truncated) to `OutputType` at construction. All arithmetic here wraps
/// (matching C++ `uint32_t` overflow, which is well-defined modular
/// arithmetic, not UB) — `wrapping_add` is required, not cosmetic: plain `+`
/// panics in a Rust debug build on the frequent overflow this formula
/// produces.
pub(crate) fn hash_combine(seed: u32, value: u32) -> u32 {
    seed ^ value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add(seed << 6)
        .wrapping_add(seed >> 2)
}

/// `hash_value(std::pair<A, B>)`: a fresh sub-hash, seed 0, first element
/// combined before the second. Source: `hash.hpp` (`hash_value` overload for
/// `std::pair`). This nested sub-hash is itself what gets fed into the outer
/// per-round accumulator via one more `hash_combine` call — not a flat
/// two-call fold into the same seed as the caller.
fn hash_pair(a: u32, b: u32) -> u32 {
    let seed = hash_combine(0, a);
    hash_combine(seed, b)
}

/// `hash_range` / `hash_value(std::vector<T>)`: sequential left fold, seed
/// starting at 0. Source: `hash.hpp`. Used for the radius-0 connectivity
/// invariant, whose component count is *not* fixed (see
/// [`connectivity_invariant`]) — a shorter or longer vector changes the
/// number of `hash_combine` calls, not just one component's value.
pub(crate) fn hash_vec(values: &[u32]) -> u32 {
    values.iter().fold(0u32, |seed, &v| hash_combine(seed, v))
}

/// `MorganAtomInvGenerator::getConnectivityInvariants`, radius-0 invariant.
/// Source: `Code/GraphMol/Fingerprints/FingerprintUtil.cpp`. Builds
/// `[atomicNum, totalDegree, totalNumHs, formalCharge, deltaMass]` as
/// `uint32_t` (charge and mass-delta are two's-complement-wrapped, not
/// clamped — `(v as i32) as u32`, matching `int` → `uint32_t` implicit
/// conversion) and appends a **literal trailing `1` only when the atom is in
/// a ring** — an out-of-ring atom hashes a 5-element vector, an in-ring atom
/// a 6-element one. This is the one component chematic's existing
/// `EcfpInvariantMode::RdkitMorgan` byte layout gets structurally wrong for
/// bit-exactness purposes: that path always writes a ring-membership byte
/// (0 or 1), which is a fixed-length encoding fundamentally incompatible
/// with RDKit's variable-length-vector-then-hash approach, even though both
/// encode the same *information*. `total_degree`/`total_h_count`/
/// `isotope_delta` values are reused verbatim from `ecfp.rs` (already
/// verified RDKit-equivalent by partition, per `EcfpInvariantMode::RdkitMorgan`'s
/// own doc comment) — only the assembly and hash function differ here.
fn connectivity_invariant(mol: &Molecule, idx: AtomIdx, ring_set: &RingSet) -> u32 {
    let atom = mol.atom(idx);
    let mut components: SmallVec<[u32; 6]> = SmallVec::new();
    components.push(atom.element.atomic_number() as u32);
    components.push(rdkit_total_degree(mol, idx) as u32);
    components.push(rdkit_total_h_count(mol, idx) as u32);
    components.push((atom.charge as i32) as u32);
    components.push(rdkit_isotope_delta(mol, idx) as u32);
    if ring_set.contains_atom(idx) {
        components.push(1);
    }
    hash_vec(&components)
}

/// `MorganBondInvGenerator::getBondInvariants`, `useBondTypes=true,
/// useChirality=false` branch only (this diagnostic's pinned oracle config —
/// see the module docs). Source: `FingerprintUtil.cpp`. Value is
/// `static_cast<int32_t>(bond->getBondType())` reinterpreted as `uint32_t`,
/// i.e. RDKit's own `Bond::BondType` C++ enum ordinal (`Code/GraphMol/Bond.h`,
/// commit `8afba32e...`, verbatim: `UNSPECIFIED=0, SINGLE, DOUBLE, TRIPLE,
/// QUADRUPLE, QUINTUPLE, HEXTUPLE, ONEANDAHALF, TWOANDAHALF, THREEANDAHALF,
/// FOURANDAHALF, FIVEANDAHALF, AROMATIC, IONIC, HYDROGEN, THREECENTER,
/// DATIVEONE, DATIVE, DATIVEL, DATIVER, OTHER, ZERO` — only `UNSPECIFIED`
/// has an explicit value, the rest auto-increment): SINGLE=1, DOUBLE=2,
/// TRIPLE=3, QUADRUPLE=4, AROMATIC=**12**, DATIVE=17, ZERO=21 — every value
/// used below independently source-quoted this session, none guessed.
/// **Not** chematic's own `bond_type_int` (`ecfp.rs`), whose AROMATIC=4
/// collides with RDKit's own QUADRUPLE value — reusing it here would
/// silently misclassify every aromatic bond.
///
/// Chematic's `BondOrder::Query*` variants are SMARTS-pattern-only bond
/// kinds with no RDKit `BondType` counterpart at all (not merely
/// unverified — RDKit's `BondType` enum has no query concept; SMARTS bond
/// queries are a distinct RDKit type) and cannot appear in a SMILES-parsed
/// molecule (this diagnostic's entire corpus) — mapped to `None`. The
/// production `rdkit_morgan_ecfp4` path (which reuses this function) turns
/// `None` into an explicit `RdkitMorganError::UnsupportedBondOrder`, never
/// an implicit/guessed mapping.
pub(crate) fn checked_bond_invariant(order: BondOrder) -> Option<u32> {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => Some(1),
        BondOrder::Double => Some(2),
        BondOrder::Triple => Some(3),
        BondOrder::Quadruple => Some(4),
        BondOrder::Aromatic => Some(12),
        BondOrder::Dative => Some(17),
        BondOrder::Zero => Some(21),
        BondOrder::QueryAny
        | BondOrder::QuerySingleOrDouble
        | BondOrder::QuerySingleOrAromatic
        | BondOrder::QueryDoubleOrAromatic => None,
    }
}

/// Diagnostic-only wrapper: unsupported bond orders map to `u32::MAX`, a
/// deliberately out-of-band placeholder, so the diagnostic can still trace
/// every bond rather than abort — no production caller uses this; the
/// production path uses [`checked_bond_invariant`] and surfaces an explicit
/// `Err` instead.
fn bond_invariant(order: BondOrder) -> u32 {
    checked_bond_invariant(order).unwrap_or(u32::MAX)
}

/// One computed `(atom, radius)` combination from [`rdkit_morgan_raw_trace`].
///
/// `raw_identifier_full`/`raw_identifier_default` are tracked as
/// **independent** fields, not one shared value plus two emitted booleans,
/// because they genuinely can differ: RDKit's `includeRedundantEnvironments`
/// isn't a filter applied after one shared computation, it changes which
/// atoms are still *alive* to contribute their own fresh invariant at each
/// round (a suppressed atom's neighbors freeze/zero its contribution one
/// and two rounds after its death respectively -- see the module docs'
/// grace period note), so an atom that survives under `full` but was
/// suppressed under `default` can legitimately compute a **different** raw
/// identifier at the same `(atom, radius)` between the two lifecycles. An
/// earlier, buggier version of this function shared one `dead` array across
/// both lifecycles; caught by the full-corpus comparator (not by hand
/// verification, which only checked value equality on entries present in
/// both, not entry *count*) -- ethane's atom 1 silently lost its radius-2
/// `full` entry because it had already been suppressed under `default`
/// semantics one round earlier. Fixed by running two fully independent
/// passes ([`expand_one_pass`]) and merging by key.
///
/// `None` in either field means that lifecycle never computed this atom at
/// this radius: past radius 0, this only happens for a degree-0 atom (in
/// `raw_identifier_full`) or a suppressed non-winner (in
/// `raw_identifier_default` -- RDKit's own `default`-mode `bitInfoMap`
/// likewise never records losing candidates, only the winner, so there is
/// no oracle ground truth for a "losing" identifier to track here either).
#[derive(Debug, Clone, Copy)]
pub struct RdkitMorganRawTraceEntry {
    pub atom_idx: u32,
    pub radius: u32,
    pub raw_identifier_full: Option<u32>,
    pub raw_identifier_default: Option<u32>,
}

/// One independent lifecycle simulation (either RDKit's
/// `includeRedundantEnvironments=true` "full", or its default suppressed
/// mode), returning `(atom_idx, radius) -> raw_identifier` only for
/// combinations this lifecycle actually computes/emits. Structurally
/// mirrors [`crate::morgan_environment`]'s proven-correct round loop
/// (degree-0 death, cumulative per-atom `BondSet`, cross-round `seen` set,
/// one-round invariant grace period) -- only the invariant/hash computation
/// differs; see the module docs for why that duplication is deliberate.
pub(crate) fn expand_one_pass(
    mol: &Molecule,
    ring_set: &RingSet,
    bond_invariants: &[u32],
    max_radius: u32,
    suppress: bool,
) -> FxHashMap<(u32, u32), u32> {
    expand_one_pass_with_chirality(mol, ring_set, bond_invariants, max_radius, suppress, None)
}

/// Morgan expansion with RDKit's opt-in tetrahedral chirality re-fold.
///
/// RDKit adds the chiral contribution once an atom's current neighborhood
/// proves it has four distinct single-bond ligand invariants. Once proven,
/// the contribution is retained in subsequent rounds through the ordinary
/// invariant hash and is added again by RDKit's per-round loop.
pub(crate) fn expand_one_pass_with_chirality(
    mol: &Molecule,
    ring_set: &RingSet,
    bond_invariants: &[u32],
    max_radius: u32,
    suppress: bool,
    cip_codes: Option<&FxHashMap<AtomIdx, CipCode>>,
) -> FxHashMap<(u32, u32), u32> {
    let n = mol.atom_count();
    let mut out: FxHashMap<(u32, u32), u32> = FxHashMap::default();
    if n == 0 {
        return out;
    }
    let bond_count = mol.bond_count();

    let mut current_invariants: Vec<u32> = Vec::with_capacity(n);
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let id = connectivity_invariant(mol, idx, ring_set);
        current_invariants.push(id);
        out.insert((i as u32, 0), id);
    }

    let mut dead = vec![false; n];
    let mut atom_neighborhoods: Vec<BondSet> = (0..n).map(|_| BondSet::empty(bond_count)).collect();
    let mut seen: FxHashSet<BondSet> = FxHashSet::default();
    let mut chiral_atoms = vec![false; n];

    for layer in 0..max_radius {
        let mut next_invariants = vec![0u32; n];
        let mut round_atom_neighborhoods = atom_neighborhoods.clone();
        let mut groups: FxHashMap<BondSet, Vec<(u32, u32)>> = FxHashMap::default();

        for i in 0..n {
            if dead[i] {
                continue;
            }
            let idx = AtomIdx(i as u32);
            let neighbors: Vec<(AtomIdx, BondIdx)> = mol.neighbors(idx).collect();
            if neighbors.is_empty() {
                dead[i] = true;
                continue;
            }
            let mut bond_env = BondSet::empty(bond_count);
            let mut pairs: SmallVec<[(u32, u32); 6]> = SmallVec::new();
            for (nb_idx, bond_idx) in &neighbors {
                bond_env.set(bond_idx.0);
                bond_env.union_with(&atom_neighborhoods[nb_idx.0 as usize]);
                pairs.push((
                    bond_invariants[bond_idx.0 as usize],
                    current_invariants[nb_idx.0 as usize],
                ));
            }
            pairs.sort_unstable();

            let mut invar = layer;
            invar = hash_combine(invar, current_invariants[i]);
            let mut looks_chiral =
                cip_codes.is_some() && mol.atom(idx).chirality != chematic_core::Chirality::None;
            let mut previous_neighbor: Option<u32> = None;
            for &(bond_inv, nb_inv) in &pairs {
                invar = hash_combine(invar, hash_pair(bond_inv, nb_inv));
                if looks_chiral {
                    if bond_inv != 1 || previous_neighbor == Some(nb_inv) {
                        looks_chiral = false;
                    }
                    previous_neighbor = Some(nb_inv);
                }
            }
            if looks_chiral {
                chiral_atoms[i] = true;
                let code = cip_codes
                    .and_then(|codes| codes.get(&idx).copied())
                    .map(chiral_code)
                    .unwrap_or(1);
                invar = hash_combine(invar, code);
            }

            next_invariants[i] = invar;
            round_atom_neighborhoods[i] = bond_env.clone();
            groups.entry(bond_env).or_default().push((invar, i as u32));
        }

        for (bond_env, mut members) in groups {
            if !suppress {
                for &(invariant, atom_idx) in &members {
                    out.insert((atom_idx, layer + 1), invariant);
                }
            } else if seen.contains(&bond_env) {
                for &(_, atom_idx) in &members {
                    dead[atom_idx as usize] = true;
                }
            } else {
                members.sort_unstable();
                let (winner_invariant, winner_idx) = members[0];
                out.insert((winner_idx, layer + 1), winner_invariant);
                seen.insert(bond_env);
                for &(_, atom_idx) in &members[1..] {
                    dead[atom_idx as usize] = true;
                }
            }
        }

        current_invariants = next_invariants;
        atom_neighborhoods = round_atom_neighborhoods;
    }

    out
}

fn chiral_code(code: CipCode) -> u32 {
    match code {
        CipCode::R => 3,
        CipCode::S => 2,
        _ => 1,
    }
}

/// Every `(atom, radius)` combination chematic's port computes, up to
/// `max_radius`, under both of RDKit's `includeRedundantEnvironments`
/// lifecycles at once -- see [`RdkitMorganRawTraceEntry`] for why they're
/// tracked as two independent optional fields rather than one shared value.
pub fn rdkit_morgan_raw_trace(mol: &Molecule, max_radius: u32) -> Vec<RdkitMorganRawTraceEntry> {
    let max_radius = max_radius.min(MAX_ECFP_RADIUS);
    if mol.atom_count() == 0 {
        return Vec::new();
    }

    let ring_set = chematic_perception::find_sssr(mol);
    let bond_count = mol.bond_count();
    let bond_invariants: Vec<u32> = (0..bond_count)
        .map(|b| bond_invariant(mol.bond(BondIdx(b as u32)).order))
        .collect();

    let full = expand_one_pass(mol, &ring_set, &bond_invariants, max_radius, false);
    let default = expand_one_pass(mol, &ring_set, &bond_invariants, max_radius, true);

    let mut keys: Vec<(u32, u32)> = full.keys().chain(default.keys()).copied().collect();
    keys.sort_unstable();
    keys.dedup();

    keys.into_iter()
        .map(|(atom_idx, radius)| RdkitMorganRawTraceEntry {
            atom_idx,
            radius,
            raw_identifier_full: full.get(&(atom_idx, radius)).copied(),
            raw_identifier_default: default.get(&(atom_idx, radius)).copied(),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn hash_combine_matches_boost_formula_by_construction() {
        // seed=0 combine(0, v) == v.wrapping_add(0x9e3779b9) (0<<6 and 0>>2 vanish, 0^x==x).
        assert_eq!(hash_combine(0, 42), 42u32.wrapping_add(0x9e37_79b9));
    }

    #[test]
    fn hash_pair_is_nested_not_flat() {
        // hash_pair(a,b) must differ from a single hash_combine(a,b) starting
        // from a nonzero seed -- it always restarts its own seed at 0.
        let a = 123u32;
        let b = 456u32;
        let nested = hash_pair(a, b);
        let flat_from_a = hash_combine(a, b);
        assert_ne!(nested, flat_from_a);
        assert_eq!(nested, hash_combine(hash_combine(0, a), b));
    }

    #[test]
    fn connectivity_invariant_ring_membership_changes_vector_length_not_just_value() {
        // Benzene carbon (in-ring) vs. an acyclic sp2 carbon of otherwise
        // identical [atomicNum,totalDegree,totalNumHs,charge,deltaMass] must
        // NOT differ by a simple XOR-with-a-constant relationship, since the
        // real difference is one extra hash_combine call, not one changed
        // component value in a fixed-length vector.
        let ring_components: SmallVec<[u32; 6]> = SmallVec::from_slice(&[6, 3, 1, 0, 0, 1]);
        let acyclic_components: SmallVec<[u32; 6]> = SmallVec::from_slice(&[6, 3, 1, 0, 0]);
        assert_ne!(hash_vec(&ring_components), hash_vec(&acyclic_components));
        // And the acyclic hash must equal hashing the same 5 components alone
        // (proving the ring bit is appended, not substituted in place).
        assert_eq!(hash_vec(&acyclic_components), hash_vec(&[6, 3, 1, 0, 0]));
    }

    #[test]
    fn bond_invariant_aromatic_is_twelve_not_chematic_own_four() {
        assert_eq!(bond_invariant(BondOrder::Single), 1);
        assert_eq!(bond_invariant(BondOrder::Double), 2);
        assert_eq!(bond_invariant(BondOrder::Triple), 3);
        assert_eq!(bond_invariant(BondOrder::Aromatic), 12);
    }

    #[test]
    fn radius_zero_trace_covers_every_atom_exactly_once() {
        let mol = parse("c1ccccc1").unwrap();
        let trace = rdkit_morgan_raw_trace(&mol, 2);
        let radius0: Vec<_> = trace.iter().filter(|e| e.radius == 0).collect();
        assert_eq!(radius0.len(), 6);
        for e in &radius0 {
            assert!(e.raw_identifier_full.is_some() && e.raw_identifier_default.is_some());
        }
    }

    #[test]
    fn degree_zero_atom_never_appears_past_radius_zero() {
        let mol = parse("[Cl-]").unwrap();
        let trace = rdkit_morgan_raw_trace(&mol, 2);
        assert_eq!(trace.len(), 1);
        assert_eq!(trace[0].radius, 0);
    }

    #[test]
    fn ethane_suppresses_one_of_two_symmetric_atoms_at_radius_one() {
        let mol = parse("CC").unwrap();
        let trace = rdkit_morgan_raw_trace(&mol, 2);
        let radius1_emitted_default = trace
            .iter()
            .filter(|e| e.radius == 1 && e.raw_identifier_default.is_some())
            .count();
        assert_eq!(radius1_emitted_default, 1);
        let radius1_computed = trace.iter().filter(|e| e.radius == 1).count();
        assert_eq!(
            radius1_computed, 2,
            "both atoms must still be *computed* under the full lifecycle"
        );
    }

    /// Regression pin for the bug this session's full-corpus comparator run
    /// caught (see [`RdkitMorganRawTraceEntry`]'s doc comment): ethane's
    /// atom 1 loses the radius-1 representative tie under `default`
    /// suppression, but must still independently compute a real radius-2
    /// `full`-lifecycle identifier — matching real RDKit's
    /// `includeRedundantEnvironments=true` output, which keeps computing a
    /// suppressed atom every round (only degree-0 death is unconditional).
    #[test]
    fn suppressed_atom_still_computes_full_lifecycle_entries_at_later_radii() {
        let mol = parse("CC").unwrap();
        let trace = rdkit_morgan_raw_trace(&mol, 2);
        let atom1_radius2 = trace
            .iter()
            .find(|e| e.atom_idx == 1 && e.radius == 2)
            .expect("atom 1 must have a radius-2 entry even though it lost the radius-1 tie");
        assert!(
            atom1_radius2.raw_identifier_full.is_some(),
            "atom 1's radius-2 full-lifecycle identifier must not be dropped"
        );
        assert!(
            atom1_radius2.raw_identifier_default.is_none(),
            "atom 1 stays suppressed under default lifecycle at radius 2"
        );
    }
}

#[cfg(test)]
mod rdkit_ground_truth_radius0 {
    //! Radius-0 invariants hand-verified against real RDKit 2026.03.3
    //! (`rdMolDescriptors.GetConnectivityInvariants(mol, includeRingMembership=True)`),
    //! not derived from this port itself — pinned as a permanent regression
    //! fixture, not just a one-time sanity check. All 8 values matched
    //! bit-for-bit on first run (2026-07-20), covering: plain/aromatic/
    //! Kekule-vs-aromatic-pyridine ring membership, a degree-zero anion, an
    //! isotope whose mass delta truncates to zero, and a formal-charge atom.
    use super::*;
    use chematic_smiles::parse;

    fn radius0(smi: &str) -> Vec<u32> {
        let mol = parse(smi).unwrap();
        let ring_set = chematic_perception::find_sssr(&mol);
        (0..mol.atom_count())
            .map(|i| connectivity_invariant(&mol, AtomIdx(i as u32), &ring_set))
            .collect()
    }

    #[test]
    fn matches_real_rdkit_getconnectivityinvariants() {
        assert_eq!(radius0("C"), vec![2246733040]);
        assert_eq!(radius0("CC"), vec![2246728737, 2246728737]);
        assert_eq!(radius0("c1ccccc1"), vec![3218693969; 6]);
        assert_eq!(
            radius0("c1ccncc1"),
            vec![
                3218693969, 3218693969, 3218693969, 2041434490, 3218693969, 3218693969
            ]
        );
        // Kekule-spelled pyridine must match the aromatic spelling exactly —
        // radius-0 depends only on atom-level properties + ring membership,
        // not on aromatic-bond flags or the specific Kekule structure.
        assert_eq!(
            radius0("C1=CC=NC=C1"),
            vec![
                3218693969, 3218693969, 3218693969, 2041434490, 3218693969, 3218693969
            ]
        );
        assert_eq!(radius0("[Cl-]"), vec![3855292234]);
        // Carbon-13's mass delta (13.00335 - 12.011 = 0.99235) truncates
        // toward zero to 0 -- same invariant as unlabeled carbon.
        assert_eq!(radius0("[13CH4]"), vec![2246733040]);
        assert_eq!(
            radius0("CC(=O)[O-]"),
            vec![2246728737, 2246699815, 864942730, 864942795]
        );
    }

    /// Full radius 0-2 trace for `c1ccncc1`, hand-verified against real
    /// RDKit's `GetSparseFingerprint(includeRedundantEnvironments=True)`
    /// `AdditionalOutput.GetBitInfoMap()`, inverted to `(atom, radius) ->
    /// raw_id & 0xFFFFFFFF` (2026-07-20). Every one of 18 `(atom, radius)`
    /// identifiers matched bit-for-bit on first run.
    #[test]
    fn matches_real_rdkit_full_trace_aromatic_pyridine() {
        let mol = parse("c1ccncc1").unwrap();
        let trace = rdkit_morgan_raw_trace(&mol, 2);
        let mut got: Vec<((u32, u32), u32)> = trace
            .iter()
            .filter_map(|e| {
                e.raw_identifier_full
                    .map(|rid| ((e.atom_idx, e.radius), rid))
            })
            .collect();
        got.sort_unstable();
        let expected = vec![
            ((0, 0), 3218693969),
            ((0, 1), 98513984),
            ((0, 2), 2763854213),
            ((1, 0), 3218693969),
            ((1, 1), 98513984),
            ((1, 2), 1207774339),
            ((2, 0), 3218693969),
            ((2, 1), 3776905034),
            ((2, 2), 1821698485),
            ((3, 0), 2041434490),
            ((3, 1), 3118255683),
            ((3, 2), 1343371647),
            ((4, 0), 3218693969),
            ((4, 1), 3776905034),
            ((4, 2), 1821698485),
            ((5, 0), 3218693969),
            ((5, 1), 98513984),
            ((5, 2), 1207774339),
        ];
        assert_eq!(got, expected);
    }

    /// The same molecule, Kekule-spelled — must produce the *identical*
    /// trace as the aromatic spelling above, but only once
    /// `apply_aromaticity` runs first (see the module docs' caller
    /// contract). Without it, ring bonds hash as literal Single/Double
    /// instead of RDKit's Aromatic=12 and every radius >= 1 identifier
    /// diverges — hand-confirmed this session, pinned here so a future
    /// change can't silently drop the precondition.
    #[test]
    fn kekule_input_matches_aromatic_input_only_after_apply_aromaticity() {
        let kekule = parse("C1=CC=NC=C1").unwrap();
        let aromatic = parse("c1ccncc1").unwrap();

        let kekule_raw = rdkit_morgan_raw_trace(&kekule, 2);
        let aromatic_trace = rdkit_morgan_raw_trace(&aromatic, 2);
        assert_ne!(
            sorted_ids(&kekule_raw),
            sorted_ids(&aromatic_trace),
            "un-aromatized Kekule input must NOT match by construction (documents the \
             precondition, doesn't just assume it)"
        );

        let kekule_perceived = chematic_perception::apply_aromaticity(&kekule);
        let kekule_perceived_trace = rdkit_morgan_raw_trace(&kekule_perceived, 2);
        assert_eq!(
            sorted_ids(&kekule_perceived_trace),
            sorted_ids(&aromatic_trace)
        );
    }

    fn sorted_ids(trace: &[RdkitMorganRawTraceEntry]) -> Vec<((u32, u32), u32)> {
        let mut ids: Vec<((u32, u32), u32)> = trace
            .iter()
            .filter_map(|e| {
                e.raw_identifier_full
                    .map(|rid| ((e.atom_idx, e.radius), rid))
            })
            .collect();
        ids.sort_unstable();
        ids
    }
}
