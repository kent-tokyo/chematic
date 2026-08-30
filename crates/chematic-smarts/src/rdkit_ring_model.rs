//! An auxiliary, SMARTS-predicate-only ring model approximating RDKit's
//! default `RingInfo` for `[RN]` (ring-*count*) matching.
//!
//! **Why this exists.** RDKit's default sanitization
//! (`Chem.MolFromSmiles(..., sanitize=True)`) does not stop at a minimal
//! SSSR: `MolOps::symmetrizeSSSR` (RDKit source,
//! `Code/GraphMol/FindRings.cpp:996-1093`, pinned commit
//! `8afba32ec539dcb2369bc84549d802aca3f7eb39`) adds "extra" rings back into
//! `RingInfo` whenever an extra candidate (a) is the same size as some basis
//! ring, (b) shares at least one bond with that basis ring, and (c) does not
//! drop any bond the basis ring is the *sole* provider of (see
//! `symmetrizeSSSR`'s inline comment at that file, "may miss extra rings
//! that would need to swap two (or three...) rings", and the
//! `bondCounts`/`replacesAllUniqueBonds` logic immediately below it). RDKit's
//! own candidates come from *its own* SSSR search's rejected duplicate D2
//! candidates (`findSSSRforDupCands`, same file, line ~283); this module
//! does not have access to that intermediate state (chematic's `find_sssr`
//! doesn't expose it, and this module must not change `find_sssr`'s public
//! output/behavior per this crate's scope). Instead it re-derives candidate
//! "extra" rings from chematic's own already-computed SSSR basis, by
//! enumerating all root-centered shortest rings and re-applying RDKit's own
//! substitution rule (b)/(c) above. This keeps the candidate pool bounded to
//! D2-like shortest cycles rather than admitting arbitrary longer simple
//! cycles, while remaining an approximation until the rejected Horton
//! candidates are exposed. See `docs/rdkit_compat.md`'s "SMARTS-R2" section
//! for the measured over-/under-generation split relative to a live RDKit
//! oracle.
//!
//! **What this changes, and what it deliberately does not.** A graph-theory
//! fact bounds the blast radius tightly: an edge lies on *some* basis cycle
//! of a graph if and only if it lies on *any* cycle at all (a non-bridge
//! edge), independent of which basis is chosen — because a cycle's edge set
//! is a member of the cycle space over GF(2), so if an edge appeared in zero
//! basis cycles it could not appear in any GF(2) combination of them either.
//! Consequently `[R]`/`[R0]` (ring *membership*, boolean), `[x]`/`[xN]`
//! (ring-*bond*-count) and ring-bond `@`/`!@` are **provably invariant** to
//! which valid SSSR basis is chosen — adding RDKit's "extra" rings can never
//! change any of those three predicates' verdicts, only `[RN]` (N ≥ 1, exact
//! ring *count*) can move, because an atom can gain membership in an
//! additional same-size ring beyond its raw-SSSR count. `[rN]`/`[kN]`
//! (min-ring-size / any-ring-of-size-N) *could* in principle move too (an
//! extra ring's atom set differs from the basis ring it substitutes for),
//! but both already measure at ~100%/99.98% against RDKit using chematic's
//! plain SSSR alone (`docs/rdkit_compat.md`'s SMARTS-R1 section) — the
//! opt-in matcher (`crate::rdkit_parity_match`) deliberately leaves them
//! wired to the plain SSSR, unchanged, and only routes `AtomPrimitive::RingCount`
//! through this module's model.
//!
//! **Termination.** Simple-cycle enumeration is depth-bounded by the basis's
//! largest ring size (not by atom count or ring count), so cost scales with
//! `atoms × max_degree^max_ring_size`, not with the size of the ring system.
//! A hard candidate-count cap (independent defense-in-depth, not the primary
//! bound) turns any unexpectedly-adversarial topology into
//! [`RdkitParityError::RingModelBudgetExceeded`] instead of a hang — see
//! [`build_rdkit_parity_ring_model`].

use rustc_hash::{FxHashMap, FxHashSet};

use chematic_core::{AtomIdx, BondIdx, Molecule};
use chematic_perception::{
    RingSet, find_smallest_rings_bfs, find_smallest_rings_bfs_with_blocked_bonds,
    find_smallest_rings_bfs_with_rdkit_tree, select_rdkit_d2_roots, trim_ring_bonds,
};

/// Typed error for chematic-smarts's opt-in RDKit-parity matching mode.
///
/// Never silently falls back to a partial or non-parity result — see each
/// variant's doc comment for what it guards.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RdkitParityError {
    /// The RDKit-parity ring-count model's bounded simple-cycle search hit
    /// its candidate-count cap before finishing. Rather than silently
    /// returning a ring-count model built from a truncated candidate search
    /// (which could under-report `[RN]` matches while *appearing* to be a
    /// full RDKit-parity result), the whole call fails closed.
    ///
    /// This is a resource/complexity bound on this module's own
    /// candidate-generation approximation of RDKit's `symmetrizeSSSR` — it
    /// is not the shared VF2 [`crate::MatchOutcome::BudgetExhausted`]
    /// (VF2 state-space search budget), which is a separate, orthogonal
    /// budget covering the actual subgraph-isomorphism search.
    RingModelBudgetExceeded {
        /// Number of candidate cycles examined before the cap was hit.
        candidates_examined: usize,
        /// The configured cap ([`RdkitRingModelBudget::max_candidates`]).
        cap: usize,
    },
    /// RDKit-parity aromaticity preprocessing
    /// (`chematic_perception::apply_aromaticity_rdkit_parity_experimental`)
    /// was requested (`RdkitParityConfig::use_rdkit_parity_aromaticity`) and
    /// failed — most commonly because the target molecule's Kekulé
    /// structure could not be derived, a known limitation of that engine on
    /// a small class of bridgehead-heteroatom-fused rings (see that
    /// function's own doc comment). Never silently substituted with the
    /// default Hückel engine's aromatic flags, since the two are not
    /// guaranteed to agree.
    Aromaticity(chematic_perception::AromaticityError),
}

impl std::fmt::Display for RdkitParityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RdkitParityError::RingModelBudgetExceeded {
                candidates_examined,
                cap,
            } => write!(
                f,
                "RDKit-parity ring model: candidate cycle search exceeded budget \
                 ({candidates_examined} candidates examined, cap {cap}); \
                 no ring-count model was produced"
            ),
            RdkitParityError::Aromaticity(e) => {
                write!(f, "RDKit-parity aromaticity preprocessing failed: {e}")
            }
        }
    }
}

impl std::error::Error for RdkitParityError {}

/// Resource bound for [`build_rdkit_parity_ring_model`]'s simple-cycle
/// search. See the module doc comment's "Termination" section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RdkitRingModelBudget {
    /// Maximum number of candidate cycles examined (closed simple paths
    /// tested against the substitution rule) before failing closed with
    /// [`RdkitParityError::RingModelBudgetExceeded`].
    pub max_candidates: usize,
}

impl Default for RdkitRingModelBudget {
    fn default() -> Self {
        // Generous for any realistic drug-like or cage/cluster molecule: a
        // depth-bounded DFS over a few hundred atoms with typical organic
        // valence (degree <= 4) and a ring size <= ~8 stays several orders
        // of magnitude below this, per the module doc's cost model.
        Self {
            max_candidates: 2_000_000,
        }
    }
}

/// RDKit-parity ring-*count* model: an auxiliary, atom-indexed table
/// approximating `RingInfo::numAtomRings` under RDKit's default
/// `symmetrizeSSSR`, built from chematic's own SSSR without modifying it.
///
/// Only backs `AtomPrimitive::RingCount` (`[RN]`, N >= 1) in
/// [`crate::rdkit_parity_match`] — see the module doc comment for why every
/// other ring-shaped SMARTS predicate is left on plain SSSR.
#[derive(Debug, Clone)]
pub struct RdkitParityRingModel {
    ring_count_by_atom: FxHashMap<AtomIdx, u8>,
    /// Extra rings accepted beyond the raw SSSR basis (for diagnostics/tests).
    extra_ring_count: usize,
}

impl RdkitParityRingModel {
    /// Number of accepted rings (SSSR basis + extras) containing `atom`.
    pub fn ring_count(&self, atom: AtomIdx) -> u8 {
        self.ring_count_by_atom.get(&atom).copied().unwrap_or(0)
    }

    /// Number of "extra" (beyond the raw SSSR basis) rings this model
    /// accepted. Zero for the overwhelming majority of molecules — see the
    /// module doc comment; only certain highly symmetric cage/cluster
    /// topologies (cubane, prismane, adamantane, bicyclo[2.2.2]octane-type
    /// bridged cages) are expected to produce any.
    pub fn extra_ring_count(&self) -> usize {
        self.extra_ring_count
    }
}

/// Build an [`RdkitParityRingModel`] for `mol`, given its already-computed
/// `sssr` (never recomputed here — this module never calls `find_sssr`
/// itself, it only consumes a caller-supplied [`RingSet`], matching the
/// existing `find_matches_with_rings`-style call convention).
pub fn build_rdkit_parity_ring_model(
    mol: &Molecule,
    sssr: &RingSet,
    budget: &RdkitRingModelBudget,
) -> Result<RdkitParityRingModel, RdkitParityError> {
    let base_rings = sssr.rings();
    if base_rings.is_empty() {
        return Ok(RdkitParityRingModel {
            ring_count_by_atom: FxHashMap::default(),
            extra_ring_count: 0,
        });
    }

    let base_bond_sets: Vec<FxHashSet<BondIdx>> =
        base_rings.iter().map(|r| ring_bond_set(mol, r)).collect();
    let base_bond_key: FxHashSet<Vec<u32>> = base_bond_sets.iter().map(bond_set_key).collect();

    // How many base rings each bond appears in -- needed for the
    // "does the candidate drop a uniquely-provided bond" check.
    let mut bond_ring_count: FxHashMap<BondIdx, u32> = FxHashMap::default();
    for bs in &base_bond_sets {
        for &b in bs {
            *bond_ring_count.entry(b).or_insert(0) += 1;
        }
    }

    // Base ring indices grouped by size (only sizes that appear in the SSSR
    // basis are ever eligible substitution targets, per RDKit's rule).
    let mut by_size: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for (i, r) in base_rings.iter().enumerate() {
        by_size.entry(r.len()).or_default().push(i);
    }
    let mut seen_extra: FxHashSet<Vec<u32>> = FxHashSet::default();
    let mut extra_rings: Vec<Vec<AtomIdx>> = Vec::new();
    let mut candidates_examined = 0usize;

    // RDKit's duplicate-D2 pool is rooted and shortest-ring based. Consume
    // the same bounded primitive used by the perception crate instead of
    // enumerating every simple cycle up to the largest SSSR size; the latter
    // admits non-D2 cycles and over-produces extras in macrocycles.
    let d2_roots = select_rdkit_d2_roots(mol);
    // RDKit's Figueras search starts with D2 nodes.  Highly symmetric cages
    // such as cubane have no D2 nodes and are handled by the later D3-style
    // fallback, for which considering every root is the useful approximation.
    let has_d2_roots = !d2_roots.is_empty();
    let roots: Box<dyn Iterator<Item = AtomIdx>> = if has_d2_roots {
        Box::new(d2_roots.iter().copied())
    } else {
        Box::new((0..mol.atom_count()).map(|raw| AtomIdx(raw as u32)))
    };
    let mut accept_candidate = |candidate: Vec<AtomIdx>| -> bool {
        let len = candidate.len();
        let Some(ring_idxs) = by_size.get(&len) else {
            return false;
        };
        let cand_set = ring_bond_set(mol, &candidate);
        let key = bond_set_key(&cand_set);
        if base_bond_key.contains(&key) || !seen_extra.insert(key) {
            return false;
        }
        for &ring_idx in ring_idxs {
            let base_set = &base_bond_sets[ring_idx];
            if base_set.iter().any(|b| cand_set.contains(b))
                && base_set.iter().all(|b| {
                    bond_ring_count.get(b).copied().unwrap_or(0) != 1 || cand_set.contains(b)
                })
            {
                extra_rings.push(candidate);
                return true;
            }
        }
        false
    };

    if !has_d2_roots {
        for root in roots {
            for candidate in find_smallest_rings_bfs(mol, root) {
                candidates_examined += 1;
                if candidates_examined > budget.max_candidates {
                    return Err(RdkitParityError::RingModelBudgetExceeded {
                        candidates_examined,
                        cap: budget.max_candidates,
                    });
                }
                accept_candidate(candidate);
            }
        }
    } else {
        // First identify rings discovered by more than one D2 node. RDKit
        // treats these as duplicate candidates and does not accept the
        // original ring again; it removes the competing D2 connections and
        // searches for the next-smallest replacement ring.
        let mut duplicate_groups: FxHashMap<Vec<u32>, (Vec<AtomIdx>, Vec<AtomIdx>)> =
            FxHashMap::default();
        let mut active_blocked = FxHashSet::default();
        for root in roots {
            let candidates = find_smallest_rings_bfs_with_blocked_bonds(mol, root, &active_blocked);
            if candidates.is_empty() {
                for (_, bond) in mol.neighbors(root) {
                    if is_ring_bond(mol, bond) {
                        active_blocked.insert(bond);
                    }
                }
                active_blocked = trim_ring_bonds(mol, &active_blocked);
                continue;
            }
            for candidate in candidates {
                candidates_examined += 1;
                if candidates_examined > budget.max_candidates {
                    return Err(RdkitParityError::RingModelBudgetExceeded {
                        candidates_examined,
                        cap: budget.max_candidates,
                    });
                }
                let key = bond_set_key(&ring_bond_set(mol, &candidate));
                let entry = duplicate_groups
                    .entry(key)
                    .or_insert_with(|| (candidate.clone(), Vec::new()));
                if !entry.1.contains(&root) {
                    entry.1.push(root);
                }
            }
        }

        for (_, (original_candidate, duplicate_roots)) in duplicate_groups {
            if duplicate_roots.len() < 2 {
                // A ring found by only one D2 node is not a duplicate-D2
                // candidate; retain the normal Figueras acceptance path.
                accept_candidate(original_candidate);
                continue;
            }
            let mut replacement_candidates = Vec::new();
            for &root in &duplicate_roots {
                let mut blocked_bonds = FxHashSet::default();
                for &other in &duplicate_roots {
                    if other == root {
                        continue;
                    }
                    for (_, bond) in mol.neighbors(other) {
                        if is_ring_bond(mol, bond) {
                            blocked_bonds.insert(bond);
                        }
                    }
                }
                for candidate in find_smallest_rings_bfs_with_rdkit_tree(mol, root, &blocked_bonds)
                {
                    candidates_examined += 1;
                    if candidates_examined > budget.max_candidates {
                        return Err(RdkitParityError::RingModelBudgetExceeded {
                            candidates_examined,
                            cap: budget.max_candidates,
                        });
                    }
                    replacement_candidates.push(candidate);
                }
            }
            if let Some(min_size) = replacement_candidates.iter().map(Vec::len).min() {
                replacement_candidates.retain(|candidate| candidate.len() == min_size);
            }
            // RDKit gathers the minimum-size replacements for the duplicate
            // group and then inserts that minimum set. Choosing the first
            // accepted candidate in a stable bond-key order avoids adding all
            // equivalent paths emitted by the shortest-path enumerator (the
            // latter was the source of three extras for one macrocycle).
            replacement_candidates
                .sort_by_key(|candidate| bond_set_key(&ring_bond_set(mol, candidate)));
            replacement_candidates
                .into_iter()
                .any(&mut accept_candidate);
        }
    }

    // Several D2 groups can describe the same connected symmetry class. Keep
    // one stable representative per bond-overlap component; retaining every
    // overlapping replacement would inflate the symmetrized count (notably
    // the issue-337 macrocycle family) while adding no new ring frontier.
    extra_rings.sort_by_key(|ring| bond_set_key(&ring_bond_set(mol, ring)));
    let mut selected_extra_sets: Vec<FxHashSet<BondIdx>> = Vec::new();
    extra_rings.retain(|ring| {
        let candidate_set = ring_bond_set(mol, ring);
        if selected_extra_sets
            .iter()
            .any(|selected| selected.iter().any(|bond| candidate_set.contains(bond)))
        {
            false
        } else {
            selected_extra_sets.push(candidate_set);
            true
        }
    });

    let mut ring_count_by_atom: FxHashMap<AtomIdx, u8> = FxHashMap::default();
    for ring in base_rings.iter().chain(extra_rings.iter()) {
        for &a in ring {
            *ring_count_by_atom.entry(a).or_insert(0) += 1;
        }
    }

    Ok(RdkitParityRingModel {
        ring_count_by_atom,
        extra_ring_count: extra_rings.len(),
    })
}

/// Bond-index set of a ring given as a cyclic atom sequence (consecutive
/// atoms are bonded — the same convention `chematic_perception::RingSet`
/// documents its own rings with).
fn ring_bond_set(mol: &Molecule, ring: &[AtomIdx]) -> FxHashSet<BondIdx> {
    let n = ring.len();
    (0..n)
        .filter_map(|i| {
            let a = ring[i];
            let b = ring[(i + 1) % n];
            mol.bond_between(a, b).map(|(bidx, _)| bidx)
        })
        .collect()
}

fn is_ring_bond(mol: &Molecule, bond: BondIdx) -> bool {
    matches!(
        mol.bond(bond).order,
        chematic_core::BondOrder::Single
            | chematic_core::BondOrder::Double
            | chematic_core::BondOrder::Triple
            | chematic_core::BondOrder::Quadruple
            | chematic_core::BondOrder::Aromatic
            | chematic_core::BondOrder::Up
            | chematic_core::BondOrder::Down
    )
}

fn bond_set_key(set: &FxHashSet<BondIdx>) -> Vec<u32> {
    let mut v: Vec<u32> = set.iter().map(|b| b.0).collect();
    v.sort_unstable();
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_perception::find_sssr;
    use chematic_smiles::parse;

    fn model_for(smiles: &str) -> (Molecule, RdkitParityRingModel) {
        let mol = parse(smiles).unwrap();
        let sssr = find_sssr(&mol);
        let model =
            build_rdkit_parity_ring_model(&mol, &sssr, &RdkitRingModelBudget::default()).unwrap();
        (mol, model)
    }

    // -- No-extra cases: plain SSSR already matches RDKit's symmetrized count --

    #[test]
    fn benzene_no_extra() {
        let (_, model) = model_for("c1ccccc1");
        assert_eq!(model.extra_ring_count(), 0);
    }

    #[test]
    fn naphthalene_no_extra() {
        let (mol, model) = model_for("c1ccc2ccccc2c1");
        assert_eq!(model.extra_ring_count(), 0);
        // Fusion atoms: ring count 2, matches RDKit ground truth.
        let fusion_count = (0..mol.atom_count())
            .filter(|&i| model.ring_count(AtomIdx(i as u32)) == 2)
            .count();
        assert_eq!(fusion_count, 2);
    }

    #[test]
    fn spiro_no_extra() {
        // spiro[4.4]nonane
        let (mol, model) = model_for("C1CCC2(C1)CCCC2");
        assert_eq!(model.extra_ring_count(), 0);
        let spiro_count = (0..mol.atom_count())
            .filter(|&i| model.ring_count(AtomIdx(i as u32)) == 2)
            .count();
        assert_eq!(spiro_count, 1);
    }

    #[test]
    fn norbornane_no_extra() {
        // Ground truth (rdkit==2026.03.3, live oracle): NumRings()==2, same as
        // chematic's raw SSSR (cycle_rank 2) -- no basis-cardinality gap here.
        let (_, model) = model_for("C1CC2CCC1C2");
        assert_eq!(model.extra_ring_count(), 0);
    }

    // -- Extra-ring cases: RDKit's symmetrizeSSSR adds one same-size alternate --

    #[test]
    fn cubane_one_extra_face() {
        // Ground truth (rdkit==2026.03.3, live oracle, scratch/rdkit_ring_ground_truth.py):
        // cycle_rank 5, RDKit NumRings()==6 (the 6th face is the GF(2) sum of
        // all 5 basis 4-rings) -- every atom ends up in exactly 3 rings.
        let (mol, model) = model_for("C12C3C4C1C5C4C3C25");
        assert_eq!(model.extra_ring_count(), 1);
        for i in 0..mol.atom_count() {
            assert_eq!(model.ring_count(AtomIdx(i as u32)), 3, "atom {i}");
        }
    }

    #[test]
    fn adamantane_one_extra_ring() {
        // Ground truth: cycle_rank 3, RDKit NumRings()==4.
        let (_, model) = model_for("C1C2CC3CC1CC(C2)C3");
        assert_eq!(model.extra_ring_count(), 1);
    }

    #[test]
    fn bicyclo_2_2_2_octane_one_extra_ring() {
        // Ground truth: cycle_rank 2, RDKit NumRings()==3 (all three
        // six-membered "faces" of the cage are captured).
        let (mol, model) = model_for("C1CC2CCC1CC2");
        assert_eq!(model.extra_ring_count(), 1);
        // The two bridgehead atoms sit in all 3 rings.
        let bridgeheads = (0..mol.atom_count())
            .filter(|&i| model.ring_count(AtomIdx(i as u32)) == 3)
            .count();
        assert_eq!(bridgeheads, 2);
    }

    // -- Adversarial highly-symmetric target: dodecahedrane --
    //
    // Ground truth (rdkit==2026.03.3, live oracle, PubChem CID 123218 SMILES
    // `C12C3C4C5C1C6C7C2C8C3C9C4C1C5C6C2C7C8C9C12`, independently curl'd from
    // PubChem's PUG REST API and cross-checked against the live RDKit
    // oracle): 20 atoms, 30 bonds, cycle_rank 11 (the minimal SSSR basis),
    // RDKit's symmetrized `NumRings()` == 12 -- the classic dodecahedron's 12
    // real pentagonal faces, one more than the minimal basis needs, with
    // *every* atom sitting in exactly 3 faces (a real vertex of a regular
    // dodecahedron touches exactly 3 pentagons). This is exactly the
    // "large symmetric cage" class the gate calls for: real, chemically
    // valid, maximally symmetric (all 20 atoms equivalent), forces this
    // module's candidate search to consider all 11 basis rings together to
    // recover the 12th face (same shape as cubane's 6th face, one order of
    // magnitude bigger).
    #[test]
    fn dodecahedrane_terminates_and_matches_rdkit_symmetrized_count() {
        let smiles = "C12C3C4C5C1C6C7C2C8C3C9C4C1C5C6C2C7C8C9C12";
        let mol = parse(smiles).unwrap();
        assert_eq!(mol.atom_count(), 20);
        let sssr = find_sssr(&mol);
        assert_eq!(sssr.ring_count(), 11, "minimal SSSR basis size");

        let budget = RdkitRingModelBudget::default();
        let model = build_rdkit_parity_ring_model(&mol, &sssr, &budget)
            .expect("must terminate within the default budget, not hang or exceed it");
        assert_eq!(
            model.extra_ring_count(),
            1,
            "RDKit finds exactly 1 extra face"
        );
        for i in 0..mol.atom_count() {
            assert_eq!(
                model.ring_count(AtomIdx(i as u32)),
                3,
                "every dodecahedrane atom sits in exactly 3 faces (atom {i})"
            );
        }
    }

    #[test]
    fn dodecahedrane_with_starved_budget_fails_closed_not_a_hang() {
        // A budget too small to find the 12th face must report the typed
        // error, never hang and never silently return the (wrong,
        // under-augmented) 11-ring model as if it were complete.
        let smiles = "C12C3C4C5C1C6C7C2C8C3C9C4C1C5C6C2C7C8C9C12";
        let mol = parse(smiles).unwrap();
        let sssr = find_sssr(&mol);
        let budget = RdkitRingModelBudget { max_candidates: 5 };
        let result = build_rdkit_parity_ring_model(&mol, &sssr, &budget);
        assert!(matches!(
            result,
            Err(RdkitParityError::RingModelBudgetExceeded { .. })
        ));
    }

    #[test]
    fn issue_337_macrocycle_keeps_one_extra_ring_representative() {
        let smiles = "c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2";
        let (_, model) = model_for(smiles);
        assert_eq!(model.extra_ring_count(), 1);
    }

    #[test]
    fn tiny_budget_reports_typed_error_not_a_hang_or_silent_partial() {
        let mol = parse("C1C2CC3CC1CC(C2)C3").unwrap(); // adamantane
        let sssr = find_sssr(&mol);
        let budget = RdkitRingModelBudget { max_candidates: 0 };
        let result = build_rdkit_parity_ring_model(&mol, &sssr, &budget);
        assert!(matches!(
            result,
            Err(RdkitParityError::RingModelBudgetExceeded { .. })
        ));
    }
}
