//! Morgan/ECFP RDKit environment-parity diagnostic PR: a per-`(atom, radius)`
//! trace of chematic's real production Morgan expansion, for an atom×radius
//! comparison against an RDKit oracle (see
//! `scripts/ecfp_rdkit_environment_parity.py`).
//!
//! This module changes no production behavior — [`morgan_trace`] calls the
//! exact same `pub(crate)` primitives [`crate::ecfp::initial_atom_id`] and
//! [`crate::ecfp::expand_atom_id`] that [`crate::ecfp::ecfp`] and
//! [`crate::ecfp::ecfp_with_bitinfo`] already use, rather than
//! reimplementing hashing (which could silently drift from production). The
//! only genuinely new logic is [`atom_ball`], a hash-independent graph BFS
//! that gives the "environment membership" ground truth on chematic's side,
//! the counterpart to RDKit's `Chem.FindAtomEnvironmentOfRadiusN`.
//!
//! Not a stable public API: this module is declared private in `lib.rs`, so
//! these items are reachable outside the crate only via the `diagnostics`
//! Cargo feature's re-export (off by default), matching the precedent in
//! `chematic-perception`'s `diagnostics` module.

#![allow(dead_code)] // only reachable via the `diagnostics` feature + this file's own tests

use chematic_core::{AtomIdx, Molecule};
use chematic_perception::find_sssr;
use rustc_hash::FxHashSet;

use crate::ecfp::{
    EcfpConfig, EcfpInvariantMode, MAX_ECFP_RADIUS, expand_atom_id, initial_atom_id,
};

/// One `(atom, radius)` environment as chematic's real production code
/// computed it.
#[derive(Debug, Clone)]
pub struct MorganTraceEntry {
    pub atom_idx: u32,
    pub radius: u32,
    /// Raw (unfolded) FNV-1a id: `initial_atom_id`'s return at radius 0,
    /// `expand_atom_id`'s return at radius >= 1.
    pub raw_environment_id: u64,
    /// `raw_environment_id % nbits` — the bit this environment would set in
    /// the folded fingerprint.
    pub folded_bit: usize,
    /// Always `true` today: chematic has no redundant-environment
    /// suppression, so every atom is emitted at every radius. Present so
    /// this trace shape can also represent RDKit's oracle side (which does
    /// suppress), without a reshape once a future PR adds real suppression.
    pub emitted: bool,
    /// Sorted atom indices within `radius` bond-hops of `atom_idx` — the
    /// hash-independent BFS ball, chematic's counterpart to RDKit's
    /// `Chem.FindAtomEnvironmentOfRadiusN`-derived atom set.
    pub atom_ball: Vec<u32>,
}

/// Trace every `(atom, radius)` environment chematic's Morgan expansion
/// computes for `mol` under `config`/`mode`, iteration 0 through
/// `config.radius`.
pub fn morgan_trace(
    mol: &Molecule,
    config: &EcfpConfig,
    mode: EcfpInvariantMode,
) -> Vec<MorganTraceEntry> {
    let n = mol.atom_count();
    let nbits = config.nbits;
    let radius = config.radius.min(MAX_ECFP_RADIUS);

    let mut trace = Vec::new();
    if n == 0 {
        return trace;
    }

    let ring_set = find_sssr(mol);

    // Iteration 0: initial atom identifiers (same call ecfp_with_bitinfo_and_mode makes).
    let mut ids: Vec<u64> = Vec::with_capacity(n);
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let id = initial_atom_id(mol, idx, &ring_set, config.use_chirality, mode);
        trace.push(MorganTraceEntry {
            atom_idx: i as u32,
            radius: 0,
            raw_environment_id: id,
            folded_bit: (id % nbits as u64) as usize,
            emitted: true,
            atom_ball: vec![i as u32],
        });
        ids.push(id);
    }

    // Iterations 1..=radius: expansion (same call ecfp_with_bitinfo_and_mode makes).
    let mut new_ids = vec![0u64; n];
    for r in 1..=radius {
        for (i, slot) in new_ids.iter_mut().enumerate() {
            let new_id = expand_atom_id(mol, i, r, &ids);
            *slot = new_id;
            trace.push(MorganTraceEntry {
                atom_idx: i as u32,
                radius: r,
                raw_environment_id: new_id,
                folded_bit: (new_id % nbits as u64) as usize,
                emitted: true,
                atom_ball: atom_ball(mol, AtomIdx(i as u32), r),
            });
        }
        core::mem::swap(&mut ids, &mut new_ids);
    }

    trace
}

/// Sorted atom indices within `radius` bond-hops of `center` — a plain graph
/// BFS, no hashing involved. RDKit's `Chem.FindAtomEnvironmentOfRadiusN`
/// ground truth is this same notion (bonds within radius hops); comparing
/// this set against it detects environment-membership divergence
/// independently of any hash-identifier divergence.
pub fn atom_ball(mol: &Molecule, center: AtomIdx, radius: u32) -> Vec<u32> {
    let mut visited: FxHashSet<u32> = FxHashSet::default();
    visited.insert(center.0);
    let mut frontier = vec![center];
    for _ in 0..radius {
        if frontier.is_empty() {
            break;
        }
        let mut next = Vec::new();
        for &atom in &frontier {
            for (nb, _bond) in mol.neighbors(atom) {
                if visited.insert(nb.0) {
                    next.push(nb);
                }
            }
        }
        frontier = next;
    }
    let mut result: Vec<u32> = visited.into_iter().collect();
    result.sort_unstable();
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecfp::{
        ecfp, ecfp_with_bitinfo, ecfp_with_bitinfo_and_mode, ecfp_with_invariant_mode,
        morgan_fp_counts,
    };
    use chematic_smiles::parse;
    use rustc_hash::FxHashMap;
    use std::collections::HashSet;

    fn molecules() -> Vec<Molecule> {
        vec![
            parse("c1ccccc1").unwrap(),              // benzene
            parse("CC(=O)Oc1ccccc1C(=O)O").unwrap(), // aspirin
            parse("Cc1ccccc1").unwrap(),             // toluene
        ]
    }

    fn folded_bits(trace: &[MorganTraceEntry]) -> HashSet<usize> {
        trace.iter().map(|e| e.folded_bit).collect()
    }

    fn bitinfo_from_trace(trace: &[MorganTraceEntry]) -> FxHashMap<usize, Vec<(u32, u32)>> {
        let mut info: FxHashMap<usize, Vec<(u32, u32)>> = FxHashMap::default();
        for e in trace {
            info.entry(e.folded_bit)
                .or_default()
                .push((e.atom_idx, e.radius));
        }
        for v in info.values_mut() {
            v.sort_unstable();
        }
        info
    }

    #[test]
    fn trace_matches_production_chematic_mode() {
        let config = EcfpConfig::default();
        for mol in molecules() {
            let trace = morgan_trace(&mol, &config, EcfpInvariantMode::Chematic);
            assert!(trace.iter().all(|e| e.emitted), "chematic never suppresses");

            let fp = ecfp(&mol, &config);
            let fp_bits: HashSet<usize> = (0..config.nbits).filter(|&i| fp.get(i)).collect();
            assert_eq!(
                folded_bits(&trace),
                fp_bits,
                "trace folded-bit set must match ecfp()'s real output"
            );

            let (_, mut info) = ecfp_with_bitinfo(&mol, &config);
            for v in info.values_mut() {
                v.sort_unstable();
            }
            assert_eq!(
                bitinfo_from_trace(&trace),
                info,
                "trace bitinfo grouping must match ecfp_with_bitinfo()'s real output"
            );

            let mut trace_counts: FxHashMap<u64, u32> = FxHashMap::default();
            for e in &trace {
                *trace_counts.entry(e.raw_environment_id).or_insert(0) += 1;
            }
            assert_eq!(
                trace_counts,
                morgan_fp_counts(&mol, config.radius),
                "trace raw-id counts must match morgan_fp_counts()'s real output"
            );
        }
    }

    #[test]
    fn trace_matches_production_rdkit_morgan_mode() {
        let config = EcfpConfig::default();
        for mol in molecules() {
            let trace = morgan_trace(&mol, &config, EcfpInvariantMode::RdkitMorgan);

            let fp = ecfp_with_invariant_mode(&mol, &config, EcfpInvariantMode::RdkitMorgan);
            let fp_bits: HashSet<usize> = (0..config.nbits).filter(|&i| fp.get(i)).collect();
            assert_eq!(folded_bits(&trace), fp_bits);

            let (_, mut info) =
                ecfp_with_bitinfo_and_mode(&mol, &config, EcfpInvariantMode::RdkitMorgan);
            for v in info.values_mut() {
                v.sort_unstable();
            }
            assert_eq!(bitinfo_from_trace(&trace), info);
        }
    }

    #[test]
    fn atom_ball_radius0_is_self() {
        let mol = parse("CCO").unwrap();
        assert_eq!(atom_ball(&mol, AtomIdx(0), 0), vec![0]);
    }

    #[test]
    fn atom_ball_grows_with_radius() {
        let mol = parse("CCCCC").unwrap(); // pentane, atoms 0..4 in a chain
        assert_eq!(atom_ball(&mol, AtomIdx(0), 1), vec![0, 1]);
        assert_eq!(atom_ball(&mol, AtomIdx(0), 2), vec![0, 1, 2]);
    }

    #[test]
    fn atom_ball_stops_growing_past_molecule_diameter() {
        // Ethane: atom 0's ball is already the whole (2-atom) molecule at
        // radius 1 -- radius 2 must not panic or grow past that.
        let mol = parse("CC").unwrap();
        assert_eq!(atom_ball(&mol, AtomIdx(0), 1), vec![0, 1]);
        assert_eq!(atom_ball(&mol, AtomIdx(0), 2), vec![0, 1]);
    }

    #[test]
    fn empty_molecule_trace_is_empty() {
        use chematic_core::MoleculeBuilder;
        let mol = MoleculeBuilder::new().build();
        let trace = morgan_trace(&mol, &EcfpConfig::default(), EcfpInvariantMode::Chematic);
        assert!(trace.is_empty());
    }
}
