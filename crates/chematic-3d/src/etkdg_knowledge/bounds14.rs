//! Macrocycle 1-4 bound adjustments (Wave 2 spec §6).
//!
//! A **pure API** Coordinator can call when constructing bounds later --
//! this PR does **not** apply these adjustments to any real bounds matrix
//! (`dg_fft::build_bound_matrix` / `distance_geometry_v2.rs`, neither of
//! which this PR edits). Evaluated only via this module's own tests and the
//! `torsion_knowledge_v2_gap_check` example's Arm D (which applies the
//! adjustments to a **copy** of the harness's own working bounds to simulate
//! what a future integration would do, never to the production embedder).
//!
//! Translated/adapted from RDKit's real macrocycle-1-4 algorithm:
//! - `Code/Geometry/Utils.h`'s `compute14DistCis`/`compute14DistTrans`
//!   (law-of-cosines-in-a-rotated-plane construction, Copyright (C)
//!   2004-2006 Rational Discovery LLC, BSD-3-Clause) -- translated near-
//!   verbatim as [`dist14_cis`]/[`dist14_trans`] below.
//! - `Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp`'s
//!   `_setMacrocycleTwoInSameRing14Bounds`/`_setMacrocycleAllInSameRing14Bounds`
//!   (amide/ester bonds pinned to one configuration; every other 1-4 pair
//!   relaxed to a wide cis-or-trans band) -- the core decision this module
//!   reproduces, re-expressed against `chematic_core::Molecule` and this
//!   crate's own `classify::BondClassification::amide_like` (rather than
//!   RDKit's own bond-type-triple check, which finds the same physical
//!   case with less duplicated code -- a documented simplification, not a
//!   different rule).
//! - `minMacrocycleRingSize = 9` (`BoundsMatrixBuilder.cpp:36`), identical to
//!   this crate's own [`super::classify::MACROCYCLE_MIN`].
//!
//! "Old" bounds represent what a naive, macrocycle-**unaware** 1-4 rule
//! would set: pinned to the single trans-configuration distance
//! (`dist14_trans`) with a small tolerance, mirroring how a generic
//! (non-macrocycle) 1-4 rule commits to one computed value rather than a
//! band. "New" bounds are this module's proposed macrocycle-aware
//! adjustment: for amide/ester-like central bonds, still pinned to a single
//! configuration but the **cis** one (matching real amide-bond rigidity,
//! `_setMacrocycleTwoInSameRing14Bounds`'s amide/ester branch); for every
//! other 1-4 pair, relaxed to the full `[min(cis,trans), max(cis,trans)]`
//! band ("assume anything is possible", RDKit's own inline comment in the
//! generic branch of `_setMacrocycleAllInSameRing14Bounds`).

use chematic_core::{AtomIdx, Molecule};

use super::classify::{MACROCYCLE_MIN, RingMembershipIndex, classify_bond};
use super::types::{TorsionKnowledgeConfig, TorsionKnowledgeError, TorsionKnowledgeSource};
use crate::dg_fft::{ideal_bond_angle, ideal_bond_length};

/// Small fixed tolerance (Å) around a pinned single-configuration 1-4
/// distance -- same order of magnitude as RDKit's own `GEN_DIST_TOL`
/// constant (`BoundsMatrixBuilder.cpp`), not independently re-derived.
const PIN_TOLERANCE_ANGSTROM: f64 = 0.05;

/// One proposed adjustment to a single atom pair's `[lower, upper]` bound.
#[derive(Clone, Debug)]
pub struct PairBoundAdjustment {
    pub atom_pair: (AtomIdx, AtomIdx),
    pub old_lower: f64,
    pub new_lower: f64,
    pub old_upper: f64,
    pub new_upper: f64,
    pub rule_id: String,
    pub source: TorsionKnowledgeSource,
    pub ring_size: usize,
    pub reason: String,
}

/// `d1`, `d2`, `d3`: 1-2, 2-3, 3-4 bond lengths (Å). `ang12`, `ang23`: the
/// bond angles at atoms 2 and 3 (radians). Returns the 1-4 (atom1-atom4)
/// distance assuming a **cis** (0-degree torsion) configuration around the
/// 2-3 bond. Translated from `RDGeom::compute14DistCis`
/// (`Code/Geometry/Utils.h`).
pub fn dist14_cis(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() - d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}

/// Same as [`dist14_cis`] but for a **trans** (180-degree torsion)
/// configuration. Translated from `RDGeom::compute14DistTrans`.
pub fn dist14_trans(d1: f64, d2: f64, d3: f64, ang12: f64, ang23: f64) -> f64 {
    let dx = d2 - d3 * ang23.cos() - d1 * ang12.cos();
    let dy = d3 * ang23.sin() + d1 * ang12.sin();
    (dx * dx + dy * dy).sqrt()
}

/// Propose macrocycle-aware 1-4 bound relaxations for every eligible 1-4
/// atom pair (atom1-atom4, connected through a central B-C bond that sits
/// in a ring of size >= [`MACROCYCLE_MIN`]) in `mol`.
///
/// Returns `Ok(vec![])` (a true no-op) when
/// `config.use_macrocycle_14_bounds` is false -- never silently computes
/// adjustments a caller didn't ask for.
///
/// Every returned pair is a genuine 1-4 (never a 1-2 or 1-3, checked via
/// `mol.bond_between`): `atom1` and `atom4` are not directly bonded, and
/// neither is `atom1`-`atom3` nor `atom2`-`atom4` bonded (the same
/// fused-ring false-1-4 guard RDKit's own comment attributes to
/// "sf.net bug 2835784").
pub fn macrocycle_14_bound_adjustments(
    mol: &Molecule,
    config: &TorsionKnowledgeConfig,
) -> Result<Vec<PairBoundAdjustment>, TorsionKnowledgeError> {
    let mut out = Vec::new();
    if !config.use_macrocycle_14_bounds {
        return Ok(out);
    }

    let ring_index = RingMembershipIndex::build(mol);

    for (_, bond2) in mol.bonds() {
        let (atom2, atom3) = (bond2.atom1, bond2.atom2);
        if mol.atom(atom2).wildcard || mol.atom(atom3).wildcard {
            return Err(TorsionKnowledgeError::InvalidTopology);
        }
        if mol.neighbors(atom2).count() < 2 || mol.neighbors(atom3).count() < 2 {
            continue; // no 1-4 pair possible through this bond
        }

        let bc = classify_bond(mol, &ring_index, atom2, atom3);
        let ring_size = match &bc.ring {
            super::classify::RingMembership::NotInRing => continue,
            super::classify::RingMembership::SmallRing(_) => continue,
            super::classify::RingMembership::Macrocycle(n) => *n,
            super::classify::RingMembership::FusedOrBridged { chosen_size } => *chosen_size,
        };
        if ring_size < MACROCYCLE_MIN {
            continue;
        }

        for (atom1, _) in mol.neighbors(atom2) {
            if atom1 == atom3 {
                continue;
            }
            for (atom4, _) in mol.neighbors(atom3) {
                if atom4 == atom2 || atom4 == atom1 {
                    continue;
                }
                // Genuine-1-4 guard: reject if atom1/atom4 are directly
                // bonded (that would make this a 1-3 or ring-closure short
                // circuit, not a real 1-4), or if the "diagonal" bonds
                // atom1-atom3 / atom2-atom4 exist (fused-ring false-1-4,
                // per RDKit's own cited bug).
                if mol.bond_between(atom1, atom4).is_some()
                    || mol.bond_between(atom1, atom3).is_some()
                    || mol.bond_between(atom2, atom4).is_some()
                {
                    continue;
                }

                let bl1 = ideal_bond_length(mol, atom1, atom2);
                let bl2 = ideal_bond_length(mol, atom2, atom3);
                let bl3 = ideal_bond_length(mol, atom3, atom4);
                let ba12 = ideal_bond_angle(mol, atom2);
                let ba23 = ideal_bond_angle(mol, atom3);

                let d_cis = dist14_cis(bl1, bl2, bl3, ba12, ba23);
                let d_trans = dist14_trans(bl1, bl2, bl3, ba12, ba23);

                let old_lower = d_trans - PIN_TOLERANCE_ANGSTROM;
                let old_upper = d_trans + PIN_TOLERANCE_ANGSTROM;

                let (new_lower, new_upper, rule_id, reason) = if bc.amide_like {
                    (
                        d_cis - PIN_TOLERANCE_ANGSTROM,
                        d_cis + PIN_TOLERANCE_ANGSTROM,
                        "macrocycle_14:amide_ester_pinned",
                        format!(
                            "amide/ester-like central bond in a {ring_size}-membered ring: 1-4 pair pinned to the single cis configuration ({d_cis:.3} Å +/- {PIN_TOLERANCE_ANGSTROM}), not relaxed to a band"
                        ),
                    )
                } else {
                    let lo = d_cis.min(d_trans);
                    let hi = d_cis.max(d_trans);
                    (
                        lo,
                        hi,
                        "macrocycle_14:relaxed_band",
                        format!(
                            "generic 1-4 pair in a {ring_size}-membered macrocycle: bounds relaxed to the full cis-or-trans band [{lo:.3}, {hi:.3}] Å instead of the naive single-trans-configuration pin"
                        ),
                    )
                };

                out.push(PairBoundAdjustment {
                    atom_pair: (atom1, atom4),
                    old_lower,
                    new_lower,
                    old_upper,
                    new_upper,
                    rule_id: rule_id.to_string(),
                    source: TorsionKnowledgeSource::MacrocycleAdaptation,
                    ring_size,
                    reason,
                });
            }
        }
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn disabled_flag_is_true_no_op() {
        let mol = parse("C1CCCCCCCCCCC1").unwrap(); // cyclododecane
        let config = TorsionKnowledgeConfig::default();
        let out = macrocycle_14_bound_adjustments(&mol, &config).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn cyclododecane_produces_relaxed_band_pairs() {
        let mol = parse("C1CCCCCCCCCCC1").unwrap();
        let config = TorsionKnowledgeConfig {
            use_macrocycle_14_bounds: true,
            ..TorsionKnowledgeConfig::default()
        };
        let out = macrocycle_14_bound_adjustments(&mol, &config).unwrap();
        assert!(!out.is_empty());
        for adj in &out {
            assert!(adj.new_lower <= adj.new_upper, "{adj:?}");
            assert!(adj.new_lower.is_finite() && adj.new_upper.is_finite());
            assert_eq!(adj.ring_size, 12);
        }
    }

    #[test]
    fn small_ring_produces_no_macrocycle_adjustments() {
        let mol = parse("C1CCCCC1").unwrap(); // cyclohexane, not a macrocycle
        let config = TorsionKnowledgeConfig {
            use_macrocycle_14_bounds: true,
            ..TorsionKnowledgeConfig::default()
        };
        let out = macrocycle_14_bound_adjustments(&mol, &config).unwrap();
        assert!(out.is_empty(), "{out:?}");
    }

    #[test]
    fn every_pair_is_a_genuine_14_not_12_or_13() {
        let mol = parse("C1CCCCCCCCCCC1").unwrap();
        let config = TorsionKnowledgeConfig {
            use_macrocycle_14_bounds: true,
            ..TorsionKnowledgeConfig::default()
        };
        let out = macrocycle_14_bound_adjustments(&mol, &config).unwrap();
        for adj in &out {
            let (a, b) = adj.atom_pair;
            assert!(
                mol.bond_between(a, b).is_none(),
                "1-4 pair must not be directly bonded: {adj:?}"
            );
            assert_ne!(a, b);
        }
    }

    #[test]
    fn macrolactam_bond_is_pinned_not_relaxed() {
        // A 12-membered ring lactam: one amide bond inside the ring.
        let mol = parse("O=C1CCCCCCCCCCN1").unwrap();
        let config = TorsionKnowledgeConfig {
            use_macrocycle_14_bounds: true,
            ..TorsionKnowledgeConfig::default()
        };
        let out = macrocycle_14_bound_adjustments(&mol, &config).unwrap();
        assert!(
            out.iter()
                .any(|a| a.rule_id == "macrocycle_14:amide_ester_pinned"),
            "{out:?}"
        );
    }
}
