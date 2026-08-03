//! Macrocycle 1-4 bound adjustments (Wave 2 spec §6).
//!
//! A **pure API** — this module (PR #191, "Wave 2") does not itself apply these
//! adjustments to any real bounds matrix; it was evaluated only via this module's
//! own tests and the `torsion_knowledge_v2_gap_check` example's Arm D (which applied
//! the adjustments to a **copy** of the harness's own working bounds to simulate a
//! future integration, never the production embedder).
//!
//! **Update (Wave 2 → Wave 3 Coordinator integration, `pipeline_v2.rs`):** this
//! function's output IS now wired into the real embedder in production, via a
//! minimal `pub(crate)` hook added to `distance_geometry_v2.rs`
//! (`embed_distance_geometry_v2_with_adjustments`) that `pipeline_v2::embed_pipeline_v2`
//! calls at stage 4, after converting each [`PairBoundAdjustment`] into that hook's
//! own `DistanceBoundAdjustment` type. This module itself is unchanged by that
//! integration (still a pure, standalone API) — only the paragraph above, which
//! predates the integration, needed correcting so it doesn't read as still true.
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
//! configuration, but WHICH one (cis or trans) depends on whether atom1 and
//! atom4 play the "same role" relative to the ring -- both continue the
//! macrocycle, or both are exocyclic substituents (TRANS, the dominant
//! trans-amide macrolactam conformation) -- or "different roles" -- one
//! continues the ring, the other is exocyclic (CIS, forced by the same
//! planar amide geometry). Verified against a live RDKit oracle (issue
//! found while surveying RDKit's open issues, analogous to RDKit #9266):
//! unconditionally pinning every combinatorial 1-4 pair through a
//! *tertiary* amide/ester bond to cis (this module's original behavior) is
//! a geometrically unsatisfiable constraint set whenever there's more than
//! one atom1/atom4 combination, since the real planar geometry always
//! splits the combinations into 2 cis + 2 trans, never 4-cis; for every
//! other (non-amide-like) 1-4 pair, relaxed to the full
//! `[min(cis,trans), max(cis,trans)]` band ("assume anything is possible",
//! RDKit's own inline comment in the generic branch of
//! `_setMacrocycleAllInSameRing14Bounds`).

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
                let lo = d_cis.min(d_trans);
                let hi = d_cis.max(d_trans);
                let relaxed_reason = |detail: &str| {
                    format!(
                        "{detail}: bounds relaxed to the full cis-or-trans band [{lo:.3}, {hi:.3}] Å instead of the naive single-trans-configuration pin"
                    )
                };

                // For a tertiary amide/ester central bond, atom1 and atom4
                // each have (up to) two candidates: one continues the
                // macrocycle ring, the other is an exocyclic substituent
                // (e.g. an N-methyl, or the carbonyl =O itself, which is
                // exocyclic relative to the ring). Verified against a live
                // RDKit oracle (embedded conformers across 20 seeds, plus
                // RDKit's own `_checkMacrocycleAllInSameRingAmideEster14` /
                // `BoundsMatrixBuilder.cpp`): the two ring-continuation
                // atoms are TRANS to each other (the dominant trans-amide
                // macrolactam conformation, which keeps the macrocycle
                // extended), and by the same planar-amide geometry, the two
                // EXOCYCLIC atoms are also TRANS to each other -- only a
                // ring-continuation/exocyclic CROSS pair is CIS. A single
                // "same role -> trans, different role -> cis" rule, not
                // four independently memorized cases: pinning all four
                // combinations to cis (the original, unconditional
                // behavior) was a geometrically unsatisfiable constraint
                // set for any tertiary amide/ester. Applies just as much to
                // a SECONDARY (N-H) amide's 2 real combinations (one N-side
                // candidate x two C-side candidates) -- both a cis and a
                // trans value are still needed there, so "blanket cis is
                // harmless for secondary amides" (an earlier, incorrect
                // claim about this same code) does not hold either.
                //
                // "Continues the ring" must mean the SAME macrocycle as the
                // central bond, not merely "is in some ring of any kind" --
                // ring SIZE alone can't tell an N-cyclohexyl substituent's
                // own 6-ring apart from a coincidentally same-sized
                // macrocycle segment. Resolved via ring IDENTITY
                // (`ring_ids_for`, index into the real ring list) instead of
                // `ring_sizes_for`: a side bond only "continues the ring" if
                // its own ring-ID set intersects the central bond's
                // eligible-macrocycle ring-ID set.
                let central_ring_ids: Vec<usize> = ring_index
                    .ring_ids_for(atom2, atom3)
                    .iter()
                    .copied()
                    .filter(|&id| ring_index.rings()[id].len() >= MACROCYCLE_MIN)
                    .collect();

                let (new_lower, new_upper, rule_id, reason) = match (
                    bc.amide_like,
                    central_ring_ids.as_slice(),
                ) {
                    (true, [central_ring_id]) => {
                        let continues_central_macrocycle = |x: AtomIdx, y: AtomIdx| {
                            ring_index.ring_ids_for(x, y).contains(central_ring_id)
                        };
                        let atom1_continues_ring = continues_central_macrocycle(atom2, atom1);
                        let atom4_continues_ring = continues_central_macrocycle(atom3, atom4);
                        let same_role = atom1_continues_ring == atom4_continues_ring;
                        let (d, config_name) = if same_role {
                            (d_trans, "trans")
                        } else {
                            (d_cis, "cis")
                        };
                        (
                            d - PIN_TOLERANCE_ANGSTROM,
                            d + PIN_TOLERANCE_ANGSTROM,
                            "macrocycle_14:amide_ester_pinned",
                            format!(
                                "amide/ester-like central bond in a {ring_size}-membered ring: 1-4 pair pinned to the single {config_name} configuration ({d:.3} Å +/- {PIN_TOLERANCE_ANGSTROM}) -- {} ring-continuation ({}/{}), not relaxed to a band",
                                if same_role {
                                    "both atoms share the same"
                                } else {
                                    "atoms differ in"
                                },
                                atom1_continues_ring,
                                atom4_continues_ring,
                            ),
                        )
                    }
                    (true, ids) => {
                        // Fail closed: the central bond sits in more than
                        // one eligible macrocycle (a fused/bridged system),
                        // so which ring's path the cis/trans split should
                        // follow is not well-defined from ring membership
                        // alone. Abstain from the pin rather than guess.
                        (
                            lo,
                            hi,
                            "macrocycle_14:relaxed_band",
                            relaxed_reason(&format!(
                                "amide/ester-like central bond belongs to {} distinct eligible \
                                 macrocycles ({ring_size}-membered chosen as smallest) -- role \
                                 (ring-continuation vs exocyclic) is ambiguous, abstained from \
                                 pinning",
                                ids.len()
                            )),
                        )
                    }
                    (false, _) => (
                        lo,
                        hi,
                        "macrocycle_14:relaxed_band",
                        relaxed_reason(&format!(
                            "generic 1-4 pair in a {ring_size}-membered macrocycle"
                        )),
                    ),
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

    /// Issue found while surveying RDKit's open issues (analogous to RDKit
    /// #9266, macrocyclic tertiary amides embedding twisted): a TERTIARY
    /// amide/ester central bond has up to 4 combinatorial (atom1, atom4)
    /// 1-4 pairs (2 candidates per side). Pinning all 4 to the same cis
    /// configuration (the previous, unconditional behavior) is geometrically
    /// unsatisfiable -- verified against a live RDKit oracle (embedded
    /// conformers across 20 seeds, converged unambiguously): the real planar
    /// amide geometry always splits the 4 combinations 2-cis + 2-trans,
    /// specifically along "does this atom continue the macrocycle ring, or
    /// is it an exocyclic substituent" -- same role (both ring-continuation,
    /// or both exocyclic) -> trans; different role -> cis.
    #[test]
    fn tertiary_amide_1_4_pairs_split_cis_and_trans_by_ring_role() {
        // N-methyl 13-membered macrolactam: amide N has a ring-continuation
        // neighbor (the preceding ring carbon) and an exocyclic neighbor
        // (the N-methyl); the carbonyl C has a ring-continuation neighbor
        // (the following ring carbon) and an exocyclic neighbor (=O).
        let mol = parse("O=C1CCCCCCCCCCN1C").unwrap();
        // Atom indices per this exact SMILES's left-to-right parse order:
        // 0=O, 1=C(carbonyl), 2..11=ring CH2 x10, 12=N, 13=C(methyl).
        let o_idx = AtomIdx(0);
        let carbonyl_c_idx = AtomIdx(1);
        let ring_c_idx = AtomIdx(2); // carbonyl C's ring-continuation neighbor
        let ring_n_idx = AtomIdx(11); // amide N's ring-continuation neighbor
        let amide_n_idx = AtomIdx(12);
        let methyl_idx = AtomIdx(13); // amide N's exocyclic neighbor

        let config = TorsionKnowledgeConfig {
            use_macrocycle_14_bounds: true,
            ..TorsionKnowledgeConfig::default()
        };
        let out = macrocycle_14_bound_adjustments(&mol, &config).unwrap();

        let bl2 = ideal_bond_length(&mol, amide_n_idx, carbonyl_c_idx);
        let ba_n = ideal_bond_angle(&mol, amide_n_idx);
        let ba_c = ideal_bond_angle(&mol, carbonyl_c_idx);
        let expected_config_for = |a1: AtomIdx, a4: AtomIdx, same_role: bool| {
            let bl1 = ideal_bond_length(&mol, a1, amide_n_idx);
            let bl3 = ideal_bond_length(&mol, carbonyl_c_idx, a4);
            if same_role {
                dist14_trans(bl1, bl2, bl3, ba_n, ba_c)
            } else {
                dist14_cis(bl1, bl2, bl3, ba_n, ba_c)
            }
        };

        // (ring, ring) -> same role -> trans.
        let ring_ring = expected_bound(&out, ring_n_idx, ring_c_idx);
        assert_close(
            ring_ring,
            expected_config_for(ring_n_idx, ring_c_idx, true),
            "ring_N--ring_C (both ring-continuation) must be pinned TRANS",
        );

        // (exocyclic, exocyclic) -> same role -> trans.
        let methyl_o = expected_bound(&out, methyl_idx, o_idx);
        assert_close(
            methyl_o,
            expected_config_for(methyl_idx, o_idx, true),
            "methyl_N--O (both exocyclic) must be pinned TRANS",
        );

        // (ring, exocyclic) -> different role -> cis.
        let ring_o = expected_bound(&out, ring_n_idx, o_idx);
        assert_close(
            ring_o,
            expected_config_for(ring_n_idx, o_idx, false),
            "ring_N--O (ring-continuation vs exocyclic) must be pinned CIS",
        );

        // (exocyclic, ring) -> different role -> cis.
        let methyl_ring_c = expected_bound(&out, methyl_idx, ring_c_idx);
        assert_close(
            methyl_ring_c,
            expected_config_for(methyl_idx, ring_c_idx, false),
            "methyl_N--ring_C (exocyclic vs ring-continuation) must be pinned CIS",
        );

        // The two same-role pairs and the two different-role pairs must NOT
        // collapse to the same pinned distance -- this is the literal
        // defect being fixed (previously all 4 pinned to the same cis
        // value).
        assert!(
            (ring_ring - ring_o).abs() > 2.0 * PIN_TOLERANCE_ANGSTROM,
            "trans-pinned and cis-pinned distances must be clearly distinct: \
             ring_ring={ring_ring:.3} ring_o={ring_o:.3}"
        );
    }

    /// Helper for the test above: find the pinned midpoint distance
    /// `(new_lower + new_upper) / 2` for a given unordered atom pair.
    fn expected_bound(out: &[PairBoundAdjustment], a: AtomIdx, b: AtomIdx) -> f64 {
        let adj = out
            .iter()
            .find(|adj| {
                let (x, y) = adj.atom_pair;
                (x == a && y == b) || (x == b && y == a)
            })
            .unwrap_or_else(|| panic!("no adjustment found for pair ({a:?}, {b:?}): {out:?}"));
        assert_eq!(
            adj.rule_id, "macrocycle_14:amide_ester_pinned",
            "expected the amide-pinned rule for {a:?}--{b:?}: {adj:?}"
        );
        (adj.new_lower + adj.new_upper) / 2.0
    }

    fn assert_close(actual: f64, expected: f64, msg: &str) {
        assert!(
            (actual - expected).abs() < 1e-6,
            "{msg}: actual={actual:.4} expected={expected:.4}"
        );
    }

    #[test]
    fn cyclododecane_adjustments_unaffected_by_amide_role_logic() {
        // No amide/ester bond at all -- the ring-role diagonal-split logic
        // must never be reached, so every adjustment must still be the
        // plain relaxed_band rule with its original [min(cis,trans),
        // max(cis,trans)] bounds, unchanged by this fix.
        let mol = parse("C1CCCCCCCCCCC1").unwrap();
        let config = TorsionKnowledgeConfig {
            use_macrocycle_14_bounds: true,
            ..TorsionKnowledgeConfig::default()
        };
        let out = macrocycle_14_bound_adjustments(&mol, &config).unwrap();
        assert!(!out.is_empty());
        for adj in &out {
            assert_eq!(adj.rule_id, "macrocycle_14:relaxed_band", "{adj:?}");
        }
    }

    /// Issue found in review: the previous role check
    /// (`!ring_sizes_for(...).is_empty()`) asks "is this bond in ANY ring,"
    /// not "does this bond continue the SAME macrocycle as the central
    /// bond" -- ring SIZE alone can't tell a coincidentally same-sized,
    /// unrelated ring apart from a genuine continuation of the central
    /// bond's own macrocycle. A plain N-cyclohexyl/N-piperidyl substituent
    /// turns out NOT to trigger this in practice (verified empirically,
    /// see below) since a simple pendant ring's attachment bond isn't
    /// itself a ring bond at all -- but a genuinely shared amide bond (the
    /// central N-C bond itself belonging to two independent macrocycles, a
    /// "theta graph" between the two atoms) does, and ring size cannot
    /// distinguish the ambiguous case from the unambiguous one either way.
    #[test]
    fn ncyclohexyl_macrolactam_substituent_attachment_is_not_a_ring_bond() {
        // Empirical check backing the note above: a plain N-cyclohexyl
        // substituent's own attachment bond (amide N -- cyclohexyl ipso
        // carbon) is NOT itself a ring bond (only the 6 internal cyclohexyl
        // C-C bonds are), so both the old ring-size check and the new
        // ring-ID check agree here -- this is a no-regression sanity check,
        // not a reproduction of the reviewer's stated bug.
        let mol = parse("O=C1CCCCCCCCCCN1C1CCCCC1").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        // atom 12 = amide N, atom 13 = cyclohexyl ipso carbon (SMILES parse order).
        assert!(
            rings.ring_ids_for(AtomIdx(12), AtomIdx(13)).is_empty(),
            "the N-to-cyclohexyl attachment bond must not itself be a ring bond"
        );

        let config = TorsionKnowledgeConfig {
            use_macrocycle_14_bounds: true,
            ..TorsionKnowledgeConfig::default()
        };
        let out = macrocycle_14_bound_adjustments(&mol, &config).unwrap();
        assert!(
            out.iter()
                .any(|a| a.rule_id == "macrocycle_14:amide_ester_pinned"),
            "the amide pin must still fire for a simple N-cyclohexyl substituent: {out:?}"
        );
    }

    /// The genuine ring-ID-vs-ring-size divergence: an amide bond that is
    /// itself the shared edge of two independent, both macrocycle-sized
    /// rings (a "theta graph" between the amide N and the carbonyl C).
    /// Ring SIZE alone can't even ask "which ring" here -- there's no
    /// single answer for "does atom1 continue THE macrocycle," since there
    /// are two. Correct behavior is to abstain from the pin (fail closed)
    /// rather than guess which ring's role assignment to trust.
    #[test]
    fn amide_bond_shared_by_two_macrocycles_abstains_from_pin() {
        let mol = parse("O=C1C(CCCCCCC2)CCCCCCCCN12").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        // atom 1 = carbonyl C, atom 18 = amide N (SMILES parse order);
        // verified via direct probe: this bond belongs to both a
        // 10-membered and an 11-membered ring in the augmented ring set.
        let ids = rings.ring_ids_for(AtomIdx(1), AtomIdx(18));
        assert_eq!(
            ids.len(),
            2,
            "fixture must genuinely put the amide bond in 2 eligible macrocycles: {ids:?}"
        );

        let config = TorsionKnowledgeConfig {
            use_macrocycle_14_bounds: true,
            ..TorsionKnowledgeConfig::default()
        };
        let out = macrocycle_14_bound_adjustments(&mol, &config).unwrap();
        let amide_bond_pairs: Vec<_> = out
            .iter()
            .filter(|a| {
                let (x, y) = a.atom_pair;
                // A 1-4 pair through the (1,18) central bond has atom1 or
                // atom4 adjacent to atom 1 or atom 18 -- simplest robust
                // check: at least one endpoint is a direct neighbor of
                // atom 1 or atom 18.
                mol.bond_between(x, AtomIdx(1)).is_some()
                    || mol.bond_between(x, AtomIdx(18)).is_some()
                    || mol.bond_between(y, AtomIdx(1)).is_some()
                    || mol.bond_between(y, AtomIdx(18)).is_some()
            })
            .collect();
        assert!(!amide_bond_pairs.is_empty(), "{out:?}");
        for adj in &amide_bond_pairs {
            assert_eq!(
                adj.rule_id, "macrocycle_14:relaxed_band",
                "an amide bond shared by 2 eligible macrocycles must abstain from \
                 pinning (ambiguous role), not guess: {adj:?}"
            );
        }
    }

    /// A real, ordinary fused bicyclic macrolactam (bridge far from the
    /// amide bond) must NOT trip the ambiguity-abstain path above -- the
    /// amide bond here belongs to exactly one eligible macrocycle even
    /// though the molecule as a whole is bicyclic. Guards against the
    /// abstain logic over-firing on the exact molecule class (fused/bridged
    /// macrocyclic systems) this whole fix is about.
    #[test]
    fn ordinary_fused_bicyclic_macrolactam_amide_bond_still_pins() {
        let mol = parse("O=C1CCC2CCCCC2CCCN1C").unwrap();
        let rings = RingMembershipIndex::build(&mol);
        // atom 1 = carbonyl C, atom 13 = amide N (SMILES parse order);
        // verified via direct probe: exactly one eligible macrocycle (the
        // 9-membered ring) contains this bond, despite the fused cyclohexane.
        let ids = rings.ring_ids_for(AtomIdx(1), AtomIdx(13));
        assert_eq!(
            ids.len(),
            1,
            "an ordinary fused bicyclic macrolactam's amide bond must belong to \
             exactly one eligible macrocycle: {ids:?}"
        );

        let config = TorsionKnowledgeConfig {
            use_macrocycle_14_bounds: true,
            ..TorsionKnowledgeConfig::default()
        };
        let out = macrocycle_14_bound_adjustments(&mol, &config).unwrap();
        assert!(
            out.iter()
                .any(|a| a.rule_id == "macrocycle_14:amide_ester_pinned"),
            "the amide pin must still fire for an ordinary fused bicyclic macrolactam: {out:?}"
        );
    }

    /// Corrected claim: the original PR body said blanket-cis pinning was
    /// "harmless" for a secondary (N-H) amide. It is not -- even a
    /// secondary amide has 2 real (atom1, atom4) combinations (1 N-side
    /// candidate x 2 C-side candidates: ring-continuation and the carbonyl
    /// =O), and by the same planar-amide geometry those 2 combinations
    /// still split into one cis and one trans pair, not a single value.
    #[test]
    fn secondary_macrolactam_still_needs_both_cis_and_trans() {
        // O=C1CCCCCCCCCCN1: atom0=O, atom1=carbonyl C, atom2..11=10 ring
        // CH2's, atom12=N-H (no substituent). N's only non-carbonyl
        // neighbor is atom11 (ring-continuation); carbonylC's two
        // candidates are atom2 (ring-continuation) and atom0 (=O, exocyclic).
        let mol = parse("O=C1CCCCCCCCCCN1").unwrap();
        let o_idx = AtomIdx(0);
        let carbonyl_c_idx = AtomIdx(1);
        let ring_c_idx = AtomIdx(2);
        let ring_n_idx = AtomIdx(11);
        let amide_n_idx = AtomIdx(12);

        let config = TorsionKnowledgeConfig {
            use_macrocycle_14_bounds: true,
            ..TorsionKnowledgeConfig::default()
        };
        let out = macrocycle_14_bound_adjustments(&mol, &config).unwrap();
        assert_eq!(
            out.iter()
                .filter(|a| a.rule_id == "macrocycle_14:amide_ester_pinned")
                .count(),
            2,
            "a secondary amide has exactly 2 combinatorial 1-4 pairs: {out:?}"
        );

        let bl2 = ideal_bond_length(&mol, amide_n_idx, carbonyl_c_idx);
        let ba_n = ideal_bond_angle(&mol, amide_n_idx);
        let ba_c = ideal_bond_angle(&mol, carbonyl_c_idx);
        let expected_config_for = |a1: AtomIdx, a4: AtomIdx, same_role: bool| {
            let bl1 = ideal_bond_length(&mol, a1, amide_n_idx);
            let bl3 = ideal_bond_length(&mol, carbonyl_c_idx, a4);
            if same_role {
                dist14_trans(bl1, bl2, bl3, ba_n, ba_c)
            } else {
                dist14_cis(bl1, bl2, bl3, ba_n, ba_c)
            }
        };

        // (ring, ring) -> same role -> trans.
        let ring_ring = expected_bound(&out, ring_n_idx, ring_c_idx);
        assert_close(
            ring_ring,
            expected_config_for(ring_n_idx, ring_c_idx, true),
            "ring_N--ring_C (both ring-continuation) must be pinned TRANS",
        );

        // (ring, =O) -> different role -> cis.
        let ring_o = expected_bound(&out, ring_n_idx, o_idx);
        assert_close(
            ring_o,
            expected_config_for(ring_n_idx, o_idx, false),
            "ring_N--O (ring-continuation vs exocyclic) must be pinned CIS",
        );

        // The corrected claim this test pins: even a secondary (N-H) amide
        // needs BOTH a cis and a trans value, not one blanket value --
        // contradicting an earlier, incorrect "harmless for secondary
        // amides" claim about this same code.
        assert!(
            (ring_ring - ring_o).abs() > 2.0 * PIN_TOLERANCE_ANGSTROM,
            "trans-pinned and cis-pinned distances must be clearly distinct: \
             ring_ring={ring_ring:.3} ring_o={ring_o:.3}"
        );
    }
}
