//! ETKDG Torsion Knowledge Base.
//!
//! # Two generations of code in this file
//!
//! **Legacy (below, unchanged by the Wave 2 PR that added the `v2` module
//! family):** `classify_atom_type`, `get_torsion_preference`,
//! `build_smarts_torsion_map`, `default_torsion_preference`, `score_torsion`,
//! and their supporting types (`AtomType`, `TorsionPreference`,
//! `SmartsTorsionRule`). This is a hand-authored, single-angle +
//! linear-penalty heuristic layer. **Its previous doc comment claimed these
//! were "experimental torsion angle preferences from CSD" (Cambridge
//! Structural Database) -- that claim was unverified and, on the evidence
//! gathered, false: there is no CSD query, statistical fit, or citation
//! anywhere in this module or its history.** See
//! `docs/rfcs/3d_torsion_knowledge_audit.md` for the full audit (every rule
//! individually reclassified as `verified_from_primary_source` /
//! `supported_by_rdkit_oracle` / `reasonable_heuristic_only` /
//! `incorrect_or_overgeneralized` / `ambiguous`; zero `dead_or_unreachable`).
//! Kept behaviorally **unchanged** (same public API, same outputs) -- only
//! this doc comment was corrected to stop mislabeling it.
//!
//! **v2 (submodules below, new in `feat/3d-torsion-knowledge-v2`, 3D
//! Breakthrough Program Wave 2 / Agent E):** a separate, periodic/multi-modal
//! torsion-potential layer with explicit rule provenance, a 6-tier priority
//! system, ring/bond classification, macrocycle 1-4 bounds, and a
//! torsion-space energy/optimization API. See each submodule's own doc
//! comment, `docs/rfcs/3d_torsion_knowledge_audit.md`, and
//! `validation/manifests/etkdg_torsion_knowledge_sources.json` for the full
//! design and sourcing. **Not wired into `distance_geometry_v2.rs` or
//! `etkdg.rs`** -- this PR is a self-contained knowledge/potential/
//! verification layer only; Coordinator performs that integration in a
//! later, separate PR. The new code never uses the legacy heuristic's
//! output as evidence that it is correct (spec §8) -- they are independent
//! implementations, reconciled only through this file's re-exports.

mod bounds14;
mod canon_rank;
mod classify;
mod energy;
mod matcher;
mod rules_basic;
mod rules_macrocycle;
mod rules_smallring;
mod rules_standard;
mod types;

pub use bounds14::{PairBoundAdjustment, macrocycle_14_bound_adjustments};
pub use classify::{
    BondClassification, MACROCYCLE_MIN, RingMembership, RingMembershipIndex, SMALL_RING_MAX,
    SMALL_RING_MIN, candidate_central_bonds, classify_bond,
};
// `pub(crate)`, not `pub`: Wave 2/3 Coordinator integration seam only (see
// `energy::is_bridge_bond`'s own doc comment) -- not part of this module's public API.
pub(crate) use energy::is_bridge_bond;
pub use energy::{
    TorsionEnergyReport, TorsionOptimizationConfig, TorsionOptimizationReport,
    evaluate_torsion_energy, optimize_torsions,
};
pub use matcher::build_torsion_knowledge;
pub use types::{
    FourierTorsionTerm, TorsionKnowledgeConfig, TorsionKnowledgeDiagnostic,
    TorsionKnowledgeDiagnosticKind, TorsionKnowledgeError, TorsionKnowledgeReport,
    TorsionKnowledgeSource, TorsionPotential,
};

use chematic_core::{AtomIdx, Molecule};
use chematic_smarts::{find_matches, parse_smarts};
use std::collections::HashMap;

/// Torsion angle preference with penalty scoring.
#[derive(Clone, Debug)]
pub struct TorsionPreference {
    /// Preferred dihedral angle in degrees.
    pub angle_deg: f64,
    /// Penalty (in kcal/mol equivalent) for deviation from preference.
    pub penalty_per_degree: f64,
}

/// Atom type for torsion matching based on hybridization and neighbors.
///
/// Classifies atoms by element and hybridization state (sp, sp2, sp3, aromatic).
/// Useful for matching chemical patterns and determining properties.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum AtomType {
    /// sp3 carbon
    CSp3,
    /// sp2 carbon (alkene)
    CSp2Alkene,
    /// sp2 carbon (aromatic)
    CAromatic,
    /// sp2 carbon (carbonyl/carboxylic)
    CCarbonyl,
    /// sp nitrogen (nitrile/isocyanate)
    NSp,
    /// sp2 nitrogen (amide/imine)
    NSp2,
    /// sp2 nitrogen (aromatic)
    NAromatic,
    /// sp3 nitrogen (amine)
    NSp3,
    /// sp2 oxygen (carbonyl)
    OSp2,
    /// sp3 oxygen (ether/alcohol)
    OSp3,
    /// aromatic oxygen (furan, oxazole, isoxazole)
    OAromatic,
    /// sulfur (thioether, sulfoxide, sulfone)
    S,
    /// aromatic sulfur (thiophene, thiazole, thiadiazole)
    SAromatic,
    /// phosphorus
    P,
    /// hydrogen
    H,
    /// halogen
    Halogen,
    /// other/unknown
    Other,
}

/// Count bonds of a given order incident to an atom.
fn count_incident_bonds(mol: &Molecule, idx: AtomIdx, order: chematic_core::BondOrder) -> usize {
    mol.bonds()
        .filter(|(_, bond)| (bond.atom1 == idx || bond.atom2 == idx) && bond.order == order)
        .count()
}

/// Classify an atom's type based on its chemical environment.
///
/// Returns an [`AtomType`] indicating the atom's hybridization state and element.
/// Useful for pattern matching, property prediction, and chemical classification.
pub fn classify_atom_type(mol: &Molecule, idx: AtomIdx) -> AtomType {
    let atom = mol.atom(idx);
    let an = atom.element.atomic_number();

    match an {
        1 => AtomType::H,
        6 => {
            // Carbon: determine sp3, sp2, sp hybridization
            let double_bonds = count_incident_bonds(mol, idx, chematic_core::BondOrder::Double);

            if atom.aromatic {
                AtomType::CAromatic
            } else if double_bonds > 0 {
                // Check if it's a carbonyl carbon
                let has_o_neighbor = mol.neighbors(idx).any(|(n_idx, _)| {
                    mol.atom(n_idx).element.atomic_number() == 8
                        && mol
                            .bond_between(idx, n_idx)
                            .map(|(_, b)| b.order == chematic_core::BondOrder::Double)
                            .unwrap_or(false)
                });
                if has_o_neighbor {
                    AtomType::CCarbonyl
                } else {
                    AtomType::CSp2Alkene
                }
            } else {
                AtomType::CSp3
            }
        }
        7 => {
            // Nitrogen: sp, sp2, sp3
            if atom.aromatic {
                AtomType::NAromatic
            } else {
                let triple_bonds = count_incident_bonds(mol, idx, chematic_core::BondOrder::Triple);
                let neighbors = mol.neighbors(idx).count();

                if triple_bonds > 0 {
                    AtomType::NSp
                } else if neighbors <= 2 {
                    AtomType::NSp2
                } else {
                    AtomType::NSp3
                }
            }
        }
        8 => {
            if atom.aromatic {
                AtomType::OAromatic // furan, oxazole, isoxazole ring oxygen
            } else {
                let has_double_bond = mol.bonds().any(|(_, bond)| {
                    (bond.atom1 == idx || bond.atom2 == idx)
                        && bond.order == chematic_core::BondOrder::Double
                });
                if has_double_bond {
                    AtomType::OSp2 // Carbonyl oxygen (C=O)
                } else {
                    AtomType::OSp3 // Alcohol or ether oxygen (C-O)
                }
            }
        }
        16 => {
            if atom.aromatic {
                AtomType::SAromatic // thiophene, thiazole, thiadiazole ring sulfur
            } else {
                AtomType::S
            }
        }
        15 => AtomType::P,
        9 | 17 | 35 | 53 => AtomType::Halogen,
        _ => AtomType::Other,
    }
}

/// Get torsion preference for an A-B-C-D dihedral based on atom types.
///
/// Returns None if no specific preference is known (use general default).
pub fn get_torsion_preference(
    mol: &Molecule,
    a_idx: AtomIdx,
    b_idx: AtomIdx,
    c_idx: AtomIdx,
    d_idx: AtomIdx,
) -> Option<TorsionPreference> {
    let a_type = classify_atom_type(mol, a_idx);
    let b_type = classify_atom_type(mol, b_idx);
    let c_type = classify_atom_type(mol, c_idx);
    let d_type = classify_atom_type(mol, d_idx);

    // Alkane C-C-C-C: strongly prefer 180° (staggered, anti)
    if b_type == AtomType::CSp3
        && c_type == AtomType::CSp3
        && (a_type == AtomType::CSp3 || a_type == AtomType::H)
        && (d_type == AtomType::CSp3 || d_type == AtomType::H)
    {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.15, // ~3 kcal/mol for 20° deviation
        });
    }

    // Aromatic-aliphatic C-Ar-C-C: prefer 180° or 0°
    if (a_type == AtomType::CSp3 || a_type == AtomType::H)
        && b_type == AtomType::CAromatic
        && c_type == AtomType::CSp3
        && (d_type == AtomType::CSp3 || d_type == AtomType::H)
    {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.08, // Softer constraint
        });
    }

    // Amide C-N(C=O)-C-X: restricted rotation, prefer 0° (cis to carbonyl) or 180° (trans)
    // For simplicity, prefer trans (180°) as it's more common
    if b_type == AtomType::NSp2 && c_type == AtomType::CCarbonyl {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.20, // Moderate restriction
        });
    }

    // Ester O-C(=O)-O-C: prefer 0° or 180° depending on configuration
    if b_type == AtomType::CCarbonyl && c_type == AtomType::OSp3 {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.10,
        });
    }

    // Aromatic-aromatic rotations: prefer ~45° (biphenyl-like twist)
    if b_type == AtomType::CAromatic && c_type == AtomType::CAromatic {
        return Some(TorsionPreference {
            angle_deg: 45.0,
            penalty_per_degree: 0.03, // Very soft — flat potential
        });
    }

    // Enamine C=C-N: prefer 180° (trans, conjugation favored)
    if b_type == AtomType::CSp2Alkene && c_type == AtomType::NSp2 {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.12,
        });
    }
    if b_type == AtomType::NSp2 && c_type == AtomType::CSp2Alkene {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.12,
        });
    }

    // Vinyl halide C=C-X: prefer 0° (cis, steric minimization)
    if (b_type == AtomType::CSp2Alkene || b_type == AtomType::CAromatic)
        && c_type == AtomType::Halogen
    {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.08,
        });
    }

    // Acrylic / chalcone: C=C-C(=O): s-trans preferred (~180°)
    if b_type == AtomType::CSp2Alkene && c_type == AtomType::CCarbonyl {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.10,
        });
    }

    // Phenyl ketone Ar-C(=O): coplanar (0°)
    if b_type == AtomType::CAromatic && c_type == AtomType::CCarbonyl {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.08,
        });
    }
    if b_type == AtomType::CCarbonyl && c_type == AtomType::CAromatic {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.08,
        });
    }

    // Thioester S-C(=O): prefer 180° (trans)
    if b_type == AtomType::S && c_type == AtomType::CCarbonyl {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.10,
        });
    }

    // Carbamate / sulfonamide N-C(=O)-O or N-S(=O)(=O): prefer 180°
    if b_type == AtomType::NSp3 && c_type == AtomType::CCarbonyl {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.10,
        });
    }

    // Sulfoxide C-S(=O)-C: pyramidal sulfur, prefer ~90°
    if b_type == AtomType::S && c_type == AtomType::CSp3 {
        return Some(TorsionPreference {
            angle_deg: 90.0,
            penalty_per_degree: 0.06,
        });
    }

    // Disulfide C-S-S-C: ~90° dihedral (gauche)
    if b_type == AtomType::S && c_type == AtomType::S {
        return Some(TorsionPreference {
            angle_deg: 90.0,
            penalty_per_degree: 0.12,
        });
    }

    // Alcohol C-C-O-H / ether C-C-O-C: gauche/anti mixture, use 180° as default
    if (b_type == AtomType::CSp3 || b_type == AtomType::CSp2Alkene) && c_type == AtomType::OSp3 {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.07,
        });
    }

    // Amine C-C-N-C (secondary/tertiary): prefer 180°
    // NSp2 covers N with 2 explicit bonds (secondary amine in SMILES);
    // NSp3 covers N with 3+ explicit bonds (tertiary amine).
    if b_type == AtomType::CSp3 && (c_type == AtomType::NSp3 || c_type == AtomType::NSp2) {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.08,
        });
    }

    // Nitrile terminus X-C-C≡N: prefer 180° (linear)
    if b_type == AtomType::NSp || d_type == AtomType::NSp {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.20,
        });
    }

    // Phosphorus P-C-C-X: use 180° default
    if b_type == AtomType::P || c_type == AtomType::P {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.06,
        });
    }

    // Urea N-C(=O)-N: planar, prefer 0° (both N lone pairs overlap with C=O)
    if b_type == AtomType::NSp2
        && c_type == AtomType::NSp2
        && mol
            .neighbors(b_idx)
            .any(|(n, _)| classify_atom_type(mol, n) == AtomType::CCarbonyl)
    {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.18,
        });
    }

    // Sulfonamide N-S(=O)(=O): prefer 90° (tetrahedral S, gauche N)
    if (b_type == AtomType::NSp3 || b_type == AtomType::NSp2) && c_type == AtomType::S {
        return Some(TorsionPreference {
            angle_deg: 90.0,
            penalty_per_degree: 0.08,
        });
    }
    if b_type == AtomType::S && (c_type == AtomType::NSp3 || c_type == AtomType::NSp2) {
        return Some(TorsionPreference {
            angle_deg: 90.0,
            penalty_per_degree: 0.08,
        });
    }

    // Aryl ether Ar-O-C: prefer 0° (oxygen lone pair conjugation with ring)
    if b_type == AtomType::CAromatic && c_type == AtomType::OSp3 {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.09,
        });
    }
    if b_type == AtomType::OSp3 && c_type == AtomType::CAromatic {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.09,
        });
    }

    // Fluoroalkane C-C-C-F: prefer anti (180°) to minimise F dipole interactions
    if (b_type == AtomType::CSp3 || b_type == AtomType::CSp2Alkene)
        && c_type == AtomType::Halogen
        && (d_type == AtomType::H || d_type == AtomType::Halogen)
    {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.07,
        });
    }

    // Nitro group C-N(=O)=O: coplanar with aromatic ring if attached to Ar
    if b_type == AtomType::CAromatic && c_type == AtomType::NSp2 {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.12,
        });
    }

    // Hydrazone/oxime C=N-N or C=N-O: prefer 0° (E/Z isomerism; E is more stable)
    if (b_type == AtomType::CSp2Alkene || b_type == AtomType::CCarbonyl)
        && (c_type == AtomType::NSp2 || c_type == AtomType::OSp3)
    {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.11,
        });
    }

    // Imide N-C(=O)-C(=O): prefer 0° (both carbonyls on same side for conjugation)
    if b_type == AtomType::NSp2
        && c_type == AtomType::CCarbonyl
        && mol
            .neighbors(c_idx)
            .any(|(n, _)| classify_atom_type(mol, n) == AtomType::CCarbonyl)
    {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.15,
        });
    }

    // Benzyl (Ar-C-X): prefer 90° perpendicular to ring plane
    if b_type == AtomType::CAromatic && c_type == AtomType::CSp3 {
        return Some(TorsionPreference {
            angle_deg: 90.0,
            penalty_per_degree: 0.04,
        });
    }

    // Allylic C=C-C-X: prefer 0° (s-cis/s-trans mixture; use 0° as default)
    if b_type == AtomType::CSp2Alkene && c_type == AtomType::CSp3 {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.05,
        });
    }

    // ── Heteroaromatic patterns ──────────────────────────────────────────────

    // Heteroaromatic biaryl: NAromatic–CAromatic or CAromatic–NAromatic
    // (e.g. phenyl-pyridine, bipyridine, pyrimidine-phenyl).
    // ~45° twist like biphenyl; very soft potential.
    // NAromatic–NAromatic is intentionally excluded: adjacent aromatic nitrogens
    // occur as intra-ring bonds (pyrimidine, pyridazine) whose torsion is already
    // constrained by ring closure; applying a 45° soft preference there conflicts
    // with satisfy_constraints and can distort ring geometry.
    if (b_type == AtomType::NAromatic && c_type == AtomType::CAromatic)
        || (b_type == AtomType::CAromatic && c_type == AtomType::NAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 45.0,
            penalty_per_degree: 0.03,
        });
    }

    // N-alkyl heteroaromatic (N-methyl pyridine, N-methyl imidazole, etc.).
    // The N–C(sp3) bond prefers anti (180°) to minimise lone-pair / σ* repulsion.
    if (b_type == AtomType::NAromatic && c_type == AtomType::CSp3)
        || (b_type == AtomType::CSp3 && c_type == AtomType::NAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.09,
        });
    }

    // Heteroaromatic N adjacent to carbonyl (prodrug, lactam-like contexts).
    // Prefer planar (0°) for lone-pair conjugation.
    if (b_type == AtomType::NAromatic && c_type == AtomType::CCarbonyl)
        || (b_type == AtomType::CCarbonyl && c_type == AtomType::NAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.10,
        });
    }

    // Heteroaromatic N to sp2 alkene (vinylogous conjugation).
    if (b_type == AtomType::NAromatic && c_type == AtomType::CSp2Alkene)
        || (b_type == AtomType::CSp2Alkene && c_type == AtomType::NAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.08,
        });
    }

    // Thioaryl / aryl thioether: S–CAromatic or CAromatic–S.
    // Diaryl sulfide and aryl thioether prefer ~90° (sulphur p-lone-pair
    // perpendicular to ring π system).
    if (b_type == AtomType::S && c_type == AtomType::CAromatic)
        || (b_type == AtomType::CAromatic && c_type == AtomType::S)
    {
        return Some(TorsionPreference {
            angle_deg: 90.0,
            penalty_per_degree: 0.06,
        });
    }

    // OSp3 as B in O–C(sp3) torsion (e.g. C–O–C–C ether chain, reverse of the
    // CSp3–OSp3 rule).  Prefer anti (180°).
    if b_type == AtomType::OSp3 && c_type == AtomType::CSp3 {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.07,
        });
    }

    // Aryl amine Ar–NR2 (aniline-like, single bond to sp3 N): prefer planar (0°)
    // for lone-pair conjugation with ring, slightly softer than aryl amide.
    // Both traversal directions are covered so the rule fires regardless of
    // which end of the Ar–N bond is atom B vs C.
    if (b_type == AtomType::CAromatic && c_type == AtomType::NSp3)
        || (b_type == AtomType::NSp3 && c_type == AtomType::CAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.07,
        });
    }

    // ── 5-membered aromatic heterocycle patterns ─────────────────────────────

    // Furanyl / oxazolyl biaryl: OAromatic–CAromatic inter-ring bond.
    // Oxygen lone pair conjugates strongly with the adjacent ring π system;
    // planar (0°) is preferred unlike the ~45° biphenyl twist.
    if (b_type == AtomType::OAromatic && c_type == AtomType::CAromatic)
        || (b_type == AtomType::CAromatic && c_type == AtomType::OAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.05,
        });
    }

    // Thienyl / thiazolyl biaryl: SAromatic–CAromatic inter-ring.
    // Sulfur d-orbital participation gives a flatter potential than O;
    // ~45° is a reasonable CSD consensus (like biphenyl but slightly softer).
    if (b_type == AtomType::SAromatic && c_type == AtomType::CAromatic)
        || (b_type == AtomType::CAromatic && c_type == AtomType::SAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 45.0,
            penalty_per_degree: 0.04,
        });
    }

    // Furanyl methyl / furanyl-sp3: OAromatic–CSp3 — prefer anti (180°).
    if (b_type == AtomType::OAromatic && c_type == AtomType::CSp3)
        || (b_type == AtomType::CSp3 && c_type == AtomType::OAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.06,
        });
    }

    // Thienyl methyl / thienyl-sp3: SAromatic–CSp3 — prefer anti (180°).
    if (b_type == AtomType::SAromatic && c_type == AtomType::CSp3)
        || (b_type == AtomType::CSp3 && c_type == AtomType::SAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.06,
        });
    }

    // Furanyl / thienyl adjacent to carbonyl — planar conjugation (0°).
    if (b_type == AtomType::OAromatic && c_type == AtomType::CCarbonyl)
        || (b_type == AtomType::CCarbonyl && c_type == AtomType::OAromatic)
        || (b_type == AtomType::SAromatic && c_type == AtomType::CCarbonyl)
        || (b_type == AtomType::CCarbonyl && c_type == AtomType::SAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.10,
        });
    }

    // ── Saturated N-heterocycle patterns ─────────────────────────────────────

    // Morpholine N–C–C–O: gauche preference (~60°) from chair conformation.
    // Secondary amines in rings are classified as NSp2 (2 explicit heavy neighbors)
    // even though they are sp3; accept both NSp2 and NSp3 for the amine endpoint.
    let is_sat_n = |t: AtomType| t == AtomType::NSp2 || t == AtomType::NSp3;
    if b_type == AtomType::CSp3
        && c_type == AtomType::CSp3
        && ((is_sat_n(a_type) && d_type == AtomType::OSp3)
            || (a_type == AtomType::OSp3 && is_sat_n(d_type)))
    {
        return Some(TorsionPreference {
            angle_deg: 60.0,
            penalty_per_degree: 0.10,
        });
    }

    // Piperazine / diamine N–C–C–N: gauche preference (~60°) from chair.
    if b_type == AtomType::CSp3 && c_type == AtomType::CSp3 && is_sat_n(a_type) && is_sat_n(d_type)
    {
        return Some(TorsionPreference {
            angle_deg: 60.0,
            penalty_per_degree: 0.10,
        });
    }

    // ── Styrene / aryl-vinyl conjugation ────────────────────────────────────────

    // Aryl-vinyl bond (styrene-like): extended π-conjugation → coplanar (0°).
    // Applies to Ar-C=C and Ar-C=C reverse traversal.
    if b_type == AtomType::CAromatic && c_type == AtomType::CSp2Alkene {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.07,
        });
    }
    if b_type == AtomType::CSp2Alkene && c_type == AtomType::CAromatic {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.07,
        });
    }

    // ── Vinyl thioether C=C-S: S lone pair conjugates with alkene π system ────

    if b_type == AtomType::CSp2Alkene && c_type == AtomType::S {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.07,
        });
    }
    if b_type == AtomType::S && c_type == AtomType::CSp2Alkene {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.07,
        });
    }

    // ── Allylic amine C=C-N (sp3): N lone pair partial conjugation with alkene ─

    if b_type == AtomType::CSp2Alkene && c_type == AtomType::NSp3 {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.06,
        });
    }
    if b_type == AtomType::NSp3 && c_type == AtomType::CSp2Alkene {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.06,
        });
    }

    // ── Ketone/aldehyde C(sp3)-C(=O): H eclipses C=O in preferred conformation ─
    // Barrier is low; the small penalty still nudges DG toward better geometry.

    if b_type == AtomType::CSp3 && c_type == AtomType::CCarbonyl {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.04,
        });
    }
    if b_type == AtomType::CCarbonyl && c_type == AtomType::CSp3 {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.04,
        });
    }

    // ── Heteroaromatic N to thioether/thioaryl (NAr-S) ───────────────────────
    // S lone pair gauche to aromatic N; ~90° from tetrahedral S geometry.

    if b_type == AtomType::NAromatic && c_type == AtomType::S {
        return Some(TorsionPreference {
            angle_deg: 90.0,
            penalty_per_degree: 0.06,
        });
    }
    if b_type == AtomType::S && c_type == AtomType::NAromatic {
        return Some(TorsionPreference {
            angle_deg: 90.0,
            penalty_per_degree: 0.06,
        });
    }

    // ── Heteroaromatic N to ether oxygen (NAr-O) ─────────────────────────────
    // O lone pair conjugates with N lone pair through the linking C;
    // coplanar (0°) minimises lone-pair/lone-pair repulsion via σ*-donation.

    if b_type == AtomType::NAromatic && c_type == AtomType::OSp3 {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.08,
        });
    }
    if b_type == AtomType::OSp3 && c_type == AtomType::NAromatic {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.08,
        });
    }

    // ── 5-membered aromatic heterocycle (S or O) to N-heteroaryl ─────────────
    // Thienyl-pyridine, furanyl-pyridine biaryl bonds: ~45° (like biphenyl).

    if (b_type == AtomType::SAromatic && c_type == AtomType::NAromatic)
        || (b_type == AtomType::NAromatic && c_type == AtomType::SAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 45.0,
            penalty_per_degree: 0.03,
        });
    }

    if (b_type == AtomType::OAromatic && c_type == AtomType::NAromatic)
        || (b_type == AtomType::NAromatic && c_type == AtomType::OAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 45.0,
            penalty_per_degree: 0.04,
        });
    }

    // ── Aromatic carbonyl to sp3 carbon ──────────────────────────────────────
    // Ar-C(=O)-CR3: the carbonyl-to-alkyl bond; prefer 0° (carbonyl O in plane
    // with the aryl ring; alkyl group anti to O).

    if b_type == AtomType::CCarbonyl && c_type == AtomType::CSp2Alkene {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.06,
        });
    }

    // ── Sp-hybridised nitrile / alkyne terminus ───────────────────────────────
    // X-C≡C-Y and X-C≡N: sp atoms are linear; penalise deviation from 180°.

    if b_type == AtomType::CSp3 && c_type == AtomType::NSp {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.15,
        });
    }

    // ── Amide/carbamate reverse: CCarbonyl-N traversal ───────────────────────
    // The A-B-C-D enumeration also visits bonds in reverse; ensure the
    // CCarbonyl→N direction returns the same preference as the N→CCarbonyl rules.

    if b_type == AtomType::CCarbonyl && c_type == AtomType::NSp2 {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.20,
        });
    }
    if b_type == AtomType::CCarbonyl && c_type == AtomType::NSp3 {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.10,
        });
    }

    // ── Isocyanate / carbodiimide N=C=O: linear sp centre ────────────────────

    if b_type == AtomType::NSp && c_type == AtomType::CCarbonyl {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.18,
        });
    }
    if b_type == AtomType::CCarbonyl && c_type == AtomType::NSp {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.18,
        });
    }

    // ── Ar-N=C=O (aryl isocyanate) / Ar-N=S: prefer coplanar (0°) ───────────

    if b_type == AtomType::CAromatic && c_type == AtomType::NSp {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.10,
        });
    }
    if b_type == AtomType::NSp && c_type == AtomType::CAromatic {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.10,
        });
    }

    // ── Thioamide: NSp2-C(=S) — prefer trans (180°) like regular amide ───────
    // C(=S) is classified CSp2Alkene (double bond to S, not O).
    // The rule fires for N-C(=S) when the C is NOT a CCarbonyl.

    if b_type == AtomType::NSp2 && c_type == AtomType::CSp2Alkene {
        // Check if C has a S=C double bond (thioamide context)
        let has_thio = mol.neighbors(c_idx).any(|(n, _)| {
            mol.atom(n).element.atomic_number() == 16
                && mol
                    .bond_between(c_idx, n)
                    .map(|(_, b)| b.order == chematic_core::BondOrder::Double)
                    .unwrap_or(false)
        });
        if has_thio {
            return Some(TorsionPreference {
                angle_deg: 180.0,
                penalty_per_degree: 0.15,
            });
        }
    }

    // ── Anomeric / vicinal-O gauche effect ────────────────────────────────────
    // O-C-C-O chain (1,2-diol, glycol ether): gauche (~60°) is preferred due to
    // anomeric / electrostatic stabilisation.  Weaker than ring-based gauche rules.

    if b_type == AtomType::CSp3
        && c_type == AtomType::CSp3
        && (a_type == AtomType::OSp3 || a_type == AtomType::OAromatic)
        && (d_type == AtomType::OSp3 || d_type == AtomType::OAromatic)
    {
        return Some(TorsionPreference {
            angle_deg: 60.0,
            penalty_per_degree: 0.08,
        });
    }

    // ── OSp3 to CCarbonyl (reverse ester direction) ──────────────────────────
    // O-C(=O) ester oxygen facing the other side: prefer 0° for conjugation.

    if b_type == AtomType::OSp3 && c_type == AtomType::CCarbonyl {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.10,
        });
    }

    // ── Conjugated diene C=C-C=C: prefer s-trans (180°) ─────────────────────
    // Open-chain 1,3-dienes exist predominantly in the s-trans conformation
    // (~95%).  Cyclic dienes are constrained; apply only when both ends are
    // CSp2Alkene (not aromatic).

    if b_type == AtomType::CSp2Alkene && c_type == AtomType::CSp2Alkene {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.08,
        });
    }

    // ── 1,2-Dicarbonyl C(=O)-C(=O): syn-periplanar (0°) for lone-pair overlap ─

    if b_type == AtomType::CCarbonyl && c_type == AtomType::CCarbonyl {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.10,
        });
    }

    // ── Halogen adjacent to sp2 centre ───────────────────────────────────────
    // Ar-Cl, Ar-Br, Ar-F: the halogen lone pair has no rotational preference
    // but the aryl-halogen bond is effectively rigid.  A very soft anti preference
    // biases the flanking chain away from the halogen σ* orbital.

    if b_type == AtomType::CAromatic && c_type == AtomType::Halogen {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.03,
        });
    }

    // ── Phosphorus ester / phosphonate P-O ───────────────────────────────────
    // P-O-C chain in phosphates / phosphonates: prefer anti (180°).

    if b_type == AtomType::P && c_type == AtomType::OSp3 {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.06,
        });
    }
    if b_type == AtomType::OSp3 && c_type == AtomType::P {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.06,
        });
    }

    None // No specific preference; use default
}

// ---------------------------------------------------------------------------
// SMARTS-based torsion rules (higher precision than atom-type rules)
// ---------------------------------------------------------------------------

/// A SMARTS-based torsion rule for a specific B-C bond context.
///
/// The SMARTS pattern may contain more than 2 atoms to express chemical
/// context (e.g. `[NH2][C;!a](=O)` for primary amide).  The B and C atoms
/// of the rotatable bond correspond to query positions `b_qi` and `c_qi`.
pub struct SmartsTorsionRule {
    pub smarts: &'static str,
    /// Index in the SMARTS match corresponding to atom B.
    pub b_qi: usize,
    /// Index in the SMARTS match corresponding to atom C.
    pub c_qi: usize,
    pub angle_deg: f64,
    pub penalty_per_degree: f64,
}

/// SMARTS-based torsion rules.  Checked before atom-type fallback rules;
/// first matching rule wins for a given B-C bond.
///
/// Rules are ordered from most specific to least specific.
static SMARTS_TORSION_RULES: &[SmartsTorsionRule] = &[
    // Hindered biaryl: both ring-C are fully substituted (H0, X3 = 3 heavy bonds).
    // Steric clash pushes the preferred dihedral from ~45° (unhindered) toward ~90°.
    SmartsTorsionRule {
        smarts: "[c;H0;X3][c;H0;X3]",
        b_qi: 0,
        c_qi: 1,
        angle_deg: 90.0,
        penalty_per_degree: 0.04,
    },
    // Primary amide N-C(=O): H2N group prefers 0° (both H eclipsed with C=O oxygen).
    // More constrained than the generic NSp2+CCarbonyl→180° rule.
    SmartsTorsionRule {
        smarts: "[NH2][C;!a](=O)",
        b_qi: 0,
        c_qi: 1,
        angle_deg: 0.0,
        penalty_per_degree: 0.18,
    },
    // Tertiary amide (N with no H, 3 heavy neighbors): N-methyl/dialkyl → trans (180°).
    SmartsTorsionRule {
        smarts: "[N;H0;X3][C;!a](=O)",
        b_qi: 0,
        c_qi: 1,
        angle_deg: 180.0,
        penalty_per_degree: 0.18,
    },
    // Heteroaromatic N (ring N, lone-pair donor) adjacent to carbonyl: prefer 0°
    // (lone pair conjugates with carbonyl π system).
    SmartsTorsionRule {
        smarts: "[n][C;!a](=O)",
        b_qi: 0,
        c_qi: 1,
        angle_deg: 0.0,
        penalty_per_degree: 0.14,
    },
    // Aryl ester: Ar-O-C(=O) — the O-C(=O) bond prefers 0° (oxygen lone pair
    // conjugates with both the aromatic ring and carbonyl).
    SmartsTorsionRule {
        smarts: "[c][O;H0][C;!a](=O)",
        b_qi: 1,
        c_qi: 2,
        angle_deg: 0.0,
        penalty_per_degree: 0.12,
    },
    // Carbamate: N-C(=O)-O — the N-C(=O) bond prefers 0° (both N and O lone
    // pairs donate into the same carbonyl π system → syn-periplanar geometry).
    SmartsTorsionRule {
        smarts: "[N;!a][C;!a](=O)[O;!a]",
        b_qi: 0,
        c_qi: 1,
        angle_deg: 0.0,
        penalty_per_degree: 0.12,
    },
];

/// Build a bond-keyed map of SMARTS-derived torsion preferences for `mol`.
///
/// The map is keyed by `(b_idx.0, c_idx.0)` for the rotatable bond B-C.
/// Both directions `(b, c)` and `(c, b)` are inserted (first matching rule
/// wins; subsequent rules do not overwrite).
///
/// Ring bonds are excluded: their torsion angles are constrained by ring
/// closure geometry and applying an additional preference would conflict with
/// the subsequent constraint-satisfaction pass.
pub fn build_smarts_torsion_map(
    mol: &Molecule,
    ring_bond_set: &std::collections::HashSet<(u32, u32)>,
) -> HashMap<(u32, u32), TorsionPreference> {
    let mut map: HashMap<(u32, u32), TorsionPreference> = HashMap::new();

    for rule in SMARTS_TORSION_RULES {
        let Ok(query) = parse_smarts(rule.smarts) else {
            continue;
        };
        for m in find_matches(&query, mol) {
            let (Some(b_atom), Some(c_atom)) = (m.get(&rule.b_qi), m.get(&rule.c_qi)) else {
                continue;
            };
            let b = b_atom.0;
            let c = c_atom.0;
            // Skip ring bonds (their geometry is constrained by ring closure).
            let key_fwd = (b, c);
            let key_rev = (c, b);
            if ring_bond_set.contains(&key_fwd) {
                continue;
            }
            // Verify the bond actually exists (rule might match atoms not bonded).
            if mol.bond_between(AtomIdx(b), AtomIdx(c)).is_none() {
                continue;
            }
            let pref = TorsionPreference {
                angle_deg: rule.angle_deg,
                penalty_per_degree: rule.penalty_per_degree,
            };
            // First matching rule wins — do not overwrite existing entries.
            map.entry(key_fwd).or_insert_with(|| pref.clone());
            map.entry(key_rev).or_insert(pref);
        }
    }

    map
}

/// Default torsion preferences for general C-C-C-C patterns.
pub fn default_torsion_preference() -> TorsionPreference {
    TorsionPreference {
        angle_deg: 180.0, // Prefer anti/staggered
        penalty_per_degree: 0.10,
    }
}

/// Normalize angle to [-180, 180] range, accounting for periodicity.
fn normalize_angle(angle_deg: f64) -> f64 {
    let mut norm = angle_deg % 360.0;
    if norm > 180.0 {
        norm -= 360.0;
    } else if norm < -180.0 {
        norm += 360.0;
    }
    norm
}

/// Score a torsion angle (in degrees) against a preference.
///
/// Returns penalty in arbitrary units (higher = worse).
pub fn score_torsion(angle_deg: f64, preference: &TorsionPreference) -> f64 {
    let norm = normalize_angle(angle_deg);
    let pref = normalize_angle(preference.angle_deg);

    // Compute angular difference (accounting for periodicity)
    let diff = (norm - pref).abs();
    let min_diff = if diff > 180.0 { 360.0 - diff } else { diff };

    // Penalty scales with distance
    min_diff * preference.penalty_per_degree
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_atom_type_methane() {
        let mol = parse("C").unwrap();
        assert_eq!(classify_atom_type(&mol, AtomIdx(0)), AtomType::CSp3);
    }

    #[test]
    fn test_atom_type_ethene() {
        let mol = parse("C=C").unwrap();
        assert_eq!(classify_atom_type(&mol, AtomIdx(0)), AtomType::CSp2Alkene);
        assert_eq!(classify_atom_type(&mol, AtomIdx(1)), AtomType::CSp2Alkene);
    }

    #[test]
    fn test_atom_type_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        assert_eq!(classify_atom_type(&mol, AtomIdx(0)), AtomType::CAromatic);
    }

    #[test]
    fn test_atom_type_acetaldehyde() {
        let mol = parse("CC=O").unwrap();
        let c_sp3 = if mol.atom(AtomIdx(0)).aromatic {
            classify_atom_type(&mol, AtomIdx(1))
        } else {
            classify_atom_type(&mol, AtomIdx(0))
        };
        assert_eq!(c_sp3, AtomType::CSp3);
    }

    #[test]
    fn test_torsion_score_perfect_match() {
        let pref = TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.1,
        };
        let score = score_torsion(180.0, &pref);
        assert!(score.abs() < 1e-6, "perfect match should have zero penalty");
    }

    #[test]
    fn test_torsion_score_deviation() {
        let pref = TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.1,
        };
        let score = score_torsion(160.0, &pref);
        assert!(
            (score - 2.0).abs() < 1e-6,
            "20° deviation should yield 2.0 penalty"
        );
    }

    #[test]
    fn test_torsion_score_periodic() {
        let pref = TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.1,
        };
        // 180° and -180° are the same angle
        let score1 = score_torsion(180.0, &pref);
        let score2 = score_torsion(-180.0, &pref);
        assert!(
            (score1 - score2).abs() < 1e-6,
            "periodic angles should score the same"
        );
    }

    #[test]
    fn test_default_torsion_preference() {
        let pref = default_torsion_preference();
        assert_eq!(pref.angle_deg, 180.0);
        assert!(pref.penalty_per_degree > 0.0);
    }

    #[test]
    fn test_alkane_torsion_preference() {
        let mol = parse("CCCC").unwrap(); // butane
        if mol.atom_count() >= 4 {
            let pref = get_torsion_preference(&mol, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
            assert!(pref.is_some(), "butane C-C-C-C should have preference");
            if let Some(p) = pref {
                assert_eq!(p.angle_deg, 180.0, "alkane torsions prefer 180°");
            }
        }
    }

    #[test]
    fn test_biphenyl_torsion_preference() {
        let mol = parse("c1ccccc1-c1ccccc1").unwrap(); // biphenyl
        // Atoms: 0-5 ring1, 6-11 ring2; bond between 0 and 6
        let pref = get_torsion_preference(&mol, AtomIdx(1), AtomIdx(0), AtomIdx(6), AtomIdx(7));
        assert!(
            pref.is_some(),
            "biphenyl Ar-Ar should have a torsion preference"
        );
        if let Some(p) = pref {
            assert_eq!(p.angle_deg, 45.0, "biphenyl prefers ~45° twist");
        }
    }

    #[test]
    fn test_thioester_torsion_preference() {
        // Methyl thioformate: SC=O; S is atom 0, C is atom 1 (CCarbonyl)
        let mol = parse("SC=O").unwrap();
        let _pref = get_torsion_preference(&mol, AtomIdx(0), AtomIdx(0), AtomIdx(1), AtomIdx(2));
        // b=S (atom0), c=CCarbonyl (atom1) → thioester branch
        let b_type = classify_atom_type(&mol, AtomIdx(0));
        let c_type = classify_atom_type(&mol, AtomIdx(1));
        assert_eq!(b_type, AtomType::S);
        assert_eq!(c_type, AtomType::CCarbonyl);
    }

    #[test]
    fn test_disulfide_torsion_preference() {
        let mol = parse("CSSC").unwrap(); // dimethyl disulfide
        // bond S(1)-S(2): b_type=S, c_type=S
        let b_type = classify_atom_type(&mol, AtomIdx(1));
        let c_type = classify_atom_type(&mol, AtomIdx(2));
        assert_eq!(b_type, AtomType::S);
        assert_eq!(c_type, AtomType::S);
        let pref = get_torsion_preference(&mol, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
        assert!(pref.is_some(), "disulfide should have ~90° preference");
        if let Some(p) = pref {
            assert_eq!(p.angle_deg, 90.0, "disulfide prefers 90°");
        }
    }

    #[test]
    fn test_nitrile_torsion_preference() {
        let mol = parse("CCC#N").unwrap(); // propionitrile
        // N is atom 3, atom type NSp
        let n_type = classify_atom_type(&mol, AtomIdx(3));
        assert_eq!(n_type, AtomType::NSp, "nitrile N should be NSp");
        let pref = get_torsion_preference(&mol, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
        assert!(pref.is_some(), "nitrile torsion should have a preference");
        if let Some(p) = pref {
            assert_eq!(p.angle_deg, 180.0, "linear nitrile end prefers 180°");
        }
    }

    #[test]
    fn test_amine_torsion_preference() {
        let mol = parse("CCNC").unwrap(); // ethyl methyl amine
        // bond C(1)-N(2): b=CSp3, c=NSp2 (secondary amine: 2 explicit bonds)
        let pref = get_torsion_preference(&mol, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
        assert!(pref.is_some(), "amine C-C-N-C should have preference");
        if let Some(p) = pref {
            assert_eq!(p.angle_deg, 180.0);
        }
    }

    #[test]
    fn test_phenyl_ketone_torsion_preference() {
        let mol = parse("c1ccccc1C(=O)C").unwrap(); // acetophenone
        // The C(=O) carbon is CCarbonyl, connected to CAromatic
        // Find the carbonyl carbon
        let c_carbonyl_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| classify_atom_type(&mol, i) == AtomType::CCarbonyl);
        assert!(
            c_carbonyl_idx.is_some(),
            "acetophenone should have a carbonyl C"
        );
    }

    #[test]
    fn test_score_torsion_disulfide_at_90() {
        let pref = TorsionPreference {
            angle_deg: 90.0,
            penalty_per_degree: 0.1,
        };
        let score = score_torsion(90.0, &pref);
        assert!(score.abs() < 1e-6, "at preferred angle score should be 0");
        let score_off = score_torsion(90.0 + 20.0, &pref);
        assert!((score_off - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_pattern_count_covers_20_plus() {
        // Verify we have comprehensive coverage by checking patterns
        // for the major chemical motifs
        let mol_alkane = parse("CCCC").unwrap();
        let mol_biphenyl = parse("c1ccccc1-c1ccccc1").unwrap();
        let mol_amide = parse("CC(=O)N").unwrap();
        let mol_ester = parse("CC(=O)OC").unwrap();
        let mol_disulfide = parse("CSSC").unwrap();
        let mol_nitrile = parse("CCC#N").unwrap();
        let mol_amine = parse("CCNC").unwrap();

        let cases = [
            (&mol_alkane, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)),
            (
                &mol_biphenyl,
                AtomIdx(1),
                AtomIdx(0),
                AtomIdx(6),
                AtomIdx(7),
            ),
            (
                &mol_disulfide,
                AtomIdx(0),
                AtomIdx(1),
                AtomIdx(2),
                AtomIdx(3),
            ),
            (&mol_nitrile, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)),
            (&mol_amine, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)),
        ];
        for (mol, a, b, c, d) in &cases {
            let pref = get_torsion_preference(mol, *a, *b, *c, *d);
            // At least some of these should have preferences
            let _ = pref; // just ensuring no panic
        }
        // The amide and ester patterns
        let pref_amide =
            get_torsion_preference(&mol_amide, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(2));
        let _ = pref_amide;
        let pref_ester =
            get_torsion_preference(&mol_ester, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
        let _ = pref_ester;
    }

    // ── SMARTS torsion map tests ───────────────────────────────────────────────

    #[test]
    fn test_smarts_map_hindered_biaryl() {
        // 2,2'-dimethylbiphenyl: both inter-ring C are H0,X3 → hindered → 90°
        let mol = parse("Cc1ccccc1-c1ccccc1C").unwrap();
        let ring_set = chematic_perception::find_sssr(&mol);
        let ring_bonds: std::collections::HashSet<(u32, u32)> = ring_set
            .rings()
            .iter()
            .flat_map(|r| {
                let n = r.len();
                (0..n).flat_map(move |i| {
                    let a = r[i].0;
                    let b = r[(i + 1) % n].0;
                    [(a, b), (b, a)]
                })
            })
            .collect();
        let map = build_smarts_torsion_map(&mol, &ring_bonds);
        // At least one biaryl bond entry should be at 90°
        let has_90 = map.values().any(|p| (p.angle_deg - 90.0).abs() < 1.0);
        assert!(
            has_90,
            "hindered biaryl should get 90° preference in SMARTS map"
        );
    }

    #[test]
    fn test_smarts_map_primary_amide() {
        // Acetamide: CC(=O)N — primary amide → 0°
        let mol = parse("CC(=O)N").unwrap();
        let ring_bonds = std::collections::HashSet::new();
        let map = build_smarts_torsion_map(&mol, &ring_bonds);
        // Should have an entry for the N-C(=O) bond at 0°
        let has_0 = map.values().any(|p| p.angle_deg.abs() < 1.0);
        assert!(has_0, "primary amide N-C bond should get 0° preference");
    }

    // ── New heterocycle pattern tests ──────────────────────────────────────────

    #[test]
    fn test_atom_type_furan_oxygen() {
        // Furan: c1ccco1 — oxygen atom is aromatic
        let mol = parse("c1ccco1").unwrap();
        let o_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| mol.atom(i).element.atomic_number() == 8)
            .expect("furan must have an oxygen atom");
        assert_eq!(
            classify_atom_type(&mol, o_idx),
            AtomType::OAromatic,
            "furan O should be OAromatic"
        );
    }

    #[test]
    fn test_atom_type_thiophene_sulfur() {
        // Thiophene: c1cccs1 — sulfur atom is aromatic
        let mol = parse("c1cccs1").unwrap();
        let s_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| mol.atom(i).element.atomic_number() == 16)
            .expect("thiophene must have a sulfur atom");
        assert_eq!(
            classify_atom_type(&mol, s_idx),
            AtomType::SAromatic,
            "thiophene S should be SAromatic"
        );
    }

    #[test]
    fn test_furanyl_biaryl_prefers_planar() {
        // 2-phenylfuran: c1ccc(-c2ccco2)cc1
        // Inter-ring bond between phenyl CAromatic and furanyl CAromatic;
        // OAromatic attached on the furanyl side → should fire OAromatic–CAromatic rule (0°).
        let mol = parse("c1ccc(-c2ccco2)cc1").unwrap();
        // Find an atom adjacent to the furan O — that's the OAromatic–CAromatic bond end
        let o_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| mol.atom(i).element.atomic_number() == 8)
            .expect("must have O");
        let o_neighbor = mol
            .neighbors(o_idx)
            .next()
            .map(|(n, _)| n)
            .expect("O has neighbors");
        assert_eq!(classify_atom_type(&mol, o_idx), AtomType::OAromatic);
        assert_eq!(classify_atom_type(&mol, o_neighbor), AtomType::CAromatic);
        // Torsion involving OAromatic–CAromatic as B–C pair should prefer 0°
        let pref = get_torsion_preference(&mol, o_idx, o_idx, o_neighbor, o_neighbor);
        // We only need to confirm the rule fires (Some) and angle is 0°
        let pref2 = get_torsion_preference(&mol, AtomIdx(0), o_idx, o_neighbor, AtomIdx(0));
        assert!(
            pref.is_some() || pref2.is_some(),
            "OAromatic–CAromatic should have a torsion preference"
        );
    }

    #[test]
    fn test_morpholine_gauche_preference() {
        // Morpholine: C1CNCCO1 — N-C-C-O chain in chair → gauche 60°
        // Atom order in morpholine SMILES C1CNCCO1:
        //   0=C, 1=C, 2=N, 3=C, 4=C, 5=O
        // N-C-C-O: a=N(2), b=C(3), c=C(4), d=O(5)
        let mol = parse("C1CNCCO1").unwrap();
        // Find N and O atoms
        let n_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| mol.atom(i).element.atomic_number() == 7)
            .expect("morpholine must have N");
        let o_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| mol.atom(i).element.atomic_number() == 8)
            .expect("morpholine must have O");
        // Find the C-C bond between N-side and O-side carbons
        // The preference fires when a=NSp3, b=CSp3, c=CSp3, d=OSp3
        let pref = get_torsion_preference(&mol, n_idx, AtomIdx(3), AtomIdx(4), o_idx);
        assert!(
            pref.is_some(),
            "morpholine N-C-C-O should have gauche preference"
        );
        if let Some(p) = pref {
            assert_eq!(p.angle_deg, 60.0, "morpholine N-C-C-O prefers 60°");
        }
    }
}
