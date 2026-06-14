//! ETKDG Torsion Knowledge Base — experimental torsion angle preferences from CSD.

use chematic_core::{AtomIdx, Molecule};

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
    /// sulfur
    S,
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
        .filter(|(_, bond)| {
            (bond.atom1 == idx || bond.atom2 == idx) && bond.order == order
        })
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
            // Oxygen: distinguish sp2 (carbonyl C=O) from sp3 (alcohol, ether C-O)
            let has_double_bond = mol
                .bonds()
                .any(|(_, bond)| {
                    (bond.atom1 == idx || bond.atom2 == idx) && bond.order == chematic_core::BondOrder::Double
                });
            if has_double_bond {
                AtomType::OSp2  // Carbonyl oxygen (C=O)
            } else {
                AtomType::OSp3  // Alcohol or ether oxygen (C-O)
            }
        }
        16 => AtomType::S,
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
            penalty_per_degree: 0.15,  // ~3 kcal/mol for 20° deviation
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
            penalty_per_degree: 0.08,  // Softer constraint
        });
    }

    // Amide C-N(C=O)-C-X: restricted rotation, prefer 0° (cis to carbonyl) or 180° (trans)
    // For simplicity, prefer trans (180°) as it's more common
    if b_type == AtomType::NSp2 && c_type == AtomType::CCarbonyl {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.20,  // Moderate restriction
        });
    }

    // Ester O-C(=O)-O-C: prefer 0° or 180° depending on configuration
    if b_type == AtomType::CCarbonyl && c_type == AtomType::OSp3 {
        return Some(TorsionPreference {
            angle_deg: 0.0,
            penalty_per_degree: 0.10,
        });
    }

    // Aromatic-aromatic rotations: less restricted, prefer 180°
    if b_type == AtomType::CAromatic && c_type == AtomType::CAromatic {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.05,  // Very soft constraint
        });
    }

    None  // No specific preference; use default
}

/// Default torsion preferences for general C-C-C-C patterns.
pub fn default_torsion_preference() -> TorsionPreference {
    TorsionPreference {
        angle_deg: 180.0,  // Prefer anti/staggered
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
        assert!((score - 2.0).abs() < 1e-6, "20° deviation should yield 2.0 penalty");
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
        assert!((score1 - score2).abs() < 1e-6, "periodic angles should score the same");
    }

    #[test]
    fn test_default_torsion_preference() {
        let pref = default_torsion_preference();
        assert_eq!(pref.angle_deg, 180.0);
        assert!(pref.penalty_per_degree > 0.0);
    }

    #[test]
    fn test_alkane_torsion_preference() {
        let mol = parse("CCCC").unwrap();  // butane
        if mol.atom_count() >= 4 {
            let pref = get_torsion_preference(&mol, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
            assert!(pref.is_some(), "butane C-C-C-C should have preference");
            if let Some(p) = pref {
                assert_eq!(p.angle_deg, 180.0, "alkane torsions prefer 180°");
            }
        }
    }
}
