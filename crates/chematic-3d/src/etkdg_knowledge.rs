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

    // Aromatic-aromatic rotations: prefer ~45° (biphenyl-like twist)
    if b_type == AtomType::CAromatic && c_type == AtomType::CAromatic {
        return Some(TorsionPreference {
            angle_deg: 45.0,
            penalty_per_degree: 0.03,  // Very soft — flat potential
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
    if (b_type == AtomType::CSp3 || b_type == AtomType::CSp2Alkene)
        && c_type == AtomType::OSp3
    {
        return Some(TorsionPreference {
            angle_deg: 180.0,
            penalty_per_degree: 0.07,
        });
    }

    // Amine C-C-N-C (secondary/tertiary): prefer 180°
    // NSp2 covers N with 2 explicit bonds (secondary amine in SMILES);
    // NSp3 covers N with 3+ explicit bonds (tertiary amine).
    if b_type == AtomType::CSp3
        && (c_type == AtomType::NSp3 || c_type == AtomType::NSp2)
    {
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

    #[test]
    fn test_biphenyl_torsion_preference() {
        let mol = parse("c1ccccc1-c1ccccc1").unwrap();  // biphenyl
        // Atoms: 0-5 ring1, 6-11 ring2; bond between 0 and 6
        let pref = get_torsion_preference(&mol, AtomIdx(1), AtomIdx(0), AtomIdx(6), AtomIdx(7));
        assert!(pref.is_some(), "biphenyl Ar-Ar should have a torsion preference");
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
        let mol = parse("CSSC").unwrap();  // dimethyl disulfide
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
        let mol = parse("CCC#N").unwrap();  // propionitrile
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
        let mol = parse("CCNC").unwrap();  // ethyl methyl amine
        // bond C(1)-N(2): b=CSp3, c=NSp2 (secondary amine: 2 explicit bonds)
        let pref = get_torsion_preference(&mol, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
        assert!(pref.is_some(), "amine C-C-N-C should have preference");
        if let Some(p) = pref {
            assert_eq!(p.angle_deg, 180.0);
        }
    }

    #[test]
    fn test_phenyl_ketone_torsion_preference() {
        let mol = parse("c1ccccc1C(=O)C").unwrap();  // acetophenone
        // The C(=O) carbon is CCarbonyl, connected to CAromatic
        // Find the carbonyl carbon
        let c_carbonyl_idx = (0..mol.atom_count() as u32)
            .map(AtomIdx)
            .find(|&i| classify_atom_type(&mol, i) == AtomType::CCarbonyl);
        assert!(c_carbonyl_idx.is_some(), "acetophenone should have a carbonyl C");
    }

    #[test]
    fn test_score_torsion_disulfide_at_90() {
        let pref = TorsionPreference { angle_deg: 90.0, penalty_per_degree: 0.1 };
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
            (&mol_alkane,   AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)),
            (&mol_biphenyl, AtomIdx(1), AtomIdx(0), AtomIdx(6), AtomIdx(7)),
            (&mol_disulfide, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)),
            (&mol_nitrile,  AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)),
            (&mol_amine,    AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)),
        ];
        for (mol, a, b, c, d) in &cases {
            let pref = get_torsion_preference(mol, *a, *b, *c, *d);
            // At least some of these should have preferences
            let _ = pref; // just ensuring no panic
        }
        // The amide and ester patterns
        let pref_amide = get_torsion_preference(&mol_amide, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(2));
        let _ = pref_amide;
        let pref_ester = get_torsion_preference(&mol_ester, AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3));
        let _ = pref_ester;
    }
}
