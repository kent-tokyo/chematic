//! Pharmacophore feature detection and fingerprinting.
//!
//! Pharmacophore features are abstract representations of functional groups important
//! for drug binding. This module detects: donors, acceptors, aromatic, hydrophobic,
//! positive, and negative features.

use chematic_core::{AtomIdx, Molecule};

/// Pharmacophore feature types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FeatureType {
    /// Hydrogen bond donor (N-H, O-H)
    Donor,
    /// Hydrogen bond acceptor (N, O with lone pairs)
    Acceptor,
    /// Aromatic ring
    Aromatic,
    /// Hydrophobic (aliphatic C, Br, Cl, I)
    Hydrophobic,
    /// Positive charge (N+, guanidinium)
    Positive,
    /// Negative charge (O-, S-)
    Negative,
}

/// A detected pharmacophore feature in a molecule.
#[derive(Debug, Clone)]
pub struct Feature {
    /// Type of pharmacophore feature
    pub ftype: FeatureType,
    /// Primary atom index
    pub atom: AtomIdx,
    /// Secondary atoms (e.g., all atoms in aromatic ring)
    pub neighbors: Vec<AtomIdx>,
}

/// Detect all pharmacophore features in a molecule.
pub fn detect_features(mol: &Molecule) -> Vec<Feature> {
    let mut features = Vec::new();

    for idx in 0..mol.atom_count() {
        let atom_idx = AtomIdx(idx as u32);
        let atom = mol.atom(atom_idx);

        match atom.element {
            chematic_core::Element::N => {
                // N is donor (implicit hydrogens)
                if atom.charge == 0 {
                    features.push(Feature {
                        ftype: FeatureType::Donor,
                        atom: atom_idx,
                        neighbors: vec![],
                    });
                }

                // N acceptor (lone pair)
                if !atom.aromatic {
                    features.push(Feature {
                        ftype: FeatureType::Acceptor,
                        atom: atom_idx,
                        neighbors: vec![],
                    });
                }

                // N+ positive
                if atom.charge > 0 {
                    features.push(Feature {
                        ftype: FeatureType::Positive,
                        atom: atom_idx,
                        neighbors: vec![],
                    });
                }
            }
            chematic_core::Element::O => {
                // O is donor (O-H implicit)
                features.push(Feature {
                    ftype: FeatureType::Donor,
                    atom: atom_idx,
                    neighbors: vec![],
                });

                // O acceptor (lone pair)
                features.push(Feature {
                    ftype: FeatureType::Acceptor,
                    atom: atom_idx,
                    neighbors: vec![],
                });

                // O- negative
                if atom.charge < 0 {
                    features.push(Feature {
                        ftype: FeatureType::Negative,
                        atom: atom_idx,
                        neighbors: vec![],
                    });
                }
            }
            chematic_core::Element::C => {
                // Aromatic carbon
                if atom.aromatic {
                    let ring_atoms = collect_aromatic_ring(mol, atom_idx);
                    if !ring_atoms.is_empty() {
                        features.push(Feature {
                            ftype: FeatureType::Aromatic,
                            atom: atom_idx,
                            neighbors: ring_atoms,
                        });
                    }
                }

                // Hydrophobic aliphatic C
                if !atom.aromatic {
                    features.push(Feature {
                        ftype: FeatureType::Hydrophobic,
                        atom: atom_idx,
                        neighbors: vec![],
                    });
                }
            }
            chematic_core::Element::BR | chematic_core::Element::CL | chematic_core::Element::I => {
                // Halides are hydrophobic
                features.push(Feature {
                    ftype: FeatureType::Hydrophobic,
                    atom: atom_idx,
                    neighbors: vec![],
                });
            }
            chematic_core::Element::S => {
                // Sulfur: donor (S-H) and acceptor
                let h_count = count_hydrogens(mol, atom_idx);
                if h_count > 0 {
                    features.push(Feature {
                        ftype: FeatureType::Donor,
                        atom: atom_idx,
                        neighbors: vec![],
                    });
                }
                features.push(Feature {
                    ftype: FeatureType::Acceptor,
                    atom: atom_idx,
                    neighbors: vec![],
                });
            }
            _ => {}
        }
    }

    features
}

/// Count hydrogen atoms bonded to a given atom.
fn count_hydrogens(mol: &Molecule, atom: AtomIdx) -> u8 {
    (0..mol.atom_count())
        .filter(|&i| {
            let idx = AtomIdx(i as u32);
            mol.atom(idx).element == chematic_core::Element::H
                && mol.bond_between(atom, idx).is_some()
        })
        .count() as u8
}

/// Collect all atoms in the aromatic ring containing the given atom.
fn collect_aromatic_ring(mol: &Molecule, start: AtomIdx) -> Vec<AtomIdx> {
    if !mol.atom(start).aromatic {
        return vec![];
    }

    let mut ring = vec![start];
    let mut visited = std::collections::HashSet::from([start]);
    let mut queue = std::collections::VecDeque::from([start]);

    while let Some(current) = queue.pop_front() {
        for neighbor_idx in 0..mol.atom_count() {
            let neighbor = AtomIdx(neighbor_idx as u32);
            if visited.contains(&neighbor) {
                continue;
            }

            if mol.bond_between(current, neighbor).is_some() && mol.atom(neighbor).aromatic {
                visited.insert(neighbor);
                queue.push_back(neighbor);
                ring.push(neighbor);
            }
        }
    }

    ring.sort_by_key(|a| a.0);
    ring
}

/// Generate a pharmacophore fingerprint from detected features.
/// Returns a bit vector where each bit represents a feature type at a specific atom.
pub fn features_to_bitvec(features: &[Feature]) -> u64 {
    let mut bits = 0u64;

    for (i, feature) in features.iter().enumerate() {
        if i >= 8 {
            break; // Limit to 8 features per molecule
        }

        let feature_bits = match feature.ftype {
            FeatureType::Donor => 0b001,
            FeatureType::Acceptor => 0b010,
            FeatureType::Aromatic => 0b011,
            FeatureType::Hydrophobic => 0b100,
            FeatureType::Positive => 0b101,
            FeatureType::Negative => 0b110,
        };

        bits |= (feature_bits as u64) << (i * 3);
    }

    bits
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_benzene_aromatic() {
        let mol = parse("c1ccccc1").unwrap();
        let features = detect_features(&mol);
        let aromatic_count = features
            .iter()
            .filter(|f| f.ftype == FeatureType::Aromatic)
            .count();
        assert!(aromatic_count > 0, "benzene should have aromatic feature");
    }

    #[test]
    fn test_ethanol_donor_acceptor() {
        let mol = parse("CCO").unwrap();
        let features = detect_features(&mol);
        let has_donor = features.iter().any(|f| f.ftype == FeatureType::Donor);
        let has_acceptor = features.iter().any(|f| f.ftype == FeatureType::Acceptor);
        assert!(has_donor, "ethanol should have donor (O-H)");
        assert!(has_acceptor, "ethanol should have acceptor (O)");
    }

    #[test]
    fn test_aniline_donor_aromatic() {
        let mol = parse("Nc1ccccc1").unwrap();
        let features = detect_features(&mol);
        let has_donor = features.iter().any(|f| f.ftype == FeatureType::Donor);
        let has_aromatic = features.iter().any(|f| f.ftype == FeatureType::Aromatic);
        assert!(has_donor, "aniline should have N-H donor");
        assert!(has_aromatic, "aniline should have aromatic rings");
    }
}
