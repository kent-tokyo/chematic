//! C-Series Phase 1: ERG (Extended Reduced Graph) fingerprints.
//!
//! Extended Reduced Graph fingerprints encode the topological structure of molecules
//! using reduced graphs where nodes represent functional groups and edges represent
//! the connectivity between them.
//!
//! Useful for:
//! - Functional group-based similarity searching
//! - Molecule classification by functional group composition
//! - Scaffold identification
//! - Chemotype clustering

use chematic_core::{Atom, BondOrder};
use crate::bitvec::BitVec2048;

/// Atom property encoding for ERG nodes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ErgAtomType {
    /// Aliphatic carbon
    CAliphatic = 0,
    /// Aromatic carbon
    CAromatic = 1,
    /// Nitrogen (any)
    N = 2,
    /// Oxygen (any)
    O = 3,
    /// Sulfur (any)
    S = 4,
    /// Halogen (F, Cl, Br, I)
    Halogen = 5,
    /// Other heteroatom
    Other = 6,
}

impl ErgAtomType {
    /// Get atom type from Atom and aromaticity.
    pub fn from_atom(atom: &Atom) -> Self {
        let an = atom.element.atomic_number();
        let aromatic = atom.aromatic;

        match an {
            6 => {
                if aromatic {
                    ErgAtomType::CAromatic
                } else {
                    ErgAtomType::CAliphatic
                }
            }
            7 => ErgAtomType::N,
            8 => ErgAtomType::O,
            16 => ErgAtomType::S,
            9 | 17 | 35 | 53 => ErgAtomType::Halogen,
            _ => ErgAtomType::Other,
        }
    }
}

/// Bond type encoding for ERG edges.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ErgBondType {
    /// Single bond
    Single = 0,
    /// Double bond
    Double = 1,
    /// Triple bond
    Triple = 2,
    /// Aromatic bond
    Aromatic = 3,
}

impl ErgBondType {
    /// Get bond type from BondOrder.
    pub fn from_bond(order: BondOrder) -> Self {
        match order {
            BondOrder::Single => ErgBondType::Single,
            BondOrder::Double => ErgBondType::Double,
            BondOrder::Triple => ErgBondType::Triple,
            BondOrder::Aromatic => ErgBondType::Aromatic,
            _ => ErgBondType::Single,
        }
    }
}

/// Configuration for ERG fingerprint generation.
#[derive(Clone, Debug)]
pub struct ErgConfig {
    /// Use atom counts in fingerprint
    pub use_atom_counts: bool,
    /// Use bond types
    pub use_bond_types: bool,
}

impl Default for ErgConfig {
    fn default() -> Self {
        ErgConfig {
            use_atom_counts: true,
            use_bond_types: true,
        }
    }
}

/// Extended Reduced Graph fingerprint.
#[derive(Clone, Debug)]
pub struct ErgFingerprint {
    /// Fingerprint bits (2048-bit)
    pub bits: BitVec2048,
    /// Atom type counts
    pub atom_counts: [u32; 7],
    /// Bond type counts
    pub bond_counts: [u32; 4],
}

impl ErgFingerprint {
    /// Calculate Tanimoto similarity between ERG fingerprints.
    pub fn tanimoto(&self, other: &ErgFingerprint) -> f64 {
        self.bits.tanimoto(&other.bits)
    }
}

/// Generate ERG fingerprint from a molecule.
pub fn erg(mol: &chematic_core::Molecule) -> ErgFingerprint {
    erg_with_config(mol, &ErgConfig::default())
}

/// Generate ERG fingerprint with custom configuration.
pub fn erg_with_config(
    mol: &chematic_core::Molecule,
    config: &ErgConfig,
) -> ErgFingerprint {
    let mut bits = BitVec2048::new();
    let mut atom_counts = [0u32; 7];
    let mut bond_counts = [0u32; 4];

    // Count atom types
    for (_, atom) in mol.atoms() {
        let erg_type = ErgAtomType::from_atom(atom);
        atom_counts[erg_type as usize] += 1;

        // Set bits based on atom type
        let bit_pos = (erg_type as usize) * 16;
        if bit_pos < 2048 {
            bits.set(bit_pos);
        }
    }

    // Count bond types
    for (_, bond) in mol.bonds() {
        let erg_type = ErgBondType::from_bond(bond.order);
        bond_counts[erg_type as usize] += 1;

        // Set bits based on bond type
        let bit_pos = 112 + (erg_type as usize) * 16;
        if bit_pos < 2048 {
            bits.set(bit_pos);
        }
    }

    // Encode atom counts in fingerprint
    if config.use_atom_counts {
        for (i, &count) in atom_counts.iter().enumerate() {
            for j in 0..4 {
                if ((count >> j) & 1) != 0 {
                    let bit_pos = 200 + i * 4 + j;
                    if bit_pos < 2048 {
                        bits.set(bit_pos);
                    }
                }
            }
        }
    }

    // Encode bond counts in fingerprint
    if config.use_bond_types {
        for (i, &count) in bond_counts.iter().enumerate() {
            for j in 0..4 {
                if ((count >> j) & 1) != 0 {
                    let bit_pos = 228 + i * 4 + j;
                    if bit_pos < 2048 {
                        bits.set(bit_pos);
                    }
                }
            }
        }
    }

    ErgFingerprint {
        bits,
        atom_counts,
        bond_counts,
    }
}

/// Convenience function for standard ERG generation.
pub fn erg_extended(mol: &chematic_core::Molecule) -> ErgFingerprint {
    erg(mol)
}

/// Calculate Tanimoto similarity between two molecules using ERG.
pub fn tanimoto_erg(mol1: &chematic_core::Molecule, mol2: &chematic_core::Molecule) -> f64 {
    let fp1 = erg(mol1);
    let fp2 = erg(mol2);
    fp1.tanimoto(&fp2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_erg_simple() {
        let mol = parse("CC").unwrap();
        let fp = erg(&mol);

        assert_eq!(fp.atom_counts[ErgAtomType::CAliphatic as usize], 2);
        assert!(fp.bits.popcount() > 0);
    }

    #[test]
    fn test_erg_identical() {
        let mol = parse("CC").unwrap();
        let fp1 = erg(&mol);
        let fp2 = erg(&mol);

        // Identical molecules should have identical fingerprints
        assert_eq!(fp1.bits.tanimoto(&fp2.bits), 1.0);
        assert_eq!(fp1.atom_counts, fp2.atom_counts);
    }

    #[test]
    fn test_erg_different_molecules() {
        let mol1 = parse("CC").unwrap();
        let mol2 = parse("c1ccccc1").unwrap();

        let fp1 = erg(&mol1);
        let fp2 = erg(&mol2);

        // Aliphatic vs aromatic should differ
        assert!(fp1.atom_counts[ErgAtomType::CAromatic as usize] == 0);
        assert!(fp2.atom_counts[ErgAtomType::CAromatic as usize] > 0);
    }

    #[test]
    fn test_erg_symmetry() {
        let mol1 = parse("CC").unwrap();
        let mol2 = parse("c1ccccc1").unwrap();

        let sim12 = tanimoto_erg(&mol1, &mol2);
        let sim21 = tanimoto_erg(&mol2, &mol1);

        // Similarity should be symmetric
        assert!((sim12 - sim21).abs() < 1e-10);
    }

    #[test]
    fn test_erg_heteroatom_detection() {
        let mol = parse("CCO").unwrap();
        let fp = erg(&mol);

        assert!(fp.atom_counts[ErgAtomType::O as usize] > 0);
    }

    #[test]
    fn test_erg_config() {
        let mol = parse("CC").unwrap();
        let config = ErgConfig {
            use_atom_counts: false,
            use_bond_types: true,
        };

        let fp = erg_with_config(&mol, &config);
        assert!(fp.bits.popcount() > 0);
    }

    #[test]
    fn test_erg_aromatic_vs_aliphatic() {
        let aliphatic = parse("CCCC").unwrap();
        let aromatic = parse("c1ccccc1").unwrap();

        let fp_aliphatic = erg(&aliphatic);
        let fp_aromatic = erg(&aromatic);

        // Aromatic should have aromatic carbon counts
        assert_eq!(fp_aliphatic.atom_counts[ErgAtomType::CAromatic as usize], 0);
        assert!(fp_aromatic.atom_counts[ErgAtomType::CAromatic as usize] > 0);
    }

    #[test]
    fn test_erg_bond_counting() {
        let single_bond = parse("CC").unwrap();
        let double_bond = parse("C=C").unwrap();

        let fp_single = erg(&single_bond);
        let fp_double = erg(&double_bond);

        // Different bond types
        assert!(fp_single.bond_counts[ErgBondType::Single as usize] > 0);
        assert!(fp_double.bond_counts[ErgBondType::Double as usize] > 0);
    }
}
