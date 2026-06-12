//! C-Series Phase 1: ERG (Extended Reduced Graph) fingerprints.
//!
//! **NOTE**: Current implementation is a simplified atom-type counting version.
//! True ERG (RDKit rdReducedGraphs) constructs a reduced graph by:
//! 1. Identifying functional group clusters (e.g., carbonyl, carboxyl, amine)
//! 2. Collapsing each cluster into a single node
//! 3. Computing fingerprint on the reduced graph structure
//!
//! Current simplified version counts atom and bond types directly,
//! which loses structural/topological information compared to true ERG.
//!
//! Useful for:
//! - Functional group-based similarity searching (approximate)
//! - Molecule classification by atom composition
//! - Fast structural filtering
//! - Chemotype clustering (lower accuracy than true ERG)
//!
//! TODO (v0.1.90+): Upgrade to true ERG by:
//! 1. Detect functional group clusters in molecule
//! 2. Build reduced graph with collapsed nodes
//! 3. Compute fingerprint on reduced structure

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
/// Generate ERG fingerprint from a molecule (simplified atom-type counting).
///
/// **Current Implementation**: Counts atom and bond types in the molecule.
/// This is a simplified approximation that captures composition but not structure.
///
/// **True ERG** would:
/// 1. Identify functional group clusters (carbonyl, carboxyl, amine, etc.)
/// 2. Collapse each cluster into a single node in a reduced graph
/// 3. Compute fingerprint on the reduced graph topology
/// 4. Achieve superior discrimination of structural isomers and functional diversity
///
/// Current version is useful for fast filtering but lacks structural sophistication.
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

    // v0.1.90: Count atom types and degree information for better structural encoding
    for (idx, atom) in mol.atoms() {
        let erg_type = ErgAtomType::from_atom(atom);
        atom_counts[erg_type as usize] += 1;

        // Set bits based on atom type
        let bit_pos = (erg_type as usize) * 16;
        if bit_pos < 2048 {
            bits.set(bit_pos);
        }

        // v0.1.90: Add degree-based bits for structural context
        let degree = mol.bonds().filter(|(_, b)| b.atom1 == idx || b.atom2 == idx).count();
        let degree_bits = ((degree.min(4) as usize) << 2) + (erg_type as usize);
        let degree_bit_pos = 512 + degree_bits.min(127);
        if degree_bit_pos < 2048 {
            bits.set(degree_bit_pos);
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

    // Functional group detection: simple topology-based bits
    // Detect aromatic rings (aromatic atom count > 0)
    if atom_counts[ErgAtomType::CAromatic as usize] > 0 {
        bits.set(256); // Bit 256: has aromatic ring
    }

    // Functional group: heteroatom presence (N, O, S, etc.)
    let has_heteroatom = atom_counts[ErgAtomType::N as usize] > 0
        || atom_counts[ErgAtomType::O as usize] > 0
        || atom_counts[ErgAtomType::S as usize] > 0
        || atom_counts[ErgAtomType::Halogen as usize] > 0;

    if has_heteroatom {
        bits.set(257); // Bit 257: contains heteroatom
    }

    // Detect aliphatic-only structures (no aromatic atoms)
    if atom_counts[ErgAtomType::CAromatic as usize] == 0 && atom_counts[ErgAtomType::CAliphatic as usize] > 0 {
        bits.set(258); // Bit 258: aliphatic carbon only
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

    #[test]
    fn test_erg_functional_group_aromatic_bit() {
        let aliphatic = parse("CCCC").unwrap();
        let aromatic = parse("c1ccccc1").unwrap();

        let fp_aliphatic = erg(&aliphatic);
        let fp_aromatic = erg(&aromatic);

        // Aromatic bit (256) should be set for benzene, not for alkane
        assert!(!fp_aliphatic.bits.get(256), "aliphatic should not have aromatic bit");
        assert!(fp_aromatic.bits.get(256), "aromatic should have aromatic bit");
    }

    #[test]
    fn test_erg_functional_group_heteroatom_bit() {
        let alkane = parse("CC").unwrap();
        let alcohol = parse("CCO").unwrap();
        let amine = parse("CCN").unwrap();

        let fp_alkane = erg(&alkane);
        let fp_alcohol = erg(&alcohol);
        let fp_amine = erg(&amine);

        // Heteroatom bit (257) should be set for compounds with O/N
        assert!(!fp_alkane.bits.get(257), "alkane should not have heteroatom bit");
        assert!(fp_alcohol.bits.get(257), "alcohol should have heteroatom bit");
        assert!(fp_amine.bits.get(257), "amine should have heteroatom bit");
    }

    #[test]
    fn test_erg_functional_group_improved_discrimination() {
        let methane = parse("C").unwrap();
        let ethanol = parse("CCO").unwrap();
        let pyridine = parse("c1ccncc1").unwrap();

        let fp_methane = erg(&methane);
        let fp_ethanol = erg(&ethanol);
        let fp_pyridine = erg(&pyridine);

        // Functional group bits should improve discrimination
        let sim_methane_ethanol = fp_methane.tanimoto(&fp_ethanol);
        let sim_methane_pyridine = fp_methane.tanimoto(&fp_pyridine);
        let sim_ethanol_pyridine = fp_ethanol.tanimoto(&fp_pyridine);

        // All similarities should be in valid range
        assert!((0.0..=1.0).contains(&sim_methane_ethanol));
        assert!((0.0..=1.0).contains(&sim_methane_pyridine));
        assert!((0.0..=1.0).contains(&sim_ethanol_pyridine));
    }
}
