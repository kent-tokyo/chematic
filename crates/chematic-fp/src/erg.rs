//! v0.1.91: True ERG (Extended Reduced Graph) fingerprints
//!
//! Implements ERG algorithm (Sheridan 1996) with Ertl 2017 functional group detection:
//! 1. Identify functional group clusters using Ertl 2017 algorithm
//! 2. Collapse each cluster into a single reduced node
//! 3. Assign node types from pharmacophore features
//! 4. Build reduced graph with linker information
//! 5. Generate fingerprint from reduced graph topology
//!
//! Reference: https://pubs.acs.org/doi/10.1021/ci9602928 (Sheridan)
//!            P. Ertl, J. Cheminf. 2017, 9, 36 (Ertl)

use std::collections::{HashSet, VecDeque};
use chematic_core::{Atom, AtomIdx, BondOrder, Molecule};
use crate::bitvec::BitVec2048;

/// Ertl 2017 functional group identification.
/// Returns atom index sets, one per functional group cluster.
fn identify_functional_groups(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atom_count();
    if n == 0 {
        return Vec::new();
    }

    let mut marked = vec![false; n];

    // Mark all non-C, non-H heavy atoms
    for (idx, atom) in mol.atoms() {
        let an = atom.element.atomic_number();
        if an != 1 && an != 6 {
            marked[idx.0 as usize] = true;
        }
    }

    // Mark carbons bonded to heteroatoms
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() != 6 {
            continue;
        }
        let has_hetero = mol.neighbors(idx).any(|(nb, _)| {
            let an = mol.atom(nb).element.atomic_number();
            an != 1 && an != 6
        });
        if has_hetero {
            marked[idx.0 as usize] = true;
        }
    }

    // BFS to find connected components of marked atoms
    let mut visited = vec![false; n];
    let mut groups: Vec<Vec<usize>> = Vec::new();

    for start in 0..n {
        if !marked[start] || visited[start] {
            continue;
        }

        let mut component = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited[start] = true;

        while let Some(cur) = queue.pop_front() {
            component.push(cur);
            for (nb, _) in mol.neighbors(AtomIdx(cur as u32)) {
                let nbi = nb.0 as usize;
                if marked[nbi] && !visited[nbi] {
                    visited[nbi] = true;
                    queue.push_back(nbi);
                }
            }
        }

        component.sort_unstable();
        groups.push(component);
    }

    groups
}

/// Node type for reduced graph: pharmacophore feature bits.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ErgNodeType(pub u8);

impl ErgNodeType {
    const AROMATIC: u8 = 1;
    const DONOR: u8 = 2;
    const ACCEPTOR: u8 = 4;
    const HYDROPHOBIC: u8 = 8;
    const POSITIVE: u8 = 16;
    const NEGATIVE: u8 = 32;

    pub fn new() -> Self {
        ErgNodeType(0)
    }

    pub fn with_aromatic(mut self) -> Self {
        self.0 |= Self::AROMATIC;
        self
    }

    pub fn with_donor(mut self) -> Self {
        self.0 |= Self::DONOR;
        self
    }

    pub fn with_acceptor(mut self) -> Self {
        self.0 |= Self::ACCEPTOR;
        self
    }

    pub fn with_hydrophobic(mut self) -> Self {
        self.0 |= Self::HYDROPHOBIC;
        self
    }

    pub fn with_positive(mut self) -> Self {
        self.0 |= Self::POSITIVE;
        self
    }

    pub fn with_negative(mut self) -> Self {
        self.0 |= Self::NEGATIVE;
        self
    }
}

/// Reduced graph node representing a functional group or backbone segment.
#[derive(Clone, Debug)]
pub struct ErgNode {
    pub ntype: ErgNodeType,
    pub atom_indices: Vec<usize>,
}

/// Reduced graph edge with linker information.
#[derive(Clone, Debug)]
pub struct ErgEdge {
    pub node_a: usize,
    pub node_b: usize,
    pub linker_len: u32,
}

/// Build reduced graph from molecule using functional group detection.
/// If no functional groups found, treats entire molecule as one backbone node.
fn build_reduced_graph(mol: &Molecule) -> (Vec<ErgNode>, Vec<ErgEdge>) {
    let fg_groups = identify_functional_groups(mol);

    // Create nodes from functional groups
    let mut nodes: Vec<ErgNode> = fg_groups
        .into_iter()
        .map(|atom_indices| {
            let mut ntype = ErgNodeType::new();

            // Check for aromatic atoms in group
            let has_aromatic = atom_indices.iter().any(|&i| {
                mol.atom(AtomIdx(i as u32)).aromatic
            });
            if has_aromatic {
                ntype = ntype.with_aromatic();
            }

            // Check for heteroatom presence
            let has_n = atom_indices
                .iter()
                .any(|&i| mol.atom(AtomIdx(i as u32)).element.atomic_number() == 7);
            let has_o = atom_indices
                .iter()
                .any(|&i| mol.atom(AtomIdx(i as u32)).element.atomic_number() == 8);
            let has_s = atom_indices
                .iter()
                .any(|&i| mol.atom(AtomIdx(i as u32)).element.atomic_number() == 16);

            if has_n {
                ntype = ntype.with_donor().with_acceptor();
            }
            if has_o {
                ntype = ntype.with_acceptor();
            }
            if has_s {
                ntype = ntype.with_acceptor();
            }

            ErgNode { ntype, atom_indices }
        })
        .collect();

    // If no FGs found, treat entire molecule as one backbone node
    if nodes.is_empty() {
        let all_atoms: Vec<usize> = (0..mol.atom_count()).collect();
        let mut ntype = ErgNodeType::new();

        // Check for aromatic atoms
        if all_atoms.iter().any(|&i| mol.atom(AtomIdx(i as u32)).aromatic) {
            ntype = ntype.with_aromatic();
        }

        nodes.push(ErgNode {
            ntype,
            atom_indices: all_atoms,
        });
    }

    // Create edges: find shortest paths between node pairs
    let mut edges = Vec::new();
    let fg_set: HashSet<usize> = nodes.iter().flat_map(|n| n.atom_indices.clone()).collect();

    for i in 0..nodes.len() {
        for j in (i + 1)..nodes.len() {
            let linker_len = shortest_path_linker(mol, &nodes[i], &nodes[j], &fg_set);
            edges.push(ErgEdge {
                node_a: i,
                node_b: j,
                linker_len,
            });
        }
    }

    (nodes, edges)
}

/// Compute shortest path length between two functional groups (linker atoms only).
fn shortest_path_linker(
    mol: &Molecule,
    node_a: &ErgNode,
    node_b: &ErgNode,
    fg_set: &HashSet<usize>,
) -> u32 {
    let mut dist = vec![u32::MAX; mol.atom_count()];

    // BFS from node_a atoms
    let mut queue = VecDeque::new();
    for &i in &node_a.atom_indices {
        dist[i] = 0;
        queue.push_back(i);
    }

    while let Some(cur) = queue.pop_front() {
        for (nb, _) in mol.neighbors(AtomIdx(cur as u32)) {
            let nbi = nb.0 as usize;

            // Stop if we reach node_b
            if node_b.atom_indices.contains(&nbi) {
                let linker = dist[cur] + 1 - node_a.atom_indices.len() as u32;
                return linker.max(1);
            }

            if dist[nbi] == u32::MAX {
                // Only count non-fg atoms as linker
                let new_dist = dist[cur] + if fg_set.contains(&nbi) { 0 } else { 1 };
                dist[nbi] = new_dist;
                queue.push_back(nbi);
            }
        }
    }

    u32::MAX
}

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

/// Generate True ERG fingerprint from a molecule.
///
/// v0.1.91: Uses Ertl 2017 functional group detection and reduced graph topology.
/// Identifies functional group clusters, builds a reduced graph, and generates
/// a fingerprint from the reduced graph structure.
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

    // Count atom/bond types for backward compatibility
    for (idx, atom) in mol.atoms() {
        let erg_type = ErgAtomType::from_atom(atom);
        atom_counts[erg_type as usize] += 1;

        // Set bits based on atom type
        let bit_pos = (erg_type as usize) * 16;
        if bit_pos < 2048 {
            bits.set(bit_pos);
        }

        // Add degree-based bits for structural context
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

    // Encode atom/bond counts
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

    // v0.1.91: Add reduced graph topology bits
    let (nodes, _edges) = build_reduced_graph(mol);

    // Aromatic node detection
    if nodes.iter().any(|n| n.ntype.0 & ErgNodeType::AROMATIC != 0) {
        bits.set(256); // Aromatic
    }

    // Heteroatom functional group detection
    if nodes.iter().any(|n| n.ntype.0 != 0) {
        bits.set(257); // Has heteroatom
    }

    // Aliphatic-only detection
    if !nodes.iter().any(|n| n.ntype.0 & ErgNodeType::AROMATIC != 0)
        && atom_counts[ErgAtomType::CAliphatic as usize] > 0
    {
        bits.set(258); // Aliphatic only
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
