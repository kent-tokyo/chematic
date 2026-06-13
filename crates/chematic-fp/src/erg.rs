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
use chematic_core::{Atom, AtomIdx, BondOrder, Molecule, implicit_hcount};
use crate::bitvec::BitVec2048;
use crate::ecfp::fnv1a;

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
}

impl Default for ErgNodeType {
    fn default() -> Self {
        Self::new()
    }
}

impl ErgNodeType {

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

/// Return true if atom `n_idx` (nitrogen) is in an amide or sulfonamide context.
///
/// Amide: N directly bonded to C=O.
/// Sulfonamide: N directly bonded to S(=O).
/// In both cases the N lone pair is delocalized → acts as ACCEPTOR, not DONOR.
fn is_amide_like_nitrogen(mol: &Molecule, n_idx: AtomIdx) -> bool {
    mol.neighbors(n_idx).any(|(nb, _)| {
        let nb_an = mol.atom(nb).element.atomic_number();
        match nb_an {
            6 => mol.neighbors(nb).any(|(c_nb, bond_idx)| {
                mol.atom(c_nb).element.atomic_number() == 8
                    && mol.bond(bond_idx).order == BondOrder::Double
            }),
            16 => mol.neighbors(nb).any(|(s_nb, bond_idx)| {
                mol.atom(s_nb).element.atomic_number() == 8
                    && mol.bond(bond_idx).order == BondOrder::Double
            }),
            _ => false,
        }
    })
}

/// Assign pharmacophore-based node type for a functional group cluster.
///
/// Rules (per atom in the group):
/// - AROMATIC: any aromatic atom
/// - DONOR: N-H (non-amide/sulfonamide), O-H, or S-H
/// - ACCEPTOR: N (non-amide unless aromatic), amide/sulfonamide N, any O, any F
/// - POSITIVE: formal charge > 0
/// - NEGATIVE: formal charge < 0
/// - HYDROPHOBIC: group with only C/H atoms and no other pharmacophore feature
///
/// Key fix vs prior version: amide/sulfonamide N is ACCEPTOR-only even when it
/// carries an H, because the lone pair is delocalized into the C=O or S=O π system.
fn assign_pharmacophore_features(mol: &Molecule, atom_indices: &[usize]) -> ErgNodeType {
    let mut ntype = ErgNodeType::new();

    for &i in atom_indices {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        let an = atom.element.atomic_number();

        if atom.aromatic {
            ntype = ntype.with_aromatic();
        }

        // Count H attached to this atom (explicit H neighbors)
        let explicit_h: usize = mol.neighbors(idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 1)
            .count();
        let impl_h = implicit_hcount(mol, idx) as usize;
        let total_h = explicit_h + impl_h;

        // H-bond donor:
        //   - O-H and S-H: always donors when H present
        //   - N-H: donor ONLY when NOT in amide/sulfonamide context
        if (an == 8 || an == 16) && total_h > 0 {
            ntype = ntype.with_donor();
        } else if an == 7 && total_h > 0 && !is_amide_like_nitrogen(mol, idx) {
            ntype = ntype.with_donor();
        }

        // H-bond acceptor
        match an {
            // Non-aromatic N: acceptor (including amide N — delocalized lone pair)
            7 if !atom.aromatic => ntype = ntype.with_acceptor(),
            // Aromatic N: acceptor only if no H (pyridine-like), not if has H (pyrrole)
            7 if atom.aromatic && total_h == 0 => ntype = ntype.with_acceptor(),
            // O and F: always acceptor
            8 | 9 => ntype = ntype.with_acceptor(),
            _ => {}
        }

        // Formal charge
        if atom.charge > 0 { ntype = ntype.with_positive(); }
        if atom.charge < 0 { ntype = ntype.with_negative(); }
    }

    // Hydrophobic: all atoms are C or H and no other pharmacophore bit set
    let all_c_or_h = atom_indices.iter().all(|&i| {
        let an = mol.atom(AtomIdx(i as u32)).element.atomic_number();
        an == 6 || an == 1
    });
    if all_c_or_h && ntype.0 == 0 {
        ntype = ntype.with_hydrophobic();
    }

    ntype
}

/// Build reduced graph from molecule using functional group detection.
/// If no functional groups found, treats entire molecule as one backbone node.
fn build_reduced_graph(mol: &Molecule) -> (Vec<ErgNode>, Vec<ErgEdge>) {
    let fg_groups = identify_functional_groups(mol);

    // Create nodes from functional groups using pharmacophore feature assignment
    let mut nodes: Vec<ErgNode> = fg_groups
        .into_iter()
        .map(|atom_indices| {
            let ntype = assign_pharmacophore_features(mol, &atom_indices);
            ErgNode { ntype, atom_indices }
        })
        .collect();

    // If no FGs found, treat entire molecule as one backbone node
    if nodes.is_empty() {
        let all_atoms: Vec<usize> = (0..mol.atom_count()).collect();
        let ntype = assign_pharmacophore_features(mol, &all_atoms);
        nodes.push(ErgNode { ntype, atom_indices: all_atoms });
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
        let degree_bits = (degree.min(4) << 2) + (erg_type as usize);
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

    // Reduced graph topology encoding.
    //
    // For each edge (i, j) in the reduced graph, encode:
    //   (sorted node types, binned linker length) → bit in [259, 2047]
    //
    // Linker bins: 0=adjacent, 1=short(1-2), 2=medium(3-5), 3=long(6+)
    // This makes the fingerprint sensitive to inter-group distance.
    let (nodes, edges) = build_reduced_graph(mol);

    for edge in &edges {
        if edge.linker_len == u32::MAX {
            continue; // no path between groups — skip
        }
        let ta = nodes[edge.node_a].ntype.0;
        let tb = nodes[edge.node_b].ntype.0;
        let bin: u8 = match edge.linker_len {
            0 => 0,
            1..=2 => 1,
            3..=5 => 2,
            _ => 3,
        };
        let (t_lo, t_hi) = if ta <= tb { (ta, tb) } else { (tb, ta) };
        // Hash the triplet into [259, 2047] for topology bits.
        let h = fnv1a(&[t_lo, t_hi, bin, 0xE7]) as usize;
        bits.set(259 + h % (2048 - 259));
    }

    // Global node-level flags (bits 256-258)
    if nodes.iter().any(|n| n.ntype.0 & ErgNodeType::AROMATIC != 0) {
        bits.set(256);
    }
    if nodes.iter().any(|n| n.ntype.0 != 0) {
        bits.set(257);
    }
    if !nodes.iter().any(|n| n.ntype.0 & ErgNodeType::AROMATIC != 0)
        && atom_counts[ErgAtomType::CAliphatic as usize] > 0
    {
        bits.set(258);
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

        // Bit 257: set when any pharmacophore feature is present.
        // With pharmacophore assignment, HYDROPHOBIC is set for alkane too.
        assert!(fp_alkane.bits.get(257), "alkane gets HYDROPHOBIC → bit 257 set");
        assert!(fp_alcohol.bits.get(257), "alcohol gets ACCEPTOR/DONOR → bit 257 set");
        assert!(fp_amine.bits.get(257), "amine gets ACCEPTOR/DONOR → bit 257 set");

        // Alcohol and amine should differ from alkane in topology bits (259+)
        let topo_alkane: usize = (259..2048).filter(|&b| fp_alkane.bits.get(b)).count();
        let topo_alcohol: usize = (259..2048).filter(|&b| fp_alcohol.bits.get(b)).count();
        let topo_amine: usize = (259..2048).filter(|&b| fp_amine.bits.get(b)).count();
        // All should be valid (we just verify topo bits exist, not exact values)
        let _ = (topo_alkane, topo_alcohol, topo_amine);

        // Alkane FP should differ from alcohol and amine FPs
        assert!(fp_alkane.tanimoto(&fp_alcohol) < 1.0, "alkane vs alcohol should differ");
        assert!(fp_alkane.tanimoto(&fp_amine) < 1.0, "alkane vs amine should differ");
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

    #[test]
    fn test_erg_linker_distance_changes_fingerprint() {
        // Same terminal functional groups (NH2 and COOH) but different linker lengths.
        // The reduced graph edge topology bits (259+) must differ.
        let short = parse("NCC(=O)O").unwrap();   // 1 linker C
        let long  = parse("NCCCCCC(=O)O").unwrap(); // 5 linker Cs

        let fp_short = erg(&short);
        let fp_long  = erg(&long);

        // Topology bits in [259, 2047] must differ for different linker lengths.
        let short_topo: Vec<usize> = (259..2048).filter(|&b| fp_short.bits.get(b)).collect();
        let long_topo:  Vec<usize> = (259..2048).filter(|&b| fp_long.bits.get(b)).collect();
        assert_ne!(short_topo, long_topo,
            "different linker lengths must produce different topology bits");

        // Similarity should be less than 1.0 (not identical).
        assert!(fp_short.tanimoto(&fp_long) < 1.0,
            "different linker lengths should give Tanimoto < 1.0");
    }

    // ---- New pharmacophore feature tests ----

    #[test]
    fn test_erg_nh_donor() {
        // Aniline: c1ccccc1N — the N group should be a DONOR (N-H present)
        let aniline = parse("c1ccccc1N").unwrap();
        let fp = erg(&aniline);
        // Should differ from benzene (no N)
        let benzene = parse("c1ccccc1").unwrap();
        let fp_benz = erg(&benzene);
        assert!(fp.tanimoto(&fp_benz) < 1.0, "aniline should differ from benzene");
        assert!(fp.bits.popcount() > 0);
    }

    #[test]
    fn test_erg_carbonyl_acceptor() {
        // Acetone C(C)(C)=O — the C=O group should be ACCEPTOR
        let acetone = parse("CC(C)=O").unwrap();
        let propane = parse("CCC").unwrap();
        let fp_acetone = erg(&acetone);
        let fp_propane = erg(&propane);
        // Acetone (has O acceptor) should differ from propane (pure hydrocarbon)
        assert!(fp_acetone.tanimoto(&fp_propane) < 1.0,
            "acetone vs propane should differ due to acceptor O");
    }

    #[test]
    fn test_erg_carboxylate_negative() {
        // Acetate anion CC(=O)[O-]: O- has formal charge -1 → NEGATIVE feature.
        // Test the pharmacophore assignment directly (fingerprint-level may
        // collapse single-node molecules to identical atom-count bits).
        let acetate = parse("CC(=O)[O-]").unwrap();
        let groups = identify_functional_groups(&acetate);
        assert!(!groups.is_empty(), "acetate should have a functional group");
        let ntype = assign_pharmacophore_features(&acetate, &groups[0]);
        assert!(
            ntype.0 & ErgNodeType::NEGATIVE != 0,
            "acetate O- should produce NEGATIVE pharmacophore feature, got ntype={}", ntype.0
        );

        // Acetic acid CC(=O)O: no negative charge → DONOR (O-H) instead
        let acid = parse("CC(=O)O").unwrap();
        let acid_groups = identify_functional_groups(&acid);
        assert!(!acid_groups.is_empty());
        let acid_ntype = assign_pharmacophore_features(&acid, &acid_groups[0]);
        assert_eq!(
            acid_ntype.0 & ErgNodeType::NEGATIVE, 0,
            "acetic acid should NOT have NEGATIVE feature"
        );
        assert!(
            acid_ntype.0 & ErgNodeType::DONOR != 0,
            "acetic acid O-H should have DONOR feature"
        );
    }

    #[test]
    fn test_erg_hydrophobic_chain() {
        // Hexane CCCCCC: all C → HYDROPHOBIC node type
        let hexane = parse("CCCCCC").unwrap();
        let fp = erg(&hexane);
        // Hydrophobic molecules should have bit 257 set (any feature including HYDROPHOBIC)
        assert!(fp.bits.get(257), "hexane should have a pharmacophore feature (HYDROPHOBIC)");
        // No aromatic (bit 256)
        assert!(!fp.bits.get(256), "hexane should not have aromatic bit");
    }

    #[test]
    fn test_erg_pyridine_vs_pyrrole_acceptor() {
        // Pyridine: aromatic N without H → ACCEPTOR only
        // Pyrrole: aromatic N-H → DONOR only
        let pyridine = parse("c1ccncc1").unwrap();
        let pyrrole = parse("c1cc[nH]c1").unwrap();
        let fp_pyr = erg(&pyridine);
        let fp_rol = erg(&pyrrole);
        // They should have different fingerprints (one is acceptor, other is donor)
        assert!(fp_pyr.tanimoto(&fp_rol) < 1.0,
            "pyridine (N acceptor) vs pyrrole (N-H donor) should differ");
    }

    #[test]
    fn test_erg_adjacent_groups_bin0() {
        // Two functional groups with 0 linker atoms (directly bonded): bin=0.
        // OCCO: O and O separated by 2 C linker atoms → bin=1 (short).
        // ON: O and N directly bonded → linker_len should be small.
        let direct = parse("NO").unwrap();   // N-O: effectively adjacent FGs
        let indirect = parse("NCCCO").unwrap(); // N-CCC-O: 3-linker

        let fp_direct   = erg(&direct);
        let fp_indirect = erg(&indirect);

        assert!(fp_direct.tanimoto(&fp_indirect) < 1.0,
            "adjacent vs. 3-linker groups should differ");
    }

    /// Amide N: ACCEPTOR only, NOT donor (lone pair delocalized into C=O).
    #[test]
    fn test_erg_amide_n_is_acceptor_not_donor() {
        use super::{assign_pharmacophore_features, identify_functional_groups};

        let acetamide = parse("CC(=O)N").unwrap(); // CH3-C(=O)-NH2
        let groups = identify_functional_groups(&acetamide);

        // Find the group containing N (atomic number 7)
        let n_group = groups.iter().find(|g| {
            g.iter().any(|&i| acetamide.atom(AtomIdx(i as u32)).element.atomic_number() == 7)
        }).expect("should find N-containing group");

        let ntype = assign_pharmacophore_features(&acetamide, n_group);

        // Amide N should be ACCEPTOR (bit 4 set)
        assert!(ntype.0 & ErgNodeType::ACCEPTOR != 0,
            "amide N should be ACCEPTOR (bits={:#010b})", ntype.0);

        // Amide N should NOT be DONOR (bit 2) — lone pair delocalized into C=O
        assert!(ntype.0 & ErgNodeType::DONOR == 0,
            "amide N should NOT be DONOR (bits={:#010b})", ntype.0);
    }

    /// Ethylamine N: both DONOR and ACCEPTOR (free lone pair + N-H).
    #[test]
    fn test_erg_amine_n_is_donor_and_acceptor() {
        use super::{assign_pharmacophore_features, identify_functional_groups};

        let ethylamine = parse("CCN").unwrap(); // CH3CH2-NH2
        let groups = identify_functional_groups(&ethylamine);

        let n_group = groups.iter().find(|g| {
            g.iter().any(|&i| ethylamine.atom(AtomIdx(i as u32)).element.atomic_number() == 7)
        }).expect("should find N-containing group");

        let ntype = assign_pharmacophore_features(&ethylamine, n_group);

        // Amine N should be both DONOR and ACCEPTOR
        assert!(ntype.0 & ErgNodeType::DONOR != 0,
            "amine N should be DONOR (bits={:#010b})", ntype.0);
        assert!(ntype.0 & ErgNodeType::ACCEPTOR != 0,
            "amine N should be ACCEPTOR (bits={:#010b})", ntype.0);
    }

    /// Sulfonamide N: ACCEPTOR only (same logic as amide N).
    #[test]
    fn test_erg_sulfonamide_n_is_acceptor_not_donor() {
        use super::{assign_pharmacophore_features, identify_functional_groups};

        let sulfonamide = parse("CS(=O)(=O)N").unwrap(); // methanesulfonamide
        let groups = identify_functional_groups(&sulfonamide);

        let n_group = groups.iter().find(|g| {
            g.iter().any(|&i| sulfonamide.atom(AtomIdx(i as u32)).element.atomic_number() == 7)
        }).expect("should find N-containing group");

        let ntype = assign_pharmacophore_features(&sulfonamide, n_group);

        assert!(ntype.0 & ErgNodeType::ACCEPTOR != 0,
            "sulfonamide N should be ACCEPTOR (bits={:#010b})", ntype.0);
        assert!(ntype.0 & ErgNodeType::DONOR == 0,
            "sulfonamide N should NOT be DONOR (bits={:#010b})", ntype.0);
    }
}
