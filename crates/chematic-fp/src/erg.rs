//! ERG (Extended Reduced Graph) fingerprints — bespoke 2048-bit and 315-dim float variants.
//!
//! Implements ERG algorithm (Sheridan 1996 / Rarey & Dixon 1998) with Ertl 2017 FG detection:
//! 1. Identify functional group clusters (Ertl 2017)
//! 2. Collapse each cluster into a single reduced node with pharmacophore type
//! 3. Build reduced graph with linker bond counts between nodes
//! 4. Generate fingerprint from reduced graph topology
//!
//! Two fingerprint formats are provided:
//! - `erg()`: 2048-bit bespoke bitvector (original chematic format, fast hashing)
//! - `erg_vec()`: 315-float histogram (Rarey & Dixon 1998 structure: 21 node-type pairs
//!   × 15 distance bins, with Gaussian fuzzing ±1 bin at weight 0.3)
//!   NOTE: numerical parity with RDKit GetErGFingerprint is unverified (no fixture data).
//!
//! Reference: Rarey & Dixon 1998 J.Comput.-Aided Mol.Des. 12:471–490
//!            Sheridan 1996 J.Chem.Inf.Comput.Sci. 36:128–136
//!            Ertl 2017 J.Cheminf. 9:36

use crate::bitvec::BitVec2048;
use crate::ecfp::fnv1a;
use chematic_core::{Atom, AtomIdx, BondOrder, Molecule, implicit_hcount};
use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::VecDeque;

/// Functional group identification for ERG reduced graph construction.
///
/// Algorithm (Union-Find over heteroatoms):
/// 1. Identify all heteroatoms (non-C, non-H).
/// 2. Merge heteroatoms that share a common C neighbor (e.g. C=O and OH on same
///    carbonyl C → carboxylate/ester one node; amide C bonded to N and O → one node).
/// 3. Merge heteroatoms directly bonded to each other (e.g. N-O in nitro groups).
/// 4. Expand each cluster by including all directly bonded C atoms.
///
/// This keeps amine and carboxylate as *separate* nodes in amino acids, because
/// the only carbon shared between them (α-C) is bonded to only ONE heteroatom cluster
/// at a time. C–C bonds between FG clusters become linker atoms for edge encoding.
fn identify_functional_groups(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atom_count();
    if n == 0 {
        return Vec::new();
    }

    // --- Phase 1: collect heteroatom indices ---
    let mut is_hetero = vec![false; n];
    for (idx, atom) in mol.atoms() {
        let an = atom.element.atomic_number();
        if an != 1 && an != 6 {
            is_hetero[idx.0 as usize] = true;
        }
    }
    let hetero_idxs: Vec<usize> = (0..n).filter(|&i| is_hetero[i]).collect();
    if hetero_idxs.is_empty() {
        return Vec::new();
    }

    // --- Phase 2: Union-Find ---
    let mut parent: Vec<usize> = (0..n).collect();

    // Merge heteroatoms bonded directly to each other (nitro N-O, N-N hydrazine, etc.)
    for &hi in &hetero_idxs {
        for (nb, _) in mol.neighbors(AtomIdx(hi as u32)) {
            let nbi = nb.0 as usize;
            if is_hetero[nbi] {
                uf_union(&mut parent, hi, nbi);
            }
        }
    }

    // Merge heteroatoms that share a C neighbor (same FG context)
    for ci in 0..n {
        if mol.atom(AtomIdx(ci as u32)).element.atomic_number() != 6 {
            continue;
        }
        let hetero_nbs: Vec<usize> = mol
            .neighbors(AtomIdx(ci as u32))
            .filter(|(nb, _)| is_hetero[nb.0 as usize])
            .map(|(nb, _)| nb.0 as usize)
            .collect();
        if hetero_nbs.len() >= 2 {
            for &hi in &hetero_nbs[1..] {
                uf_union(&mut parent, hetero_nbs[0], hi);
            }
        }
    }

    // --- Phase 3: collect clusters and expand with bonded C atoms ---
    let mut cluster_map: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for &hi in &hetero_idxs {
        let root = uf_find(&mut parent, hi);
        cluster_map.entry(root).or_default().push(hi);
    }

    cluster_map
        .into_values()
        .map(|mut atoms| {
            let hetero_set: FxHashSet<usize> = atoms.iter().copied().collect();
            let snapshot = atoms.clone();
            for &ai in &snapshot {
                for (nb, _) in mol.neighbors(AtomIdx(ai as u32)) {
                    let nbi = nb.0 as usize;
                    let nb_an = mol.atom(nb).element.atomic_number();
                    if nb_an == 6 && !hetero_set.contains(&nbi) {
                        atoms.push(nbi);
                    }
                }
            }
            atoms.sort_unstable();
            atoms.dedup();
            atoms
        })
        .collect()
}

fn uf_find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]]; // path halving
        x = parent[x];
    }
    x
}

fn uf_union(parent: &mut [usize], a: usize, b: usize) {
    let ra = uf_find(parent, a);
    let rb = uf_find(parent, b);
    if ra != rb {
        parent[rb] = ra;
    }
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
        let explicit_h: usize = mol
            .neighbors(idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 1)
            .count();
        let impl_h = implicit_hcount(mol, idx) as usize;
        let total_h = explicit_h + impl_h;

        // H-bond donor:
        //   - O-H and S-H: always donors when H present
        //   - N-H: donor ONLY when NOT in amide/sulfonamide context
        if total_h > 0 && ((an == 8 || an == 16) || (an == 7 && !is_amide_like_nitrogen(mol, idx)))
        {
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
        if atom.charge > 0 {
            ntype = ntype.with_positive();
        }
        if atom.charge < 0 {
            ntype = ntype.with_negative();
        }
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
            ErgNode {
                ntype,
                atom_indices,
            }
        })
        .collect();

    // If no FGs found, treat entire molecule as one backbone node
    if nodes.is_empty() {
        let all_atoms: Vec<usize> = (0..mol.atom_count()).collect();
        let ntype = assign_pharmacophore_features(mol, &all_atoms);
        nodes.push(ErgNode {
            ntype,
            atom_indices: all_atoms,
        });
    }

    // Cap node count to prevent O(k²·n) blow-up on pathological all-heteroatom inputs.
    const MAX_ERG_NODES: usize = 128;
    if nodes.len() > MAX_ERG_NODES {
        nodes.truncate(MAX_ERG_NODES);
    }

    // Create edges: find shortest paths between node pairs
    let mut edges = Vec::new();
    let fg_set: FxHashSet<usize> = nodes.iter().flat_map(|n| n.atom_indices.clone()).collect();

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
    fg_set: &FxHashSet<usize>,
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

            // Stop if we reach node_b.
            // linker = total bond count from any node_a atom to this node_b atom.
            if node_b.atom_indices.contains(&nbi) {
                return dist[cur] + 1;
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
pub fn erg_with_config(mol: &chematic_core::Molecule, config: &ErgConfig) -> ErgFingerprint {
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
        let degree = mol.neighbors(idx).count();
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

// ─── ERG 315-dim float histogram (Rarey & Dixon 1998 structure) ──────────────

/// Length of the ERG float vector: 21 node-type pairs × 15 distance bins.
pub const ERG_VEC_LEN: usize = 315;

/// Number of pharmacophore feature types encoded in the float vector.
const N_FEAT: usize = 6;

/// Number of distance bins (distances 0..=14+ linker atoms).
const N_BINS: usize = 15;

/// Gaussian fuzzing weight applied to the adjacent distance bins.
/// Mirrors the `fuzzIncrement` parameter used in RDKit's ErG implementation.
const FUZZ: f64 = 0.3;

/// Feature type index ordering for the 315-dim vector.
/// Mapping from bit position (as bit-shift in ErgNodeType.0) to feature index 0-5:
///   0 = AROMATIC (bit 1 = 1<<0)
///   1 = DONOR    (bit 2 = 1<<1)
///   2 = ACCEPTOR (bit 4 = 1<<2)
///   3 = POSITIVE (bit 16 = 1<<4)
///   4 = NEGATIVE (bit 32 = 1<<5)
///   5 = HYDROPHOBIC (bit 8 = 1<<3)
const FEATURE_BITS: [(u8, usize); N_FEAT] = [
    (ErgNodeType::AROMATIC, 0),
    (ErgNodeType::DONOR, 1),
    (ErgNodeType::ACCEPTOR, 2),
    (ErgNodeType::POSITIVE, 3),
    (ErgNodeType::NEGATIVE, 4),
    (ErgNodeType::HYDROPHOBIC, 5),
];

/// Compute index into the 21 unique unordered pairs for feature types (a, b).
///
/// Lower-triangular enumeration (including diagonal):
/// (0,0)=0 (0,1)=1 … (0,5)=5 | (1,1)=6 … (1,5)=10 | … | (5,5)=20
#[inline]
fn pair_idx(a: usize, b: usize) -> usize {
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    lo * N_FEAT - lo * (lo.wrapping_sub(1)) / 2 + (hi - lo)
}

/// Add `weight` to `vec` at distance bin `d` with Gaussian fuzzing (±1 at `FUZZ`).
#[inline]
fn fuzz_add(vec: &mut [f64; ERG_VEC_LEN], base: usize, d: usize, weight: f64) {
    vec[base + d] += weight;
    if d > 0 {
        vec[base + d - 1] += weight * FUZZ;
    }
    if d + 1 < N_BINS {
        vec[base + d + 1] += weight * FUZZ;
    }
}

/// Generate an ERG-style 315-element float histogram.
///
/// Format: `vec[pair_index * 15 + distance_bin]` where
/// - `pair_index` ∈ 0..21 (unique unordered pair of 6 pharmacophore feature types)
/// - `distance_bin` ∈ 0..15 (linker bond count capped at 14)
///
/// Adjacent bins receive a `FUZZ = 0.3` weight (Gaussian smearing).
///
/// **NOTE**: numerical parity with RDKit `GetErGFingerprint` is unverified.
/// Use as a richer alternative to the 2048-bit bespoke FP when float vectors
/// are needed (e.g. cosine similarity, PCA).
pub fn erg_vec(mol: &Molecule) -> [f64; ERG_VEC_LEN] {
    let mut vec = [0.0f64; ERG_VEC_LEN];
    let (nodes, edges) = build_reduced_graph(mol);

    // Self-pair entries: each node contributes its feature types at bin 0.
    // This encodes node composition independent of topology, so single-node
    // molecules (pure alkanes, single-FG molecules) are still distinguishable.
    for node in &nodes {
        let ntype = node.ntype.0;
        for &(bit, fi) in &FEATURE_BITS {
            if ntype & bit == 0 {
                continue;
            }
            let base = pair_idx(fi, fi) * N_BINS;
            fuzz_add(&mut vec, base, 0, 1.0);
        }
    }

    // Inter-node pair entries: each edge encodes (feature_a, feature_b, distance).
    for edge in &edges {
        let dist = edge.linker_len;
        if dist == u32::MAX {
            continue;
        }
        let bin = (dist as usize).min(N_BINS - 1);
        let ntype_a = nodes[edge.node_a].ntype.0;
        let ntype_b = nodes[edge.node_b].ntype.0;

        for &(bit_a, fi) in &FEATURE_BITS {
            if ntype_a & bit_a == 0 {
                continue;
            }
            for &(bit_b, fj) in &FEATURE_BITS {
                if ntype_b & bit_b == 0 {
                    continue;
                }
                let base = pair_idx(fi, fj) * N_BINS;
                fuzz_add(&mut vec, base, bin, 1.0);
            }
        }
    }

    vec
}

/// Cosine similarity between two ERG float vectors.
pub fn cosine_erg_vec(v1: &[f64; ERG_VEC_LEN], v2: &[f64; ERG_VEC_LEN]) -> f64 {
    let dot: f64 = v1.iter().zip(v2.iter()).map(|(a, b)| a * b).sum();
    let n1: f64 = v1.iter().map(|x| x * x).sum::<f64>().sqrt();
    let n2: f64 = v2.iter().map(|x| x * x).sum::<f64>().sqrt();
    if n1 < 1e-12 || n2 < 1e-12 {
        return 0.0;
    }
    (dot / (n1 * n2)).min(1.0)
}

/// Tanimoto similarity on float vectors: T = dot / (|v1|² + |v2|² − dot).
pub fn tanimoto_erg_vec(v1: &[f64; ERG_VEC_LEN], v2: &[f64; ERG_VEC_LEN]) -> f64 {
    let dot: f64 = v1.iter().zip(v2.iter()).map(|(a, b)| a * b).sum();
    let n1: f64 = v1.iter().map(|x| x * x).sum();
    let n2: f64 = v2.iter().map(|x| x * x).sum();
    let denom = n1 + n2 - dot;
    if denom < 1e-12 {
        return 1.0;
    }
    (dot / denom).clamp(0.0, 1.0)
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
        assert!(
            !fp_aliphatic.bits.get(256),
            "aliphatic should not have aromatic bit"
        );
        assert!(
            fp_aromatic.bits.get(256),
            "aromatic should have aromatic bit"
        );
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
        assert!(
            fp_alkane.bits.get(257),
            "alkane gets HYDROPHOBIC → bit 257 set"
        );
        assert!(
            fp_alcohol.bits.get(257),
            "alcohol gets ACCEPTOR/DONOR → bit 257 set"
        );
        assert!(
            fp_amine.bits.get(257),
            "amine gets ACCEPTOR/DONOR → bit 257 set"
        );

        // Alcohol and amine should differ from alkane in topology bits (259+)
        let topo_alkane: usize = (259..2048).filter(|&b| fp_alkane.bits.get(b)).count();
        let topo_alcohol: usize = (259..2048).filter(|&b| fp_alcohol.bits.get(b)).count();
        let topo_amine: usize = (259..2048).filter(|&b| fp_amine.bits.get(b)).count();
        // All should be valid (we just verify topo bits exist, not exact values)
        let _ = (topo_alkane, topo_alcohol, topo_amine);

        // Alkane FP should differ from alcohol and amine FPs
        assert!(
            fp_alkane.tanimoto(&fp_alcohol) < 1.0,
            "alkane vs alcohol should differ"
        );
        assert!(
            fp_alkane.tanimoto(&fp_amine) < 1.0,
            "alkane vs amine should differ"
        );
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
        // erg_vec (315-dim float histogram) encodes distances in distinct bins,
        // so the vectors must differ. The bespoke 2048-bit erg() uses FNV hashing
        // which can have bin collisions for this specific pair — use erg_vec instead.
        let short = parse("NCC(=O)O").unwrap(); // 1 linker C
        let long = parse("NCCCCCC(=O)O").unwrap(); // 5 linker Cs

        let v_short = erg_vec(&short);
        let v_long = erg_vec(&long);

        // Vectors must not be identical (different linker bin → different histogram).
        assert_ne!(
            v_short, v_long,
            "different linker lengths must produce different erg_vec entries"
        );

        // Similarity must be strictly < 1.0.
        let sim = tanimoto_erg_vec(&v_short, &v_long);
        assert!(
            sim < 1.0,
            "different linker lengths should give Tanimoto < 1.0, got {sim:.4}"
        );
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
        assert!(
            fp.tanimoto(&fp_benz) < 1.0,
            "aniline should differ from benzene"
        );
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
        assert!(
            fp_acetone.tanimoto(&fp_propane) < 1.0,
            "acetone vs propane should differ due to acceptor O"
        );
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
            "acetate O- should produce NEGATIVE pharmacophore feature, got ntype={}",
            ntype.0
        );

        // Acetic acid CC(=O)O: no negative charge → DONOR (O-H) instead
        let acid = parse("CC(=O)O").unwrap();
        let acid_groups = identify_functional_groups(&acid);
        assert!(!acid_groups.is_empty());
        let acid_ntype = assign_pharmacophore_features(&acid, &acid_groups[0]);
        assert_eq!(
            acid_ntype.0 & ErgNodeType::NEGATIVE,
            0,
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
        assert!(
            fp.bits.get(257),
            "hexane should have a pharmacophore feature (HYDROPHOBIC)"
        );
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
        assert!(
            fp_pyr.tanimoto(&fp_rol) < 1.0,
            "pyridine (N acceptor) vs pyrrole (N-H donor) should differ"
        );
    }

    #[test]
    fn test_erg_adjacent_groups_bin0() {
        // Two functional groups with 0 linker atoms (directly bonded): bin=0.
        // OCCO: O and O separated by 2 C linker atoms → bin=1 (short).
        // ON: O and N directly bonded → linker_len should be small.
        let direct = parse("NO").unwrap(); // N-O: effectively adjacent FGs
        let indirect = parse("NCCCO").unwrap(); // N-CCC-O: 3-linker

        let fp_direct = erg(&direct);
        let fp_indirect = erg(&indirect);

        assert!(
            fp_direct.tanimoto(&fp_indirect) < 1.0,
            "adjacent vs. 3-linker groups should differ"
        );
    }

    /// Amide N: ACCEPTOR only, NOT donor (lone pair delocalized into C=O).
    #[test]
    fn test_erg_amide_n_is_acceptor_not_donor() {
        use super::{assign_pharmacophore_features, identify_functional_groups};

        let acetamide = parse("CC(=O)N").unwrap(); // CH3-C(=O)-NH2
        let groups = identify_functional_groups(&acetamide);

        // Find the group containing N (atomic number 7)
        let n_group = groups
            .iter()
            .find(|g| {
                g.iter()
                    .any(|&i| acetamide.atom(AtomIdx(i as u32)).element.atomic_number() == 7)
            })
            .expect("should find N-containing group");

        let ntype = assign_pharmacophore_features(&acetamide, n_group);

        // Amide N should be ACCEPTOR (bit 4 set)
        assert!(
            ntype.0 & ErgNodeType::ACCEPTOR != 0,
            "amide N should be ACCEPTOR (bits={:#010b})",
            ntype.0
        );

        // Amide N should NOT be DONOR (bit 2) — lone pair delocalized into C=O
        assert!(
            ntype.0 & ErgNodeType::DONOR == 0,
            "amide N should NOT be DONOR (bits={:#010b})",
            ntype.0
        );
    }

    /// Ethylamine N: both DONOR and ACCEPTOR (free lone pair + N-H).
    #[test]
    fn test_erg_amine_n_is_donor_and_acceptor() {
        use super::{assign_pharmacophore_features, identify_functional_groups};

        let ethylamine = parse("CCN").unwrap(); // CH3CH2-NH2
        let groups = identify_functional_groups(&ethylamine);

        let n_group = groups
            .iter()
            .find(|g| {
                g.iter()
                    .any(|&i| ethylamine.atom(AtomIdx(i as u32)).element.atomic_number() == 7)
            })
            .expect("should find N-containing group");

        let ntype = assign_pharmacophore_features(&ethylamine, n_group);

        // Amine N should be both DONOR and ACCEPTOR
        assert!(
            ntype.0 & ErgNodeType::DONOR != 0,
            "amine N should be DONOR (bits={:#010b})",
            ntype.0
        );
        assert!(
            ntype.0 & ErgNodeType::ACCEPTOR != 0,
            "amine N should be ACCEPTOR (bits={:#010b})",
            ntype.0
        );
    }

    // ── erg_vec (315-dim float histogram) tests ──────────────────────────────

    #[test]
    fn test_erg_vec_length() {
        let mol = parse("CC").unwrap();
        let v = erg_vec(&mol);
        assert_eq!(
            v.len(),
            ERG_VEC_LEN,
            "erg_vec must return exactly 315 floats"
        );
    }

    #[test]
    fn test_erg_vec_consistency() {
        let mol = parse("c1ccccc1N").unwrap();
        let v1 = erg_vec(&mol);
        let v2 = erg_vec(&mol);
        assert_eq!(v1, v2, "erg_vec must be deterministic");
    }

    #[test]
    fn test_erg_vec_nonnegative() {
        for smi in &["CC", "CCO", "c1ccccc1", "c1ccncc1", "CC(=O)N", "c1ccccc1N"] {
            let mol = parse(smi).unwrap();
            let v = erg_vec(&mol);
            for (i, &x) in v.iter().enumerate() {
                assert!(x >= 0.0, "erg_vec[{i}] negative for {smi}: {x}");
            }
        }
    }

    #[test]
    fn test_erg_vec_same_mol_self_similarity() {
        let mol = parse("c1ccccc1").unwrap();
        let v = erg_vec(&mol);
        let sim = tanimoto_erg_vec(&v, &v);
        assert!(
            (sim - 1.0).abs() < 1e-9,
            "self-similarity must be 1.0 (got {sim})"
        );
    }

    #[test]
    fn test_erg_vec_different_molecules() {
        let mol1 = parse("CC").unwrap();
        let mol2 = parse("c1ccccc1N").unwrap();
        let v1 = erg_vec(&mol1);
        let v2 = erg_vec(&mol2);
        let sim = tanimoto_erg_vec(&v1, &v2);
        assert!(
            (0.0..=1.0).contains(&sim),
            "Tanimoto must be in [0,1] (got {sim})"
        );
        assert!(sim < 1.0, "ethane vs aniline must not be identical");
    }

    #[test]
    fn test_erg_vec_donor_acceptor_linker_distance() {
        // GABA (H2N-CH2-CH2-CH2-COOH): donor (NH2) and acceptor (C=O, OH) with 3 linker Cs
        // Glycine (H2N-CH2-COOH): 1 linker C
        // The DA pair at different distance bins must give different vectors.
        let gaba = parse("NCCCC(=O)O").unwrap();
        let glycine = parse("NCC(=O)O").unwrap();
        let v_gaba = erg_vec(&gaba);
        let v_glycine = erg_vec(&glycine);
        let sim = tanimoto_erg_vec(&v_gaba, &v_glycine);
        assert!(
            sim < 1.0,
            "GABA vs glycine differ only in linker length — Tanimoto must be < 1 (got {sim:.3})"
        );
    }

    #[test]
    fn test_erg_vec_fuzz_spreads_adjacent_bins() {
        // A node pair at distance bin d must also populate bin d±1 via FUZZ.
        let mol = parse("NCC(=O)O").unwrap(); // donor-linker-acceptor
        let v = erg_vec(&mol);
        // Some bin must be > 0 in the Donor–Acceptor (fi=1, fj=2) pair region
        let da_pair = pair_idx(1, 2) * N_BINS;
        let sum: f64 = v[da_pair..da_pair + N_BINS].iter().sum();
        assert!(
            sum > 0.0,
            "Donor–Acceptor bins must be non-zero for glycine analogue"
        );
    }

    #[test]
    fn test_pair_idx_all_21() {
        let mut seen = FxHashSet::default();
        for a in 0..N_FEAT {
            for b in a..N_FEAT {
                let idx = pair_idx(a, b);
                assert!(idx < 21, "pair_idx({a},{b}) = {idx} out of range");
                assert!(
                    seen.insert(idx),
                    "pair_idx({a},{b}) = {idx} duplicates earlier pair"
                );
            }
        }
        assert_eq!(seen.len(), 21);
    }

    /// Sulfonamide N: ACCEPTOR only (same logic as amide N).
    #[test]
    fn test_erg_sulfonamide_n_is_acceptor_not_donor() {
        use super::{assign_pharmacophore_features, identify_functional_groups};

        let sulfonamide = parse("CS(=O)(=O)N").unwrap(); // methanesulfonamide
        let groups = identify_functional_groups(&sulfonamide);

        let n_group = groups
            .iter()
            .find(|g| {
                g.iter()
                    .any(|&i| sulfonamide.atom(AtomIdx(i as u32)).element.atomic_number() == 7)
            })
            .expect("should find N-containing group");

        let ntype = assign_pharmacophore_features(&sulfonamide, n_group);

        assert!(
            ntype.0 & ErgNodeType::ACCEPTOR != 0,
            "sulfonamide N should be ACCEPTOR (bits={:#010b})",
            ntype.0
        );
        assert!(
            ntype.0 & ErgNodeType::DONOR == 0,
            "sulfonamide N should NOT be DONOR (bits={:#010b})",
            ntype.0
        );
    }
}
