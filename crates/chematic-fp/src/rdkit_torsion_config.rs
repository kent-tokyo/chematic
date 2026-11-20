//! Configurable RDKit topological torsions.

use std::fmt;

use chematic_core::{AtomIdx, Molecule};
use rustc_hash::FxHashMap;

use crate::bitvec::BitVecN;
use crate::rdkit_count_simulation::{
    COUNT_BOUNDS, ConfiguredCounts, MAX_CONFIGURED_BITS, valid_count_bounds,
};
use crate::rdkit_morgan_hash::hash_vec;
use crate::rdkit_torsion::{
    CODE_SIZE, atom_code, chirality_bits, legacy_atom_stereo, orient_torsion_codes,
    pack_torsion_codes, torsion_endpoint_code,
};

const MAX_ATOMS: usize = (u64::BITS as usize - 1) / CODE_SIZE as usize;

/// Options for the RDKit Topological Torsion generator.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitTorsionConfig {
    /// Number of atom positions, including the repeated end of a closed ring.
    pub atom_count: usize,
    /// Count only minimum-distance paths between each pair of endpoints.
    pub only_shortest_paths: bool,
    /// Apply RDKit's legacy R/S atom invariants and wider sparse encoding.
    pub include_chirality: bool,
    pub nbits: usize,
    /// Count thresholds in bit order. `None` disables count simulation.
    pub count_bounds: Option<Vec<u32>>,
}

impl Default for RdkitTorsionConfig {
    fn default() -> Self {
        Self {
            atom_count: 4,
            only_shortest_paths: false,
            include_chirality: false,
            nbits: 2048,
            count_bounds: Some(COUNT_BOUNDS.to_vec()),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitTorsionResult {
    pub fingerprint: BitVecN,
    pub folded_counts: Vec<(u32, u32)>,
    pub sparse_counts: Vec<(u64, u32)>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RdkitTorsionError {
    InvalidAtomCount { maximum: usize },
    InvalidBitWidth,
    InvalidCountBounds,
    Chirality(chematic_smiles::RdkitSmilesError),
}

impl fmt::Display for RdkitTorsionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomCount { maximum } => {
                write!(f, "Topological Torsion atom count must be between 1 and {maximum}")
            }
            Self::InvalidBitWidth => {
                write!(f, "Topological Torsion bit width must be between 1 and {MAX_CONFIGURED_BITS}")
            }
            Self::InvalidCountBounds => f.write_str(
                "Topological Torsion count simulation needs nonempty bounds shorter than nbits and at most 2147483647",
            ),
            Self::Chirality(error) => write!(f, "Topological Torsion chirality: {error}"),
        }
    }
}

impl std::error::Error for RdkitTorsionError {}

/// Return packed bits, folded counts and unhashed counts in one traversal.
/// Supports one to seven atom positions, or one to five with chirality, so
/// the sparse identifiers and the declared result size fit in 64 bits.
pub fn rdkit_torsion_config(
    mol: &Molecule,
    config: &RdkitTorsionConfig,
) -> Result<RdkitTorsionResult, RdkitTorsionError> {
    let code_size = CODE_SIZE + if config.include_chirality { 2 } else { 0 };
    let maximum = (u64::BITS as usize - 1) / code_size as usize;
    if !(1..=maximum).contains(&config.atom_count) {
        return Err(RdkitTorsionError::InvalidAtomCount { maximum });
    }
    if !(1..=MAX_CONFIGURED_BITS).contains(&config.nbits) {
        return Err(RdkitTorsionError::InvalidBitWidth);
    }
    if !valid_count_bounds(config.count_bounds.as_deref(), config.nbits) {
        return Err(RdkitTorsionError::InvalidCountBounds);
    }

    chematic_smiles::with_rdkit_model_molecule(mol, |model| {
        let chirality = legacy_atom_stereo(mol, model, config.include_chirality)
            .map_err(RdkitTorsionError::Chirality)?;
        Ok(chematic_perception::with_rdkit_parity_view(model, |view| {
            make_fingerprint(view.unwrap_or(model), config, &chirality, code_size)
        }))
    })
}

fn make_fingerprint(
    mol: &Molecule,
    config: &RdkitTorsionConfig,
    chirality: &[Option<u8>],
    code_size: u32,
) -> RdkitTorsionResult {
    let endpoint_codes: Vec<_> = mol
        .atoms()
        .map(|(atom, _)| {
            torsion_endpoint_code(
                atom_code(mol, atom, 0) | chirality_bits(chirality, atom.0 as usize),
            )
        })
        .collect();
    let mut counts = ConfiguredCounts::new(config.nbits, config.count_bounds.as_deref());
    let mut sparse = FxHashMap::<u64, u32>::default();
    visit_paths(mol, config, |path| {
        let mut storage = [0u32; MAX_ATOMS];
        let codes = &mut storage[..path.len()];
        for (index, code) in codes.iter_mut().enumerate() {
            let interior = index > 0 && index + 1 < path.len();
            *code = endpoint_codes[path[index].0 as usize] - u32::from(interior);
        }
        orient_torsion_codes(codes);
        counts.add(hash_vec(codes));
        let count = sparse
            .entry(pack_torsion_codes(codes, code_size))
            .or_default();
        *count = count.saturating_add(1);
    });
    let (fingerprint, folded_counts) = counts.finish();
    let mut sparse_counts: Vec<_> = sparse.into_iter().collect();
    sparse_counts.sort_unstable_by_key(|&(code, _)| code);
    RdkitTorsionResult {
        fingerprint,
        folded_counts,
        sparse_counts,
    }
}

fn eligible(mol: &Molecule, atom: AtomIdx) -> bool {
    let atom = mol.atom(atom);
    atom.wildcard || atom.element != chematic_core::Element::H
}

fn visit_paths(mol: &Molecule, config: &RdkitTorsionConfig, mut visit: impl FnMut(&[AtomIdx])) {
    let mut storage = [AtomIdx(0); MAX_ATOMS];
    let path = &mut storage[..config.atom_count];
    let mut distances = config
        .only_shortest_paths
        .then(|| ShortestDistances::new(mol.atom_count()));
    for (start, _) in mol.atoms() {
        path[0] = start;
        // RDKit seeds every atom for a length-one path, including retained H.
        if path.len() == 1 {
            visit(path);
        } else if eligible(mol, start) {
            let distances = distances
                .as_mut()
                .map(|distances| distances.fill_from(mol, start, path.len() - 1));
            extend_path(mol, path, 1, distances, &mut visit);
        }
    }
}

fn extend_path(
    mol: &Molecule,
    path: &mut [AtomIdx],
    filled: usize,
    distances: Option<&[u8]>,
    visit: &mut impl FnMut(&[AtomIdx]),
) {
    let start = path[0];
    let end = path[filled - 1];
    let last = filled + 1 == path.len();
    for (neighbor, _) in mol.neighbors(end) {
        if !eligible(mol, neighbor)
            || distances.is_some_and(|d| usize::from(d[neighbor.0 as usize]) < filled)
        {
            continue;
        }
        if path[..filled].contains(&neighbor) {
            // RDKit keeps the first walk of a ring's bond set. Its adjacency
            // scan starts at the smallest atom index, then the smaller neighbor.
            if last && neighbor == start && closes_first_ring_walk(&path[..filled]) {
                path[filled] = neighbor;
                visit(path);
            }
            continue;
        }
        path[filled] = neighbor;
        if last {
            if start < neighbor {
                visit(path);
            }
        } else {
            extend_path(mol, path, filled + 1, distances, visit);
        }
    }
}

fn closes_first_ring_walk(path: &[AtomIdx]) -> bool {
    path.len() >= 3 && path[1] < path[path.len() - 1] && path[1..].iter().all(|&a| a > path[0])
}

struct ShortestDistances {
    distance: Vec<u8>,
    queue: Vec<AtomIdx>,
}

impl ShortestDistances {
    fn new(atom_count: usize) -> Self {
        Self {
            distance: vec![u8::MAX; atom_count],
            queue: Vec::new(),
        }
    }

    fn fill_from(&mut self, mol: &Molecule, start: AtomIdx, max: usize) -> &[u8] {
        for atom in self.queue.drain(..) {
            self.distance[atom.0 as usize] = u8::MAX;
        }
        self.distance[start.0 as usize] = 0;
        self.queue.push(start);
        let mut head = 0;
        while let Some(&atom) = self.queue.get(head) {
            head += 1;
            let next = self.distance[atom.0 as usize] + 1;
            if usize::from(next) > max {
                continue;
            }
            for (neighbor, _) in mol.neighbors(atom) {
                let slot = &mut self.distance[neighbor.0 as usize];
                if *slot == u8::MAX {
                    *slot = next;
                    self.queue.push(neighbor);
                }
            }
        }
        &self.distance
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn closed_paths_match_rdkit_and_are_excluded_by_shortest_path_filter() {
        let mol = parse("CC1CCC1").unwrap();
        let mut config = RdkitTorsionConfig {
            atom_count: 5,
            nbits: 65,
            count_bounds: Some(vec![2, 1, 4]),
            ..Default::default()
        };
        // RDKit 2026.09.2 retains two open paths and one ring closure.
        let result = rdkit_torsion_config(&mol, &config).unwrap();
        assert_eq!(
            result.fingerprint.to_le_bytes(),
            [0, 0, 1, 0, 6, 0, 0, 0, 0]
        );
        assert_eq!(result.folded_counts, [(16, 2), (27, 1)]);
        assert_eq!(
            result.sparse_counts,
            [(2_272_046_105_120, 2), (2_340_765_581_346, 1)]
        );
        config.only_shortest_paths = true;
        let result = rdkit_torsion_config(&mol, &config).unwrap();
        assert_eq!(result.fingerprint.to_le_bytes(), [0; 9]);
        assert!(result.sparse_counts.is_empty());
    }

    #[test]
    fn encodes_chiral_paths_and_wide_sparse_ids_like_rdkit() {
        let mut config = RdkitTorsionConfig {
            atom_count: 5,
            nbits: 31,
            include_chirality: true,
            count_bounds: None,
            ..Default::default()
        };
        let mol = parse("C[C@H](F)CCCC").unwrap();
        let result = rdkit_torsion_config(&mol, &config).unwrap();
        assert_eq!(result.fingerprint.to_le_bytes(), [0, 0, 0, 72]);
        assert_eq!(
            result.sparse_counts,
            [
                (580_817_151_662_112, 1),
                (633_593_709_789_216, 1),
                (2_252_100_595_679_265, 1)
            ]
        );
        config.atom_count = 7;
        assert_eq!(
            rdkit_torsion_config(&mol, &config),
            Err(RdkitTorsionError::InvalidAtomCount { maximum: 5 })
        );

        config.include_chirality = false;
        config.only_shortest_paths = true;
        let result = rdkit_torsion_config(&parse("CCCCCCC").unwrap(), &config).unwrap();
        assert_eq!(result.fingerprint.to_le_bytes(), [0, 2, 0, 0]);
        assert_eq!(result.sparse_counts, [(577_588_855_536_893_984, 1)]);
    }
}
