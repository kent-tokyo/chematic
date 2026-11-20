//! Configured RDK generator on an indexed, RDKit-sanitized graph.

use super::*;
use crate::bitvec::BitVecN;
use crate::rdkit_count_simulation::{ConfiguredCounts, MAX_CONFIGURED_BITS, valid_count_bounds};
use rustc_hash::FxHashMap;

mod linear;

const MAX_CONFIGURED_PATH: usize = 12;

#[derive(Clone, Debug)]
pub struct RdkitRdkConfig {
    pub min_path: usize,
    pub max_path: usize,
    pub nbits: usize,
    pub bits_per_feature: usize,
    pub use_hs: bool,
    pub branched_paths: bool,
    pub use_bond_order: bool,
    /// Empty selects all roots. Repeated roots contribute repeated counts.
    pub from_atoms: Vec<AtomIdx>,
    /// Ignore every path containing any of these input atom indices.
    pub ignore_atoms: Vec<AtomIdx>,
    pub atom_invariants: Option<Vec<u32>>,
    /// None disables simulation. These thresholds affect only the bit vector.
    pub count_bounds: Option<Vec<u32>>,
    pub max_visits: usize,
    pub max_paths: usize,
}

impl Default for RdkitRdkConfig {
    fn default() -> Self {
        Self {
            min_path: MIN_PATH,
            max_path: MAX_PATH,
            nbits: FP_SIZE as usize,
            bits_per_feature: NUM_BITS_PER_FEATURE as usize,
            use_hs: true,
            branched_paths: true,
            use_bond_order: true,
            from_atoms: Vec::new(),
            ignore_atoms: Vec::new(),
            atom_invariants: None,
            count_bounds: None,
            max_visits: 1_000_000,
            max_paths: 100_000,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RdkitRdkConfigError {
    PathRange,
    FingerprintSize,
    CountBounds,
    AtomIndices,
    AtomInvariants,
    Graph,
    Budget,
    VisitLimit,
    PathLimit,
    IndexedView,
    Sanitization(String),
}

impl core::fmt::Display for RdkitRdkConfigError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::PathRange => f.write_str("RDK path range must satisfy 1 <= min_path <= max_path <= 12"),
            Self::FingerprintSize => f.write_str("RDK nbits must be 1..=1048576 and bits_per_feature must be 1..=32"),
            Self::CountBounds => f.write_str("RDK count bounds require a nonempty list shorter than nbits and signed 32-bit values"),
            Self::AtomIndices => f.write_str("RDK root or ignored atom list exceeds 4096 indices or contains an out-of-range index"),
            Self::AtomInvariants => f.write_str("RDK atom invariants must cover every atom"),
            Self::Graph => f.write_str("RDK requires a simple graph with at most 1024 atoms and 4096 bonds"),
            Self::Budget => f.write_str("RDK visit and path limits must be 1..=10000000 and 1..=1000000 respectively"),
            Self::VisitLimit => f.write_str("RDK search visit limit reached"),
            Self::PathLimit => f.write_str("RDK selected path limit reached"),
            Self::IndexedView => f.write_str("RDK sanitization changed input atom or bond correspondence"),
            Self::Sanitization(message) => write!(f, "RDK sanitization: {message}"),
        }
    }
}

impl std::error::Error for RdkitRdkConfigError {}

#[derive(Clone, Debug)]
pub struct RdkitRdkConfigResult {
    pub fingerprint: BitVecN,
    pub folded_counts: Vec<(u32, u32)>,
    pub sparse_counts: Vec<(u32, u32)>,
}

/// The configured RDK fingerprint uses the input graph's atom indices and
/// preserves explicit hydrogen atoms. Set `use_hs` to false to exclude their
/// paths. The legacy fixed-size API uses SMILES-model hydrogen suppression.
/// Counts include the additional identifiers from `bits_per_feature` and
/// are independent of count simulation. Resource errors return no result.
pub fn rdkit_rdk_config(
    mol: &Molecule,
    config: &RdkitRdkConfig,
) -> Result<RdkitRdkConfigResult, RdkitRdkConfigError> {
    validate(mol, config)?;
    chematic_perception::with_rdkit_parity_view(mol, |view| {
        let view = view.map_err(|e| RdkitRdkConfigError::Sanitization(e.to_string()))?;
        if !same_indexed_graph(mol, view) {
            return Err(RdkitRdkConfigError::IndexedView);
        }
        compute(&crate::rdkit_fp_view(view), config)
    })
}

fn validate(mol: &Molecule, c: &RdkitRdkConfig) -> Result<(), RdkitRdkConfigError> {
    use RdkitRdkConfigError::*;
    if c.min_path == 0 || c.min_path > c.max_path || c.max_path > MAX_CONFIGURED_PATH {
        return Err(PathRange);
    }
    if c.nbits == 0 || c.nbits > MAX_CONFIGURED_BITS || !(1..=32).contains(&c.bits_per_feature) {
        return Err(FingerprintSize);
    }
    if !valid_count_bounds(c.count_bounds.as_deref(), c.nbits) {
        return Err(CountBounds);
    }
    if mol.atom_count() > 1024
        || mol.bond_count() > 4096
        || mol.bonds().any(|(_, b)| b.atom1 == b.atom2)
    {
        return Err(Graph);
    }
    for indices in [&c.from_atoms, &c.ignore_atoms] {
        if indices.len() > 4096
            || indices
                .iter()
                .any(|atom| atom.0 as usize >= mol.atom_count())
        {
            return Err(AtomIndices);
        }
    }
    if c.atom_invariants
        .as_ref()
        .is_some_and(|values| values.len() < mol.atom_count())
    {
        return Err(AtomInvariants);
    }
    if !(1..=10_000_000).contains(&c.max_visits) || !(1..=1_000_000).contains(&c.max_paths) {
        return Err(Budget);
    }
    Ok(())
}

fn same_indexed_graph(source: &Molecule, view: &Molecule) -> bool {
    source.atom_count() == view.atom_count()
        && source.bond_count() == view.bond_count()
        && source.atoms().zip(view.atoms()).all(|((_, a), (_, b))| {
            a.element == b.element && a.isotope == b.isotope && a.wildcard == b.wildcard
        })
        && source.bonds().zip(view.bonds()).all(|((_, a), (_, b))| {
            (a.atom1, a.atom2) == (b.atom1, b.atom2) || (a.atom1, a.atom2) == (b.atom2, b.atom1)
        })
}

fn compute(
    mol: &Molecule,
    c: &RdkitRdkConfig,
) -> Result<RdkitRdkConfigResult, RdkitRdkConfigError> {
    let default_invariants;
    let invariants = match &c.atom_invariants {
        Some(invariants) => invariants.as_slice(),
        None => {
            default_invariants = mol
                .atoms()
                .map(|(i, a)| (crate::rdkit_atomic_num(mol, i) % 128) << 1 | u32::from(a.aromatic))
                .collect::<Vec<_>>();
            &default_invariants
        }
    };
    let bond_types: Vec<_> = mol
        .bonds()
        .map(|(_, b)| {
            if c.use_bond_order {
                rdkit_bond_type_code(b.order)
            } else {
                1
            }
        })
        .collect();
    let ends: Vec<_> = mol.bonds().map(|(_, b)| (b.atom1.0, b.atom2.0)).collect();
    let mut counts = ConfiguredCounts::new(c.nbits, c.count_bounds.as_deref());
    let mut sparse = FxHashMap::default();
    let mut budget = Budget {
        visits_left: c.max_visits,
        paths_left: c.max_paths,
    };
    let excluded = excluded_bonds(mol, c);
    let mut add_path = |path: &[u32], degrees: &[u32], atom_count| {
        let seed = path_seed_small::<MAX_CONFIGURED_PATH>(
            path,
            &ends,
            invariants,
            &bond_types,
            degrees,
            atom_count,
        );
        let mut add = |id| {
            counts.add(id);
            *sparse.entry(id).or_insert(0) += 1;
        };
        add(seed);
        if c.bits_per_feature > 1 {
            let mut mt = RdkitWeakMt19937::seeded(seed);
            for _ in 1..c.bits_per_feature {
                add(mt.next_u32() / 2);
            }
        }
    };

    if c.branched_paths {
        let mut walker = SubgraphWalker::with_mask(mol, excluded.clone(), c.max_path);
        for root in roots(c) {
            walker.forbidden.clone_from(&excluded);
            walker.try_for_each(root, 1, c.max_path, |path, degrees, atom_count| {
                budget.step()?;
                if path.len() >= c.min_path {
                    budget.select()?;
                    add_path(path, degrees, atom_count);
                }
                Ok(())
            })?;
        }
    } else {
        linear::enumerate(mol, c, &excluded, &mut budget, add_path)?;
    }
    let (fingerprint, folded_counts) = counts.finish();
    let mut sparse_counts: Vec<_> = sparse.into_iter().collect();
    sparse_counts.sort_unstable();
    Ok(RdkitRdkConfigResult {
        fingerprint,
        folded_counts,
        sparse_counts,
    })
}

fn roots(c: &RdkitRdkConfig) -> impl Iterator<Item = Option<AtomIdx>> + '_ {
    c.from_atoms
        .is_empty()
        .then_some(None)
        .into_iter()
        .chain(c.from_atoms.iter().copied().map(Some))
}

fn excluded_bonds(mol: &Molecule, c: &RdkitRdkConfig) -> Vec<bool> {
    let mut excluded_atoms = vec![false; mol.atom_count()];
    if !c.use_hs {
        for (i, _) in mol.atoms() {
            excluded_atoms[i.0 as usize] = crate::rdkit_atomic_num(mol, i) == 1;
        }
    }
    for a in &c.ignore_atoms {
        excluded_atoms[a.0 as usize] = true;
    }
    mol.bonds()
        .map(|(_, b)| {
            excluded_atoms[b.atom1.0 as usize]
                || excluded_atoms[b.atom2.0 as usize]
                || matches!(
                    b.order,
                    BondOrder::QueryAny
                        | BondOrder::QuerySingleOrDouble
                        | BondOrder::QuerySingleOrAromatic
                        | BondOrder::QueryDoubleOrAromatic
                )
        })
        .collect()
}

struct Budget {
    visits_left: usize,
    paths_left: usize,
}

impl Budget {
    fn step(&mut self) -> Result<(), RdkitRdkConfigError> {
        self.visits_left = self
            .visits_left
            .checked_sub(1)
            .ok_or(RdkitRdkConfigError::VisitLimit)?;
        Ok(())
    }

    fn select(&mut self) -> Result<(), RdkitRdkConfigError> {
        self.paths_left = self
            .paths_left
            .checked_sub(1)
            .ok_or(RdkitRdkConfigError::PathLimit)?;
        Ok(())
    }
}
