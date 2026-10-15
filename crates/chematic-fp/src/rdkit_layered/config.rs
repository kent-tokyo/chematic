//! Checked configurable Layered fingerprints, pinned to RDKit 2026.09.1.
//! RDKit-derived behavior is covered by THIRD_PARTY_NOTICES.md (BSD-3-Clause).
use super::*;
mod enumeration;
mod validation;

#[derive(Debug, Clone)]
pub struct RdkitLayeredConfig {
    pub layer_flags: u32,
    pub min_path: usize,
    pub max_path: usize,
    pub fp_size: usize,
    pub branched_paths: bool,
    /// Empty means all paths; duplicate roots retain their count multiplicity.
    pub from_atoms: Vec<AtomIdx>,
    /// Packed LSB-first mask, with zero padding above fp_size.
    pub set_only_bits: Option<Vec<u8>>,
    /// Initial u32 counters; extra trailing entries are returned unchanged.
    pub atom_counts: Option<Vec<u32>>,
    pub max_visits: u64,
    pub max_paths: usize,
}
impl Default for RdkitLayeredConfig {
    fn default() -> Self {
        Self {
            layer_flags: u32::MAX,
            min_path: 1,
            max_path: 7,
            fp_size: 2048,
            branched_paths: true,
            from_atoms: Vec::new(),
            set_only_bits: None,
            atom_counts: None,
            max_visits: 1_000_000,
            max_paths: 100_000,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitLayeredError {
    pub code: &'static str,
    pub detail: String,
}
impl core::fmt::Display for RdkitLayeredError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Layered {}: {}", self.code, self.detail)
    }
}
impl std::error::Error for RdkitLayeredError {}
fn error(code: &'static str, detail: impl Into<String>) -> RdkitLayeredError {
    RdkitLayeredError {
        code,
        detail: detail.into(),
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitLayeredResult {
    pub fp_size: usize,
    pub bits: Vec<u8>,
    pub atom_counts: Vec<u32>,
    pub selected_paths: usize,
    pub search_visits: u64,
}

/// No caller counters are mutated, and resource exhaustion returns no partial result.
/// Requires a simple indexed molecule; complex query bonds are rejected.
pub fn rdkit_layered_config(
    mol: &Molecule,
    config: &RdkitLayeredConfig,
) -> Result<RdkitLayeredResult, RdkitLayeredError> {
    validation::validate(mol, config)?;
    chematic_perception::with_rdkit_parity_view(mol, |view| {
        let view = view.map_err(|e| error("aromaticity", e.to_string()))?;
        validation::same_indices(mol, view)?;
        let fp_view = crate::rdkit_fp_view(view);
        compute(&fp_view, config)
    })
}
fn compute(
    mol: &Molecule,
    c: &RdkitLayeredConfig,
) -> Result<RdkitLayeredResult, RdkitLayeredError> {
    let prepared = features::Prepared::new(mol, true);
    let mut result = RdkitLayeredResult {
        fp_size: c.fp_size,
        bits: vec![0; c.fp_size.div_ceil(8)],
        atom_counts: c
            .atom_counts
            .clone()
            .unwrap_or_else(|| vec![0; mol.atom_count()]),
        selected_paths: 0,
        search_visits: 0,
    };
    let mut budget = enumeration::Budget::new(c);
    enumeration::enumerate(mol, c, &mut budget, |path, multiplicity| {
        let hashes = prepared.hashes(path);
        let mut accepted = false;
        for (layer, seed) in hashes.seeds.into_iter().enumerate() {
            if c.layer_flags & (1 << layer) == 0 {
                continue;
            }
            let Some(seed) = seed else {
                continue;
            };
            let bit = seed as usize % c.fp_size;
            if c.set_only_bits
                .as_ref()
                .is_some_and(|mask| mask[bit / 8] & (1 << (bit % 8)) == 0)
            {
                continue;
            }
            result.bits[bit / 8] |= 1 << (bit % 8);
            accepted = true;
        }
        if accepted {
            for atom in hashes.atoms {
                result.atom_counts[atom] =
                    result.atom_counts[atom].wrapping_add(multiplicity as u32);
            }
        }
    })?;
    result.selected_paths = budget.paths;
    result.search_visits = budget.visits;
    Ok(result)
}
