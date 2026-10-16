//! Checked custom SMARTS feature masks and native circular fingerprints.

use crate::bitvec::BitVecN;
use crate::ecfp::{bond_type_int, fnv1a};
use chematic_core::{AtomIdx, Molecule};
use chematic_smarts::{MatchConfig, QueryMolecule};

/// Invalid input, unsupported size or exhausted matching budget.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeaturePatternError(pub String);

impl std::fmt::Display for FeaturePatternError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for FeaturePatternError {}

/// Owned, ordered custom feature definitions. Pattern i sets bit i on every
/// target atom participating in any of its matches. No default patterns are
/// substituted, including for an empty set. Maximum 32 patterns, 4096 bytes
/// and 32 top-level query atoms per pattern.
pub struct FeaturePatternSet {
    sources: Vec<String>,
    queries: Vec<QueryMolecule>,
}

impl FeaturePatternSet {
    /// Compile an iterable, consuming at most 33 entries to enforce the cap.
    pub fn new<S: AsRef<str>>(
        patterns: impl IntoIterator<Item = S>,
    ) -> Result<Self, FeaturePatternError> {
        let mut sources = Vec::new();
        let mut queries = Vec::new();
        for pattern in patterns {
            if sources.len() == 32 {
                return Err(FeaturePatternError(
                    "at most 32 feature patterns are supported".into(),
                ));
            }
            let text = pattern.as_ref();
            if text.len() > 4096 {
                return Err(FeaturePatternError(
                    "feature pattern exceeds 4096 bytes".into(),
                ));
            }
            let query = chematic_smarts::parse_smarts(text)
                .map_err(|error| FeaturePatternError(error.to_string()))?;
            if query.atoms.is_empty() || query.atoms.len() > 32 {
                return Err(FeaturePatternError(
                    "feature pattern must have 1 to 32 query atoms".into(),
                ));
            }
            sources.push(text.to_owned());
            queries.push(query);
        }
        Ok(Self { sources, queries })
    }

    /// Original SMARTS in bit order.
    pub fn patterns(&self) -> &[String] {
        &self.sources
    }

    /// Full-width atom masks, preserving source atom indices. Matching ignores
    /// chirality, enforces isotopes and uses the RDKit-parity aromatic view with
    /// native SSSR SMARTS ring semantics. Failure never returns partial masks.
    /// Maximum 256 atoms / 512 bonds, and 1..=1,000,000 visits per pattern.
    pub fn invariants(
        &self,
        mol: &Molecule,
        max_visits: u64,
    ) -> Result<Vec<u32>, FeaturePatternError> {
        if mol.atom_count() > 256 || mol.bond_count() > 512 {
            return Err(FeaturePatternError(
                "feature target exceeds 256 atoms or 512 bonds".into(),
            ));
        }
        if !(1..=1_000_000).contains(&max_visits) {
            return Err(FeaturePatternError(
                "max_visits must be 1 to 1000000 per pattern".into(),
            ));
        }
        let mut masks = vec![0u32; mol.atom_count()];
        if self.queries.is_empty() || masks.is_empty() {
            return Ok(masks);
        }
        let perceived = chematic_perception::apply_aromaticity_rdkit_parity_shared(mol);
        let target = perceived
            .as_ref()
            .as_ref()
            .map_err(|error| FeaturePatternError(format!("aromatic perception failed: {error}")))?;
        let rings = chematic_perception::find_sssr(target);
        let config = MatchConfig {
            max_visit_budget: Some(max_visits),
            ..Default::default()
        };
        for (index, query) in self.queries.iter().enumerate() {
            let (matches, exhausted) = chematic_smarts::find_matches_with_rings_and_config_checked(
                query, target, &rings, &config,
            );
            if exhausted {
                return Err(FeaturePatternError(format!(
                    "feature pattern {index} exhausted visit budget"
                )));
            }
            let bit = 1u32 << index;
            for mapping in matches {
                for atom in mapping.values() {
                    masks[atom.0 as usize] |= bit;
                }
            }
        }
        Ok(masks)
    }

    /// Native FNV circular feature bits, not RDKit Morgan bit identifiers.
    /// Seeds hash all four little-endian bytes of each u32 mask; every atom is
    /// emitted at every radius, including zero-mask atoms. Radius 0..=3;
    /// widths 128, 256, 512, 1024 or 2048. Existing FCFP APIs are unchanged.
    pub fn fingerprint(
        &self,
        mol: &Molecule,
        radius: u32,
        nbits: usize,
        max_visits: u64,
    ) -> Result<BitVecN, FeaturePatternError> {
        if radius > 3 || ![128, 256, 512, 1024, 2048].contains(&nbits) {
            return Err(FeaturePatternError(
                "radius must be 0 to 3 and nbits 128/256/512/1024/2048".into(),
            ));
        }
        let masks = self.invariants(mol, max_visits)?;
        let perceived = chematic_perception::apply_aromaticity_rdkit_parity_shared(mol);
        let mol = perceived
            .as_ref()
            .as_ref()
            .map_err(|error| FeaturePatternError(format!("aromatic perception failed: {error}")))?;
        let mut ids: Vec<_> = masks
            .iter()
            .map(|mask| fnv1a(&mask.to_le_bytes()))
            .collect();
        let mut fp = BitVecN::new(nbits);
        for id in &ids {
            fp.set((*id % nbits as u64) as usize);
        }
        for round in 1..=radius {
            let mut next = Vec::with_capacity(ids.len());
            for (index, id) in ids.iter().enumerate() {
                let mut neighbors: Vec<_> = mol
                    .neighbors(AtomIdx(index as u32))
                    .map(|(atom, bond)| (bond_type_int(mol.bond(bond).order), ids[atom.0 as usize]))
                    .collect();
                neighbors.sort_unstable();
                let mut bytes = vec![round as u8];
                bytes.extend_from_slice(&id.to_le_bytes());
                for (bond, neighbor) in neighbors {
                    bytes.push(bond);
                    bytes.extend_from_slice(&neighbor.to_le_bytes());
                }
                let new_id = fnv1a(&bytes);
                fp.set((new_id % nbits as u64) as usize);
                next.push(new_id);
            }
            ids = next;
        }
        Ok(fp)
    }
}
