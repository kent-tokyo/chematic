//! R-group decomposition: split molecules into a common scaffold core and variable
//! R-group substituents.
//!
//! Given a SMARTS scaffold pattern and a list of molecules, this module identifies
//! for each molecule which atoms belong to the scaffold core and extracts the
//! substituents at each attachment point as separate R-group SMILES strings.
//!
//! # Usage
//!
//! ```
//! # use chematic_smiles::parse;
//! # use chematic_chem::rgroup::{rgroup_decompose, RGroupError};
//! let mols = vec![
//!     parse("CCc1ccccc1").unwrap(),
//!     parse("CCCc1ccccc1").unwrap(),
//!     parse("Nc1ccccc1").unwrap(),
//! ];
//! let refs: Vec<&_> = mols.iter().collect();
//! let results = rgroup_decompose("c1ccccc1", &refs).unwrap();
//! for r in results.iter().flatten() {
//!     println!("R1: {:?}", r.r_groups.get(&1));
//! }
//! ```

use std::collections::{HashMap, HashSet, VecDeque};

use chematic_core::{Atom, AtomIdx, BondIdx, BondOrder, Molecule, MoleculeBuilder};
use chematic_smarts::{MatchConfig, find_matches_with_config, parse_smarts};
use chematic_smiles::canonical_smiles;

/// Error type for R-group decomposition.
#[derive(Debug, Clone, PartialEq)]
pub enum RGroupError {
    /// The scaffold SMARTS string could not be parsed.
    InvalidSmarts(String),
}

impl std::fmt::Display for RGroupError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidSmarts(s) => write!(f, "invalid SMARTS: {s}"),
        }
    }
}

impl std::error::Error for RGroupError {}

/// Result of decomposing one molecule against a scaffold.
#[derive(Debug, Clone)]
pub struct RGroupResult {
    /// Index of the molecule in the input slice.
    pub mol_idx: usize,
    /// Canonical SMILES of the scaffold core.
    /// Attachment positions (where R-groups were cut) carry `[*]` wildcard atoms.
    pub core_smiles: String,
    /// R-groups keyed by 1-based attachment-point index.
    ///
    /// The attachment points are ordered by the scaffold-atom index (ascending)
    /// at which the cut was made.  Each value is a SMILES string where `[*]`
    /// represents the bond back to the scaffold core.
    pub r_groups: HashMap<u8, String>,
}

/// Decompose `mols` against `scaffold_smarts`.
///
/// For each molecule the algorithm:
/// 1. Finds the first SMARTS match of the scaffold.
/// 2. Treats the matched atoms as the **core**.
/// 3. Identifies every bond crossing from a core atom to a non-core atom —
///    each such bond defines one R-group attachment point.
/// 4. Extracts the non-core subgraph reachable from that bond's endpoint (BFS,
///    blocked at core atoms) and writes it as a SMILES fragment with `[*]` at
///    the attachment.
///
/// Attachment points are numbered 1, 2, … in ascending order of the core atom's
/// molecule-atom index.
///
/// Returns `Vec<Option<RGroupResult>>` parallel to `mols`:
/// - `Some(result)` — the scaffold matched and R-groups were extracted.
/// - `None` — the scaffold did not match this molecule.
///
/// # Errors
///
/// Returns `RGroupError::InvalidSmarts` if `scaffold_smarts` fails to parse.
pub fn rgroup_decompose(
    scaffold_smarts: &str,
    mols: &[&Molecule],
) -> Result<Vec<Option<RGroupResult>>, RGroupError> {
    let query =
        parse_smarts(scaffold_smarts).map_err(|e| RGroupError::InvalidSmarts(e.to_string()))?;

    let cfg = MatchConfig {
        max_matches: Some(1),
        ..Default::default()
    };

    let mut results = Vec::with_capacity(mols.len());

    for (mol_idx, &mol) in mols.iter().enumerate() {
        let matches = find_matches_with_config(&query, mol, &cfg);

        if matches.is_empty() {
            results.push(None);
            continue;
        }

        // Core atoms: molecule atom indices that were matched by the scaffold.
        // Note: find_matches returns HashMap<query_atom_idx, mol_atom_idx>.
        let core: HashSet<AtomIdx> = matches[0].values().copied().collect();

        // Find attachment points: pairs (core_atom, rgroup_root) for bonds
        // crossing the core boundary.  Collect them sorted by core atom index
        // for stable R-group numbering.
        let mut attachment_pairs: Vec<(AtomIdx, AtomIdx)> = Vec::new();
        {
            let mut seen_bonds: HashSet<(u32, u32)> = HashSet::new();
            let mut core_sorted: Vec<AtomIdx> = core.iter().copied().collect();
            core_sorted.sort_by_key(|a| a.0);

            for &ca in &core_sorted {
                for (nb, _) in mol.neighbors(ca) {
                    if !core.contains(&nb) {
                        let key = (ca.0.min(nb.0), ca.0.max(nb.0));
                        if seen_bonds.insert(key) {
                            attachment_pairs.push((ca, nb));
                        }
                    }
                }
            }
        }

        // For each attachment point, BFS the non-core subgraph and build a fragment.
        let mut r_groups: HashMap<u8, String> = HashMap::new();
        for (rg_idx, &(_core_atom, rg_root)) in attachment_pairs.iter().enumerate() {
            let rg_num = (rg_idx + 1) as u8;
            let fragment = extract_rgroup(mol, rg_root, &core);
            r_groups.insert(rg_num, canonical_smiles(&fragment));
        }

        // Build core SMILES: copy core atoms, replacing boundary bonds with [*].
        let core_mol = build_core_mol(mol, &core, &attachment_pairs);
        let core_smiles = canonical_smiles(&core_mol);

        results.push(Some(RGroupResult {
            mol_idx,
            core_smiles,
            r_groups,
        }));
    }

    Ok(results)
}

/// Extract the R-group subgraph starting at `root`, blocked by the `core` atom set.
///
/// Returns a `Molecule` containing all atoms reachable from `root` through bonds
/// that do not enter the core, plus a `[*]` wildcard stub at the attachment point.
fn extract_rgroup(mol: &Molecule, root: AtomIdx, core: &HashSet<AtomIdx>) -> Molecule {
    // BFS to collect the R-group atom set.
    let mut rg_atoms: HashSet<AtomIdx> = HashSet::new();
    let mut queue: VecDeque<AtomIdx> = VecDeque::new();
    rg_atoms.insert(root);
    queue.push_back(root);

    while let Some(cur) = queue.pop_front() {
        for (nb, _) in mol.neighbors(cur) {
            if !core.contains(&nb) && !rg_atoms.contains(&nb) {
                rg_atoms.insert(nb);
                queue.push_back(nb);
            }
        }
    }

    // Build subgraph molecule, remapping indices.
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Add root first so its index is predictable (index 0).
    let new_root = builder.add_atom(mol.atom(root).clone());
    remap.insert(root, new_root);

    let mut sorted: Vec<AtomIdx> = rg_atoms.iter().copied().filter(|&a| a != root).collect();
    sorted.sort_by_key(|a| a.0);
    for &old_idx in &sorted {
        let new_idx = builder.add_atom(mol.atom(old_idx).clone());
        remap.insert(old_idx, new_idx);
    }

    // Add bonds within the R-group.
    for bidx in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(bidx as u32));
        if let (Some(&new_a), Some(&new_b)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(new_a, new_b, bond.order);
        }
    }

    // Add a [*] wildcard attached to root representing the attachment to the core.
    let wc = builder.add_atom(Atom::wildcard());
    let _ = builder.add_bond(new_root, wc, BondOrder::Single);

    builder.build()
}

/// Build the core molecule with `[*]` stubs at each attachment point.
fn build_core_mol(
    mol: &Molecule,
    core: &HashSet<AtomIdx>,
    attachment_pairs: &[(AtomIdx, AtomIdx)],
) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut remap: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Add all core atoms (stable order).
    let mut core_sorted: Vec<AtomIdx> = core.iter().copied().collect();
    core_sorted.sort_by_key(|a| a.0);
    for &old_idx in &core_sorted {
        let new_idx = builder.add_atom(mol.atom(old_idx).clone());
        remap.insert(old_idx, new_idx);
    }

    // Add bonds between core atoms.
    for bidx in 0..mol.bond_count() {
        let bond = mol.bond(BondIdx(bidx as u32));
        if let (Some(&na), Some(&nb)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(na, nb, bond.order);
        }
    }

    // Add [*] wildcards at each attachment point.
    for &(core_atom, _rg_root) in attachment_pairs {
        if let Some(&new_ca) = remap.get(&core_atom) {
            let wc = builder.add_atom(Atom::wildcard());
            let _ = builder.add_bond(new_ca, wc, BondOrder::Single);
        }
    }

    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    fn decompose<'a>(
        scaffold: &str,
        smiles: impl IntoIterator<Item = &'a str>,
    ) -> Vec<Option<RGroupResult>> {
        let mols: Vec<Molecule> = smiles.into_iter().map(mol).collect();
        let refs: Vec<&Molecule> = mols.iter().collect();
        rgroup_decompose(scaffold, &refs).expect("decompose failed")
    }

    #[test]
    fn test_rgroup_monosubstituted_benzene() {
        // Three alkylbenzenes — scaffold = benzene ring (6 C, no attachment markers).
        let results = decompose("c1ccccc1", ["Cc1ccccc1", "CCc1ccccc1", "CCCc1ccccc1"]);
        assert_eq!(results.len(), 3);
        for (i, r) in results.iter().enumerate() {
            let r = r.as_ref().unwrap_or_else(|| panic!("mol {i} should match"));
            assert_eq!(r.mol_idx, i);
            // Must have exactly one R-group (one substituted position).
            assert_eq!(r.r_groups.len(), 1, "mol {i}: expected 1 R-group");
            let r1 = r.r_groups.get(&1).expect("R1 missing");
            // R1 must contain a wildcard marker and the appropriate alkyl chain.
            assert!(
                r1.contains('*'),
                "R1 must contain [*] attachment, got {r1:?}"
            );
        }
    }

    #[test]
    fn test_rgroup_mol_idx_preserved() {
        let results = decompose("c1ccccc1", ["Cc1ccccc1", "CCCC", "CCCc1ccccc1"]);
        // CCCC does not match benzene scaffold → None
        assert!(results[0].is_some());
        assert!(results[1].is_none());
        assert!(results[2].is_some());
        assert_eq!(results[0].as_ref().unwrap().mol_idx, 0);
        assert_eq!(results[2].as_ref().unwrap().mol_idx, 2);
    }

    #[test]
    fn test_rgroup_no_substituents() {
        // Benzene itself: scaffold matches, zero R-groups.
        let results = decompose("c1ccccc1", ["c1ccccc1"]);
        let r = results[0]
            .as_ref()
            .expect("benzene should match its own scaffold");
        assert_eq!(r.r_groups.len(), 0, "benzene has no R-groups vs itself");
    }

    #[test]
    fn test_rgroup_two_substituents() {
        // para-xylene: two methyl groups → two R-groups.
        let results = decompose("c1ccccc1", ["Cc1ccc(C)cc1"]);
        let r = results[0].as_ref().expect("should match");
        assert_eq!(r.r_groups.len(), 2, "para-xylene has two R-groups");
    }

    #[test]
    fn test_rgroup_invalid_smarts() {
        let mols = [mol("c1ccccc1")];
        let refs: Vec<&Molecule> = mols.iter().collect();
        let err = rgroup_decompose("((invalid", &refs);
        assert!(matches!(err, Err(RGroupError::InvalidSmarts(_))));
    }

    #[test]
    fn test_rgroup_different_heteroatom_substituents() {
        // Aniline and phenol vs benzene scaffold.
        let results = decompose("c1ccccc1", ["Nc1ccccc1", "Oc1ccccc1"]);
        let aniline = results[0].as_ref().expect("aniline matches");
        let phenol = results[1].as_ref().expect("phenol matches");
        assert_eq!(aniline.r_groups.len(), 1);
        assert_eq!(phenol.r_groups.len(), 1);
        // Both R1 values must start with N or O (amine / hydroxyl).
        let r1_a = aniline.r_groups.get(&1).unwrap();
        let r1_p = phenol.r_groups.get(&1).unwrap();
        assert!(
            r1_a.starts_with('N') || r1_a.contains('N'),
            "aniline R1 should contain N, got {r1_a}"
        );
        assert!(
            r1_p.starts_with('O') || r1_p.contains('O'),
            "phenol R1 should contain O, got {r1_p}"
        );
    }
}
