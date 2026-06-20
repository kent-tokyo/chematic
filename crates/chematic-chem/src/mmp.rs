//! Matched Molecular Pair (MMP) analysis.
//!
//! A matched molecular pair is two molecules that differ by exactly one
//! structural transformation at a single-bond (BRICS) cut point.
//!
//! # Algorithm
//!
//! For each molecule, every BRICS-breakable bond is severed to give a
//! **(core, substituent)** pair where the core is the larger fragment and the
//! substituent is the smaller fragment.  Both fragments include a `[*]` atom
//! marking the attachment point.
//!
//! Two molecules form a MMP when they share the same core SMILES but have
//! different substituent SMILES at that core position.
//!
//! # Limitations
//!
//! Only single BRICS cuts are considered.  Transformations that require two
//! cuts (e.g. replacing a ring) are not detected.

use std::collections::{HashMap, HashSet};

use chematic_core::{Atom, AtomIdx, BondOrder, Molecule, MoleculeBuilder};
use chematic_smiles::canonical_smiles;

use crate::brics::brics_bonds;

/// A single matched molecular pair.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmpPair {
    /// Canonical SMILES of the first molecule.
    pub mol_a: String,
    /// Canonical SMILES of the second molecule.
    pub mol_b: String,
    /// Fragment SMILES of the common core (contains `[*]`).
    pub core: String,
    /// Substituent fragment of `mol_a` (contains `[*]`).
    pub fragment_a: String,
    /// Substituent fragment of `mol_b` (contains `[*]`).
    pub fragment_b: String,
}

/// Find all matched molecular pairs in `mols`.
///
/// Returns one `MmpPair` per unique (mol_a, mol_b, core, {frag_a, frag_b}) combination.
/// Pairs are deduplicated: (A,B,core,fA,fB) and (B,A,core,fB,fA) are the same pair.
pub fn find_mmp(mols: &[&Molecule]) -> Vec<MmpPair> {
    // 1. For every molecule, collect all (core_smiles, sub_smiles) from all BRICS cuts.
    let mol_smiles: Vec<String> = mols.iter().map(|m| canonical_smiles(m)).collect();

    // index: core_smiles → Vec<(mol_idx, sub_smiles)>
    let mut index: HashMap<String, Vec<(usize, String)>> = HashMap::new();

    for (mol_idx, mol) in mols.iter().enumerate() {
        for (core_smi, sub_smi) in all_cuts(mol) {
            index.entry(core_smi).or_default().push((mol_idx, sub_smi));
        }
    }

    // 2. For each core with ≥ 2 entries, find molecule pairs with different substituents.
    // Dedup key: (min_mol_idx, max_mol_idx, core, sorted sub pair)
    let mut seen: HashSet<(usize, usize, String, String, String)> = HashSet::new();
    let mut pairs: Vec<MmpPair> = Vec::new();

    for (core_smi, entries) in &index {
        if entries.len() < 2 {
            continue;
        }
        for i in 0..entries.len() {
            for j in (i + 1)..entries.len() {
                let (mi, sub_i) = &entries[i];
                let (mj, sub_j) = &entries[j];
                if mi == mj || sub_i == sub_j {
                    continue; // same mol or same substituent → not a MMP
                }
                // Canonical dedup key: mol indices in ascending order.
                let (lo, hi, sub_lo, sub_hi) = if mi < mj {
                    (*mi, *mj, sub_i.clone(), sub_j.clone())
                } else {
                    (*mj, *mi, sub_j.clone(), sub_i.clone())
                };
                let key = (lo, hi, core_smi.clone(), sub_lo.clone(), sub_hi.clone());
                if seen.insert(key) {
                    pairs.push(MmpPair {
                        mol_a: mol_smiles[lo].clone(),
                        mol_b: mol_smiles[hi].clone(),
                        core: core_smi.clone(),
                        fragment_a: sub_lo,
                        fragment_b: sub_hi,
                    });
                }
            }
        }
    }

    pairs.sort_by(|a, b| {
        a.mol_a
            .cmp(&b.mol_a)
            .then(a.mol_b.cmp(&b.mol_b))
            .then(a.core.cmp(&b.core))
    });
    pairs
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// All (core_smiles, sub_smiles) pairs from every BRICS cut of `mol`.
/// Convention: smaller atom-count side → substituent; larger → core.
fn all_cuts(mol: &Molecule) -> Vec<(String, String)> {
    let mut result = Vec::new();
    for (a1, a2) in brics_bonds(mol) {
        let side1 = atoms_on_side(mol, a1, a2);
        let side2 = atoms_on_side(mol, a2, a1);
        let (sub, core, at_sub, at_core) = if side1.len() <= side2.len() {
            (side1, side2, a1, a2)
        } else {
            (side2, side1, a2, a1)
        };
        let core_smi = fragment_smiles(mol, &core, at_core);
        let sub_smi = fragment_smiles(mol, &sub, at_sub);
        result.push((core_smi, sub_smi));
    }
    result
}

/// BFS from `from`, never stepping through `not_via`.
fn atoms_on_side(mol: &Molecule, from: AtomIdx, not_via: AtomIdx) -> HashSet<AtomIdx> {
    let mut visited = HashSet::new();
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(from);
    while let Some(idx) = queue.pop_front() {
        if visited.contains(&idx) {
            continue;
        }
        visited.insert(idx);
        for (nb, _) in mol.neighbors(idx) {
            if nb != not_via && !visited.contains(&nb) {
                queue.push_back(nb);
            }
        }
    }
    visited
}

/// Build a sub-molecule from `side` atoms and add a `[*]` atom bonded to `attach`,
/// then return its canonical SMILES.
fn fragment_smiles(mol: &Molecule, side: &HashSet<AtomIdx>, attach: AtomIdx) -> String {
    let mut builder = MoleculeBuilder::new();
    let mut idx_map: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    // Wildcard attachment marker.
    let mut wc = Atom::new(chematic_core::Element::C);
    wc.wildcard = true;
    let wc_idx = builder.add_atom(wc);

    // Side atoms.
    for &orig in side {
        let atom = mol.atom(orig);
        let mut a = Atom::new(atom.element);
        a.charge = atom.charge;
        a.isotope = atom.isotope;
        a.aromatic = atom.aromatic;
        a.chirality = atom.chirality;
        a.hydrogen_count = atom.hydrogen_count;
        a.atom_map = atom.atom_map;
        let new_idx = builder.add_atom(a);
        idx_map.insert(orig, new_idx);
    }

    // Bond: wildcard → attachment atom.
    let _ = builder.add_bond(wc_idx, *idx_map.get(&attach).unwrap(), BondOrder::Single);

    // Intra-side bonds.
    for (_, bond) in mol.bonds() {
        if side.contains(&bond.atom1)
            && side.contains(&bond.atom2)
            && let (Some(&n1), Some(&n2)) = (idx_map.get(&bond.atom1), idx_map.get(&bond.atom2))
        {
            let _ = builder.add_bond(n1, n2, bond.order);
        }
    }

    canonical_smiles(&builder.build())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    #[test]
    fn mmp_ethylbenzene_propylbenzene() {
        // ethylbenzene ↔ propylbenzene: same benzene core, subs differ by one CH2.
        let a = mol("CCc1ccccc1");
        let b = mol("CCCc1ccccc1");
        let pairs = find_mmp(&[&a, &b]);

        // There must be exactly 1 MMP (the ring-chain cut).
        // Canonical form of monosubstituted benzene core after #14 bond-order fix.
        let matching: Vec<_> = pairs.iter().filter(|p| p.core == "c1c([*])cccc1").collect();
        assert_eq!(
            matching.len(),
            1,
            "expected 1 pair with benzene core, got: {pairs:?}"
        );

        let pair = &matching[0];
        // Oracle values: canonical SMILES updated for bond-order-aware Morgan ranks (#14).
        assert_eq!(
            pair.fragment_a, "C(C)[*]",
            "ethylbenzene substituent should be C(C)[*]: {pair:?}"
        );
        assert_eq!(
            pair.fragment_b, "[*]CCC",
            "propylbenzene substituent should be [*]CCC: {pair:?}"
        );
    }

    #[test]
    fn mmp_no_pairs_for_single_molecule() {
        let a = mol("CCc1ccccc1");
        let pairs = find_mmp(&[&a]);
        assert!(pairs.is_empty(), "single molecule has no MMP pairs");
    }

    #[test]
    fn mmp_no_pairs_when_no_brics_bonds() {
        // Benzene has no BRICS bonds → no cuts → no pairs.
        let a = mol("c1ccccc1");
        let b = mol("c1ccncc1");
        let pairs = find_mmp(&[&a, &b]);
        assert!(
            pairs.is_empty(),
            "benzene/pyridine have no BRICS bonds, expect 0 pairs: {pairs:?}"
        );
    }

    #[test]
    fn mmp_dedup_direction() {
        // (A,B) and (B,A) should NOT produce two entries.
        let a = mol("CCc1ccccc1");
        let b = mol("CCCc1ccccc1");
        let pairs = find_mmp(&[&a, &b]);
        let n = pairs.len();
        // flip order
        let pairs2 = find_mmp(&[&b, &a]);
        let n2 = pairs2.len();
        assert_eq!(n, n2, "pair count must be order-independent: {n} vs {n2}");
    }

    #[test]
    fn mmp_three_molecules_correct_count() {
        // Three molecules sharing the same benzene core: 3 MMP pairs (C(C), C(CC), C(CCC) subs).
        let a = mol("CCc1ccccc1");
        let b = mol("CCCc1ccccc1");
        let c = mol("CCCCc1ccccc1"); // butylbenzene
        let pairs = find_mmp(&[&a, &b, &c]);
        // Expected: (a,b), (a,c), (b,c) — 3 pairs at minimum.
        // Canonical form updated for bond-order-aware Morgan ranks (#14).
        let benzene_pairs: Vec<_> = pairs.iter().filter(|p| p.core == "c1c([*])cccc1").collect();
        assert_eq!(
            benzene_pairs.len(),
            3,
            "3 molecules → 3 benzene-core MMP pairs: {pairs:?}"
        );
    }
}
