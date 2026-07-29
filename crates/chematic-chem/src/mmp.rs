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

/// Build a fragment molecule from `side` atoms with a wildcard `[*]` attachment.
fn build_fragment_mol(
    mol: &Molecule,
    side: &HashSet<AtomIdx>,
    attach: AtomIdx,
) -> chematic_core::Molecule {
    let mut builder = MoleculeBuilder::new();
    let mut idx_map: HashMap<AtomIdx, AtomIdx> = HashMap::new();

    let mut wc = Atom::new(chematic_core::Element::C);
    wc.wildcard = true;
    let wc_idx = builder.add_atom(wc);

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

    let _ = builder.add_bond(wc_idx, *idx_map.get(&attach).unwrap(), BondOrder::Single);

    for (_, bond) in mol.bonds() {
        if side.contains(&bond.atom1)
            && side.contains(&bond.atom2)
            && let (Some(&n1), Some(&n2)) = (idx_map.get(&bond.atom1), idx_map.get(&bond.atom2))
        {
            let _ = builder.add_bond(n1, n2, bond.order);
        }
    }

    builder.build()
}

/// Build a sub-molecule from `side` atoms and return its canonical SMILES.
fn fragment_smiles(mol: &Molecule, side: &HashSet<AtomIdx>, attach: AtomIdx) -> String {
    canonical_smiles(&build_fragment_mol(mol, side, attach))
}

// ---------------------------------------------------------------------------
// Matched Molecular Series (MMS)
// ---------------------------------------------------------------------------

/// One member of a matched molecular series.
#[derive(Debug, Clone)]
pub struct MmsMember {
    /// Canonical SMILES of the molecule.
    pub smiles: String,
    /// Substituent fragment SMILES (contains `[*]`).
    pub fragment: String,
    /// Molecular weight of the substituent fragment (Da), used as sort key.
    pub mw: f64,
}

/// A matched molecular series — a set of ≥ 3 molecules sharing a common core.
#[derive(Debug, Clone)]
pub struct MmsSeries {
    /// Canonical SMILES of the common core (contains `[*]`).
    pub core: String,
    /// Series members sorted ascending by substituent molecular weight.
    pub members: Vec<MmsMember>,
}

/// Find all matched molecular series in `mols`.
///
/// A series groups ≥ 3 molecules that share the same structural core at a single
/// BRICS cut point, differing only in the substituent at that position.
/// Members are sorted by ascending substituent molecular weight.
pub fn find_mms(mols: &[&Molecule]) -> Vec<MmsSeries> {
    let mol_smiles: Vec<String> = mols.iter().map(|m| canonical_smiles(m)).collect();

    // Build core → Vec<(mol_idx, sub_smiles, sub_mw)>
    let mut index: HashMap<String, Vec<(usize, String, f64)>> = HashMap::new();

    for (mol_idx, mol) in mols.iter().enumerate() {
        // Deduplicate (core, sub) pairs within a single molecule
        let mut seen: HashSet<(String, String)> = HashSet::new();
        for (core_smi, sub_smi, sub_mw) in all_cuts_with_mw(mol) {
            if seen.insert((core_smi.clone(), sub_smi.clone())) {
                index
                    .entry(core_smi)
                    .or_default()
                    .push((mol_idx, sub_smi, sub_mw));
            }
        }
    }

    let mut series_list: Vec<MmsSeries> = Vec::new();

    for (core_smi, entries) in &index {
        // Require ≥ 3 distinct molecules
        let distinct_mols: HashSet<usize> = entries.iter().map(|(i, _, _)| *i).collect();
        if distinct_mols.len() < 3 {
            continue;
        }

        // Per-molecule dedup in the index-building loop already guarantees
        // unique (mol_idx, sub_smi) pairs per core, so no second dedup needed.
        let mut members: Vec<MmsMember> = entries
            .iter()
            .map(|(mol_idx, sub_smi, sub_mw)| MmsMember {
                smiles: mol_smiles[*mol_idx].clone(),
                fragment: sub_smi.clone(),
                mw: *sub_mw,
            })
            .collect();

        members.sort_by(|a, b| a.mw.partial_cmp(&b.mw).unwrap_or(std::cmp::Ordering::Equal));

        series_list.push(MmsSeries {
            core: core_smi.clone(),
            members,
        });
    }

    series_list.sort_by(|a, b| a.core.cmp(&b.core));
    series_list
}

/// Build all (core_smiles, sub_smiles, sub_mw) triples from every BRICS cut of `mol`.
fn all_cuts_with_mw(mol: &Molecule) -> Vec<(String, String, f64)> {
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
        let (sub_smi, sub_mw) = fragment_smiles_with_mw(mol, &sub, at_sub);
        result.push((core_smi, sub_smi, sub_mw));
    }
    result
}

/// Build a fragment molecule and return (canonical_smiles, molecular_weight).
fn fragment_smiles_with_mw(
    mol: &Molecule,
    side: &HashSet<AtomIdx>,
    attach: AtomIdx,
) -> (String, f64) {
    let frag_mol = build_fragment_mol(mol, side, attach);
    let mw = crate::descriptors::molecular_weight(&frag_mol);
    let smi = canonical_smiles(&frag_mol);
    (smi, mw)
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
        // Canonical form of monosubstituted benzene core updated (issue
        // #205): `initial_invariant`'s explicit/implicit-H-count
        // unification changed which atom the canonical DFS starts from for
        // this ring -- same core molecule, same wildcard attachment point.
        let matching: Vec<_> = pairs.iter().filter(|p| p.core == "c1cc([*])ccc1").collect();
        assert_eq!(
            matching.len(),
            1,
            "expected 1 pair with benzene core, got: {pairs:?}"
        );

        let pair = &matching[0];
        // Oracle values updated for the same canonicalization change.
        assert_eq!(
            pair.fragment_a, "[*]CC",
            "ethylbenzene substituent should be [*]CC: {pair:?}"
        );
        assert_eq!(
            pair.fragment_b, "C(C)C[*]",
            "propylbenzene substituent should be C(C)C[*]: {pair:?}"
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
        // Canonical form updated for the individualize-refine ranking rewrite (Round 10).
        let benzene_pairs: Vec<_> = pairs.iter().filter(|p| p.core == "c1cc([*])ccc1").collect();
        assert_eq!(
            benzene_pairs.len(),
            3,
            "3 molecules → 3 benzene-core MMP pairs: {pairs:?}"
        );
    }

    // -----------------------------------------------------------------------
    // MMS tests
    // -----------------------------------------------------------------------

    #[test]
    fn mms_three_alkylbenzenes_one_series() {
        // ethyl-, propyl-, butylbenzene share a benzene core → one MmsSeries
        let a = mol("CCc1ccccc1");
        let b = mol("CCCc1ccccc1");
        let c = mol("CCCCc1ccccc1");
        let series = find_mms(&[&a, &b, &c]);
        let benzene_series: Vec<_> = series
            .iter()
            .filter(|s| s.core == "c1cc([*])ccc1")
            .collect();
        assert_eq!(
            benzene_series.len(),
            1,
            "expected 1 MMS series with benzene core: {series:?}"
        );
        let s = &benzene_series[0];
        assert_eq!(s.members.len(), 3, "series should have 3 members");
        // Members sorted by ascending substituent MW
        let mws: Vec<f64> = s.members.iter().map(|m| m.mw).collect();
        assert!(
            mws.windows(2).all(|w| w[0] <= w[1] + 1e-6),
            "members not sorted by MW: {mws:?}"
        );
    }

    #[test]
    fn mms_two_molecules_no_series() {
        // Only 2 molecules → never enough for a series (need ≥3)
        let a = mol("CCc1ccccc1");
        let b = mol("CCCc1ccccc1");
        let series = find_mms(&[&a, &b]);
        assert!(
            series.is_empty(),
            "2 molecules cannot form a MMS (need ≥3): {series:?}"
        );
    }

    #[test]
    fn mms_member_mw_excludes_wildcard() {
        // Fragment C(C)[*] (ethyl substituent) = C2H5 = 29.06 Da.
        // If the wildcard atom (Element::C, wildcard=true) were included by MW,
        // the result would be ~41 Da (+12 Da bug).
        let a = mol("CCc1ccccc1");
        let b = mol("CCCc1ccccc1");
        let c = mol("CCCCc1ccccc1");
        let series = find_mms(&[&a, &b, &c]);
        let s = series
            .iter()
            .find(|s| s.core == "c1cc([*])ccc1")
            .expect("benzene-core series");
        // Members sorted ascending by MW: ethyl < propyl < butyl
        let ethyl_mw = s.members[0].mw;
        // C2H5 = 2×12.011 + 5×1.008 = 29.062; +12 bug → ~41
        assert!(
            (ethyl_mw - 29.06).abs() < 0.5,
            "ethyl MW should be ~29.06, got {ethyl_mw} (wildcard atom leaked?)"
        );
        // Confirm propyl is larger but not by an extra +12
        let propyl_mw = s.members[1].mw;
        assert!(propyl_mw > ethyl_mw);
        // Propyl C3H7 = 43.09; butyl C4H9 = 57.12 — all < 70
        assert!(s.members[2].mw < 70.0, "butyl MW should be < 70");
    }

    #[test]
    fn mms_no_brics_bonds_no_series() {
        // Benzene has no BRICS bonds → no cuts → no series
        let a = mol("c1ccccc1");
        let b = mol("c1ccncc1");
        let c = mol("c1cccnc1");
        let series = find_mms(&[&a, &b, &c]);
        assert!(
            series.is_empty(),
            "aromatic rings without BRICS bonds → no series: {series:?}"
        );
    }
}
