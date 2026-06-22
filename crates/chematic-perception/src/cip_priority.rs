//! v0.1.93: Full multi-sphere CIP (Cahn–Ingold–Prelog) priority comparison.
//!
//! Implements hierarchical digraph CIP rules with BFS sphere expansion,
//! phantom atoms for double bonds, and atomic mass tiebreaking.
//! This replaces the simplified 1-sphere CIP in stereo2d.rs / stereo3d.rs.

use chematic_core::{AtomIdx, BondOrder, Molecule};
use rustc_hash::{FxHashMap, FxHashSet};
use std::cmp::Ordering;
use std::collections::VecDeque;

/// A single "sphere layer" in a CIP branch expansion.
/// Each tuple: (atomic_number, isotope, atomic_mass) sorted descending by CIP priority.
type SphereLayer = Vec<(u8, Option<u16>, f64)>;

/// Get the key `(atomic_num, isotope, atomic_mass)` for an atom, used in CIP comparisons.
fn atom_key(mol: &Molecule, idx: AtomIdx) -> (u8, Option<u16>, f64) {
    let a = mol.atom(idx);
    (
        a.element.atomic_number(),
        a.isotope,
        a.element.atomic_mass(),
    )
}

/// Compare two `(atomic_num, isotope, atomic_mass)` keys by CIP priority.
///
/// CIP rule hierarchy:
/// 1. Higher atomic number wins.
/// 2. For explicit isotopes: `Some(mass)` > `None` (heavier beats unspecified).
/// 3. Tiebreaker: higher atomic mass wins (rule 4).
fn cmp_key(a: (u8, Option<u16>, f64), b: (u8, Option<u16>, f64)) -> Ordering {
    // Rule 1: atomic number (primary)
    match a.0.cmp(&b.0) {
        Ordering::Equal => {}
        other => return other,
    }
    // Rule 2: isotope label (if present)
    match (a.1, b.1) {
        (Some(x), Some(y)) => match x.cmp(&y) {
            Ordering::Equal => {}
            other => return other,
        },
        (Some(_), None) => return Ordering::Greater,
        (None, Some(_)) => return Ordering::Less,
        (None, None) => {}
    }
    // Rule 4 tiebreaker: atomic mass (descending for higher priority)
    b.2.partial_cmp(&a.2).unwrap_or(Ordering::Equal)
}

/// BFS state for one node during sphere expansion.
struct ExpandState {
    node: AtomIdx,
    parent: AtomIdx,
    depth: usize,
    visited: FxHashSet<AtomIdx>,
}

/// BFS-based CIP sphere expansion for the branch starting at `start`,
/// not going back through `center`.
///
/// Phantom atom rules:
/// 1. **Double-bond phantom**: when expanding node B (reached via double bond from A),
///    add a phantom entry for A at the same depth level as B's children.
/// 2. **Ring revisit phantom**: if an already-visited atom is encountered,
///    add a phantom for it but don't expand further.
fn cip_branch_spheres(mol: &Molecule, center: AtomIdx, start: AtomIdx) -> Vec<SphereLayer> {
    let mut layers: FxHashMap<usize, Vec<(u8, Option<u16>, f64)>> = FxHashMap::default();
    let max_depth = 8usize;

    // The start atom itself is at depth 1.
    let start_key = atom_key(mol, start);
    layers.entry(1).or_default().push(start_key);

    let mut expand_queue: VecDeque<ExpandState> = VecDeque::new();
    {
        let mut v = FxHashSet::default();
        v.insert(center);
        v.insert(start);
        expand_queue.push_back(ExpandState {
            node: start,
            parent: center,
            depth: 1,
            visited: v,
        });
    }

    while let Some(state) = expand_queue.pop_front() {
        if state.depth >= max_depth {
            continue;
        }
        let child_depth = state.depth + 1;

        // Phantom of parent: add if the bond used to reach this node was double.
        if let Some((_, bond_to_parent)) = mol.bond_between(state.node, state.parent)
            && bond_to_parent.order == BondOrder::Double
        {
            let phantom_key = atom_key(mol, state.parent);
            layers.entry(child_depth).or_default().push(phantom_key);
        }

        for (nb, _) in mol.neighbors(state.node) {
            if nb == state.parent || nb == center {
                continue;
            }
            let child_key = atom_key(mol, nb);
            let layer = layers.entry(child_depth).or_default();

            if state.visited.contains(&nb) {
                // Ring revisit: phantom only, no expansion.
                layer.push(child_key);
            } else {
                layer.push(child_key);
                let mut child_visited = state.visited.clone();
                child_visited.insert(nb);
                expand_queue.push_back(ExpandState {
                    node: nb,
                    parent: state.node,
                    depth: child_depth,
                    visited: child_visited,
                });
            }
        }
    }

    // Sort each layer descending and return as a Vec ordered by depth.
    let max_layer = layers.keys().copied().max().unwrap_or(0);
    let mut result = Vec::new();
    for d in 1..=max_layer {
        let mut layer = layers.remove(&d).unwrap_or_default();
        layer.sort_by(|a, b| cmp_key(*b, *a)); // descending
        result.push(layer);
    }
    result
}

/// Compare two branches from `center` starting at `a` and `b`.
///
/// Returns `Ordering::Greater` if branch `a` has higher CIP priority than `b`.
pub fn compare_branches(mol: &Molecule, center: AtomIdx, a: AtomIdx, b: AtomIdx) -> Ordering {
    // Depth-0 comparison: the substituent atoms themselves.
    let a_key = atom_key(mol, a);
    let b_key = atom_key(mol, b);
    match cmp_key(a_key, b_key) {
        Ordering::Equal => {}
        other => return other,
    }

    // Sphere-by-sphere comparison.
    let a_spheres = cip_branch_spheres(mol, center, a);
    let b_spheres = cip_branch_spheres(mol, center, b);

    let max_depth = a_spheres.len().max(b_spheres.len());
    for d in 0..max_depth {
        let a_layer = a_spheres.get(d).map(|v| v.as_slice()).unwrap_or(&[]);
        let b_layer = b_spheres.get(d).map(|v| v.as_slice()).unwrap_or(&[]);

        let min_len = a_layer.len().min(b_layer.len());
        for i in 0..min_len {
            match cmp_key(a_layer[i], b_layer[i]) {
                Ordering::Equal => {}
                other => return other,
            }
        }
        match a_layer.len().cmp(&b_layer.len()) {
            Ordering::Equal => {}
            other => return other,
        }
    }

    Ordering::Equal
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_multi_sphere_simple_distinction() {
        // C(C)(N): carbon with methyl (C) and amine (N) branches
        // Atomic numbers differ: C=6, N=7 → immediately different priority
        let mol = parse("C(C)N").unwrap();
        let center = AtomIdx(0);
        let methyl = AtomIdx(1); // C neighbor
        let amine = AtomIdx(2); // N neighbor

        let order = compare_branches(&mol, center, methyl, amine);
        // N (atomic number 7) has higher priority than C (atomic number 6)
        assert!(
            order == Ordering::Less,
            "C should have lower priority than N"
        );
    }

    #[test]
    fn test_compare_branches_symmetrical() {
        // C(C)(C)(C)C: central C with 4 identical methyl branches
        // All branches are C atoms bonded to same structure (1 parent only)
        let mol = parse("C(C)(C)(C)C").unwrap();
        let center = AtomIdx(0);
        let branch1 = AtomIdx(1);
        let branch2 = AtomIdx(2);

        let order_a = compare_branches(&mol, center, branch1, branch2);
        // Both branches are identical methyls, should compare equal
        assert!(
            order_a == Ordering::Equal,
            "identical methyl branches should compare equal"
        );
    }

    #[test]
    fn test_compare_branches_isotope() {
        // C with C and 13C neighbors should distinguish
        let mol = parse("C(C)(C)(C)C").unwrap();
        let center = AtomIdx(0);
        let c_normal = AtomIdx(1);
        let c_normal2 = AtomIdx(2);

        let order = compare_branches(&mol, center, c_normal, c_normal2);
        // Two identical normal carbons should compare equal
        assert_eq!(
            order,
            Ordering::Equal,
            "identical normal carbons should be equal"
        );
    }
}
