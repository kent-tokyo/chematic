//! Kekulization: assign alternating single/double bonds to aromatic systems.
//!
//! Algorithm:
//! 1. Collect all aromatic atoms and aromatic bonds.
//! 2. Build an adjacency structure restricted to aromatic bonds.
//! 3. Find a maximum matching on the aromatic subgraph using
//!    augmenting paths (DFS-based, O(V*E)).
//! 4. Matched edges → Double bond; unmatched edges → Single bond.
//! 5. Return the updated bond order map, or an error if no valid
//!    Kekulé form exists (i.e. some atom that *must* be doubled is unmatched).
//!
//! Scope: handles all standard aromatic systems (benzene, naphthalene,
//! pyridine, pyrrole, furan, thiophene, indole, imidazole, …).
//! Edmond's blossom algorithm for exotic odd-member cycles is a future todo.

use std::collections::{HashMap, HashSet};

use crate::bond::BondOrder;
use crate::molecule::{AtomIdx, BondIdx, Molecule};

/// Error returned when no valid Kekulé form can be found.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KekuleError {
    pub detail: String,
}

impl core::fmt::Display for KekuleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "kekulization failed: {}", self.detail)
    }
}

impl std::error::Error for KekuleError {}

/// Result of kekulization: a map from BondIdx to the new BondOrder.
///
/// Bonds NOT in this map are unchanged (non-aromatic bonds).
pub type KekuleResult = HashMap<BondIdx, BondOrder>;

/// Kekulize a molecule that contains aromatic bonds.
///
/// Returns a map of aromatic bond indices to their new (Single or Double) orders.
/// All non-aromatic bonds are unchanged and not included in the result.
///
/// If the molecule has no aromatic bonds the result is empty (success, no-op).
pub fn kekulize(mol: &Molecule) -> Result<KekuleResult, KekuleError> {
    // Collect aromatic bonds and the atoms they touch.
    let mut aromatic_bonds: Vec<BondIdx> = Vec::new();
    let mut aromatic_atoms: HashSet<AtomIdx> = HashSet::new();

    for (bidx, bond) in mol.bonds() {
        if bond.order == BondOrder::Aromatic {
            aromatic_bonds.push(bidx);
            aromatic_atoms.insert(bond.atom1);
            aromatic_atoms.insert(bond.atom2);
        }
    }

    if aromatic_bonds.is_empty() {
        return Ok(HashMap::new());
    }

    // Determine which aromatic atoms *must* be in a double bond.
    // An aromatic atom must be double-bonded if it has no lone pair to donate.
    // Heuristic: if the atom has no explicit/implicit H and is carbon or nitrogen-imine,
    // it must be matched.
    // For simplicity: atoms that must_match = those with no spare lone pair =
    //   carbon (C), nitrogen-imine (N in pyridine — has no H when in 6-membered ring).
    // Atoms that can be unmatched (lone-pair donors): O, S, N with H (pyrrole N).
    //
    // Practical rule used here: an atom can be *unmatched* iff it has an explicit
    // H count > 0 OR it is O or S (lone-pair donors).
    // All others must appear in the matching.
    let must_match: HashSet<AtomIdx> = aromatic_atoms
        .iter()
        .copied()
        .filter(|&idx| atom_must_be_matched(mol, idx))
        .collect();

    // Build adjacency list restricted to aromatic bonds BETWEEN must-match atoms.
    //
    // Lone-pair donors (O, S, [nH]) contribute their pi electrons via the lone pair,
    // NOT via a double bond.  Including them in the matching adjacency causes the
    // augmenting-path algorithm to assign double bonds to e.g. [nH]=C in pyrrole or
    // indole, which is chemically wrong and produces incorrect implicit-H counts.
    //
    // Only bonds where BOTH endpoints are in must_match are valid double-bond
    // candidates.  Lone-pair donors remain in aromatic_atoms (so their aromatic bonds
    // become Single in the result) but are excluded from the matching graph.
    let mut adj: HashMap<AtomIdx, Vec<(AtomIdx, BondIdx)>> = HashMap::new();
    for &bidx in &aromatic_bonds {
        let bond = mol.bond(bidx);
        if must_match.contains(&bond.atom1) && must_match.contains(&bond.atom2) {
            adj.entry(bond.atom1).or_default().push((bond.atom2, bidx));
            adj.entry(bond.atom2).or_default().push((bond.atom1, bidx));
        }
    }

    // Run maximum matching via augmenting paths.
    let mut matching: HashMap<AtomIdx, AtomIdx> = HashMap::new(); // atom -> matched_partner

    // Process must-match atoms in a deterministic order (by index) for reproducibility.
    // Non-must-match atoms (lone-pair donors) are skipped — they never initiate
    // augmenting paths and are never placed in the matching.
    let mut sorted_atoms: Vec<AtomIdx> = must_match.iter().copied().collect();
    sorted_atoms.sort();

    // Pass 1: ascending order (primary).
    run_matching_pass(&sorted_atoms, &adj, &mut matching);

    // Pass 2 (fallback): descending order — avoids order-dependent dead-ends.
    //
    // A greedy ascending pass can get stuck on certain ring topologies: the first
    // matched edge blocks an augmenting path that a different starting order would
    // find.  Reversing the order is O(V·E) overhead but resolves many such cases
    // without requiring the full Edmonds blossom algorithm.
    if must_match.iter().any(|&idx| !matching.contains_key(&idx)) {
        matching.clear();
        let mut rev = sorted_atoms.clone();
        rev.reverse();
        run_matching_pass(&rev, &adj, &mut matching);
    }

    // Verify that all must_match atoms are matched.
    for &idx in &must_match {
        if !matching.contains_key(&idx) {
            return Err(KekuleError {
                detail: format!(
                    "atom {} ({}) cannot be assigned a double bond",
                    idx.0,
                    mol.atom(idx).element.symbol()
                ),
            });
        }
    }

    // Build the result map: matched aromatic bonds → Double, rest → Single.
    let mut double_bonds: HashSet<BondIdx> = HashSet::new();
    for (&atom, &partner) in &matching {
        if atom >= partner {
            continue; // each pair shows up twice in `matching`; visit one orientation
        }
        if let Some((bidx, _)) = mol.bond_between(atom, partner)
            && mol.bond(bidx).order == BondOrder::Aromatic
        {
            double_bonds.insert(bidx);
        }
    }

    let result: KekuleResult = aromatic_bonds
        .iter()
        .map(|&bidx| {
            let order = if double_bonds.contains(&bidx) {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            (bidx, order)
        })
        .collect();
    Ok(result)
}

/// Apply a Kekulé result to a molecule, returning a new Molecule with updated bond orders.
///
/// Aromatic flags on atoms are *not* cleared (the molecule retains the aromaticity
/// annotation; only bond orders change).
pub fn apply_kekule(mol: &Molecule, kekule: &KekuleResult) -> Molecule {
    use crate::molecule::MoleculeBuilder;

    let mut builder = MoleculeBuilder::new();

    // Re-add all atoms.
    for (_, atom) in mol.atoms() {
        builder.add_atom(atom.clone());
    }

    // Re-add bonds, substituting updated orders for aromatic bonds.
    for (bidx, bond) in mol.bonds() {
        let order = kekule.get(&bidx).copied().unwrap_or(bond.order);
        builder
            .add_bond(bond.atom1, bond.atom2, order)
            .expect("duplicate bond during apply_kekule");
    }

    builder.build()
}

/// Attempt to find an augmenting path starting from `start` and update `matching`.
///
/// Uses iterative BFS with parent-pointer path reconstruction to avoid stack
/// overflow on large aromatic systems (wasm32 default stack is ~1 MB).
/// `visited` must already contain `start` (prevents root re-entry in odd cycles).
///
/// Returns true if an augmenting path was found and the matching was updated.
fn augment(
    start: AtomIdx,
    adj: &HashMap<AtomIdx, Vec<(AtomIdx, BondIdx)>>,
    matching: &mut HashMap<AtomIdx, AtomIdx>,
    visited: &mut HashSet<AtomIdx>,
) -> bool {
    // parent[u] = v means "u was reached from v in the BFS tree"
    let mut parent: HashMap<AtomIdx, AtomIdx> = HashMap::new();
    let mut queue: std::collections::VecDeque<AtomIdx> = std::collections::VecDeque::new();
    queue.push_back(start);

    'bfs: while let Some(v) = queue.pop_front() {
        let Some(neighbors) = adj.get(&v) else {
            continue;
        };
        for &(u, _) in neighbors {
            if !visited.insert(u) {
                continue;
            }
            parent.insert(u, v);

            match matching.get(&u).copied() {
                None => {
                    // Found a free vertex — trace back through parent pointers
                    // and flip every edge along the augmenting path.
                    let mut cur = u;
                    loop {
                        let prev = parent[&cur];
                        let prev_old_match = matching.get(&prev).copied();
                        matching.insert(prev, cur);
                        matching.insert(cur, prev);
                        match prev_old_match {
                            None | Some(_) if prev == start => break,
                            Some(m) => cur = m,
                            None => break,
                        }
                    }
                    break 'bfs;
                }
                Some(partner) => {
                    // u is matched to partner — explore further from partner.
                    if visited.insert(partner) {
                        parent.insert(partner, u);
                        queue.push_back(partner);
                    }
                }
            }
        }
    }

    // Augmentation succeeded iff start is now matched.
    matching.contains_key(&start)
}

/// Run a single augmenting-path pass over `atoms` and update `matching`.
///
/// For each unmatched atom in `atoms`, attempts to find an augmenting path using
/// the BFS-based `augment()` function. Extracted so that `kekulize()` can try
/// multiple orderings (ascending → descending) without code duplication.
fn run_matching_pass(
    atoms: &[AtomIdx],
    adj: &HashMap<AtomIdx, Vec<(AtomIdx, BondIdx)>>,
    matching: &mut HashMap<AtomIdx, AtomIdx>,
) {
    for &start in atoms {
        if matching.contains_key(&start) {
            continue;
        }
        let mut visited: HashSet<AtomIdx> = HashSet::new();
        visited.insert(start);
        augment(start, adj, matching, &mut visited);
    }
}

/// Determine whether an aromatic atom *must* appear in the matching
/// (i.e. requires a double bond for a valid Kekulé form).
///
/// An atom can be unmatched if it contributes a lone pair to the aromatic system:
/// - O (furan-type oxygen)
/// - S (thiophene-type sulfur)
/// - N with an H (pyrrole-type nitrogen: [nH])
/// - Se, As aromatic analogs
/// - Any aromatic atom that already has an exocyclic double bond (e.g. the
///   carbonyl carbon in coumarin/warfarin `c=O` fused into an aromatic ring).
///   Such an atom's pi contribution comes from conjugation with the exocyclic
///   bond; no additional ring double bond is needed or possible.
///
/// Everything else (C, N without H like pyridine) must be matched.
fn atom_must_be_matched(mol: &Molecule, idx: AtomIdx) -> bool {
    let atom = mol.atom(idx);
    match atom.element.atomic_number() {
        // O, S, Se always donate a lone pair → don't need a double bond.
        8 | 16 | 34 => false,
        // Boron aromatic (rare) — can donate lone pair.
        5 => false,
        // N with explicit H ([nH]) is a lone-pair donor.
        7 if matches!(atom.hydrogen_count, Some(h) if h > 0) => false,
        // Anionic aromatic N ([n-]) also donates its lone pair; the extra electron
        // occupies the lone-pair slot rather than a ring π bond (same as [nH]).
        7 if atom.charge < 0 => false,
        // Neutral N with a non-aromatic substituent (e.g. N-methyl in caffeine) is a lone-pair
        // donor: the substituent "replaces" the H, and the N contributes its lone pair to
        // aromaticity rather than a π bond (same as [nH] in pyrrole).
        // Charged N (pyridinium [n+], N-oxide) must still be matched even with a substituent.
        7 if atom.charge == 0
            && mol
                .neighbors(idx)
                .any(|(_, bidx)| mol.bond(bidx).order != BondOrder::Aromatic) => false,
        // Bare aromatic N with only aromatic bonds (pyridine-type): must be matched.
        7 => true,
        // Any anionic aromatic atom (e.g. cyclopentadienyl [cH-]) donates its lone pair.
        _ if atom.charge < 0 => false,
        // Any atom (typically C) that already has an exocyclic π bond (e.g. `c(=O)` in
        // a heterocyclic carbonyl) cannot also carry a ring double bond — its π-slot is
        // occupied by the exocyclic bond. Such atoms contribute via conjugation, like
        // lone-pair donors, and should not be forced into the matching.
        _ if mol.neighbors(idx).any(|(_, bidx)| {
            let o = mol.bond(bidx).order;
            o == BondOrder::Double || o == BondOrder::Triple
        }) => false,
        // All other atoms must appear in the matching.
        _ => true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::element::Element;
    use crate::molecule::MoleculeBuilder;

    /// Build benzene as a fully aromatic SMILES-style molecule.
    fn benzene() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// Build pyridine (5 aromatic C + 1 aromatic N, ring).
    fn pyridine() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::aromatic(Element::C));
        let c2 = b.add_atom(Atom::aromatic(Element::C));
        let c3 = b.add_atom(Atom::aromatic(Element::C));
        let n = b.add_atom(Atom::aromatic(Element::N));
        let c4 = b.add_atom(Atom::aromatic(Element::C));
        let c5 = b.add_atom(Atom::aromatic(Element::C));
        let atoms = [c1, c2, c3, n, c4, c5];
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// Build furan (4 aromatic C + 1 aromatic O, ring).
    fn furan() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::aromatic(Element::O));
        let c1 = b.add_atom(Atom::aromatic(Element::C));
        let c2 = b.add_atom(Atom::aromatic(Element::C));
        let c3 = b.add_atom(Atom::aromatic(Element::C));
        let c4 = b.add_atom(Atom::aromatic(Element::C));
        let atoms = [o, c1, c2, c3, c4];
        for i in 0..5 {
            b.add_bond(atoms[i], atoms[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// Build pyrrole ([nH] + 4 aromatic C, ring).
    fn pyrrole() -> Molecule {
        let mut b = MoleculeBuilder::new();
        // N with 1 H — bracket-style (hydrogen_count = Some(1))
        let mut n_atom = Atom::aromatic(Element::N);
        n_atom.hydrogen_count = Some(1);
        let n = b.add_atom(n_atom);
        let c1 = b.add_atom(Atom::aromatic(Element::C));
        let c2 = b.add_atom(Atom::aromatic(Element::C));
        let c3 = b.add_atom(Atom::aromatic(Element::C));
        let c4 = b.add_atom(Atom::aromatic(Element::C));
        let atoms = [n, c1, c2, c3, c4];
        for i in 0..5 {
            b.add_bond(atoms[i], atoms[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    #[test]
    fn test_kekulize_benzene() {
        let mol = benzene();
        let result = kekulize(&mol).expect("benzene kekulization failed");
        assert_eq!(result.len(), 6); // all 6 aromatic bonds assigned

        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        let singles = result.values().filter(|&&o| o == BondOrder::Single).count();
        assert_eq!(doubles, 3, "benzene must have 3 double bonds");
        assert_eq!(singles, 3, "benzene must have 3 single bonds");
    }

    #[test]
    fn test_kekulize_pyridine() {
        let mol = pyridine();
        let result = kekulize(&mol).expect("pyridine kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 3, "pyridine must have 3 double bonds");
    }

    #[test]
    fn test_kekulize_furan() {
        let mol = furan();
        let result = kekulize(&mol).expect("furan kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 2, "furan must have 2 double bonds");
    }

    #[test]
    fn test_kekulize_pyrrole() {
        let mol = pyrrole();
        let result = kekulize(&mol).expect("pyrrole kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 2, "pyrrole must have 2 double bonds");
    }

    #[test]
    fn test_kekulize_naphthalene() {
        // 10 aromatic C, 11 aromatic bonds (fused bicyclic)
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..10)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        // Ring 1: 0-1-2-3-4-9
        let ring1 = [0, 1, 2, 3, 4, 9];
        for i in 0..ring1.len() {
            b.add_bond(
                atoms[ring1[i]],
                atoms[ring1[(i + 1) % ring1.len()]],
                BondOrder::Aromatic,
            )
            .unwrap();
        }
        // Ring 2: 4-5-6-7-8-9 (shares bond 4-9)
        let ring2 = [4, 5, 6, 7, 8, 9];
        for i in 0..ring2.len() {
            let a = atoms[ring2[i]];
            let bb = atoms[ring2[(i + 1) % ring2.len()]];
            // Skip already-added bond (4-9)
            if mol_has_no_bond_yet(&b, a, bb) {
                b.add_bond(a, bb, BondOrder::Aromatic).unwrap();
            }
        }
        let mol = b.build();
        let result = kekulize(&mol).expect("naphthalene kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 5, "naphthalene must have 5 double bonds");
    }

    #[test]
    fn test_apply_kekule() {
        let mol = benzene();
        let kekule = kekulize(&mol).unwrap();
        let kekule_mol = apply_kekule(&mol, &kekule);

        // After applying, no aromatic bonds should remain.
        for (_, bond) in kekule_mol.bonds() {
            assert_ne!(
                bond.order,
                BondOrder::Aromatic,
                "apply_kekule should remove all aromatic bonds"
            );
        }
    }

    #[test]
    fn test_no_aromatic_bonds_noop() {
        // A molecule with no aromatic bonds should return empty result.
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let mol = b.build();
        let result = kekulize(&mol).unwrap();
        assert!(result.is_empty());
    }

    // Helper: check that the builder does not yet have a bond between a and b.
    fn mol_has_no_bond_yet(b: &MoleculeBuilder, a: AtomIdx, bb: AtomIdx) -> bool {
        for (_, partner) in b.atom_neighbors(a) {
            if partner == bb {
                return false;
            }
        }
        true
    }

    // B6 coverage: exotic ring systems that require correct matching on non-bipartite graphs.

    /// Azulene: fused 5+7 bicyclic, 10 aromatic C, 11 bonds.
    /// Valid Kekulé form exists → must yield 5 double bonds.
    #[test]
    fn test_kekulize_azulene() {
        let mut b = MoleculeBuilder::new();
        let a: Vec<_> = (0..10)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        // 5-ring: 0-1-2-3-4-0
        for i in 0..5 {
            b.add_bond(a[i], a[(i + 1) % 5], BondOrder::Aromatic).unwrap();
        }
        // 7-ring extras (shares bond 0-4): 4-5-6-7-8-9-0
        for (x, y) in [(4usize, 5usize), (5, 6), (6, 7), (7, 8), (8, 9), (9, 0)] {
            if mol_has_no_bond_yet(&b, a[x], a[y]) {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
        }
        let mol = b.build();
        let result = kekulize(&mol).expect("azulene kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 5, "azulene needs 5 double bonds");
    }

    /// Acenaphthylene: naphthalene + fused cyclopentadiene bridge, 12 aromatic C.
    #[test]
    fn test_kekulize_acenaphthylene() {
        // Connectivity: naphthalene core (atoms 0-9) + bridge atoms 10,11
        // Naphthalene: ring1=0-1-2-3-4-9, ring2=4-5-6-7-8-9
        // Bridge: 0-11-10-1 (the 5-ring across the 1,8 positions of naphthalene)
        let mut b = MoleculeBuilder::new();
        let a: Vec<_> = (0..12)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        // Naphthalene ring 1: 0-1-2-3-4-9-0
        for (x, y) in [(0,1),(1,2),(2,3),(3,4),(4,9),(9,0)] {
            b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
        }
        // Naphthalene ring 2: 4-5-6-7-8-9 (4-9 shared)
        for (x, y) in [(4,5),(5,6),(6,7),(7,8),(8,9)] {
            b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
        }
        // Bridge: 0-11-10-1
        for (x, y) in [(0,11),(11,10),(10,1)] {
            b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
        }
        let mol = b.build();
        let result = kekulize(&mol).expect("acenaphthylene kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 6, "acenaphthylene needs 6 double bonds");
    }

    // ---- Edge-case tests added v0.1.100 ----
    //
    // Build complex PAH manually to test edge cases in the matching algorithm.

    /// Biphenylene: two benzene rings connected by a cyclobutadiene bridge.
    ///
    /// Topology: Ring A (0–5), Ring B (6–11), bridge bonds (0–11) and (5–6).
    /// The 4-membered ring 0-5-6-11-0 is the challenging part.
    fn biphenylene() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let a: Vec<_> = (0..12)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        // Ring A (6): 0-1-2-3-4-5-0
        for i in 0..6 {
            b.add_bond(a[i], a[(i + 1) % 6], BondOrder::Aromatic).unwrap();
        }
        // Ring B (6): 6-7-8-9-10-11-6
        for i in 0..6 {
            b.add_bond(a[6 + i], a[6 + (i + 1) % 6], BondOrder::Aromatic).unwrap();
        }
        // 4-membered bridge: closes ring 0-5-6-11-0
        b.add_bond(a[5], a[6], BondOrder::Aromatic).unwrap();
        b.add_bond(a[0], a[11], BondOrder::Aromatic).unwrap();
        b.build()
    }

    /// Naphtho-4-ring: three fused 6-membered rings sharing a common bond sequence.
    /// Linear anthracene topology (0-1-2-3-4-5, 5-6-7-8-9-4, 6-10-11-12-13-7).
    fn anthracene() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let a: Vec<_> = (0..14)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        // Ring A: 0-1-2-3-4-5-0
        for i in 0..6 {
            b.add_bond(a[i], a[(i + 1) % 6], BondOrder::Aromatic).unwrap();
        }
        // Ring B: 5-4-9-8-7-6-5 (shares bond 4-5 with Ring A)
        for (x, y) in [(4,9),(9,8),(8,7),(7,6),(6,5)] {
            b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
        }
        // Ring C: 6-7-13-12-11-10-6 (shares bond 6-7 with Ring B)
        for (x, y) in [(7,13),(13,12),(12,11),(11,10),(10,6)] {
            b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
        }
        b.build()
    }

    #[test]
    fn test_kekulize_biphenylene() {
        // Biphenylene: 4-membered cyclobutadiene bridge between two benzenes.
        // 12 aromatic C → 6 double bonds.
        let mol = biphenylene();
        let result = kekulize(&mol).expect("biphenylene kekulization should succeed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 6, "biphenylene needs 6 double bonds");
    }

    #[test]
    fn test_kekulize_anthracene() {
        // Anthracene (C14H10): 3 linearly fused 6-membered rings.  7 double bonds.
        let mol = anthracene();
        let result = kekulize(&mol).expect("anthracene kekulization should succeed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 7, "anthracene needs 7 double bonds");
    }

    #[test]
    fn test_kekulize_biphenylene_double_bond_count() {
        // Cross-check: all 14 bonds should be assigned single or double.
        let mol = biphenylene();
        let result = kekulize(&mol).expect("biphenylene kekulization");
        let singles = result.values().filter(|&&o| o == BondOrder::Single).count();
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(singles + doubles, 14, "biphenylene has 14 aromatic bonds");
    }

    #[test]
    fn test_kekulize_large_fused_6rings() {
        // 4 fused 6-membered rings (pyrene-like) in a simple linear topology.
        // Simulates a large all-even-ring PAH that the BFS should handle easily.
        let mol = {
            let mut b = MoleculeBuilder::new();
            let a: Vec<_> = (0..16)
                .map(|_| b.add_atom(Atom::aromatic(Element::C)))
                .collect();
            // Ring 0-1-2-3-4-5
            for i in 0..6 { b.add_bond(a[i], a[(i+1)%6], BondOrder::Aromatic).unwrap(); }
            // Ring 5-4-9-8-7-6
            for (x,y) in [(4,9),(9,8),(8,7),(7,6),(6,5)] {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
            // Ring 1-2-11-10-13-12
            for (x,y) in [(2,11),(11,10),(10,13),(13,12),(12,1)] {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
            // Ring 6-7-15-14-11-2 (closed)
            for (x,y) in [(7,15),(15,14),(14,11)] {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
            b.build()
        };
        let result = kekulize(&mol).expect("4-ring PAH kekulization should succeed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert!(doubles >= 6, "4-ring PAH needs at least 6 double bonds, got {doubles}");
    }

    #[test]
    fn test_kekulize_deterministic() {
        // kekulize() is deterministic: rebuilding the same molecule gives the same count.
        let mol1 = biphenylene();
        let mol2 = biphenylene();
        let r1 = kekulize(&mol1).expect("pass1");
        let r2 = kekulize(&mol2).expect("pass2");
        assert_eq!(
            r1.values().filter(|&&o| o == BondOrder::Double).count(),
            r2.values().filter(|&&o| o == BondOrder::Double).count(),
            "kekulization must be deterministic"
        );
    }

    /// Fluoranthene: 16 aromatic C in a fused 5+6+6+6 system.
    #[test]
    fn test_kekulize_fluoranthene() {
        // Simplified: 3 fused 6-rings + 1 fused 5-ring sharing atoms
        // Build as corannulene-like: two naphthalene units bridged by a 5-ring
        // Atoms 0-15 (16 aromatic C), total bonds = 19 (for fluoranthene)
        // Use a simplified topology that gives the right ring count
        // 6-ring A: 0-1-2-3-4-5
        // 6-ring B: 0-5-6-7-8-9
        // 6-ring C: 2-3-10-11-12-13
        // 5-ring D: 0-9-14-15-1 (bridge)
        let mut b = MoleculeBuilder::new();
        let a: Vec<_> = (0..16)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        // Ring A (6-ring): 0-1-2-3-4-5-0
        for i in 0..6 { b.add_bond(a[i], a[(i+1)%6], BondOrder::Aromatic).unwrap(); }
        // Ring B (6-ring): 0-5-6-7-8-9-0
        for (x,y) in [(5,6),(6,7),(7,8),(8,9),(9,0)] {
            if mol_has_no_bond_yet(&b, a[x], a[y]) {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
        }
        // Ring C (6-ring): 1-2-10-11-12-13-1
        for (x,y) in [(2,10),(10,11),(11,12),(12,13),(13,1)] {
            if mol_has_no_bond_yet(&b, a[x], a[y]) {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
        }
        // Ring D (5-ring): 9-8-14-15-13-9
        for (x,y) in [(8,14),(14,15),(15,13),(13,9)] {
            if mol_has_no_bond_yet(&b, a[x], a[y]) {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
        }
        let mol = b.build();
        let result = kekulize(&mol).expect("fluoranthene-like kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 8, "fluoranthene-like structure needs 8 double bonds");
    }
}
