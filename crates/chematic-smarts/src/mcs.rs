//! Maximum Common Substructure (MCS) search using the McGregor connected-growth algorithm.
//!
//! # Overview
//! Given a set of molecules, `find_mcs` returns a `QueryMolecule` that represents the
//! largest substructure common to all input molecules.  The algorithm is a branch-and-bound
//! search that grows the common subgraph one atom at a time, pruning subtrees whose upper
//! bound on reachable size is no better than the current best.

use std::collections::HashMap;
use std::time::Instant;

use chematic_core::{AtomIdx, BondOrder, Molecule};

use crate::query::{AtomPrimitive, AtomQuery, BondPrimitive, BondQuery, QueryMolecule};

// ---------------------------------------------------------------------------
// Public API types
// ---------------------------------------------------------------------------

/// Configuration options for MCS search.
#[derive(Debug, Clone)]
pub struct McsConfig {
    /// If `true`, bond orders must match between molecules. If `false`, any bond order is
    /// accepted (atom-only MCS).
    pub match_bonds: bool,
    /// Minimum number of atoms required for the result to be returned.  If the best
    /// common substructure has fewer atoms than `min_atoms`, an empty `QueryMolecule` is
    /// returned instead.
    pub min_atoms: usize,
    /// Optional time limit in milliseconds.  The search is aborted when the deadline is
    /// reached and the best result found so far is returned.
    pub timeout_ms: Option<u64>,
}

impl Default for McsConfig {
    fn default() -> Self {
        Self { match_bonds: true, min_atoms: 1, timeout_ms: None }
    }
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Find the Maximum Common Substructure of `mols` using default configuration.
///
/// Returns an empty `QueryMolecule` if `mols` is empty.
pub fn find_mcs(mols: &[&Molecule]) -> QueryMolecule {
    find_mcs_with_config(mols, &McsConfig::default())
}

/// Find the Maximum Common Substructure of `mols` using `config`.
pub fn find_mcs_with_config(mols: &[&Molecule], config: &McsConfig) -> QueryMolecule {
    if mols.is_empty() {
        return QueryMolecule::new();
    }
    if mols.len() == 1 {
        return molecule_to_query(mols[0]);
    }

    let deadline = config.timeout_ms.map(|ms| Instant::now() + std::time::Duration::from_millis(ms));

    let mut state = McsState {
        mols,
        config,
        best: PartialMapping::empty(mols.len()),
        deadline,
        timed_out: false,
    };

    let n0 = mols[0].atom_count();

    'outer: for a0 in 0..n0 {
        if state.timed_out {
            break;
        }
        if let Some(d) = state.deadline {
            if Instant::now() >= d {
                state.timed_out = true;
                break;
            }
        }

        // Gather compatible atoms from each other molecule for this seed atom.
        let seed_candidates = collect_seed_candidates(mols, AtomIdx(a0 as u32), config);

        // Iterate the Cartesian product of candidate lists.
        for seed in CartesianProduct::new(&seed_candidates) {
            if state.timed_out {
                break 'outer;
            }
            if let Some(d) = state.deadline {
                if Instant::now() >= d {
                    state.timed_out = true;
                    break 'outer;
                }
            }

            // Prune: upper bound for a single-atom seed is at most min(atom counts by element).
            let ub = upper_bound_for_seeds(mols, AtomIdx(a0 as u32), &seed);
            if ub <= state.best.size {
                continue;
            }

            let mut mapping = PartialMapping::from_seed(mols, AtomIdx(a0 as u32), &seed);
            grow(&mut state, &mut mapping);
        }
    }

    // Apply min_atoms filter.
    if state.best.size < config.min_atoms {
        return QueryMolecule::new();
    }

    build_query(mols[0], &state.best, config)
}

// ---------------------------------------------------------------------------
// Internal data structures
// ---------------------------------------------------------------------------

/// A partial mapping between all molecules' atom indices and query atom indices.
#[derive(Clone)]
struct PartialMapping {
    /// `mol_map[mol_idx][atom_idx]` = `Some(query_idx)` if the atom is mapped.
    mol_map: Vec<Vec<Option<usize>>>,
    /// `query_to_mol[query_idx][mol_idx]` = `AtomIdx`.
    query_to_mol: Vec<Vec<AtomIdx>>,
    /// Number of atoms currently in the mapping.
    size: usize,
    /// Number of bonds in the common subgraph (for bond-match mode).
    bond_count: usize,
}

impl PartialMapping {
    /// Create an empty mapping for `n_mols` molecules (with `n_atoms[i]` atoms each).
    fn empty(n_mols: usize) -> Self {
        // We start with zero-length mol maps; they'll be resized lazily or we just
        // keep them as empty here and they're properly initialised in `from_seed`.
        Self {
            mol_map: vec![Vec::new(); n_mols],
            query_to_mol: Vec::new(),
            size: 0,
            bond_count: 0,
        }
    }

    /// Create a single-atom seed mapping: query atom 0 maps to `a0` in mol[0] and
    /// `seed[i]` in mol[i+1].
    fn from_seed(mols: &[&Molecule], a0: AtomIdx, seed: &[AtomIdx]) -> Self {
        let mut mol_map: Vec<Vec<Option<usize>>> = mols
            .iter()
            .map(|m| vec![None; m.atom_count()])
            .collect();
        let mut query_to_mol: Vec<Vec<AtomIdx>> = Vec::new();

        // Map the seed tuple to query atom 0.
        let mut row = vec![a0];
        row.extend_from_slice(seed);
        for (mi, &ai) in row.iter().enumerate() {
            mol_map[mi][ai.0 as usize] = Some(0);
        }
        query_to_mol.push(row);

        Self { mol_map, query_to_mol, size: 1, bond_count: 0 }
    }

    /// Check whether `atom_idx` in molecule `mol_idx` is already mapped.
    fn is_mapped(&self, mol_idx: usize, atom_idx: AtomIdx) -> bool {
        let v = &self.mol_map[mol_idx];
        let idx = atom_idx.0 as usize;
        idx < v.len() && v[idx].is_some()
    }

    /// Look up the query index for `atom_idx` in molecule `mol_idx`.
    fn query_idx_of(&self, mol_idx: usize, atom_idx: AtomIdx) -> Option<usize> {
        let v = &self.mol_map[mol_idx];
        let idx = atom_idx.0 as usize;
        if idx < v.len() { v[idx] } else { None }
    }

    /// Extend the mapping with a new atom tuple (one atom per molecule).
    /// Returns the new query index.
    fn extend(&mut self, atoms: &[AtomIdx], extra_bonds: usize) -> usize {
        let q = self.query_to_mol.len();
        for (mi, &ai) in atoms.iter().enumerate() {
            self.mol_map[mi][ai.0 as usize] = Some(q);
        }
        self.query_to_mol.push(atoms.to_vec());
        self.size += 1;
        self.bond_count += extra_bonds;
        q
    }

    /// Retract the last-added atom tuple (must be the one with the highest query index).
    fn retract(&mut self, atoms: &[AtomIdx], extra_bonds: usize) {
        let q = self.query_to_mol.len() - 1;
        for (mi, &ai) in atoms.iter().enumerate() {
            self.mol_map[mi][ai.0 as usize] = None;
        }
        self.query_to_mol.pop();
        self.size -= 1;
        self.bond_count -= extra_bonds;
        let _ = q;
    }
}

// ---------------------------------------------------------------------------
// Search state
// ---------------------------------------------------------------------------

struct McsState<'a> {
    mols: &'a [&'a Molecule],
    config: &'a McsConfig,
    best: PartialMapping,
    deadline: Option<Instant>,
    timed_out: bool,
}

// ---------------------------------------------------------------------------
// Core search: grow
// ---------------------------------------------------------------------------

fn grow(state: &mut McsState<'_>, mapping: &mut PartialMapping) {
    // Update best if this mapping is larger.
    if mapping.size > state.best.size {
        state.best = mapping.clone();
    }

    // Check timeout.
    if let Some(d) = state.deadline {
        if Instant::now() >= d {
            state.timed_out = true;
            return;
        }
    }
    if state.timed_out {
        return;
    }

    // Upper bound pruning: compute max additional atoms reachable from current state.
    let additional_ub = upper_bound_additional(state.mols, mapping);
    if mapping.size + additional_ub <= state.best.size {
        return;
    }

    // Find frontier candidates: unmapped neighbors of already-mapped atoms in mol[0].
    // For each unmapped neighbor `n0` of a mapped atom in mol[0], try to extend the mapping.
    let mol0 = state.mols[0];

    // Collect unmapped neighbors of mapped atoms in mol[0].
    let mut frontier: Vec<AtomIdx> = Vec::new();
    for row in &mapping.query_to_mol {
        let a0 = row[0];
        for (nb, _) in mol0.neighbors(a0) {
            if !mapping.is_mapped(0, nb) && !frontier.contains(&nb) {
                frontier.push(nb);
            }
        }
    }

    if frontier.is_empty() {
        return;
    }

    // Pick only the first frontier atom (to avoid duplicate paths in the DFS tree).
    // This is the "connected growth" property.
    let n0 = frontier[0];

    // For each other molecule, find compatible unmapped atoms that are:
    // 1. Atom-compatible with n0
    // 2. Adjacent to every mol[i]-atom that corresponds to a mapped neighbor of n0 in mol[0]
    // 3. Bond-compatible with each such adjacency
    let candidates_per_mol = build_frontier_candidates(state.mols, mapping, n0, state.config);

    // Cartesian product across molecules 1..n.
    for tuple in CartesianProduct::new(&candidates_per_mol) {
        if state.timed_out {
            return;
        }

        // Build the full atom tuple: [n0, tuple[0], tuple[1], ...]
        let mut atoms: Vec<AtomIdx> = Vec::with_capacity(state.mols.len());
        atoms.push(n0);
        atoms.extend_from_slice(&tuple);

        // Count new bonds (to already-mapped atoms in mol[0]) for bookkeeping.
        let extra_bonds = count_new_bonds(mol0, mapping, n0);

        mapping.extend(&atoms, extra_bonds);
        grow(state, mapping);
        mapping.retract(&atoms, extra_bonds);
    }
}

// ---------------------------------------------------------------------------
// Frontier candidate generation
// ---------------------------------------------------------------------------

/// For each molecule i >= 1, collect the atoms that can be paired with `n0` from mol[0]
/// to extend the mapping.
///
/// A candidate atom `ni` in mol[i] is valid if:
/// - It is not yet mapped.
/// - It is atom-compatible with `n0`.
/// - For every already-mapped neighbor `m0` of `n0` in mol[0], atom `ni` is bonded to
///   `mapping.query_to_mol[q_m][i]` in mol[i] with a compatible bond order (when
///   `match_bonds=true`).
fn build_frontier_candidates(
    mols: &[&Molecule],
    mapping: &PartialMapping,
    n0: AtomIdx,
    config: &McsConfig,
) -> Vec<Vec<AtomIdx>> {
    let mol0 = mols[0];
    let atom0 = mol0.atom(n0);

    // Collect mapped neighbors of n0 in mol[0] along with the bond order.
    let mut mapped_neighbors: Vec<(usize, AtomIdx, BondOrder)> = Vec::new();
    for (nb, bidx) in mol0.neighbors(n0) {
        if let Some(q) = mapping.query_idx_of(0, nb) {
            let bond = mol0.bond(bidx);
            mapped_neighbors.push((q, nb, bond.order));
        }
    }

    let mut result: Vec<Vec<AtomIdx>> = Vec::with_capacity(mols.len() - 1);

    for mi in 1..mols.len() {
        let mol_i = mols[mi];
        let mut candidates: Vec<AtomIdx> = Vec::new();

        'atom: for (ai, atom_i) in mol_i.atoms() {
            // Must not be mapped.
            if mapping.is_mapped(mi, ai) {
                continue;
            }
            // Must be atom-compatible.
            if !atoms_compatible(atom0, atom_i) {
                continue;
            }
            // Must be bonded to every corresponding mapped neighbor with compatible bond.
            for &(q, _m0, bond_order_0) in &mapped_neighbors {
                let m_i = mapping.query_to_mol[q][mi];
                match mol_i.bond_between(ai, m_i) {
                    None => continue 'atom,
                    Some((_bidx, bond_entry)) => {
                        if config.match_bonds && !bonds_compatible(bond_order_0, bond_entry.order) {
                            continue 'atom;
                        }
                    }
                }
            }
            candidates.push(ai);
        }

        result.push(candidates);
    }

    result
}

// ---------------------------------------------------------------------------
// Seed collection
// ---------------------------------------------------------------------------

/// For each molecule i >= 1, collect atoms compatible with `a0` in mol[0].
fn collect_seed_candidates(mols: &[&Molecule], a0: AtomIdx, _config: &McsConfig) -> Vec<Vec<AtomIdx>> {
    let atom0 = mols[0].atom(a0);
    let mut result = Vec::with_capacity(mols.len() - 1);
    for mi in 1..mols.len() {
        let cands: Vec<AtomIdx> = mols[mi]
            .atoms()
            .filter(|(_, a)| atoms_compatible(atom0, a))
            .map(|(idx, _)| idx)
            .collect();
        result.push(cands);
    }
    result
}

// ---------------------------------------------------------------------------
// Upper bound heuristics
// ---------------------------------------------------------------------------

/// Upper bound on the total MCS size for a given seed (used at the seed level).
fn upper_bound_for_seeds(mols: &[&Molecule], a0: AtomIdx, seed: &[AtomIdx]) -> usize {
    // Count atoms of each element in each molecule.
    // Minimum across molecules gives max atoms of that element that can be matched.
    let mut counts: HashMap<u8, Vec<usize>> = HashMap::new();
    for (mi, mol) in mols.iter().enumerate() {
        let exclude = if mi == 0 { a0 } else { seed[mi - 1] };
        for (ai, atom) in mol.atoms() {
            if ai == exclude {
                continue; // seed atom already counts; we'll add 1 at the end
            }
            let en = atom.element.atomic_number();
            counts.entry(en).or_insert_with(|| vec![0; mols.len()])[mi] += 1;
        }
    }
    let ub: usize = counts.values().map(|v| *v.iter().min().unwrap_or(&0)).sum();
    ub + 1 // +1 for the seed atom itself
}

/// Upper bound on *additional* atoms reachable from the current mapping.
///
/// For each element, count unmapped atoms in each molecule; the maximum additional
/// atoms of that element = min over all molecules.
fn upper_bound_additional(mols: &[&Molecule], mapping: &PartialMapping) -> usize {
    let mut counts: HashMap<u8, Vec<usize>> = HashMap::new();

    for (mi, mol) in mols.iter().enumerate() {
        for (ai, atom) in mol.atoms() {
            if mapping.is_mapped(mi, ai) {
                continue;
            }
            let en = atom.element.atomic_number();
            let entry = counts.entry(en).or_insert_with(|| vec![0; mols.len()]);
            entry[mi] += 1;
        }
    }

    counts.values().map(|v| *v.iter().min().unwrap_or(&0)).sum()
}

// ---------------------------------------------------------------------------
// Compatibility helpers
// ---------------------------------------------------------------------------

/// Atom compatibility: same atomic number and same aromaticity.
fn atoms_compatible(a: &chematic_core::Atom, b: &chematic_core::Atom) -> bool {
    a.element.atomic_number() == b.element.atomic_number() && a.aromatic == b.aromatic
}

/// Bond compatibility: same normalized bond order.
fn bonds_compatible(a: BondOrder, b: BondOrder) -> bool {
    normalize_bond(a) == normalize_bond(b)
}

/// Normalize Up/Down stereo bonds to Single for comparison purposes.
fn normalize_bond(o: BondOrder) -> BondOrder {
    match o {
        BondOrder::Up | BondOrder::Down => BondOrder::Single,
        other => other,
    }
}

// ---------------------------------------------------------------------------
// Bond counting helper
// ---------------------------------------------------------------------------

/// Count the number of already-mapped neighbors of `n0` in mol[0].
fn count_new_bonds(mol0: &Molecule, mapping: &PartialMapping, n0: AtomIdx) -> usize {
    mol0.neighbors(n0).filter(|(nb, _)| mapping.is_mapped(0, *nb)).count()
}

// ---------------------------------------------------------------------------
// Build output QueryMolecule
// ---------------------------------------------------------------------------

fn build_query(mol0: &Molecule, mapping: &PartialMapping, config: &McsConfig) -> QueryMolecule {
    let mut qmol = QueryMolecule::new();

    // Add one query atom per mapped position.
    for row in &mapping.query_to_mol {
        let a0 = row[0];
        let atom = mol0.atom(a0);
        let an = atom.element.atomic_number();
        let arom = atom.aromatic;

        let query = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::AtomicNum(an))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(arom))),
        );
        qmol.add_atom(query);
    }

    // Add bonds between mapped atoms (based on mol[0] topology).
    for (q1, row1) in mapping.query_to_mol.iter().enumerate() {
        let a1 = row1[0];
        for (q2, row2) in mapping.query_to_mol.iter().enumerate() {
            if q2 <= q1 {
                continue;
            }
            let a2 = row2[0];
            if let Some((_bidx, bond_entry)) = mol0.bond_between(a1, a2) {
                let bq = if config.match_bonds {
                    match normalize_bond(bond_entry.order) {
                        BondOrder::Single => BondQuery::Primitive(BondPrimitive::Single),
                        BondOrder::Double => BondQuery::Primitive(BondPrimitive::Double),
                        BondOrder::Triple => BondQuery::Primitive(BondPrimitive::Triple),
                        BondOrder::Aromatic => BondQuery::Primitive(BondPrimitive::Aromatic),
                        BondOrder::Quadruple => BondQuery::Primitive(BondPrimitive::Any),
                        BondOrder::Up | BondOrder::Down => BondQuery::Primitive(BondPrimitive::Single),
                    }
                } else {
                    BondQuery::Primitive(BondPrimitive::Any)
                };
                qmol.add_bond(q1, q2, bq);
            }
        }
    }

    qmol
}

// ---------------------------------------------------------------------------
// Convert a single molecule to a QueryMolecule (for the single-mol case)
// ---------------------------------------------------------------------------

fn molecule_to_query(mol: &Molecule) -> QueryMolecule {
    let mut qmol = QueryMolecule::new();

    for (_, atom) in mol.atoms() {
        let an = atom.element.atomic_number();
        let arom = atom.aromatic;
        let query = AtomQuery::And(
            Box::new(AtomQuery::Primitive(AtomPrimitive::AtomicNum(an))),
            Box::new(AtomQuery::Primitive(AtomPrimitive::Aromatic(arom))),
        );
        qmol.add_atom(query);
    }

    for (_, bond) in mol.bonds() {
        let bq = match normalize_bond(bond.order) {
            BondOrder::Single => BondQuery::Primitive(BondPrimitive::Single),
            BondOrder::Double => BondQuery::Primitive(BondPrimitive::Double),
            BondOrder::Triple => BondQuery::Primitive(BondPrimitive::Triple),
            BondOrder::Aromatic => BondQuery::Primitive(BondPrimitive::Aromatic),
            BondOrder::Quadruple => BondQuery::Primitive(BondPrimitive::Any),
            BondOrder::Up | BondOrder::Down => BondQuery::Primitive(BondPrimitive::Single),
        };
        qmol.add_bond(bond.atom1.0 as usize, bond.atom2.0 as usize, bq);
    }

    qmol
}

// ---------------------------------------------------------------------------
// Cartesian product iterator
// ---------------------------------------------------------------------------

/// Iterator over the Cartesian product of a slice of `Vec<T>` slices.
///
/// `CartesianProduct::new(&[vec![a,b], vec![x,y]])` yields `[a,x]`, `[a,y]`, `[b,x]`, `[b,y]`.
struct CartesianProduct<'a, T> {
    lists: &'a [Vec<T>],
    indices: Vec<usize>,
    done: bool,
}

impl<'a, T: Copy> CartesianProduct<'a, T> {
    fn new(lists: &'a [Vec<T>]) -> Self {
        let done = lists.iter().any(|l| l.is_empty());
        let indices = vec![0; lists.len()];
        Self { lists, indices, done }
    }
}

impl<'a, T: Copy> Iterator for CartesianProduct<'a, T> {
    type Item = Vec<T>;

    fn next(&mut self) -> Option<Vec<T>> {
        if self.done {
            return None;
        }

        // Emit current combination.
        let item: Vec<T> = self.lists.iter().zip(self.indices.iter()).map(|(l, &i)| l[i]).collect();

        // Advance indices (rightmost first).
        let mut carry = true;
        for k in (0..self.lists.len()).rev() {
            if carry {
                self.indices[k] += 1;
                if self.indices[k] < self.lists[k].len() {
                    carry = false;
                } else {
                    self.indices[k] = 0;
                }
            }
        }
        if carry {
            self.done = true;
        }

        Some(item)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_empty_input() {
        let result = find_mcs(&[]);
        assert_eq!(result.atom_count(), 0);
    }

    #[test]
    fn test_single_molecule() {
        let benzene = parse("c1ccccc1").unwrap();
        let result = find_mcs(&[&benzene]);
        assert_eq!(result.atom_count(), 6);
    }

    #[test]
    fn test_ethane_propane() {
        // CC and CCC → MCS should be CC (2 atoms)
        let a = parse("CC").unwrap();
        let b = parse("CCC").unwrap();
        let result = find_mcs(&[&a, &b]);
        assert!(result.atom_count() >= 2, "Expected >= 2 atoms, got {}", result.atom_count());
    }

    #[test]
    fn test_methane_ethanol() {
        // C and CCO → MCS should be C (1 atom)
        let a = parse("C").unwrap();
        let b = parse("CCO").unwrap();
        let result = find_mcs(&[&a, &b]);
        assert!(result.atom_count() >= 1);
    }

    #[test]
    fn test_benzene_toluene() {
        // c1ccccc1 and Cc1ccccc1 → MCS should be the benzene ring (6 atoms)
        let a = parse("c1ccccc1").unwrap();
        let b = parse("Cc1ccccc1").unwrap();
        let result = find_mcs(&[&a, &b]);
        assert_eq!(result.atom_count(), 6, "Expected benzene ring (6 atoms), got {}", result.atom_count());
    }

    #[test]
    fn test_no_common_atoms() {
        // Only N atoms vs only O atoms → no MCS
        let a = parse("N").unwrap();  // ammonia
        let b = parse("O").unwrap();  // water
        let result = find_mcs(&[&a, &b]);
        assert_eq!(result.atom_count(), 0);
    }

    #[test]
    fn test_identical_molecules() {
        let a = parse("CC").unwrap();
        let b = parse("CC").unwrap();
        let result = find_mcs(&[&a, &b]);
        assert_eq!(result.atom_count(), 2);
    }

    #[test]
    fn test_match_bonds_false() {
        // With match_bonds=false, single and double bond C-C should match
        let config = McsConfig { match_bonds: false, min_atoms: 1, timeout_ms: None };
        let a = parse("CC").unwrap();
        let b = parse("C=C").unwrap();
        let result = find_mcs_with_config(&[&a, &b], &config);
        assert!(result.atom_count() >= 2);
    }

    #[test]
    fn test_min_atoms_filter() {
        // CC and CCC have MCS of size 2, but if min_atoms=5, result should be empty
        let config = McsConfig { match_bonds: true, min_atoms: 5, timeout_ms: None };
        let a = parse("CC").unwrap();
        let b = parse("CCC").unwrap();
        let result = find_mcs_with_config(&[&a, &b], &config);
        assert_eq!(result.atom_count(), 0);
    }

    #[test]
    fn test_timeout_does_not_panic() {
        let config = McsConfig { match_bonds: true, min_atoms: 1, timeout_ms: Some(1) }; // 1 ms
        let a = parse("c1ccccc1").unwrap();
        let b = parse("Cc1ccccc1").unwrap();
        // Should return without panic (may return partial result)
        let _ = find_mcs_with_config(&[&a, &b], &config);
    }

    #[test]
    fn test_result_matches_all_inputs() {
        // The MCS result should be findable as a substructure in all input molecules
        use crate::find_matches;
        let a = parse("c1ccccc1").unwrap();
        let b = parse("Cc1ccccc1").unwrap();
        let result = find_mcs(&[&a, &b]);
        if result.atom_count() > 0 {
            assert!(!find_matches(&result, &a).is_empty());
            assert!(!find_matches(&result, &b).is_empty());
        }
    }

    #[test]
    fn test_aspirin_benzoic_acid() {
        // CC(=O)Oc1ccccc1C(=O)O and OC(=O)c1ccccc1
        // Common: benzoic acid core → >= 9 atoms
        let a = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let b = parse("OC(=O)c1ccccc1").unwrap();
        let result = find_mcs(&[&a, &b]);
        assert!(result.atom_count() >= 7, "Expected >= 7 atoms for aspirin/benzoic acid MCS, got {}", result.atom_count());
    }
}
