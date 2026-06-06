//! Canonical SMILES generation via the Morgan (extended connectivity) algorithm.
//!
//! A canonical SMILES is a unique string representation of a molecule:
//! two molecules that are identical (same graph, same atom properties) will
//! always produce the same canonical SMILES string.
//!
//! Algorithm:
//! 1. Assign initial invariants to each atom (atomic number, degree, charge, …).
//! 2. Iteratively update ranks using Morgan-style neighbor aggregation until
//!    the number of distinct ranks stabilises.
//! 3. Use the resulting ranks to impose a canonical DFS traversal order.
//!    Critically, both the ring-closure discovery DFS and the write DFS
//!    use the *same* canonical traversal order so the output is stable.
//! 4. Tie-breaking when two atoms have equal Morgan rank is resolved by
//!    atomic_number → isotope → charge → aromaticity → degree.
//!
//! Reference: Weininger, D. (1988) J. Chem. Inf. Comput. Sci. 28, 31-36.

use std::collections::{HashMap, HashSet};

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};

/// Return the atom indices sorted into canonical (Morgan-rank) order.
///
/// The returned `Vec<usize>` lists atom positions (0-based) in the order they
/// would be encountered during a canonical DFS write.  Atoms with higher
/// Morgan rank appear earlier.  This is the same ordering `canonical_smiles`
/// uses internally.
///
/// Useful for normalizing atom-indexed property arrays to a canonical order.
pub fn canonical_atom_order(mol: &Molecule) -> Vec<usize> {
    let ranks = morgan_ranks(mol);
    let mut order: Vec<usize> = (0..mol.atom_count()).collect();
    // Sort descending by rank (highest rank first, as in canonical DFS).
    order.sort_unstable_by(|&a, &b| ranks[b].cmp(&ranks[a]));
    order
}

/// Return `true` if atoms `a` and `b` are topologically equivalent (symmetric).
///
/// Two atoms are considered equivalent when they have the same Morgan rank —
/// meaning no graph-based feature (element, charge, degree, neighbour
/// environment, …) can distinguish them.
///
/// # Example
/// All six carbons of benzene are equivalent; the two carbons of ethane are
/// equivalent; the two oxygens of acetic acid are **not** (different degree).
/// Assign a symmetry class number to every atom.
///
/// Atoms with the same class number are topologically equivalent (symmetric).
/// Class numbers are consecutive integers starting at 0, ordered by increasing
/// Morgan rank (lowest rank = class 0).
///
/// # Example
/// Benzene returns `[0,0,0,0,0,0]` (all 6 carbons equivalent).
/// Toluene returns `[0,1,1,1,1,1,2]` (methyl-C, ring-Cs, ipso-C).
pub fn equivalent_atom_classes(mol: &Molecule) -> Vec<usize> {
    let ranks = morgan_ranks(mol);
    // Sort unique rank values to assign stable class numbers.
    let mut unique: Vec<u64> = ranks.clone();
    unique.sort_unstable();
    unique.dedup();
    ranks.iter().map(|r| unique.partition_point(|&u| u < *r)).collect()
}

pub fn are_atoms_equivalent(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> bool {
    let ranks = morgan_ranks(mol);
    let ia = a.0 as usize;
    let ib = b.0 as usize;
    if ia >= ranks.len() || ib >= ranks.len() {
        return false;
    }
    ranks[ia] == ranks[ib]
}

/// Return the canonical SMILES for a molecule.
///
/// For molecules with no atoms, returns an empty string.
/// Disconnected fragments (multiple components) are joined with `.`.
pub fn canonical_smiles(mol: &Molecule) -> String {
    if mol.atom_count() == 0 {
        return String::new();
    }

    let ranks = morgan_ranks(mol);
    CanonicalWriter::new(mol, &ranks).write_all()
}

/// Compute Morgan (extended connectivity) ranks for all atoms.
///
/// Returns a vector of normalised ordinal ranks (0-based, gap-free)
/// indexed by atom position (same order as `mol.atoms()`).
pub fn morgan_ranks(mol: &Molecule) -> Vec<u64> {
    let n = mol.atom_count();

    let mut ranks: Vec<u64> = (0..n)
        .map(|i| initial_invariant(mol, AtomIdx(i as u32)))
        .collect();

    let max_iter = n + 2;
    for _ in 0..max_iter {
        let old_distinct = count_distinct(&ranks);

        let new_ranks: Vec<u64> = (0..n)
            .map(|i| {
                let idx = AtomIdx(i as u32);
                let mut neighbor_ranks: Vec<u64> = mol
                    .neighbors(idx)
                    .map(|(nb, _)| ranks[nb.0 as usize])
                    .collect();
                neighbor_ranks.sort_unstable();
                fnv_hash_sequence(ranks[i], &neighbor_ranks)
            })
            .collect();

        let new_distinct = count_distinct(&new_ranks);
        ranks = new_ranks;

        if new_distinct <= old_distinct {
            break;
        }
    }

    normalize_ranks(&ranks)
}

/// Initial per-atom invariant packed into a u64.
fn initial_invariant(mol: &Molecule, idx: AtomIdx) -> u64 {
    let atom = mol.atom(idx);

    if atom.wildcard {
        return 0;
    }

    let an     = atom.element.atomic_number() as u64;
    let degree = mol.degree(idx) as u64;
    let charge = (atom.charge as i64 + 128) as u64;
    let iso    = atom.isotope.unwrap_or(0) as u64;
    let arom   = atom.aromatic as u64;
    let h_flag = atom.hydrogen_count.map(|h| h as u64 + 1).unwrap_or(0);

    (an     << 56)
        | (degree << 48)
        | (charge << 40)
        | (iso    << 24)
        | (h_flag << 16)
        | (arom   <<  8)
}

fn fnv_hash_sequence(base: u64, values: &[u64]) -> u64 {
    const FNV_PRIME:  u64 = 0x0000_0100_0000_01B3;
    const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    let mut h = FNV_OFFSET ^ base.wrapping_mul(FNV_PRIME);
    for &v in values {
        h ^= v;
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
}

fn count_distinct(ranks: &[u64]) -> usize {
    let mut seen: Vec<u64> = ranks.to_vec();
    seen.sort_unstable();
    seen.dedup();
    seen.len()
}

fn normalize_ranks(ranks: &[u64]) -> Vec<u64> {
    let mut sorted: Vec<(u64, usize)> = ranks
        .iter()
        .copied()
        .enumerate()
        .map(|(i, v)| (v, i))
        .collect();
    sorted.sort_unstable_by_key(|&(v, _)| v);

    let mut result = vec![0u64; ranks.len()];
    let mut current_rank: u64 = 0;
    let mut prev_val = sorted[0].0;

    for (val, idx) in sorted {
        if val != prev_val {
            current_rank += 1;
            prev_val = val;
        }
        result[idx] = current_rank;
    }

    result
}

struct CanonicalWriter<'a> {
    mol: &'a Molecule,
    ranks: &'a [u64],
    written: Vec<bool>,
    ring_bonds: HashSet<BondIdx>,
    atom_ring_nums: HashMap<AtomIdx, Vec<(u8, BondOrder)>>,
    next_ring: u8,
    out: String,
}

impl<'a> CanonicalWriter<'a> {
    fn new(mol: &'a Molecule, ranks: &'a [u64]) -> Self {
        let n = mol.atom_count();
        Self {
            mol,
            ranks,
            written: vec![false; n],
            ring_bonds: HashSet::new(),
            atom_ring_nums: HashMap::new(),
            next_ring: 1,
            out: String::new(),
        }
    }

    fn write_all(mut self) -> String {
        // Phase 1: discover ring-closure back-edges using the SAME canonical DFS
        // order that the writer will use. This ensures ring-closure numbers are
        // stable across re-parses.
        self.find_ring_closures();

        // Phase 2: canonical DFS serialization.
        let starts = self.canonical_atom_list();
        let mut first = true;
        for start in starts {
            if self.written[start.0 as usize] { continue; }
            if !first { self.out.push('.'); }
            first = false;
            self.write_chain(start, None, None);
        }

        self.out
    }

    /// Return all atoms sorted in canonical order: highest rank first, ties
    /// broken by chemical properties invariant across re-parses.
    fn canonical_atom_list(&self) -> Vec<AtomIdx> {
        let mut atoms: Vec<AtomIdx> = (0..self.mol.atom_count())
            .map(|i| AtomIdx(i as u32))
            .collect();
        atoms.sort_by(|&a, &b| self.canonical_cmp(b, a)); // descending
        atoms
    }

    /// Canonical ordering comparator (ascending; negate for descending).
    /// Tie-breaking uses chemical properties only (not atom indices),
    /// so the order is invariant between runs on chemically identical molecules.
    fn canonical_cmp(&self, a: AtomIdx, b: AtomIdx) -> std::cmp::Ordering {
        let ra = self.ranks[a.0 as usize];
        let rb = self.ranks[b.0 as usize];
        if ra != rb { return ra.cmp(&rb); }

        let atom_a = self.mol.atom(a);
        let atom_b = self.mol.atom(b);

        // Break ties with: atomic_number → isotope → charge → aromatic → degree
        atom_a.element.atomic_number().cmp(&atom_b.element.atomic_number())
            .then_with(|| atom_a.isotope.unwrap_or(0).cmp(&atom_b.isotope.unwrap_or(0)))
            .then_with(|| atom_a.charge.cmp(&atom_b.charge))
            .then_with(|| (atom_a.aromatic as u8).cmp(&(atom_b.aromatic as u8)))
            .then_with(|| self.mol.degree(a).cmp(&self.mol.degree(b)))
    }

    /// Discover back-edges by running the same canonical DFS as the writer.
    /// Using identical traversal order ensures ring-closure numbers are stable.
    fn find_ring_closures(&mut self) {
        let n = self.mol.atom_count();
        let mut visited  = vec![false; n];
        let mut in_stack = vec![false; n];

        // Iterate in canonical order (same as write_all).
        let starts = self.canonical_atom_list();
        for start in starts {
            if !visited[start.0 as usize] {
                self.dfs_mark(start, None, &mut visited, &mut in_stack);
            }
        }
    }

    fn dfs_mark(
        &mut self,
        atom: AtomIdx,
        from_bond: Option<BondIdx>,
        visited: &mut Vec<bool>,
        in_stack: &mut Vec<bool>,
    ) {
        visited[atom.0 as usize] = true;
        in_stack[atom.0 as usize] = true;

        let mut neighbors: Vec<(AtomIdx, BondIdx)> = self.mol.neighbors(atom).collect();
        self.sort_neighbors_canonical(&mut neighbors);

        for (neighbor, bidx) in neighbors {
            if Some(bidx) == from_bond { continue; }
            if self.ring_bonds.contains(&bidx) { continue; }

            if !visited[neighbor.0 as usize] {
                self.dfs_mark(neighbor, Some(bidx), visited, in_stack);
            } else if in_stack[neighbor.0 as usize] {
                self.ring_bonds.insert(bidx);
                let rn = self.next_ring;
                self.next_ring += 1;
                let bond = self.mol.bond(bidx);
                // Direction seen from `neighbor` (the open atom) going toward `atom`.
                let order_at_open = match bond.order {
                    BondOrder::Up   => if bond.atom1 == neighbor { BondOrder::Up }   else { BondOrder::Down },
                    BondOrder::Down => if bond.atom1 == neighbor { BondOrder::Down } else { BondOrder::Up },
                    other => other,
                };
                // Suppress stereo at the close atom to avoid conflicting ring-closure chars.
                let order_at_close = match bond.order {
                    BondOrder::Up | BondOrder::Down => BondOrder::Single,
                    other => other,
                };
                self.atom_ring_nums.entry(neighbor).or_default().push((rn, order_at_open));
                self.atom_ring_nums.entry(atom).or_default().push((rn, order_at_close));
            }
        }

        in_stack[atom.0 as usize] = false;
    }

    fn write_chain(
        &mut self,
        atom: AtomIdx,
        from_atom: Option<AtomIdx>,
        incoming_bond: Option<BondOrder>,
    ) {
        self.written[atom.0 as usize] = true;

        if let Some(bond) = incoming_bond {
            self.out.push(bond.smiles_char());
        }

        self.emit_atom(atom);

        // Ring-closure digits.
        if let Some(rings) = self.atom_ring_nums.remove(&atom) {
            for (rn, bond_order) in rings {
                let atom_arom = self.mol.atom(atom).aromatic;
                if !(bond_order == BondOrder::Aromatic && atom_arom)
                    && bond_order != BondOrder::Single
                {
                    self.out.push(bond_order.smiles_char());
                }
                if rn >= 10 {
                    self.out.push('%');
                    self.out.push(char::from_digit((rn / 10) as u32, 10).unwrap());
                    self.out.push(char::from_digit((rn % 10) as u32, 10).unwrap());
                } else {
                    self.out.push(char::from_digit(rn as u32, 10).unwrap());
                }
            }
        }

        // Tree-edge children, sorted canonically.
        let mut children: Vec<(AtomIdx, BondOrder)> = self.mol
            .neighbors(atom)
            .filter(|(nb, bidx)| {
                Some(*nb) != from_atom
                    && !self.written[nb.0 as usize]
                    && !self.ring_bonds.contains(bidx)
            })
            .map(|(nb, bidx)| {
                let bond = self.mol.bond(bidx);
                // Direction seen from `atom` going toward `nb`.
                let order = match bond.order {
                    BondOrder::Up   => if bond.atom1 == atom { BondOrder::Up }   else { BondOrder::Down },
                    BondOrder::Down => if bond.atom1 == atom { BondOrder::Down } else { BondOrder::Up },
                    other => other,
                };
                (nb, order)
            })
            .collect();

        // Sort children by canonical rank (ascending → highest rank = main chain).
        children.sort_by(|&(a, _), &(b, _)| self.canonical_cmp(a, b));

        let n = children.len();
        for (i, (child, bond_order)) in children.into_iter().enumerate() {
            let is_last = i == n - 1;
            let parent_arom = self.mol.atom(atom).aromatic;
            let child_arom  = self.mol.atom(child).aromatic;
            let implicit = match bond_order {
                BondOrder::Single   => !(parent_arom && child_arom),
                BondOrder::Aromatic => parent_arom && child_arom,
                _ => false,
            };
            let written_bond = if implicit { None } else { Some(bond_order) };

            if !is_last {
                self.out.push('(');
                self.write_chain(child, Some(atom), written_bond);
                self.out.push(')');
            } else {
                self.write_chain(child, Some(atom), written_bond);
            }
        }
    }

    /// Sort a neighbor list in canonical order (for consistent DFS traversal).
    fn sort_neighbors_canonical(&self, neighbors: &mut Vec<(AtomIdx, BondIdx)>) {
        neighbors.sort_by(|&(a, _), &(b, _)| self.canonical_cmp(b, a)); // descending
    }

    fn emit_atom(&mut self, idx: AtomIdx) {
        let atom = self.mol.atom(idx);

        if atom.wildcard {
            self.out.push_str("[*]");
            return;
        }

        let needs_bracket = atom.isotope.is_some()
            || atom.charge != 0
            || atom.hydrogen_count.is_some()
            || !atom.element.is_organic_subset()
            || atom.atom_map.is_some();

        if needs_bracket {
            self.out.push('[');
            if let Some(iso) = atom.isotope {
                self.out.push_str(&iso.to_string());
            }
            let sym = if atom.aromatic {
                atom.element.symbol().to_lowercase()
            } else {
                atom.element.symbol().to_string()
            };
            self.out.push_str(&sym);

            match atom.chirality {
                chematic_core::Chirality::CounterClockwise => self.out.push('@'),
                chematic_core::Chirality::Clockwise        => self.out.push_str("@@"),
                chematic_core::Chirality::None             => {}
            }

            if let Some(h) = atom.hydrogen_count {
                if h > 0 {
                    self.out.push('H');
                    if h > 1 { self.out.push_str(&h.to_string()); }
                }
            }

            match atom.charge {
                0 => {}
                1 => self.out.push('+'),
                -1 => self.out.push('-'),
                c if c > 0 => self.out.push_str(&format!("+{c}")),
                c => self.out.push_str(&c.to_string()),
            }

            if let Some(m) = atom.atom_map {
                self.out.push(':');
                self.out.push_str(&m.to_string());
            }

            self.out.push(']');
        } else if atom.aromatic {
            self.out.push_str(&atom.element.symbol().to_lowercase());
        } else {
            self.out.push_str(atom.element.symbol());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    /// Canonical SMILES must be stable: applying it twice gives the same result.
    fn is_stable(smiles: &str) -> bool {
        let mol1 = parse(smiles).expect(smiles);
        let c1 = canonical_smiles(&mol1);
        assert!(!c1.is_empty(), "canonical_smiles returned empty for '{smiles}'");
        let mol2 = parse(&c1).unwrap_or_else(|e| {
            panic!("canonical SMILES '{c1}' is not parseable: {e}")
        });
        let c2 = canonical_smiles(&mol2);
        c1 == c2
    }

    /// Two SMILES representing the same molecule must give the same canonical form.
    fn same_canonical(a: &str, b: &str) -> bool {
        let mol_a = parse(a).expect(a);
        let mol_b = parse(b).expect(b);
        canonical_smiles(&mol_a) == canonical_smiles(&mol_b)
    }

    #[test]
    fn test_methane_stable()      { assert!(is_stable("C")); }
    #[test]
    fn test_ethane_stable()       { assert!(is_stable("CC")); }
    #[test]
    fn test_ethanol_stable()      { assert!(is_stable("CCO")); }
    #[test]
    fn test_acetic_acid_stable()  { assert!(is_stable("CC(=O)O")); }
    #[test]
    fn test_benzene_stable()      { assert!(is_stable("c1ccccc1")); }
    #[test]
    fn test_pyridine_stable()     { assert!(is_stable("c1ccncc1")); }
    #[test]
    fn test_naphthalene_stable()  { assert!(is_stable("c1ccc2ccccc2c1")); }
    #[test]
    fn test_aspirin_stable()      { assert!(is_stable("CC(=O)Oc1ccccc1C(=O)O")); }
    #[test]
    fn test_caffeine_stable()     { assert!(is_stable("Cn1cnc2c1c(=O)n(c(=O)n2C)C")); }

    #[test]
    fn test_ethanol_same_from_different_starts() {
        assert!(same_canonical("CCO", "OCC"));
    }

    #[test]
    fn test_isobutane_same_canonical() {
        // CC(C)C and C(C)(C)C are the same molecule.
        assert!(same_canonical("CC(C)C", "C(C)(C)C"));
    }

    #[test]
    fn test_wildcard_roundtrip() {
        let mol = parse("[*]CC").unwrap();
        let c = canonical_smiles(&mol);
        assert!(!c.is_empty());
        let mol2 = parse(&c).unwrap();
        assert_eq!(mol.atom_count(), mol2.atom_count());
        assert!(is_stable("[*]CC"));
    }

    #[test]
    fn test_disconnected_stable() {
        assert!(is_stable("[Na+].[Cl-]"));
    }

    // E/Z stereo bond direction tests.
    #[test]
    fn test_ez_e_stable() { assert!(is_stable("C/C=C/C")); }
    #[test]
    fn test_ez_z_stable() { assert!(is_stable("C/C=C\\C")); }
    #[test]
    fn test_ez_fluoro_e_stable() { assert!(is_stable("F/C=C/Cl")); }
    #[test]
    fn test_ez_fluoro_z_stable() { assert!(is_stable("F/C=C\\Cl")); }
    #[test]
    fn test_ez_e_ne_z() {
        // E and Z isomers of 1-fluoro-2-chloroethylene must yield different canonical forms.
        let mol_e = parse("F/C=C/Cl").unwrap();
        let mol_z = parse("F/C=C\\Cl").unwrap();
        assert_ne!(canonical_smiles(&mol_e), canonical_smiles(&mol_z));
    }
}
