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

use chematic_core::{AtomIdx, BondIdx, BondOrder, Chirality, Molecule, STEREO_H_SENTINEL};

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
    ranks
        .iter()
        .map(|r| unique.partition_point(|&u| u < *r))
        .collect()
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
                // Include bond order in the neighbor contribution so that atoms
                // bonded via different bond types (e.g. O= vs O-H in acetic acid)
                // receive distinct Morgan ranks even when neighbor atom ranks are
                // otherwise identical.
                let mut neighbor_contributions: Vec<u64> = mol
                    .neighbors(idx)
                    .map(|(nb, bidx)| {
                        let bond_val = bond_order_value(mol.bond(bidx).order);
                        fnv_hash_sequence(ranks[nb.0 as usize], &[bond_val])
                    })
                    .collect();
                neighbor_contributions.sort_unstable();
                fnv_hash_sequence(ranks[i], &neighbor_contributions)
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

    let an = atom.element.atomic_number() as u64;
    let degree = mol.degree(idx) as u64;
    let charge = (atom.charge as i64 + 128) as u64;
    let iso = atom.isotope.unwrap_or(0) as u64;
    let arom = atom.aromatic as u64;
    let h_flag = atom.hydrogen_count.map(|h| h as u64 + 1).unwrap_or(0);

    (an << 56) | (degree << 48) | (charge << 40) | (iso << 24) | (h_flag << 16) | (arom << 8)
}

/// Map a BondOrder to a stable integer for use in Morgan rank hashing.
fn bond_order_value(order: BondOrder) -> u64 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
        BondOrder::Quadruple => 5,
        _ => 0,
    }
}

fn fnv_hash_sequence(base: u64, values: &[u64]) -> u64 {
    const FNV_PRIME: u64 = 0x0000_0100_0000_01B3;
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
    /// (ring_num, bond_order, ring_partner_atom)
    atom_ring_nums: HashMap<AtomIdx, Vec<(u32, BondOrder, AtomIdx)>>,
    next_ring: u32,
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
            if self.written[start.0 as usize] {
                continue;
            }
            if !first {
                self.out.push('.');
            }
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
        if ra != rb {
            return ra.cmp(&rb);
        }

        let atom_a = self.mol.atom(a);
        let atom_b = self.mol.atom(b);

        // Break ties with: atomic_number → isotope → charge → aromatic → degree
        atom_a
            .element
            .atomic_number()
            .cmp(&atom_b.element.atomic_number())
            .then_with(|| {
                atom_a
                    .isotope
                    .unwrap_or(0)
                    .cmp(&atom_b.isotope.unwrap_or(0))
            })
            .then_with(|| atom_a.charge.cmp(&atom_b.charge))
            .then_with(|| (atom_a.aromatic as u8).cmp(&(atom_b.aromatic as u8)))
            .then_with(|| self.mol.degree(a).cmp(&self.mol.degree(b)))
    }

    /// Discover back-edges by running the same canonical DFS as the writer.
    /// Using identical traversal order ensures ring-closure numbers are stable.
    fn find_ring_closures(&mut self) {
        let n = self.mol.atom_count();
        let mut visited = vec![false; n];
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
            if Some(bidx) == from_bond {
                continue;
            }
            if self.ring_bonds.contains(&bidx) {
                continue;
            }

            if !visited[neighbor.0 as usize] {
                self.dfs_mark(neighbor, Some(bidx), visited, in_stack);
            } else if in_stack[neighbor.0 as usize] {
                self.ring_bonds.insert(bidx);
                let rn = self.next_ring;
                self.next_ring += 1;
                let bond = self.mol.bond(bidx);
                // A ring bond forced to Aromatic (e.g. adjacent to an
                // exocyclic C=N) may carry its true E/Z direction stashed
                // separately rather than in `order` itself — consult that
                // first so the direction still reaches the writer.
                let stashed_direction = self.mol.bond_direction(bidx);
                let effective_order = stashed_direction.unwrap_or(bond.order);
                // Direction seen from `neighbor` (the open atom) going toward `atom`.
                let order_at_open = match effective_order {
                    BondOrder::Up => {
                        if bond.atom1 == neighbor {
                            BondOrder::Up
                        } else {
                            BondOrder::Down
                        }
                    }
                    BondOrder::Down => {
                        if bond.atom1 == neighbor {
                            BondOrder::Down
                        } else {
                            BondOrder::Up
                        }
                    }
                    other => other,
                };
                // Suppress stereo at the close atom to avoid conflicting ring-closure
                // chars, falling back to the bond's real order — Aromatic for a
                // stashed direction (implicit ring bond, no char), Single for a
                // genuine directional single bond (existing behavior).
                let order_at_close = match effective_order {
                    BondOrder::Up | BondOrder::Down => {
                        if stashed_direction.is_some() {
                            bond.order
                        } else {
                            BondOrder::Single
                        }
                    }
                    other => other,
                };
                self.atom_ring_nums
                    .entry(neighbor)
                    .or_default()
                    .push((rn, order_at_open, atom)); // partner = close atom
                self.atom_ring_nums
                    .entry(atom)
                    .or_default()
                    .push((rn, order_at_close, neighbor)); // partner = open atom
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

        // Compute parity-corrected chirality before ring data is consumed.
        let corrected_chirality = self.corrected_chirality(atom, from_atom);
        self.emit_atom(atom, corrected_chirality);

        // Ring-closure digits.
        if let Some(rings) = self.atom_ring_nums.remove(&atom) {
            for (rn, bond_order, _partner) in rings {
                let atom_arom = self.mol.atom(atom).aromatic;
                if !(bond_order == BondOrder::Aromatic && atom_arom)
                    && bond_order != BondOrder::Single
                {
                    self.out.push(bond_order.smiles_char());
                }
                // SMILES ring-closure numbers are limited to 1–99.
                // Molecules needing ≥ 100 simultaneous open ring closures are
                // exotic beyond any known organic chemistry; skip extras rather
                // than panic from `char::from_digit` overflow.
                if rn > 99 {
                    continue;
                }
                if rn >= 10 {
                    self.out.push('%');
                    self.out.push(char::from_digit(rn / 10, 10).unwrap());
                    self.out.push(char::from_digit(rn % 10, 10).unwrap());
                } else {
                    self.out.push(char::from_digit(rn, 10).unwrap());
                }
            }
        }

        // Tree-edge children, sorted canonically.
        let mut children: Vec<(AtomIdx, BondOrder)> = self
            .mol
            .neighbors(atom)
            .filter(|(nb, bidx)| {
                Some(*nb) != from_atom
                    && !self.written[nb.0 as usize]
                    && !self.ring_bonds.contains(bidx)
            })
            .map(|(nb, bidx)| {
                let bond = self.mol.bond(bidx);
                // See the ring-closure site above: a stashed direction takes
                // priority over `order` itself (e.g. an aromatic ring bond
                // that flanks an exocyclic C=N).
                let effective_order = self.mol.bond_direction(bidx).unwrap_or(bond.order);
                // Direction seen from `atom` going toward `nb`.
                let order = match effective_order {
                    BondOrder::Up => {
                        if bond.atom1 == atom {
                            BondOrder::Up
                        } else {
                            BondOrder::Down
                        }
                    }
                    BondOrder::Down => {
                        if bond.atom1 == atom {
                            BondOrder::Down
                        } else {
                            BondOrder::Up
                        }
                    }
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
            let child_arom = self.mol.atom(child).aromatic;
            let implicit = match bond_order {
                BondOrder::Single => !(parent_arom && child_arom),
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
    fn sort_neighbors_canonical(&self, neighbors: &mut [(AtomIdx, BondIdx)]) {
        neighbors.sort_by(|&(a, _), &(b, _)| self.canonical_cmp(b, a)); // descending
    }

    fn emit_atom(&mut self, idx: AtomIdx, chirality: Chirality) {
        let atom = self.mol.atom(idx);

        if atom.wildcard {
            self.out.push_str("[*]");
            return;
        }

        let needs_bracket = atom.isotope.is_some()
            || atom.charge != 0
            || atom.hydrogen_count.is_some()
            || !atom.element.is_organic_subset()
            || atom.atom_map.is_some()
            || chirality != Chirality::None;

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

            match chirality {
                Chirality::CounterClockwise => self.out.push('@'),
                Chirality::Clockwise => self.out.push_str("@@"),
                Chirality::None => {}
            }

            if let Some(h) = atom.hydrogen_count
                && h > 0
            {
                self.out.push('H');
                if h > 1 {
                    self.out.push_str(&h.to_string());
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

    /// Compute the parity-corrected chirality for `atom` when it is written
    /// with `from_atom` as the predecessor in the canonical DFS.
    ///
    /// Returns the stored chirality unchanged when no stereo neighbor order is
    /// recorded (e.g. programmatically constructed molecules).
    fn corrected_chirality(&self, atom: AtomIdx, from_atom: Option<AtomIdx>) -> Chirality {
        let stored = self.mol.atom(atom).chirality;
        if stored == Chirality::None {
            return Chirality::None;
        }

        let Some(original) = self.mol.stereo_neighbor_order(atom) else {
            return stored; // no parse-time data → return as-is
        };

        let atom_data = self.mol.atom(atom);
        let has_h = atom_data.hydrogen_count.is_some_and(|h| h > 0);

        // Build canonical neighbor sequence in SMILES output order:
        // 1. from_atom   (or H_SENTINEL if root and has bracket H)
        // 2. bracket H   (only when from_atom is Some and has_h)
        // 3. ring-closure partners in ring-number order
        // 4. children in ascending canonical rank (branches first, main chain last)
        let mut canonical: Vec<u32> = Vec::with_capacity(original.len());

        match from_atom {
            Some(prev) => {
                canonical.push(prev.0);
                if has_h {
                    canonical.push(STEREO_H_SENTINEL);
                }
            }
            None => {
                if has_h {
                    canonical.push(STEREO_H_SENTINEL);
                }
            }
        }

        if let Some(rings) = self.atom_ring_nums.get(&atom) {
            for &(_, _, partner) in rings {
                canonical.push(partner.0);
            }
        }

        let mut children: Vec<AtomIdx> = self
            .mol
            .neighbors(atom)
            .filter(|(nb, bidx)| {
                Some(*nb) != from_atom
                    && !self.written[nb.0 as usize]
                    && !self.ring_bonds.contains(bidx)
            })
            .map(|(nb, _)| nb)
            .collect();
        children.sort_by(|&a, &b| self.canonical_cmp(a, b)); // ascending rank
        for child in children {
            canonical.push(child.0);
        }

        if canonical.len() != original.len() {
            return stored; // size mismatch → fallback
        }

        if permutation_is_odd(original, &canonical) {
            match stored {
                Chirality::CounterClockwise => Chirality::Clockwise,
                Chirality::Clockwise => Chirality::CounterClockwise,
                Chirality::None => Chirality::None,
            }
        } else {
            stored
        }
    }
}

/// Return `true` if the permutation mapping `original` order to `canonical` order
/// has odd parity (i.e. requires an odd number of transpositions).
///
/// Both slices must contain the same multiset of `u32` values.
fn permutation_is_odd(original: &[u32], canonical: &[u32]) -> bool {
    let n = original.len();
    let mut pos: HashMap<u32, usize> = HashMap::with_capacity(n);
    for (i, &v) in original.iter().enumerate() {
        pos.insert(v, i);
    }
    // perm[i] = position in `original` of the element at `canonical[i]`
    let perm: Vec<usize> = canonical
        .iter()
        .map(|v| *pos.get(v).unwrap_or(&0))
        .collect();

    // Count cycles in the permutation; parity = (n - #cycles) % 2
    let mut visited = vec![false; n];
    let mut num_cycles = 0usize;
    for start in 0..n {
        if !visited[start] {
            num_cycles += 1;
            let mut j = start;
            while !visited[j] {
                visited[j] = true;
                j = perm[j];
            }
        }
    }
    (n - num_cycles) % 2 == 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    /// Canonical SMILES must be stable: applying it twice gives the same result.
    fn is_stable(smiles: &str) -> bool {
        let mol1 = parse(smiles).expect(smiles);
        let c1 = canonical_smiles(&mol1);
        assert!(
            !c1.is_empty(),
            "canonical_smiles returned empty for '{smiles}'"
        );
        let mol2 =
            parse(&c1).unwrap_or_else(|e| panic!("canonical SMILES '{c1}' is not parseable: {e}"));
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
    fn test_methane_stable() {
        assert!(is_stable("C"));
    }
    #[test]
    fn test_ethane_stable() {
        assert!(is_stable("CC"));
    }
    #[test]
    fn test_ethanol_stable() {
        assert!(is_stable("CCO"));
    }
    #[test]
    fn test_acetic_acid_stable() {
        assert!(is_stable("CC(=O)O"));
    }
    #[test]
    fn test_benzene_stable() {
        assert!(is_stable("c1ccccc1"));
    }
    #[test]
    fn test_pyridine_stable() {
        assert!(is_stable("c1ccncc1"));
    }
    #[test]
    fn test_naphthalene_stable() {
        assert!(is_stable("c1ccc2ccccc2c1"));
    }
    #[test]
    fn test_aspirin_stable() {
        assert!(is_stable("CC(=O)Oc1ccccc1C(=O)O"));
    }
    #[test]
    fn test_caffeine_stable() {
        assert!(is_stable("Cn1cnc2c1c(=O)n(c(=O)n2C)C"));
    }

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
    fn test_ez_e_stable() {
        assert!(is_stable("C/C=C/C"));
    }
    #[test]
    fn test_ez_z_stable() {
        assert!(is_stable("C/C=C\\C"));
    }
    #[test]
    fn test_ez_fluoro_e_stable() {
        assert!(is_stable("F/C=C/Cl"));
    }
    #[test]
    fn test_ez_fluoro_z_stable() {
        assert!(is_stable("F/C=C\\Cl"));
    }
    #[test]
    fn test_ez_e_ne_z() {
        // E and Z isomers of 1-fluoro-2-chloroethylene must yield different canonical forms.
        let mol_e = parse("F/C=C/Cl").unwrap();
        let mol_z = parse("F/C=C\\Cl").unwrap();
        assert_ne!(canonical_smiles(&mol_e), canonical_smiles(&mol_z));
    }

    // ── Tetrahedral stereo parity tests ─────────────────────────────────────

    #[test]
    fn test_tetrahedral_stable_no_from_atom() {
        // Bracket-H form at start of fragment — no from-atom.
        assert!(is_stable("[C@@H](F)(Cl)Br"));
        assert!(is_stable("[C@H](F)(Cl)Br"));
    }

    #[test]
    fn test_tetrahedral_stable_with_from_atom() {
        // L-alanine: chiral atom has a from-atom (N).
        assert!(is_stable("N[C@@H](C)C(=O)O"));
        assert!(is_stable("N[C@H](C)C(=O)O"));
    }

    #[test]
    fn test_enantiomers_differ() {
        // R and S configurations must give distinct canonical SMILES.
        assert!(!same_canonical("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"));
        assert!(!same_canonical("[C@@H](F)(Cl)Br", "[C@H](F)(Cl)Br"));
    }

    #[test]
    fn test_tetrahedral_same_from_different_starts() {
        // L-alanine from N vs from methyl — odd permutation, parity correction required.
        // RDKit: N[C@@H](C)C(=O)O and C[C@H](N)C(=O)O both → C[C@H](N)C(=O)O.
        assert!(same_canonical("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O"));
        // D-alanine must differ from L-alanine.
        assert!(!same_canonical("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"));
    }

    #[test]
    fn test_rdkit_agreement_alanine() {
        // Pairs where the Morgan ranks distinguish all atoms unambiguously.
        // N[C@@H](C)C(=O)O and C[C@H](N)C(=O)O: same L-alanine (RDKit agrees).
        assert!(same_canonical("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O"));
        // Enantiomers must differ (RDKit: C[C@@H](N)C(=O)O for D-alanine).
        assert!(!same_canonical("N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O"));
        // Stability: L-alanine canonical is self-stable.
        assert!(is_stable("N[C@@H](C)C(=O)O"));
        assert!(is_stable("C[C@H](N)C(=O)O"));
    }

    #[test]
    fn test_tetrahedral_all_heavy_substituents_stable() {
        // Chiral centre with no bracket H (all four heavy substituents).
        assert!(is_stable("[C@](F)(Cl)(Br)I"));
        assert!(is_stable("[C@@](F)(Cl)(Br)I"));
    }

    #[test]
    fn test_tetrahedral_all_heavy_enantiomers_differ() {
        assert!(!same_canonical("[C@](F)(Cl)(Br)I", "[C@@](F)(Cl)(Br)I"));
    }

    #[test]
    fn test_ring_stereocentre_stable() {
        // Chiral atom inside a ring — tests ring-closure partner resolution.
        assert!(is_stable("[C@@H]1CCCC1F"));
        assert!(is_stable("[C@H]1CCCC1F"));
    }

    #[test]
    fn test_ring_stereocentre_enantiomers_differ() {
        assert!(!same_canonical("[C@@H]1CCCC1F", "[C@H]1CCCC1F"));
    }

    #[test]
    fn test_chirality_from_different_entry_points() {
        // Same chiral molecule, two SMILES with different traversal order.
        // F[C@@H](Cl)Br  ≡  Cl[C@H](F)Br  (same S-configuration, just written
        // from different entry atoms — verified by signed-tetrahedral-volume).
        // Their canonical SMILES must be identical.
        let c1 = canonical_smiles(&parse("F[C@@H](Cl)Br").unwrap());
        let c2 = canonical_smiles(&parse("Cl[C@H](F)Br").unwrap());
        assert_eq!(c1, c2, "same molecule from different starts should match");

        // Cross-check: the enantiomer gives a different canonical form.
        let c3 = canonical_smiles(&parse("F[C@H](Cl)Br").unwrap());
        assert_ne!(c1, c3, "enantiomers must differ");
    }

    // ── Bond-order canonicality tests (#14 fix) ──────────────────────────

    #[test]
    fn test_acetic_acid_canonical_same_from_different_starts() {
        // Bug #14: both oxygens in acetic acid had the same Morgan rank because
        // the refinement loop omitted bond orders.  After the fix, O= (double)
        // and O-H (single) get distinct ranks regardless of atom insertion order.
        assert!(same_canonical("CC(=O)O", "OC(C)=O"));
        assert!(same_canonical("CC(=O)O", "O=C(O)C"));
        assert!(same_canonical("CC(=O)O", "C(C)(=O)O"));
    }

    #[test]
    fn test_oxygens_in_acetic_acid_not_equivalent() {
        // The two oxygens (O= vs O-H) are chemically distinct and must receive
        // different Morgan symmetry classes.
        let mol = parse("CC(=O)O").unwrap();
        let classes = equivalent_atom_classes(&mol);
        let o_classes: Vec<usize> = mol
            .atoms()
            .filter(|(_, a)| a.element.atomic_number() == 8)
            .map(|(i, _)| classes[i.0 as usize])
            .collect();
        assert_eq!(o_classes.len(), 2);
        assert_ne!(
            o_classes[0], o_classes[1],
            "O= and O-H must be in different symmetry classes"
        );
    }

    #[test]
    fn test_formic_acid_canonical_consistent() {
        // OC=O and O=CO — same formic acid, should canonicalize identically.
        assert!(same_canonical("OC=O", "O=CO"));
    }

    // ── RDKit PR #9066: conjugated E/Z round-trip ────────────────────────────

    #[test]
    fn conjugated_double_bond_ez_round_trip() {
        // RDKit PR #9066: removeRedundantBondDirSpecs() could strip bond directions
        // on conjugated double bonds, losing E/Z stereo.  Chematic does not apply
        // aggressive direction removal, but this test guards against regressions.
        for smi in &[
            r"F/C=C/C=C/Cl", // all-E conjugated diene
            r"F/C=C\C=C\Cl", // E then Z
            r"F/C=C/C=C\Cl", // E then inverted-Z
        ] {
            let mol = parse(smi).unwrap_or_else(|e| panic!("parse {smi}: {e:?}"));
            let out = canonical_smiles(&mol);
            let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse {out}: {e:?}"));
            let out2 = canonical_smiles(&mol2);
            assert_eq!(
                out, out2,
                "conjugated E/Z must be stable after two rounds: {smi} → {out} → {out2}"
            );
        }
    }

    // ── Round 10: ring-closure directional bond flip ─────────────────────────
    //
    // A directional marker (`/`, `\`) is read "toward" the ring digit from
    // wherever it's written. At the ring-OPENING occurrence that's already
    // the open->close direction; at the CLOSING occurrence it's close->open
    // (the opposite traversal direction over the same physical bond) and must
    // be flipped before use (parser.rs `close_or_open_ring`). Before the fix,
    // the closing-side marker was stored raw/unflipped, which silently
    // produced a *different stereoisomer* whenever a random SMILES spelling
    // routed a conjugated system's connecting single bond through a
    // ring-closure digit instead of a plain adjacent chain bond -- confirmed
    // via a corpus-wide worst-of-10 sweep (RDKit-checked structural
    // correctness, not just self-stability/idempotency, which this bug class
    // passed trivially since it was deterministic-but-wrong on each input).

    #[test]
    fn ring_closure_direction_flip_real_world_repro() {
        // Real molecule found via corpus sweep. `variant` is an RDKit
        // doRandom=True re-spelling of the exact same molecule as `orig`,
        // routing the diene's connecting single bond through ring-closure
        // digit "1" instead of a plain chain bond. Before the parser fix,
        // chematic silently emitted a different (RDKit-confirmed
        // non-equivalent) stereoisomer for `variant`.
        let orig = r"CC1CCOC(=O)/C=C/C=C\C(=O)O[C@@H]2C[C@H]3O[C@@H]4C[C@@H](C)C(=O)C[C@]4(COC(=O)C1O)[C@]2(C)C31CO1";
        let variant = r"C1=C\C(=O)O[C@@H]2C[C@H]3O[C@H]4[C@@]([C@@]2(C32OC2)C)(CC(=O)[C@H](C)C4)COC(=O)C(O)C(C)CCOC(=O)/C=C/1";
        assert!(
            same_canonical(orig, variant),
            "ring-closure-routed diene must canonicalize identically to the \
             chain-form spelling of the same molecule"
        );
    }

    #[test]
    fn ring_closure_direction_minimal_ez_agreement() {
        // Minimal case isolating the same mechanism: a ring-closure bond
        // (distinct from the exocyclic C=C double bond itself) whose
        // directional markers are specified at BOTH the opening and closing
        // occurrences of the ring digit. Per the flip rule, opposite raw
        // symbols (one `/`, one `\`) describe one consistent bond and must
        // parse successfully; same-symbol at both ends is the conflicting
        // case (unchanged by this fix -- only Up/Down are flipped, so a
        // same-vs-different Double/Single conflict, e.g. "C=1CC-1", is
        // unaffected).
        let mol = parse(r"F/C=C/1CCCC\1").unwrap_or_else(|e| panic!("{e:?}"));
        let out = canonical_smiles(&mol);
        let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse {out}: {e:?}"));
        assert_eq!(
            canonical_smiles(&mol2),
            out,
            "ring-closure E/Z with opposite-symbol agreement must round-trip stably"
        );

        // Same-symbol at both ends of a ring-closure directional bond is now
        // (correctly) the conflicting combination.
        assert!(matches!(
            parse(r"F/C=C/1CCCC/1"),
            Err(crate::error::SmilesError::ConflictingRingBond { ring_num: 1, .. })
        ));
    }

    // ── Round 10: ring-digit reuse racing PendingRing resolution ─────────────
    //
    // A stereocenter that OPENS a ring whose partner closes INSIDE the
    // stereocenter's own branch subtree (e.g. `[C@]1(...[C@H]1...)`) has its
    // own stereo record still unfinalized at the moment of that first
    // closure -- the immediate-resolution fast path in `close_or_open_ring`
    // only patches already-finalized records, so this case falls through to
    // the end-of-parse fallback. Before the fix, that fallback resolved by
    // raw ring DIGIT via `ring_close_partners: HashMap<u8, AtomIdx>` -- if the
    // same digit was reused later for an unrelated ring (e.g. a trailing
    // phenyl `c1ccccc1`), the later reuse's closer silently overwrote the
    // earlier, still-pending resolution, corrupting the stereocenter's
    // neighbor order with a foreign atom index (confirmed: the wrong index
    // pointed at an aromatic carbon in the unrelated trailing ring, not
    // anywhere near the stereocenter). Fixed by keying resolution on a
    // per-occurrence slot id (`next_ring_slot`) that is never reused,
    // regardless of how many times the same ring digit is.

    #[test]
    fn ring_digit_reuse_inside_stereocenter_branch_real_world_repro() {
        // Real molecule found via corpus sweep. `variant` is an RDKit
        // doRandom=True re-spelling of the exact same molecule as `orig`,
        // where the stereocenter's ring-1 partner closes inside its own
        // branch AND ring digit 1 is reused later for a trailing phenyl.
        let orig = r"COc1ccc2c3c1OC1[C@H](O)[C@](CO)(CCCCCc4ccccc4)CC4C(C2)N(C)CCC341";
        let variant = r"C([C@@]1(CC2C34CCN(C2Cc2ccc(c(c24)OC3[C@@H]1O)OC)C)CO)CCCCc1ccccc1";
        assert!(
            same_canonical(orig, variant),
            "ring-digit reuse must not corrupt a stereocenter whose own ring \
             partner closes inside its branch"
        );
    }

    #[test]
    fn ring_digit_reuse_inside_stereocenter_branch_minimal() {
        // Minimal case matching the real repro's precondition exactly:
        // `[C@H]1` opens ring 1 and its partner closes INSIDE its own first
        // branch `(CC1)` -- i.e. before the parser ever advances to a new
        // *chain* atom for atom0, so atom0's stereo record is still
        // unfinalized at the moment of that closure (the immediate-resolution
        // fast path in `close_or_open_ring` cannot catch it; only the
        // end-of-parse fallback does). Ring digit 1 is then reused by an
        // unrelated, disconnected fragment. Before the fix, the fallback
        // resolved by raw digit and the benzene's ring closure silently
        // stole atom0's still-pending resolution.
        let smi = r"[C@H]1(CC1)Cl.c1ccccc1";
        let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e:?}"));
        let out = canonical_smiles(&mol);
        let mol2 = parse(&out).unwrap_or_else(|e| panic!("re-parse {out}: {e:?}"));
        assert_eq!(
            canonical_smiles(&mol2),
            out,
            "stereocenter with in-branch ring closure + later digit reuse must be stable"
        );
    }

    // ── Allene cumulated double bond stereo ──────────────────────────────────

    #[test]
    fn allene_stereo_two_enantiomers_differ() {
        // F[C@@H]=[C]=[C@H]Cl and F[C@H]=[C]=[C@@H]Cl must produce different canonical SMILES.
        let mol_r = parse("F[C@@H]=[C]=[C@H]Cl").unwrap();
        let mol_s = parse("F[C@H]=[C]=[C@@H]Cl").unwrap();
        let smi_r = canonical_smiles(&mol_r);
        let smi_s = canonical_smiles(&mol_s);
        assert_ne!(
            smi_r, smi_s,
            "allene enantiomers must produce different canonical SMILES: {smi_r}"
        );
    }

    #[test]
    fn allene_stereo_round_trip_stable() {
        for smi in &["F[C@@H]=[C]=[C@H]Cl", "F[C@H]=[C]=[C@@H]Cl"] {
            let mol = parse(smi).unwrap();
            let out = canonical_smiles(&mol);
            let mol2 = parse(&out).unwrap();
            let out2 = canonical_smiles(&mol2);
            assert_eq!(
                out, out2,
                "allene stereo must be stable: {smi} -> {out} -> {out2}"
            );
        }
    }

    // ── RDKit PR #8957: fused-ring stereo round-trip ────────────────────────

    #[test]
    fn ring_stereo_stable_in_fused_system() {
        // RDKit PR #8957: "modern stereo" perception inverted R/S in fused polycyclic
        // systems with multiple stereocenters.  The canonical SMILES form may use @
        // instead of @@ (both are valid encodings of the same stereoisomer depending
        // on traversal order), so the invariant is round-trip stability, not literal
        // @@ count.
        let smi = r"CC[C@@]1(C)C[C@@](CC)(c2ccccc2)CCO1";
        let mol = parse(smi).expect("fused ring stereo mol");
        let out = canonical_smiles(&mol);
        // Round-trip must be stable: canonical(canonical(x)) == canonical(x).
        let mol2 = parse(&out).expect("canonical re-parse");
        let out2 = canonical_smiles(&mol2);
        assert_eq!(
            out, out2,
            "fused ring stereo must be stable after canonical round-trip"
        );
        // The canonical SMILES must still contain at least 2 stereocenters (@/@@ count ≥ 2).
        let stereo_count = out.matches('@').count();
        assert!(
            stereo_count >= 2,
            "both stereocenters must be encoded (got {stereo_count}): {out}"
        );
    }

    // ── E/Z directional-bond canonical stability (issue: Sprint 8) ──────────
    //
    // The canonical writer emits `/`,`\` directional bonds with traversal-direction
    // correction but no separate "normalization" pass. These tests lock in that the
    // direction choice is already deterministic and idempotent for stable skeletons,
    // so a future writer change cannot silently regress E/Z output. (The residual
    // canonical_diff idempotency failures are large fused-polycyclic atom-ranking
    // non-convergence, not a `/`,`\` direction bug — see docs/rdkit_compat.md.)

    /// E/Z parity of the first C=C/C=N double bond: `Some(true)` = E (opposite
    /// outward directions), `Some(false)` = Z, `None` = no specified geometry.
    fn double_bond_is_e(smiles: &str) -> Option<bool> {
        let mol = parse(smiles).unwrap();
        let (a1, a2) = mol
            .bonds()
            .find(|(_, b)| b.order == BondOrder::Double)
            .map(|(_, b)| (b.atom1, b.atom2))?;
        let outward = |end: AtomIdx, other: AtomIdx| -> Option<bool> {
            for (nb, bidx) in mol.neighbors(end) {
                if nb == other {
                    continue;
                }
                let b = mol.bond(bidx);
                match b.order {
                    // `Up` means "up" along atom1→atom2; flip when `end` is atom2.
                    BondOrder::Up => return Some(b.atom1 == end),
                    BondOrder::Down => return Some(b.atom1 != end),
                    _ => {}
                }
            }
            None
        };
        let sa = outward(a1, a2)?;
        let sb = outward(a2, a1)?;
        Some(sa != sb)
    }

    const EZ_STABLE_CORPUS: &[&str] = &[
        "C/C=C/C",     // (E)-2-butene
        "C/C=C\\C",    // (Z)-2-butene
        "F/C=C/F",     // (E)-1,2-difluoroethene
        "F/C=C\\F",    // (Z)
        "CC/C=C/CC",   // (E)-3-hexene
        "CC/C=C\\CC",  // (Z)-3-hexene
        "C/C=C/C=C/C", // (2E,4E)-hexadiene
        "Cl/C=C/Br",
        "C/C=C/c1ccccc1", // (E)-propenylbenzene
        "C/C(F)=C(\\F)C",
    ];

    #[test]
    fn ez_canonical_smiles_is_idempotent() {
        for s in EZ_STABLE_CORPUS {
            assert!(
                is_stable(s),
                "E/Z canonical SMILES must be idempotent for {s}"
            );
        }
    }

    #[test]
    fn ez_geometry_preserved_through_canonicalization() {
        for s in EZ_STABLE_CORPUS {
            let want = double_bond_is_e(s)
                .unwrap_or_else(|| panic!("input {s} must have specified geometry"));
            let canon = canonical_smiles(&parse(s).unwrap());
            let got = double_bond_is_e(&canon)
                .unwrap_or_else(|| panic!("canonical {canon} dropped geometry from {s}"));
            assert_eq!(got, want, "E/Z geometry changed: {s} -> {canon}");
        }
    }

    #[test]
    fn ez_e_and_z_differ_for_each_skeleton() {
        // Each E form must canonicalize differently from its Z form.
        for (e, z) in [
            ("C/C=C/C", "C/C=C\\C"),
            ("F/C=C/F", "F/C=C\\F"),
            ("CC/C=C/CC", "CC/C=C\\CC"),
        ] {
            assert_ne!(
                canonical_smiles(&parse(e).unwrap()),
                canonical_smiles(&parse(z).unwrap()),
                "E and Z must produce different canonical SMILES ({e} vs {z})"
            );
        }
    }

    // ── Fused-aromatic canonical idempotency (Sprint 9) ─────────────────────
    //
    // Lock in the fused aromatics that DO round-trip consistently. The residual
    // ~1.6% canonical idempotency failures on large fused polycyclics are caused
    // by aromaticity-perception round-trip inconsistency (a molecule vs the
    // re-parse of its own canonical SMILES can disagree on which bonds are
    // aromatic — e.g. 16 vs 17 on a fluorene-type linkage — which shifts Morgan
    // ranks). That is an aromaticity/parser-core issue, not a canonical-ranking
    // bug; see docs/rdkit_compat.md. These cases are stable and guarded here.

    #[test]
    fn fused_aromatic_canonical_is_idempotent() {
        for s in [
            "c1ccc2ccccc2c1",         // naphthalene
            "c1ccc2ncccc2c1",         // quinoline
            "c1ccc2c(c1)cc[nH]2",     // indole
            "c1ccc2cc3ccccc3cc2c1",   // anthracene
            "c1ccc2[nH]c3ccccc3c2c1", // carbazole
            "c1ccc2c(c1)oc1ccccc12",  // dibenzofuran
        ] {
            assert!(
                is_stable(s),
                "fused-aromatic canonical SMILES must be idempotent for {s}"
            );
        }
    }
}
