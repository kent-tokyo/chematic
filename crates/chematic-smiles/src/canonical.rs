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

use crate::writer::{emit_bracket_hydrogens, suppress_standalone_wedge};

/// Return the atom indices sorted into canonical (Morgan-rank) order.
///
/// The returned `Vec<usize>` lists atom positions (0-based) in the order they
/// would be encountered during a canonical DFS write.  Atoms with higher
/// Morgan rank appear earlier.  This is the same ordering `canonical_smiles`
/// uses internally: raw `morgan_ranks` ties are resolved via the same
/// individualize-refine + lexicographically-smallest-string selection, not
/// left as an input-order-dependent plateau.
///
/// Useful for normalizing atom-indexed property arrays to a canonical order.
pub fn canonical_atom_order(mol: &Molecule) -> Vec<usize> {
    let n = mol.atom_count();
    if n == 0 {
        return Vec::new();
    }
    let ranks = winning_individualized_ranks(mol);
    let mut order: Vec<usize> = (0..n).collect();
    // Sort descending by rank (highest rank first, as in canonical DFS).
    order.sort_unstable_by(|&a, &b| ranks[b].cmp(&ranks[a]));
    order
}

/// Resolve `morgan_ranks` ties via individualize-refine and return the fully
/// discrete per-atom ranks of whichever branch produces the
/// lexicographically smallest canonical SMILES -- shared by
/// `canonical_smiles` and `canonical_atom_order` so both use the identical
/// tie-break, instead of `canonical_atom_order` silently falling back to raw
/// (tie-break-free) `morgan_ranks`.
fn winning_individualized_ranks(mol: &Molecule) -> Vec<u64> {
    let plateaued = morgan_ranks(mol);
    let mut budget = MAX_INDIVIDUALIZE_BRANCHES;
    let branches = enumerate_discrete_ranks(mol, plateaued, &mut budget);
    branches
        .into_iter()
        .min_by_key(|ranks| CanonicalWriter::new(mol, ranks).write_all())
        .unwrap_or_default()
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
///
/// Atom ordering is fully discretized before writing (individualize-refine,
/// see `enumerate_discrete_ranks`): when the plain Morgan refinement in
/// [`morgan_ranks`] plateaus with genuine (non-automorphism) ties still
/// present, every possible resolution is tried and the lexicographically
/// smallest resulting string is returned. This makes the output invariant to
/// which input atom ordering/spelling the molecule was parsed from, not just
/// idempotent under repeated self-canonicalization.
pub fn canonical_smiles(mol: &Molecule) -> String {
    if mol.atom_count() == 0 {
        return String::new();
    }

    let ranks = winning_individualized_ranks(mol);
    CanonicalWriter::new(mol, &ranks).write_all()
}

/// Compute Morgan (extended connectivity) ranks for all atoms.
///
/// Returns a vector of normalised ordinal ranks (0-based, gap-free)
/// indexed by atom position (same order as `mol.atoms()`). This is pure
/// neighbor-hash refinement to a fixpoint -- it does NOT individualize
/// remaining ties, so atoms in the same non-trivial automorphism orbit (or
/// in a refinement cell that merely *contains* an orbit) keep equal ranks.
/// That is the correct, useful notion of "rank" for topological-symmetry
/// queries (see [`equivalent_atom_classes`], [`are_atoms_equivalent`]).
///
/// [`canonical_smiles`] does NOT use this directly for atom ordering when
/// ties remain -- see `enumerate_discrete_ranks` for the individualize-refine
/// step that resolves ties before writing.
pub fn morgan_ranks(mol: &Molecule) -> Vec<u64> {
    let n = mol.atom_count();
    let initial: Vec<u64> = (0..n)
        .map(|i| initial_invariant(mol, AtomIdx(i as u32)))
        .collect();
    refine_ranks(mol, initial)
}

/// Refine `ranks` (any starting coloring, not necessarily the initial
/// invariant) via neighbor-hash iteration until the number of distinct
/// classes stops increasing. Used both for the plain (tie-preserving) ranks
/// in [`morgan_ranks`] and, with a perturbed starting coloring, as the
/// "refine" half of individualize-refine in `enumerate_discrete_ranks`.
fn refine_ranks(mol: &Molecule, mut ranks: Vec<u64>) -> Vec<u64> {
    let n = ranks.len();
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

/// Safety cap on the number of discrete rank assignments
/// `enumerate_discrete_ranks` will explore. It exists to guarantee
/// termination on pathologically symmetric inputs (fullerene fragments,
/// deep dendrimers) where the principled fix is automorphism-aware branch
/// pruning (nauty-style), not attempted here. Once exhausted, the
/// remaining ties in that branch fall back to `canonical_cmp`'s finite
/// tie-break chain (deterministic, but not guaranteed order-independent).
///
/// This IS hit in practice, corrected 2026-07-12 after a claim of "never
/// hit" here went unverified: measured on 5,000 real ChEMBL-derived
/// molecules, 3 (0.06%) exceeded this cap, needing up to 168,219 branches
/// (16.8x). All three are real drug-synthesis intermediates with multiple
/// Boc/pivaloyl tert-butyl protecting groups, each an independent 3-way
/// symmetric orbit that multiplies combinatorially across the molecule.
/// For all three, truncation at 10,000 was confirmed (against an
/// unbounded run, and separately against 32 independent re-spellings) to
/// still find the correct lexicographically-smallest winner -- the
/// exhausted cells in these cases are true automorphism orbits (every
/// individualization within the cell writes the same string, so the
/// blowup is redundant duplicates, not competing candidates). This is not
/// a guarantee for the general case: the failure mode this cap is meant
/// to bound is a cell that merely *contains* an orbit (genuinely different
/// candidates truncated away), which no real molecule in this corpus
/// happened to exercise. Raising the constant is not a principled fix (the
/// observed distribution has a cliff -- p99.9 is ~4,922, but the 3
/// offenders need 74k-168k -- so no fixed multiple closes the gap); the
/// real fix is the orbit-aware pruning mentioned above.
const MAX_INDIVIDUALIZE_BRANCHES: usize = 10_000;

/// Individualize atom `atom_idx` within its current rank class: insert a new
/// rank strictly between its class and the next-higher class, so a
/// subsequent refinement pass can propagate the distinction through the rest
/// of the graph. `ranks` must be gap-free ordinals (as produced by
/// `refine_ranks`/`normalize_ranks`).
fn individualize(ranks: &[u64], atom_idx: usize) -> Vec<u64> {
    let v = ranks[atom_idx];
    ranks
        .iter()
        .enumerate()
        .map(|(i, &r)| {
            if i == atom_idx {
                v + 1
            } else if r > v {
                r + 1
            } else {
                r
            }
        })
        .collect()
}

/// Enumerate every discrete (all-singleton) rank assignment reachable from
/// `ranks` (already refined to a fixpoint) via individualize-refine
/// branching.
///
/// Refinement cells are always unions of automorphism orbits (a standard
/// 1-WL / equitable-partition fact: refinement can never split an orbit).
/// So: if a cell IS an orbit, every choice of which atom to individualize
/// yields an automorphic result -- the resulting SMILES strings are
/// identical, so exploring all of them is correct but redundant. If a cell
/// properly CONTAINS an orbit, no order-independent rule can select a single
/// representative (if one existed, refinement would already have used it as
/// an invariant and the cell would not be tied) -- the only order-independent
/// resolution is to try every atom in the cell and let the caller take the
/// lexicographically smallest resulting string. Cell SELECTION (which
/// non-singleton class to branch on next) is itself order-independent:
/// always the lowest-ranked non-singleton cell.
fn enumerate_discrete_ranks(mol: &Molecule, ranks: Vec<u64>, budget: &mut usize) -> Vec<Vec<u64>> {
    let mut by_rank: Vec<Vec<usize>> = Vec::new();
    for (i, &r) in ranks.iter().enumerate() {
        let r = r as usize;
        if by_rank.len() <= r {
            by_rank.resize(r + 1, Vec::new());
        }
        by_rank[r].push(i);
    }

    // `by_rank` is indexed by ordinal rank value, so the first multi-member
    // entry found is the lowest-ranked non-singleton cell.
    let Some(members) = by_rank.iter().find(|m| m.len() > 1) else {
        return vec![ranks];
    };

    let mut results = Vec::new();
    for &atom_idx in members {
        if *budget == 0 {
            break;
        }
        *budget -= 1;
        let individualized = individualize(&ranks, atom_idx);
        let re_refined = refine_ranks(mol, individualized);
        results.extend(enumerate_discrete_ranks(mol, re_refined, budget));
    }
    results
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
    /// (ring_num, bond_order, ring_partner_atom, physical_bond)
    atom_ring_nums: HashMap<AtomIdx, Vec<(u32, BondOrder, AtomIdx, BondIdx)>>,
    next_ring: u32,
    out: String,
    /// Union-find groups of directional (`/`/`\`) bonds that jointly encode
    /// one connected E/Z system — flipping every member preserves geometry,
    /// flipping a subset does not. Keyed/rooted by `BondIdx`.
    ez_group: HashMap<BondIdx, BondIdx>,
    /// Groups whose first-encountered bond (in write order) came out `Down`;
    /// every remaining bond in the group is flipped so the first directional
    /// bond of each system is always `/`, regardless of input spelling.
    ez_flip: HashMap<BondIdx, bool>,
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
            ez_group: HashMap::new(),
            ez_flip: HashMap::new(),
        }
    }

    /// Union all directional single bonds flanking each stereo double bond
    /// into one group per connected E/Z system (order-independent — depends
    /// only on molecule topology, not on canonical ranks or write order).
    fn build_ez_groups(&mut self) {
        fn is_directional(mol: &Molecule, bidx: BondIdx, order: BondOrder) -> bool {
            matches!(order, BondOrder::Up | BondOrder::Down) || mol.bond_direction(bidx).is_some()
        }

        fn find(group: &mut HashMap<BondIdx, BondIdx>, x: BondIdx) -> BondIdx {
            let parent = *group.get(&x).unwrap_or(&x);
            if parent == x {
                x
            } else {
                let root = find(group, parent);
                group.insert(x, root);
                root
            }
        }

        fn union(group: &mut HashMap<BondIdx, BondIdx>, a: BondIdx, b: BondIdx) {
            let ra = find(group, a);
            let rb = find(group, b);
            if ra != rb {
                group.insert(ra, rb);
            }
        }

        for bidx in 0..self.mol.bond_count() {
            let bidx = BondIdx(bidx as u32);
            let bond = self.mol.bond(bidx);
            if bond.order != BondOrder::Double {
                continue;
            }
            let mut side_bonds = Vec::new();
            for endpoint in [bond.atom1, bond.atom2] {
                for (_, nb_bidx) in self.mol.neighbors(endpoint) {
                    if nb_bidx == bidx {
                        continue;
                    }
                    let nb_order = self.mol.bond(nb_bidx).order;
                    if is_directional(self.mol, nb_bidx, nb_order) {
                        side_bonds.push(nb_bidx);
                    }
                }
            }
            let Some(&first) = side_bonds.first() else {
                continue;
            };
            self.ez_group.entry(first).or_insert(first);
            for &b in &side_bonds[1..] {
                self.ez_group.entry(b).or_insert(b);
                union(&mut self.ez_group, first, b);
            }
        }
    }

    /// Normalize a directional bond order so the first occurrence of each
    /// E/Z system in canonical write order is always `Up` (`/`); every other
    /// bond in the system is flipped consistently to preserve geometry.
    fn normalize_ez(&mut self, bidx: BondIdx, order: BondOrder) -> BondOrder {
        if !matches!(order, BondOrder::Up | BondOrder::Down) {
            return order;
        }
        let root = {
            let mut x = *self.ez_group.get(&bidx).unwrap_or(&bidx);
            while let Some(&p) = self.ez_group.get(&x) {
                if p == x {
                    break;
                }
                x = p;
            }
            x
        };
        let flip = *self.ez_flip.entry(root).or_insert(order == BondOrder::Down);
        if flip {
            match order {
                BondOrder::Up => BondOrder::Down,
                BondOrder::Down => BondOrder::Up,
                other => other,
            }
        } else {
            order
        }
    }

    fn write_all(mut self) -> String {
        // Phase 0: group directional bonds into connected E/Z systems
        // (topology-only, independent of canonical order).
        self.build_ez_groups();

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
                self.atom_ring_nums.entry(neighbor).or_default().push((
                    rn,
                    order_at_open,
                    atom,
                    bidx,
                )); // partner = close atom
                self.atom_ring_nums.entry(atom).or_default().push((
                    rn,
                    order_at_close,
                    neighbor,
                    bidx,
                )); // partner = open atom
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
            for (rn, bond_order, _partner, bidx) in rings {
                let bond_order = self.normalize_ez(bidx, bond_order);
                let bond_order = suppress_standalone_wedge(self.mol, bidx, bond_order);
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
        let mut children: Vec<(AtomIdx, BondIdx, BondOrder)> = self
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
                (nb, bidx, order)
            })
            .collect();

        // Sort children by canonical rank (ascending → highest rank = main chain).
        children.sort_by(|&(a, ..), &(b, ..)| self.canonical_cmp(a, b));

        let n = children.len();
        for (i, (child, bidx, bond_order)) in children.into_iter().enumerate() {
            // Normalized here (not in the map above) so the flip decision is
            // made in true left-to-right write order: this atom's earlier
            // (lower-rank) children have already fully recursed by the time
            // a later sibling's direction is decided.
            let bond_order = self.normalize_ez(bidx, bond_order);
            let bond_order = suppress_standalone_wedge(self.mol, bidx, bond_order);
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

            emit_bracket_hydrogens(&mut self.out, self.mol, idx);

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
            for &(_, _, partner, _) in rings {
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

    /// Build a copy of `mol` with atoms reordered by `perm` (perm[new_idx] = old_idx).
    /// Bonds are remapped to the new indices; stereo/direction metadata is
    /// intentionally dropped since this helper only exists to test whether the
    /// *skeleton* rank partition is invariant under atom relabeling.
    fn permute_molecule(mol: &Molecule, perm: &[usize]) -> Molecule {
        let mut old_to_new = vec![0u32; perm.len()];
        for (new_idx, &old_idx) in perm.iter().enumerate() {
            old_to_new[old_idx] = new_idx as u32;
        }
        let mut builder = chematic_core::MoleculeBuilder::new();
        for &old_idx in perm {
            builder.add_atom(mol.atom(AtomIdx(old_idx as u32)).clone());
        }
        for (_, bond) in mol.bonds() {
            let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
            let b = AtomIdx(old_to_new[bond.atom2.0 as usize]);
            let _ = builder.add_bond(a, b, bond.order);
        }
        builder.build()
    }

    /// Relabel a rank vector into "first-seen order" group ids, so partitions
    /// can be compared structurally (which atoms are grouped together)
    /// independent of the actual numeric rank values assigned.
    fn partition_key(ranks: &[u64]) -> Vec<usize> {
        let mut seen: Vec<u64> = Vec::new();
        ranks
            .iter()
            .map(|&r| match seen.iter().position(|&s| s == r) {
                Some(pos) => pos,
                None => {
                    seen.push(r);
                    seen.len() - 1
                }
            })
            .collect()
    }

    /// Sanity check demanded before any individualize-refine rewrite: the
    /// refinement-to-plateau partition (which atoms share a rank, not the raw
    /// numeric values) must be invariant under input atom permutation. If this
    /// fails, the root cause is a non-invariant initial invariant / refinement
    /// step itself (over-splitting orbits), not a missing individualize step
    /// (under-splitting ties) -- a completely different bug to fix.
    #[test]
    fn morgan_ranks_partition_is_permutation_invariant() {
        let corpus = [
            "O=C(NCc1cccnc1)NC[C@H]1CCC[C@H](OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)[C@@H]1c1ccccc1",
            "c1ccccc1",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "O=C1CCC(=O)N1",
            "c1ccc2ccc3ccccc3c2c1",
            "O=C(NCc1cccnc1)NCC1CCCC(OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)C1c1ccccc1",
        ];

        for smi in corpus {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e}"));
            let n = mol.atom_count();
            let part_orig = partition_key(&morgan_ranks(&mol));

            // A few deterministic, non-identity permutations (no RNG dependency).
            let perms: Vec<Vec<usize>> = vec![
                (0..n).rev().collect(),
                {
                    let mut p: Vec<usize> = (0..n).collect();
                    if n > 2 {
                        p.rotate_left(n / 3 + 1);
                    }
                    p
                },
                {
                    let mut p: Vec<usize> = (0..n).rev().collect();
                    if n > 3 {
                        p.swap(1, n - 2);
                        p.rotate_right(2);
                    }
                    p
                },
            ];

            for perm in perms {
                let permuted = permute_molecule(&mol, &perm);
                let part_perm = partition_key(&morgan_ranks(&permuted));

                // inverse: where did old atom `old_idx` land in the permuted molecule?
                let mut new_of_old = vec![0usize; n];
                for (new_idx, &old_idx) in perm.iter().enumerate() {
                    new_of_old[old_idx] = new_idx;
                }

                for i in 0..n {
                    for j in (i + 1)..n {
                        let same_orig = part_orig[i] == part_orig[j];
                        let same_perm = part_perm[new_of_old[i]] == part_perm[new_of_old[j]];
                        assert_eq!(
                            same_orig, same_perm,
                            "partition not permutation-invariant for '{smi}': \
                             atoms {i},{j} (perm {perm:?})"
                        );
                    }
                }
            }
        }
    }

    /// `canonical_atom_order` must be permutation-invariant: relabeling the
    /// same molecule's atoms (different parse order) must not change WHICH
    /// symmetry class of atom appears 1st, 2nd, 3rd, ... in the returned
    /// order. Unlike `canonical_smiles`, `canonical_atom_order` does not run
    /// individualize-refine -- it sorts raw `morgan_ranks` with no tie-break,
    /// so this is expected to fail on any molecule with a genuine
    /// (non-singleton) rank tie. This is a diagnostic probe for that gap, not
    /// an already-passing invariant.
    #[test]
    fn canonical_atom_order_permutation_invariance_probe() {
        let corpus = [
            "O=C(NCc1cccnc1)NC[C@H]1CCC[C@H](OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)[C@@H]1c1ccccc1",
            "c1ccccc1",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "O=C1CCC(=O)N1",
            "c1ccc2ccc3ccccc3c2c1",
            "O=C(NCc1cccnc1)NCC1CCCC(OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)C1c1ccccc1",
        ];

        let mut bad = 0;
        let mut total = 0;
        for smi in corpus {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e}"));
            let n = mol.atom_count();
            // Ground-truth atom-class labels, in the ORIGINAL molecule's index
            // space only -- never relabeled independently for the permuted
            // copy, or two arbitrary first-seen-order numberings would be
            // compared against each other and any mismatch would be a
            // methodology artifact, not a real instability (see the
            // project's canonicalization-tie-break-theory / measurement-
            // harness-controls lessons).
            let part_orig = partition_key(&morgan_ranks(&mol));
            let order_orig = canonical_atom_order(&mol);
            let profile_orig: Vec<usize> = order_orig.iter().map(|&i| part_orig[i]).collect();

            let perms: Vec<Vec<usize>> = vec![(0..n).rev().collect(), {
                let mut p: Vec<usize> = (0..n).collect();
                if n > 2 {
                    p.rotate_left(n / 3 + 1);
                }
                p
            }];

            for perm in perms {
                total += 1;
                // perm[new_idx] = old_idx (see permute_molecule's contract).
                let permuted = permute_molecule(&mol, &perm);
                let order_perm = canonical_atom_order(&permuted);
                // Map each returned NEW index back to the class of the
                // corresponding OLD atom, via `part_orig` -- the same
                // ground-truth labeling used for profile_orig.
                let profile_perm: Vec<usize> = order_perm
                    .iter()
                    .map(|&new_i| part_orig[perm[new_i]])
                    .collect();

                if profile_orig != profile_perm {
                    bad += 1;
                    eprintln!(
                        "canonical_atom_order NOT permutation-invariant for '{smi}' (perm {perm:?}): \
                         {profile_orig:?} != {profile_perm:?}"
                    );
                }
            }
        }
        eprintln!("canonical_atom_order instability: {bad}/{total} permutation trials");
        assert_eq!(
            bad, 0,
            "{bad}/{total} permutation trials were unstable -- see stderr"
        );
    }

    /// Direct probe (no permutation needed): does `canonical_atom_order`'s
    /// naive `morgan_ranks`-only sort ever disagree with the FULLY
    /// individualized/resolved rank order that `canonical_smiles` actually
    /// verified-correct output is built from? Comparing against
    /// same-partition-class labels (as the permutation probe above does) is
    /// blind to intra-class reordering among genuinely symmetric
    /// (automorphism-equivalent) atoms, where any order is harmless -- this
    /// test instead reconstructs the winning individualized branch directly,
    /// so it also catches disagreement WITHIN a Morgan-rank tie that is not
    /// a true automorphism (exactly the class of bug individualize-refine
    /// was built to fix for `canonical_smiles`, see Round 10-12 history).
    #[test]
    fn canonical_atom_order_matches_individualized_ranks_probe() {
        let corpus = [
            "O=C(NCc1cccnc1)NC[C@H]1CCC[C@H](OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)[C@@H]1c1ccccc1",
            "c1ccccc1",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            "O=C1CCC(=O)N1",
            "c1ccc2ccc3ccccc3c2c1",
            "O=C(NCc1cccnc1)NCC1CCCC(OCc2cc(C(F)(F)F)cc(C(F)(F)F)c2)C1c1ccccc1",
            // Extra cases picked for symmetry that Weisfeiler-Leman-style
            // refinement is known to struggle with (fused/bridged systems).
            "C1CC2CCC1CC2",
            "C1CC2CC1CC2",
            "c1ccc(-c2ccccc2)cc1",
            "OC1CCC(O)CC1",
            "C12CC3CC(CC(C3)C1)C2",
        ];

        let mut needed_individualization = 0;
        let mut mismatched = 0;
        for smi in corpus {
            let mol = parse(smi).unwrap_or_else(|e| panic!("{smi}: {e}"));
            let plateaued = morgan_ranks(&mol);
            let mut budget = MAX_INDIVIDUALIZE_BRANCHES;
            let branches = enumerate_discrete_ranks(&mol, plateaued, &mut budget);
            if branches.len() > 1 {
                needed_individualization += 1;
            }
            let winning_ranks = branches
                .into_iter()
                .min_by_key(|ranks| CanonicalWriter::new(&mol, ranks).write_all())
                .expect("at least one branch");

            let n = mol.atom_count();
            let mut winning_order: Vec<usize> = (0..n).collect();
            winning_order.sort_by(|&a, &b| winning_ranks[b].cmp(&winning_ranks[a]));

            let naive_order = canonical_atom_order(&mol);
            if naive_order != winning_order {
                mismatched += 1;
                eprintln!(
                    "canonical_atom_order disagrees with resolved canonical order for '{smi}': \
                     naive={naive_order:?} resolved={winning_order:?}"
                );
            }
        }
        eprintln!(
            "{needed_individualization}/{} molecules needed individualization; \
             {mismatched}/{} disagreed with canonical_atom_order",
            corpus.len(),
            corpus.len()
        );
        assert_eq!(
            mismatched, 0,
            "canonical_atom_order must match the individualized/resolved order -- see stderr"
        );
    }

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

    // ── Round 12: simple E/Z direction normalization ────────────────────────
    //
    // `/N=N/` and `\N=N\` are two equally valid SMILES spellings of the same
    // geometry (flipping every directional bond of one connected E/Z system
    // preserves meaning). Before this fix, the writer just propagated
    // whichever direction the parser happened to read, so two spellings of
    // the same molecule could canonicalize to two different strings. The fix
    // normalizes each connected E/Z system so its first directional bond (in
    // canonical write order) is always `/`.

    #[test]
    fn ez_simple_bond_direction_normalized_azo() {
        assert!(
            same_canonical("CN(C)/N=N/c1ccccc1", r"CN(C)\N=N\c1ccccc1",),
            "isolated E/Z double bond must canonicalize identically regardless \
             of which of the two equally-valid slash spellings was parsed"
        );
    }

    #[test]
    fn ez_simple_bond_direction_normalized_symmetric() {
        assert!(same_canonical("F/C=C/F", r"F\C=C\F"));
        assert!(same_canonical("C(/F)=C/F", r"C(\F)=C\F"));
    }

    // Bracket atoms with `hydrogen_count: None` must still get their implicit
    // hydrogens written by the canonical writer too — same bug/fixtures as
    // writer::tests::test_bracket_implicit_h_*.
    use chematic_core::{Atom, Element};

    #[test]
    fn canonical_bracket_implicit_h_ammonium_charge_only() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut n = Atom::new(Element::N);
        n.charge = 1;
        b.add_atom(n);
        assert_eq!(canonical_smiles(&b.build()), "[NH4+]");
    }

    #[test]
    fn canonical_bracket_implicit_h_isotope_only() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut c = Atom::new(Element::C);
        c.isotope = Some(13);
        b.add_atom(c);
        assert_eq!(canonical_smiles(&b.build()), "[13CH4]");
    }

    #[test]
    fn canonical_bracket_implicit_h_atom_map_only() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut c0 = Atom::new(Element::C);
        c0.atom_map = Some(7);
        let c0 = b.add_atom(c0);
        let c1 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        assert_eq!(canonical_smiles(&b.build()), "C[CH3:7]");
    }

    #[test]
    fn canonical_bracket_implicit_h_isotope_and_atom_map() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut c0 = Atom::new(Element::C);
        c0.isotope = Some(13);
        c0.atom_map = Some(7);
        let c0 = b.add_atom(c0);
        let c1 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        assert_eq!(canonical_smiles(&b.build()), "C[13CH3:7]");
    }

    #[test]
    fn canonical_bracket_implicit_h_hydroxide_charge_only() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let mut o = Atom::new(Element::O);
        o.charge = -1;
        b.add_atom(o);
        assert_eq!(canonical_smiles(&b.build()), "[OH-]");
    }

    // ── Standalone wedge/hash bond must not emit a meaningless SMILES
    // directional token in canonical_smiles either (same fix as
    // `writer::tests`, applied to the canonical writer's own emission
    // sites) ──────────────────────────────────────────────────────────────

    /// Minimal repro from docs/stereo2d_reader_integration_rfc.md §3: a
    /// wedge bond with no adjacent double bond must not become `/`/`\`.
    #[test]
    fn canonical_standalone_solid_wedge_not_written_as_slash() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, br, BondOrder::Up).unwrap();
        let mol = b.build();

        let out = canonical_smiles(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone solid wedge (no adjacent double bond) must not be \
             written as a directional token: got '{out}'"
        );
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Up);

        let mol2 = crate::parser::parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), mol.atom_count());
        assert_eq!(mol2.bond_count(), mol.bond_count());
    }

    /// Same shape but with a hash wedge (`BondOrder::Down`).
    #[test]
    fn canonical_standalone_hash_wedge_not_written_as_backslash() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, br, BondOrder::Down).unwrap();
        let mol = b.build();

        let out = canonical_smiles(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone hash wedge (no adjacent double bond) must not be \
             written as a directional token: got '{out}'"
        );
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Down);

        let mol2 = crate::parser::parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), mol.atom_count());
        assert_eq!(mol2.bond_count(), mol.bond_count());
    }

    /// A wedge bond landing on a ring-closure edge must be suppressed
    /// identically in the canonical writer's ring-closure emission site.
    #[test]
    fn canonical_standalone_wedge_on_ring_closure_bond() {
        let mut b = chematic_core::MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        b.add_bond(c1, c3, BondOrder::Up).unwrap();
        let mol = b.build();

        let out = canonical_smiles(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone wedge on a ring-closure bond must not be written as \
             a directional token: got '{out}'"
        );
        let mol2 = crate::parser::parse(&out).unwrap_or_else(|e| panic!("re-parse '{out}': {e}"));
        assert_eq!(mol2.atom_count(), 3);
        assert_eq!(mol2.bond_count(), 3);
    }

    /// The fix must not disturb `Atom.chirality`/`Molecule::stereo_neighbor_order`
    /// while suppressing a standalone wedge bond's spurious token on the
    /// same stereocenter -- mirrors `writer::tests::
    /// test_standalone_wedge_does_not_disturb_stereocenter_chirality` for
    /// the canonical writer's own (rank-dependent) chirality-correction path.
    #[test]
    fn canonical_standalone_wedge_does_not_disturb_stereocenter_chirality() {
        let mut center = Atom::new(Element::C);
        center.chirality = Chirality::Clockwise;
        // Force bracket notation (h=0 still forces `[..]` without printing
        // an `H` token) so the chirality symbol is actually emitted.
        center.hydrogen_count = Some(0);
        let mut b = chematic_core::MoleculeBuilder::new();
        let c = b.add_atom(center);
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let iodine = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        let wedge_bond = b.add_bond(c, iodine, BondOrder::Up).unwrap();
        b.set_stereo_neighbor_order(c, vec![f.0, cl.0, br.0, iodine.0]);
        let mol = b.build();

        let out = canonical_smiles(&mol);
        assert!(
            !out.contains('/') && !out.contains('\\'),
            "standalone wedge on a stereocenter's substituent must still be \
             suppressed: got '{out}'"
        );
        assert!(
            out.contains('@'),
            "the stereocenter's own chirality symbol must still be emitted \
             alongside the (now-suppressed) wedge bond: got '{out}'"
        );

        assert_eq!(mol.atom(c).chirality, Chirality::Clockwise);
        assert_eq!(mol.bond(wedge_bond).order, BondOrder::Up);
        assert_eq!(
            mol.stereo_neighbor_order(c),
            Some([f.0, cl.0, br.0, iodine.0].as_slice())
        );
    }

    /// Regression guard: genuine E/Z directional markers must still be
    /// emitted by canonical_smiles (already covered by `ez_e_stable` etc.
    /// above; this pins the specific double-bond-adjacency mechanism the
    /// standalone-wedge fix now gates on).
    #[test]
    fn canonical_genuine_ez_directional_bond_still_written() {
        let mol = parse("F/C=C/F").unwrap();
        let out = canonical_smiles(&mol);
        assert!(
            out.contains('/') || out.contains('\\'),
            "genuine double-bond-adjacent directional marker must survive \
             canonical_smiles(): got '{out}'"
        );
    }

    /// The pre-existing "E/Z direction stashed on an aromatic bond" fixtures
    /// (`parser::tests::direction_stash_*`) all involve a real exocyclic
    /// double bond and must keep round-tripping unaffected by this fix --
    /// exercised end-to-end here via `canonical_smiles` directly (those
    /// tests live in parser.rs and already re-run on every `cargo test`).
    #[test]
    fn canonical_aromatic_stash_with_real_double_bond_still_emits_direction() {
        let mol = parse(r"N=c1\c(O)c(O)\c1=N").unwrap();
        let out = canonical_smiles(&mol);
        assert!(
            out.contains('/') || out.contains('\\'),
            "genuine exocyclic-double-bond-adjacent aromatic stash must \
             still emit a directional token: got '{out}'"
        );
    }
}
