//! Kekulization: assign alternating single/double bonds to aromatic systems.
//!
//! Algorithm (4 passes):
//!
//! 1. Collect all aromatic atoms and bonds.
//! 2. Determine the must-match set (atoms that need a double bond).
//! 3. Find a maximum matching:
//!    - Pass 1: BFS augmenting paths, ascending atom order.
//!    - Pass 2: BFS augmenting paths, descending order (fallback).
//!    - Pass 3: Bridgehead-N exclusion (lone-pair donors at ring junctions).
//!    - Pass 4: Edmonds' blossom for non-bipartite aromatic subgraphs.
//! 4. Matched edges → Double; unmatched → Single.

use std::collections::{HashMap, HashSet, VecDeque};

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
    kekulize_with(mol, |idx| atom_must_be_matched(mol, idx))
}

/// [`kekulize`] with the caller's rule for which aromatic atoms must take a
/// double bond (an aromatic bond's atom for which `must_match` is false only
/// gets single bonds). Reaction products use it to apply RDKit's candidate
/// rule to atoms a template gave a non-aromatic double bond.
pub fn kekulize_with(
    mol: &Molecule,
    must_match_atom: impl Fn(AtomIdx) -> bool,
) -> Result<KekuleResult, KekuleError> {
    let n = mol.atom_count();
    let (aromatic_bonds, aromatic_atoms) = aromatic_subgraph(mol);
    if aromatic_bonds.is_empty() {
        return Ok(HashMap::new());
    }
    let (must_match, sorted_atoms) = required_aromatic_atoms(&aromatic_atoms, must_match_atom);
    let adj = matching_adjacency(mol, &aromatic_bonds, &must_match);
    let mut matcher = Matcher::new(n);
    let ring_fusion_pnictogen = |idx: AtomIdx| {
        let atom = mol.atom(idx);
        matches!(atom.element.atomic_number(), 7 | 15)
            && atom.charge == 0
            && mol
                .neighbors(idx)
                .filter(|&(_, b)| mol.bond(b).order == BondOrder::Aromatic)
                .count()
                >= 3
    };
    if sorted_atoms.iter().any(|&idx| ring_fusion_pnictogen(idx))
        && let Some(result) = matching_without(
            mol,
            &aromatic_bonds,
            &must_match,
            &sorted_atoms,
            &mut matcher,
            ring_fusion_pnictogen,
        )
    {
        return Ok(result);
    }
    run_matching_both_orders(&mut matcher, &sorted_atoms, &adj);
    if !matcher.all_matched(&sorted_atoms) {
        let is_bridgehead_n = |idx: AtomIdx| {
            mol.atom(idx).element.atomic_number() == 7 && adj[idx.0 as usize].len() >= 3
        };
        if sorted_atoms.iter().any(|&idx| is_bridgehead_n(idx))
            && let Some(result) = matching_without(
                mol,
                &aromatic_bonds,
                &must_match,
                &sorted_atoms,
                &mut matcher,
                is_bridgehead_n,
            )
        {
            return Ok(result);
        }
    }
    if !matcher.all_matched(&sorted_atoms) {
        apply_blossom_matching(&mut matcher, &sorted_atoms, &adj, n);
    }
    if let Some(&idx) = sorted_atoms
        .iter()
        .find(|&&idx| matcher.mate[idx.0 as usize] == NO_MATE)
    {
        return Err(KekuleError {
            detail: format!(
                "atom {} ({}) cannot be assigned a double bond",
                idx.0,
                mol.atom(idx).element.symbol()
            ),
        });
    }

    Ok(kekule_result_from_mates(
        &aromatic_bonds,
        mol,
        &matcher.mate,
    ))
}

fn aromatic_subgraph(mol: &Molecule) -> (Vec<BondIdx>, Vec<bool>) {
    let mut bonds = Vec::new();
    let mut atoms = vec![false; mol.atom_count()];
    for (index, bond) in mol.bonds() {
        if bond.order == BondOrder::Aromatic {
            bonds.push(index);
            atoms[bond.atom1.0 as usize] = true;
            atoms[bond.atom2.0 as usize] = true;
        }
    }
    (bonds, atoms)
}

fn required_aromatic_atoms(
    aromatic: &[bool],
    predicate: impl Fn(AtomIdx) -> bool,
) -> (Vec<bool>, Vec<AtomIdx>) {
    let mut required = vec![false; aromatic.len()];
    let mut atoms = Vec::new();
    for (index, &is_aromatic) in aromatic.iter().enumerate() {
        let atom = AtomIdx(index as u32);
        if is_aromatic && predicate(atom) {
            required[index] = true;
            atoms.push(atom);
        }
    }
    (required, atoms)
}

fn run_matching_both_orders(matcher: &mut Matcher, atoms: &[AtomIdx], adjacency: &MatchAdj) {
    matcher.clear();
    matcher.run_pass(atoms.iter().copied(), adjacency);
    if !matcher.all_matched(atoms) {
        matcher.clear();
        matcher.run_pass(atoms.iter().rev().copied(), adjacency);
    }
}

fn matching_without(
    mol: &Molecule,
    aromatic_bonds: &[BondIdx],
    required: &[bool],
    atoms: &[AtomIdx],
    matcher: &mut Matcher,
    exclude: impl Fn(AtomIdx) -> bool,
) -> Option<KekuleResult> {
    let mut reduced = required.to_vec();
    let remaining: Vec<AtomIdx> = atoms
        .iter()
        .copied()
        .filter(|&atom| {
            let excluded = exclude(atom);
            if excluded {
                reduced[atom.0 as usize] = false;
            }
            !excluded
        })
        .collect();
    let adjacency = matching_adjacency(mol, aromatic_bonds, &reduced);
    run_matching_both_orders(matcher, &remaining, &adjacency);
    matcher
        .all_matched(&remaining)
        .then(|| kekule_result_from_mates(aromatic_bonds, mol, &matcher.mate))
}

fn apply_blossom_matching(
    matcher: &mut Matcher,
    atoms: &[AtomIdx],
    adjacency: &MatchAdj,
    atom_count: usize,
) {
    let mut compact_index = vec![usize::MAX; atom_count];
    for (index, &atom) in atoms.iter().enumerate() {
        compact_index[atom.0 as usize] = index;
    }
    let compact_adjacency: Vec<Vec<usize>> = atoms
        .iter()
        .map(|&atom| {
            adjacency[atom.0 as usize]
                .iter()
                .filter_map(|(neighbor, _)| {
                    let index = compact_index[neighbor.0 as usize];
                    (index != usize::MAX).then_some(index)
                })
                .collect()
        })
        .collect();
    matcher.clear();
    for (index, &mate) in blossom_max_matching(atoms.len(), &compact_adjacency)
        .iter()
        .enumerate()
    {
        if mate != usize::MAX {
            matcher.mate[atoms[index].0 as usize] = atoms[mate].0;
        }
    }
}

/// Per-atom adjacency over aromatic bonds joining two `eligible` atoms, each
/// list in aromatic-bond order (atom1's entry, then atom2's).
fn matching_adjacency(mol: &Molecule, aromatic_bonds: &[BondIdx], eligible: &[bool]) -> MatchAdj {
    let n = mol.atom_count();
    let joins = |bidx: BondIdx| {
        let bond = mol.bond(bidx);
        (eligible[bond.atom1.0 as usize] && eligible[bond.atom2.0 as usize])
            .then_some((bond.atom1, bond.atom2))
    };
    // Compressed rows, each atom's entries in the order the bonds are listed
    // (the order a per-atom `Vec` push would give), in two allocations.
    let mut start = vec![0u32; n + 1];
    for &bidx in aromatic_bonds {
        if let Some((a, b)) = joins(bidx) {
            start[a.0 as usize + 1] += 1;
            start[b.0 as usize + 1] += 1;
        }
    }
    for i in 0..n {
        start[i + 1] += start[i];
    }
    let mut fill = start.clone();
    let mut entries = vec![(AtomIdx(0), BondIdx(0)); start[n] as usize];
    for &bidx in aromatic_bonds {
        if let Some((a, b)) = joins(bidx) {
            entries[fill[a.0 as usize] as usize] = (b, bidx);
            fill[a.0 as usize] += 1;
            entries[fill[b.0 as usize] as usize] = (a, bidx);
            fill[b.0 as usize] += 1;
        }
    }
    MatchAdj { start, entries }
}

/// Per-atom matching adjacency in compressed rows; `adj[atom]` is that
/// atom's `(neighbour, bond)` slice.
struct MatchAdj {
    start: Vec<u32>,
    entries: Vec<(AtomIdx, BondIdx)>,
}

impl std::ops::Index<usize> for MatchAdj {
    type Output = [(AtomIdx, BondIdx)];
    fn index(&self, atom: usize) -> &Self::Output {
        &self.entries[self.start[atom] as usize..self.start[atom + 1] as usize]
    }
}

/// [`build_kekule_result`] from a per-atom mate array.
fn kekule_result_from_mates(
    aromatic_bonds: &[BondIdx],
    mol: &Molecule,
    mate: &[u32],
) -> KekuleResult {
    let mut double = vec![false; mol.bond_count()];
    for (atom, &partner) in mate.iter().enumerate() {
        if partner == NO_MATE || atom as u32 >= partner {
            continue;
        }
        if let Some((bidx, _)) = mol.bond_between(AtomIdx(atom as u32), AtomIdx(partner))
            && mol.bond(bidx).order == BondOrder::Aromatic
        {
            double[bidx.0 as usize] = true;
        }
    }
    aromatic_bonds
        .iter()
        .map(|&bidx| {
            let order = if double[bidx.0 as usize] {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            (bidx, order)
        })
        .collect()
}

/// Build the KekuleResult map from the current matching.
pub fn build_kekule_result(
    aromatic_bonds: &[BondIdx],
    mol: &Molecule,
    matching: &HashMap<AtomIdx, AtomIdx>,
) -> KekuleResult {
    let mut double_bonds: HashSet<BondIdx> = HashSet::new();
    for (&atom, &partner) in matching {
        if atom >= partner {
            continue;
        }
        if let Some((bidx, _)) = mol.bond_between(atom, partner)
            && mol.bond(bidx).order == BondOrder::Aromatic
        {
            double_bonds.insert(bidx);
        }
    }
    aromatic_bonds
        .iter()
        .map(|&bidx| {
            let order = if double_bonds.contains(&bidx) {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            (bidx, order)
        })
        .collect()
}

/// Apply a Kekulé result to a molecule, returning a new Molecule with updated bond orders.
///
/// Aromatic flags on atoms are *not* cleared (the molecule retains the aromaticity
/// annotation; only bond orders change).
///
/// Preserves every stereo side-channel (`stereo_groups`, `stereo_neighbor_order`,
/// `bond_directions`) unchanged. This matters even though `atom.chirality` itself is
/// copied verbatim with each `Atom`: `@`/`@@` is a relative descriptor over the SMILES
/// neighbor-encounter order recorded in `stereo_neighbor_order`, not an absolute one.
/// Losing that order and re-deriving `@`/`@@` output from a *different* neighbor
/// ordering downstream (e.g. canonical atom order) can serialize the *same* tetrahedral
/// center as the wrong configuration -- silently, with no panic and no error, since
/// `chirality` alone still round-trips as "some chirality is set". Atom/bond indices are
/// unchanged (nothing is skipped, added, or reordered), the same precondition the
/// `copy_*_from` methods themselves require.
pub fn apply_kekule(mol: &Molecule, kekule: &KekuleResult) -> Molecule {
    use crate::molecule::MoleculeBuilder;

    if kekule.is_empty() {
        return mol.clone();
    }

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

    builder.copy_stereo_groups_from(mol);
    builder.copy_stereo_from(mol);
    builder.copy_atom_tags_from(mol);
    builder.copy_bond_directions_from(mol);

    builder.build()
}

/// No partner in [`Matcher::mate`].
const NO_MATE: u32 = u32::MAX;

/// Greedy augmenting-path matcher over per-atom adjacency lists.
struct Matcher {
    /// `mate[a]` = matched partner of atom `a`, or [`NO_MATE`].
    mate: Vec<u32>,
    /// BFS scratch: `seen[a] == epoch` marks atoms visited in the current search.
    seen: Vec<u32>,
    epoch: u32,
    parent: Vec<u32>,
    queue: Vec<u32>,
}

impl Matcher {
    fn new(n: usize) -> Self {
        Self {
            mate: vec![NO_MATE; n],
            seen: vec![0; n],
            epoch: 0,
            parent: vec![NO_MATE; n],
            queue: Vec::new(),
        }
    }

    fn clear(&mut self) {
        self.mate.iter_mut().for_each(|m| *m = NO_MATE);
    }

    fn all_matched(&self, atoms: &[AtomIdx]) -> bool {
        atoms.iter().all(|a| self.mate[a.0 as usize] != NO_MATE)
    }

    /// Run a single augmenting-path pass: for each unmatched atom (in the
    /// given order), try to find an augmenting path from it.
    fn run_pass(&mut self, atoms: impl Iterator<Item = AtomIdx>, adj: &MatchAdj) {
        for start in atoms {
            if self.mate[start.0 as usize] != NO_MATE {
                continue;
            }
            self.augment(start.0, adj);
        }
    }

    /// Attempt to find an augmenting path starting from `start` and update
    /// the matching. Iterative BFS with parent-pointer path reconstruction
    /// (no recursion: wasm32's default stack is ~1 MB). Returns whether
    /// `start` is matched afterwards.
    fn augment(&mut self, start: u32, adj: &MatchAdj) -> bool {
        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.seen.iter_mut().for_each(|s| *s = 0);
            self.epoch = 1;
        }
        let epoch = self.epoch;
        self.seen[start as usize] = epoch;
        self.queue.clear();
        self.queue.push(start);
        let mut head = 0;

        'bfs: while head < self.queue.len() {
            let v = self.queue[head];
            head += 1;
            for &(u, _) in &adj[v as usize] {
                let u = u.0;
                if self.seen[u as usize] == epoch {
                    continue;
                }
                self.seen[u as usize] = epoch;
                self.parent[u as usize] = v;

                let partner = self.mate[u as usize];
                if partner == NO_MATE {
                    // Found a free vertex — trace back through parent pointers
                    // and flip every edge along the augmenting path.
                    let mut cur = u;
                    loop {
                        let prev = self.parent[cur as usize];
                        let prev_old_match = self.mate[prev as usize];
                        self.mate[prev as usize] = cur;
                        self.mate[cur as usize] = prev;
                        if prev == start || prev_old_match == NO_MATE {
                            break;
                        }
                        cur = prev_old_match;
                    }
                    break 'bfs;
                }
                // u is matched to partner — explore further from partner.
                if self.seen[partner as usize] != epoch {
                    self.seen[partner as usize] = epoch;
                    self.parent[partner as usize] = u;
                    self.queue.push(partner);
                }
            }
        }

        // Augmentation succeeded iff start is now matched.
        self.mate[start as usize] != NO_MATE
    }
}

// ─── Edmonds' blossom maximum matching ───────────────────────────────────────
//
// General (non-bipartite) maximum matching via Gabow's blossom formulation.
// Vertices are integers 0..n; the returned `mate[i] = j` if matched, else NONE.
//
// `NONE` sentinel: usize::MAX — safe because `n` is always < 32 k for drug-like
// molecules (InChI library limit) and we never index at NONE.

const NONE: usize = usize::MAX;

/// Find a maximum matching in a general graph (Edmonds' blossom, O(n²m)).
fn blossom_max_matching(n: usize, adj: &[Vec<usize>]) -> Vec<usize> {
    let mut mate = vec![NONE; n];
    for v in 0..n {
        if mate[v] == NONE {
            blossom_augment(v, n, adj, &mut mate);
        }
    }
    mate
}

/// Attempt to augment the matching from free vertex `root`.
fn blossom_augment(root: usize, n: usize, adj: &[Vec<usize>], mate: &mut [usize]) {
    // base[v]: representative of the blossom containing v.
    let mut base: Vec<usize> = (0..n).collect();
    // parent[v]: predecessor of v on the augmenting path (NONE = unlabeled).
    let mut parent: Vec<usize> = vec![NONE; n];
    // is_outer[v]: true if v is an outer (even-level) vertex in the BFS forest.
    let mut is_outer: Vec<bool> = vec![false; n];

    is_outer[root] = true;
    let mut queue: VecDeque<usize> = VecDeque::new();
    queue.push_back(root);

    'bfs: while let Some(v) = queue.pop_front() {
        for &w in &adj[v] {
            if base[v] == base[w] {
                continue;
            } // same blossom
            if mate[v] == w {
                continue;
            } // already-matched edge, skip

            if is_outer[w] {
                // Both v and w are outer → odd cycle (blossom).
                let b = blossom_lca(v, w, &base, &parent, mate, n);
                blossom_mark_path(
                    v,
                    b,
                    w,
                    &mut base,
                    &mut parent,
                    &mut is_outer,
                    &mut queue,
                    mate,
                    n,
                );
                blossom_mark_path(
                    w,
                    b,
                    v,
                    &mut base,
                    &mut parent,
                    &mut is_outer,
                    &mut queue,
                    mate,
                    n,
                );
            } else if parent[w] == NONE {
                // w is unlabeled.
                parent[w] = v;
                if mate[w] == NONE {
                    // Augmenting path ends at w.  Trace parent[] and flip matching.
                    let mut cur = w;
                    while cur != NONE {
                        let prev = parent[cur];
                        let prev_old = mate[prev];
                        mate[cur] = prev;
                        mate[prev] = cur;
                        cur = prev_old;
                    }
                    break 'bfs;
                }
                // w is matched — add mate[w] as the next outer vertex.
                let u = mate[w];
                if !is_outer[u] {
                    is_outer[u] = true;
                    parent[u] = w;
                    queue.push_back(u);
                }
            }
            // else: w is inner (labeled but not outer) → skip.
        }
    }
}

/// Lowest common ancestor of `a` and `b` in the alternating BFS tree.
///
/// Traces both paths toward `root` (following outer→matched-inner→outer chains)
/// and returns the first vertex visited by both traces.
fn blossom_lca(
    mut a: usize,
    mut b: usize,
    base: &[usize],
    parent: &[usize],
    mate: &[usize],
    n: usize,
) -> usize {
    let mut visited = vec![false; n];
    loop {
        a = base[a];
        visited[a] = true;
        if mate[a] == NONE {
            break;
        } // reached a free vertex (root or its base)
        a = parent[mate[a]]; // hop: outer→matched inner→outer parent
    }
    loop {
        b = base[b];
        if visited[b] {
            return b;
        } // first vertex seen by both traces
        b = parent[mate[b]];
    }
}

/// Walk from `x` toward blossom base `b`, updating `base[]`, `parent[]`,
/// and promoting inner vertices to outer so the BFS can traverse the blossom.
#[allow(clippy::too_many_arguments)]
fn blossom_mark_path(
    mut x: usize,
    b: usize,
    child: usize,
    base: &mut [usize],
    parent: &mut [usize],
    is_outer: &mut [bool],
    queue: &mut VecDeque<usize>,
    mate: &[usize],
    n: usize,
) {
    let mut ch = child;
    while base[x] != b {
        let bx = base[x];
        let bmx = base[mate[x]];
        // Merge blossom: all vertices in bx or bmx become part of base b.
        for slot in base.iter_mut().take(n) {
            if *slot == bx || *slot == bmx {
                *slot = b;
            }
        }
        // Update augmenting-path parent pointers inside the blossom.
        parent[x] = ch;
        // Promote mate[x] to outer so the BFS can continue through the blossom.
        let mx = mate[x];
        if !is_outer[mx] {
            is_outer[mx] = true;
            queue.push_back(mx);
        }
        ch = mx;
        x = parent[mx];
    }
}

// ─────────────────────────────────────────────────────────────────────────────

/// Determine whether an aromatic atom *must* appear in the matching
/// (i.e. requires a double bond for a valid Kekulé form).
///
/// Charge-aware per RDKit's own kekulization (`markDbondCands`,
/// `Code/GraphMol/Kekulize.cpp` at the pinned commit `8afba32ec539dcb2369bc84549d802aca3f7eb39`):
/// RDKit computes a per-atom target valence `dv = defaultValence(atomicNum) + chrg`
/// (formal charge added directly, no sign flip — `isEarlyAtom` is `false` for
/// C/N/O/P/S/Se/Te, `Code/GraphMol/Atom.cpp`) *except* for carbon, which gets an
/// explicit sign flip (`chrg = -chrg` when `atomicNum == 6 && chrg > 0`, their own
/// comment: "special case for carbon - see GitHub #539"). An atom is a double-bond
/// candidate (`dBndCands[allAtm] = 1`) iff its current bond-order sum `sbo` (ring
/// bonds + H count, each counted as 1) is exactly `dv - 1`; i.e. it needs exactly one
/// more bond order to reach `dv`. This flips the "donor" classification for a charged
/// heteroatom (charge *raises* `dv`, e.g. `[nH+]`/`[o+]` become one bond order short —
/// a candidate) but *lowers* it for a cationic carbon (a `+1` carbon's `dv` drops from
/// 4 to 3, which its existing bonds already satisfy — not a candidate). Verified
/// empirically against `rdkit==2026.03.3` (pinned to the same commit) for every case
/// below; see `docs/rfcs/kekulize_charge_aware_rdkit_parity.md`.
///
/// An atom can be unmatched (a lone-pair donor) if:
/// - O, S, Se, Te (furan/thiophene/selenophene/tellurophene-type chalcogen) — but only
///   while neutral or anionic; a positively-charged chalcogen (pyrylium's `[o+]`) has
///   consumed the lone pair and needs a double bond instead, like pyridine-type N.
/// - N or P with an H ([nH] pyrrole-type nitrogen, [pH] phosphole-type phosphorus) —
///   but only while neutral or anionic; protonation (`[nH+]`, e.g. pyridinium) consumes
///   the lone pair the same way, so the atom needs a double bond instead.
/// - Any anionic aromatic atom (cyclopentadienyl-anion-type) donates its lone pair.
/// - Any aromatic atom that already has an exocyclic double bond (e.g. the
///   carbonyl carbon in coumarin/warfarin `c=O` fused into an aromatic ring).
///   Such an atom's pi contribution comes from conjugation with the exocyclic
///   bond; no additional ring double bond is needed or possible.
///
/// A cationic aromatic carbon ([cH+], e.g. tropylium) is a symmetric case: an empty
/// p-orbital electron *acceptor*, like aromatic B, but — unlike B — needs no double
/// bond of its own (see the carbon sign-flip above).
///
/// Everything else (C, N/P without H like pyridine/phosphole's ring carbons, aromatic
/// B) must be matched.
pub fn atom_must_be_matched(mol: &Molecule, idx: AtomIdx) -> bool {
    let atom = mol.atom(idx);
    match atom.element.atomic_number() {
        // O, S, Se, Te donate a lone pair → don't need a double bond — *unless* a
        // positive charge (pyrylium's [o+]) has consumed it, in which case the atom
        // needs a double bond exactly like pyridine-type N (falls through to the
        // catch-all `_ => true` below).
        8 | 16 | 34 | 52 if atom.charge <= 0 => false,
        // A cation with a third bond (`[s+]([O-])` of a 1,2,5-thiadiazole
        // S-oxide) already has RDKit's target valence 3 and keeps a lone pair;
        // only the two-bonded one (thiopyrylium, pyrylium) needs a double bond.
        8 | 16 | 34 | 52 if atom.charge == 1 => {
            let bonds: i32 = mol
                .neighbors(idx)
                .map(|(_, b)| match mol.bond(b).order {
                    BondOrder::Aromatic => 1,
                    order => i32::from(order.order_int()),
                })
                .sum();
            bonds + i32::from(atom.hydrogen_count.unwrap_or(0)) == 2
        }
        // Aromatic B contributes an empty p orbital (electron acceptor), not a lone pair.
        // It must appear in a double bond in the Kekulé form, just like aromatic C.
        5 => true,
        // N or P with explicit H ([nH], [pH]) is a lone-pair donor — *unless* protonated
        // (a positive charge consumes the lone pair, e.g. pyridinium's [nH+], which then
        // needs a double bond exactly like neutral pyridine's bare N).
        7 | 15
            if atom.charge <= 0
                && (matches!(atom.hydrogen_count, Some(h) if h > 0)
                    || mol
                        .neighbors(idx)
                        .any(|(neighbor, _)| mol.atom(neighbor).element.atomic_number() == 1)) =>
        {
            false
        }
        // Anionic aromatic N/P ([n-]) also donates its lone pair; the extra electron
        // occupies the lone-pair slot rather than a ring π bond (same as [nH]).
        7 | 15 if atom.charge < 0 => false,
        // Neutral N/P with a non-aromatic substituent (e.g. N-methyl in caffeine) is a lone-pair
        // donor: the substituent "replaces" the H, and the atom contributes its lone pair to
        // aromaticity rather than a π bond (same as [nH] in pyrrole).
        // Charged N/P (pyridinium [n+], N-oxide) must still be matched even with a substituent.
        7 | 15
            if atom.charge == 0
                && mol.neighbors(idx).any(|(neighbor, bidx)| {
                    mol.atom(neighbor).element.atomic_number() != 1
                        && mol.bond(bidx).order != BondOrder::Aromatic
                }) =>
        {
            false
        }
        // A cation one bond order short of RDKit's target valence 4 is a
        // candidate (pyridinium `[nH+]`, `C[n+]1ccccc1`); one that already
        // reaches it (`C[n+]1(C)cccc1`, `[nH2+]1cccc1`) is not, as in RDKit.
        7 | 15 if atom.charge == 1 => {
            let bonds: i32 = mol
                .neighbors(idx)
                .map(|(_, b)| match mol.bond(b).order {
                    BondOrder::Aromatic => 1,
                    order => i32::from(order.order_int()),
                })
                .sum();
            // Bonds to explicit H atoms are already counted above.
            bonds + i32::from(atom.hydrogen_count.unwrap_or(0)) == 3
        }
        // Bare aromatic N/P with only aromatic bonds (pyridine-type): must be
        // matched.
        7 | 15 => true,
        // A cationic aromatic carbon ([cH+], e.g. tropylium) has an empty p orbital: an
        // electron acceptor like aromatic B, but *without* B's obligatory double bond.
        // RDKit's own charge-sign flip (see doc comment above) drops a +1 carbon's
        // target valence from 4 to 3, which its existing ring bonds + H already
        // satisfy — no double bond needed or possible.
        // A charged carbon without that H (`[c+]1ccccc1`, `[c-]1ccccc1`: a
        // ring bond short of the target valence 3) is a candidate, as in
        // RDKit's `sbo == dv - 1`; `[cH-]` (cyclopentadienyl) is not.
        6 if atom.charge != 0 => {
            let bonds: i32 = mol
                .neighbors(idx)
                .map(|(_, b)| match mol.bond(b).order {
                    BondOrder::Aromatic => 1,
                    order => i32::from(order.order_int()),
                })
                .sum();
            let h = i32::from(atom.hydrogen_count.unwrap_or(0));
            bonds + h == 3 - i32::from(atom.charge).abs()
        }
        // Any anionic aromatic atom (e.g. cyclopentadienyl [cH-]) donates its lone pair.
        _ if atom.charge < 0 => false,
        // Any atom (typically C) that already has an exocyclic π bond (e.g. `c(=O)` in
        // a heterocyclic carbonyl) cannot also carry a ring double bond — its π-slot is
        // occupied by the exocyclic bond. Such atoms contribute via conjugation, like
        // lone-pair donors, and should not be forced into the matching.
        _ if mol.neighbors(idx).any(|(_, bidx)| {
            let o = mol.bond(bidx).order;
            o == BondOrder::Double || o == BondOrder::Triple
        }) =>
        {
            false
        }
        // All other atoms must appear in the matching.
        _ => true,
    }
}

#[cfg(test)]
include!("kekulization_tests.rs");
