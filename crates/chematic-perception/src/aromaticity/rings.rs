// ---------------------------------------------------------------------------
// Ring augmentation (XOR sub-rings)
// ---------------------------------------------------------------------------

/// Return the sorted set of bond indices that form `ring`.
fn ring_bond_set(mol: &Molecule, ring: &[AtomIdx]) -> Vec<BondIdx> {
    let n = ring.len();
    let mut bonds: Vec<BondIdx> = (0..n)
        .filter_map(|i| {
            let a = ring[i];
            let b = ring[(i + 1) % n];
            mol.bond_between(a, b).map(|(bidx, _)| bidx)
        })
        .collect();
    bonds.sort();
    bonds
}

/// Sorted symmetric difference of two sorted slices.
fn bond_sym_diff(a: &[BondIdx], b: &[BondIdx]) -> Vec<BondIdx> {
    let mut result: Vec<BondIdx> = Vec::new();
    let mut i = 0;
    let mut j = 0;
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                result.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                result.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    result.extend_from_slice(&a[i..]);
    result.extend_from_slice(&b[j..]);
    result
}

/// Reconstruct an ordered atom sequence from a set of bond indices forming a simple cycle.
/// Returns `None` if the bonds do not form a valid simple cycle.
fn ring_atoms_from_bond_set(mol: &Molecule, bonds: &[BondIdx]) -> Option<Vec<AtomIdx>> {
    if bonds.is_empty() {
        return None;
    }
    let mut adj: FxHashMap<AtomIdx, [Option<AtomIdx>; 2]> = FxHashMap::default();
    for &bidx in bonds {
        let bond = mol.bond(bidx);
        for (a, b) in [(bond.atom1, bond.atom2), (bond.atom2, bond.atom1)] {
            let e = adj.entry(a).or_insert([None; 2]);
            if e[0].is_none() {
                e[0] = Some(b);
            } else if e[1].is_none() {
                e[1] = Some(b);
            } else {
                return None; // degree > 2 — not a simple ring
            }
        }
    }
    // All atoms must have exactly 2 neighbours.
    if adj.values().any(|e| e[1].is_none()) {
        return None;
    }
    let start = *adj.keys().next()?;
    let mut path = vec![start];
    let mut prev = start;
    let mut current = adj[&start][0]?;
    while current != start {
        path.push(current);
        let [n0, n1] = adj[&current];
        let next = if n0 == Some(prev) { n1? } else { n0? };
        prev = current;
        current = next;
    }
    if path.len() != bonds.len() {
        return None;
    }
    Some(path)
}

/// Augment the SSSR ring list with smaller XOR sub-rings found by pairwise GF(2)
/// differences between SSSR rings that share atoms.
///
/// The standard SSSR algorithm sometimes stores a large fundamental cycle rather
/// than its smaller GF(2)-reduced equivalent (e.g. the 5-ring of indolizine is
/// the XOR of the 6-ring and the 9-ring the algorithm reports).
/// This augmentation adds such missing smaller rings so that aromaticity
/// perception works on the correct smallest rings without modifying the SSSR.
///
/// The returned `Vec` starts with all SSSR rings in their original order; any
/// additional sub-rings derived by GF(2) pairwise XOR follow.  The function
/// only adds a ring if it is strictly smaller than *both* parents, ensuring
/// that envelope rings (e.g. the 10-membered perimeter of naphthalene) are
/// never introduced.
pub fn augmented_ring_set(mol: &Molecule, sssr_rings: &[Vec<AtomIdx>]) -> Vec<Vec<AtomIdx>> {
    let mut rings: Vec<Vec<AtomIdx>> = sssr_rings.to_vec();

    // Track which atom-sets we already have (as sorted atom lists).
    let mut known: FxHashSet<Vec<AtomIdx>> = sssr_rings
        .iter()
        .map(|r| {
            let mut s = r.clone();
            s.sort();
            s
        })
        .collect();

    // Iterative pairwise XOR until convergence.
    //
    // A single pass only finds rings that are the XOR of two SSSR rings.
    // Iterating also finds rings that require XOR of 3+ SSSR rings
    // (e.g. the inner hexagon of coronene, or sub-rings in multi-step
    // fused PAHs where the SSSR chose large perimeter cycles).
    // Termination is guaranteed because each new ring is strictly smaller
    // than both of its parents, so ring size can only decrease.
    loop {
        let mut changed = false;
        let n = rings.len();
        let bond_sets: Vec<Vec<BondIdx>> = rings.iter().map(|r| ring_bond_set(mol, r)).collect();

        for i in 0..n {
            for j in (i + 1)..n {
                // Only consider pairs that share atoms (fused rings).
                let shares_atom = rings[i].iter().any(|a| rings[j].contains(a));
                if !shares_atom {
                    continue;
                }
                let xor_bonds = bond_sym_diff(&bond_sets[i], &bond_sets[j]);
                if xor_bonds.is_empty() {
                    continue;
                }
                // Only interesting if the XOR ring is not larger than the larger
                // parent.  Using max() recovers cases where SSSR chose a large
                // cycle (e.g. 10-ring macro vs 6-ring benzene twin).
                // Using `>` (not `>=`) also allows same-size XOR rings, which
                // handles bridged bicyclics (e.g. tropane or dioxolane spirocycles)
                // where both parent rings are 6-membered and the missing bridge
                // ring is also 6-membered.  Termination is still guaranteed:
                // the `known` set prevents duplicates, and a finite molecule has
                // finitely many valid cycles.
                if xor_bonds.len() > rings[i].len().max(rings[j].len()) {
                    continue;
                }
                if let Some(new_ring) = ring_atoms_from_bond_set(mol, &xor_bonds) {
                    let mut key = new_ring.clone();
                    key.sort();
                    if known.insert(key) {
                        rings.push(new_ring);
                        changed = true;
                    }
                }
            }
        }

        // 3-ring XOR: catches small rings that require XOR of 3 SSSR rings
        // when no intermediate 2-ring XOR produces a valid smaller ring.
        for i in 0..n {
            for j in (i + 1)..n {
                let shares_ij = rings[i].iter().any(|a| rings[j].contains(a));
                if !shares_ij {
                    continue;
                }
                let xor_ij = bond_sym_diff(&bond_sets[i], &bond_sets[j]);
                if xor_ij.is_empty() {
                    continue;
                }
                for k in (j + 1)..n {
                    let shares_k = rings[k]
                        .iter()
                        .any(|a| rings[i].contains(a) || rings[j].contains(a));
                    if !shares_k {
                        continue;
                    }
                    let xor_ijk = bond_sym_diff(&xor_ij, &bond_sets[k]);
                    let max_size = rings[i].len().max(rings[j].len()).max(rings[k].len());
                    if xor_ijk.is_empty() || xor_ijk.len() > max_size {
                        continue;
                    }
                    if let Some(new_ring) = ring_atoms_from_bond_set(mol, &xor_ijk) {
                        let mut key = new_ring.clone();
                        key.sort();
                        if known.insert(key) {
                            rings.push(new_ring);
                            changed = true;
                        }
                    }
                }
            }
        }

        if !changed {
            break;
        }
    }

    rings
}

/// Shared inner: SSSR → augmented_ring_set → strip_envelope_rings, no aromaticity filter.
fn all_ring_list_inner(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    let sssr = crate::sssr::find_sssr(mol);
    let aug = augmented_ring_set(mol, sssr.rings());
    if aug.len() <= 1 {
        return aug;
    }
    let bond_sets: Vec<Vec<BondIdx>> = aug.iter().map(|r| ring_bond_set(mol, r)).collect();
    let mut is_envelope = vec![false; aug.len()];
    strip_envelope_rings(&aug, &bond_sets, &mut is_envelope);
    aug.into_iter()
        .zip(is_envelope)
        .filter(|(_, e)| !e)
        .map(|(r, _)| r)
        .collect()
}

/// Return all rings after augmented-ring-set expansion and envelope stripping.
///
/// Same pipeline as [`aromatic_ring_list`] but with no aromaticity filter — useful
/// for aliphatic/saturated ring counting and bridgehead detection where SSSR
/// envelope rings cause over-counting.
pub fn all_ring_list(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    all_ring_list_inner(mol)
}

/// True when all ring bonds between ring atoms are `BondOrder::Aromatic`.
///
/// Rings written with aromatic-SMILES notation but containing an explicit single
/// bond (`c-n`, `nc-2`, etc.) are NOT truly aromatic.  RDKit canonicalises such
/// SMILES with lowercase atoms and a `-` bond, which the parser stores as
/// `BondOrder::Single` between two aromatic-flagged atoms.  Returning `false`
/// here lets callers exclude them from the aromatic ring count.
pub fn ring_bonds_all_aromatic(mol: &Molecule, ring: &[AtomIdx]) -> bool {
    let n = ring.len();
    (0..n).all(|i| {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        mol.bond_between(a, b)
            .map(|(bidx, _)| mol.bond(bidx).order == BondOrder::Aromatic)
            .unwrap_or(true)
    })
}

/// Return the de-duplicated list of aromatic rings after augmented-ring-set expansion
/// and envelope stripping.  Useful for filtering (e.g. counting only aromatic heterocycles).
pub fn aromatic_ring_list(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    let mol_with_arom;
    let has_aromatic_atoms = mol.atoms().any(|(_, a)| a.aromatic);
    let mol = if !has_aromatic_atoms {
        mol_with_arom = apply_aromaticity(mol);
        &mol_with_arom
    } else {
        mol
    };
    aromatic_ring_list_preperceived(mol)
}

/// Return aromatic rings from a molecule whose aromaticity flags are already
/// authoritative, without running another perception model.
///
/// This is used by compatibility profiles that deliberately apply a narrower
/// aromaticity model which can validly produce zero aromatic atoms. Calling
/// [`aromatic_ring_list`] in that case would otherwise interpret the absence
/// of flags as "not perceived yet" and replace the chosen model.
pub fn aromatic_ring_list_preperceived(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    let mut rings = all_ring_list_inner(mol);
    append_small_chordless_aromatic_cycles(mol, &mut rings);
    rings
        .into_iter()
        .filter(|ring| {
            ring.iter().all(|&idx| mol.atom(idx).aromatic) && ring_bonds_all_aromatic(mol, ring)
        })
        .collect()
}

/// Add bounded, chordless aromatic cycles that a minimum cycle basis can omit.
///
/// A minimum basis is sufficient for many ring operations, but it can choose a
/// five-membered cage face instead of a symmetry-equivalent six-membered
/// aromatic face.  RDKit's aromatic-ring count considers the latter as well.
/// Enumerating every simple cycle would be unbounded on polycyclic graphs, so
/// this helper only considers induced aromatic cycles of size 3..=6 and fails
/// closed when the candidate cap is exceeded.  Atom-set deduplication keeps the
/// result independent of traversal direction.
fn append_small_chordless_aromatic_cycles(mol: &Molecule, rings: &mut Vec<Vec<AtomIdx>>) {
    const MAX_CANDIDATES: usize = 4096;
    const MAX_SIZE: usize = 6;
    let mut known: FxHashSet<Vec<u32>> = rings
        .iter()
        .filter(|ring| ring.len() <= MAX_SIZE)
        .map(|ring| {
            let mut key: Vec<u32> = ring.iter().map(|atom| atom.0).collect();
            key.sort_unstable();
            key
        })
        .collect();
    let mut candidates = Vec::new();
    for start in 0..mol.atom_count() {
        let start = AtomIdx(start as u32);
        if !mol.atom(start).aromatic {
            continue;
        }
        let mut search = ChordlessCycleSearch::new(mol, start, MAX_SIZE, MAX_CANDIDATES);
        search.collect(start);
        candidates.extend(search.cycles);
        if candidates.len() > MAX_CANDIDATES {
            return;
        }
    }
    for ring in candidates {
        let mut key: Vec<u32> = ring.iter().map(|atom| atom.0).collect();
        key.sort_unstable();
        if known.insert(key) {
            rings.push(ring);
        }
    }
}

struct ChordlessCycleSearch<'a> {
    mol: &'a Molecule,
    start: AtomIdx,
    path: Vec<AtomIdx>,
    visited: FxHashSet<AtomIdx>,
    cycles: Vec<Vec<AtomIdx>>,
    max_size: usize,
    max_candidates: usize,
}

impl<'a> ChordlessCycleSearch<'a> {
    fn new(mol: &'a Molecule, start: AtomIdx, max_size: usize, max_candidates: usize) -> Self {
        let mut visited = FxHashSet::default();
        visited.insert(start);
        Self {
            mol,
            start,
            path: vec![start],
            visited,
            cycles: Vec::new(),
            max_size,
            max_candidates,
        }
    }

    fn collect(&mut self, current: AtomIdx) {
        if self.cycles.len() > self.max_candidates {
            return;
        }
        for (next, bond) in self.mol.neighbors(current) {
            if self.mol.bond(bond).order != BondOrder::Aromatic || !self.mol.atom(next).aromatic {
                continue;
            }
            if next == self.start {
                if (3..=self.max_size).contains(&self.path.len())
                    && chordless_cycle(self.mol, &self.path)
                {
                    self.cycles.push(self.path.clone());
                }
                continue;
            }
            if next.0 < self.start.0
                || self.visited.contains(&next)
                || self.path.len() >= self.max_size
            {
                continue;
            }
            self.visited.insert(next);
            self.path.push(next);
            self.collect(next);
            self.path.pop();
            self.visited.remove(&next);
        }
    }
}

fn chordless_cycle(mol: &Molecule, cycle: &[AtomIdx]) -> bool {
    for i in 0..cycle.len() {
        for j in (i + 1)..cycle.len() {
            let adjacent = j == i + 1 || (i == 0 && j + 1 == cycle.len());
            if !adjacent && mol.bond_between(cycle[i], cycle[j]).is_some() {
                return false;
            }
        }
    }
    true
}

/// Mark which rings in `aromatic` are GF(2) sums (bond-XOR) of 2–4 smaller rings.
fn strip_envelope_rings(
    aromatic: &[Vec<AtomIdx>],
    bond_sets: &[Vec<BondIdx>],
    is_envelope: &mut [bool],
) {
    for (index, ring) in aromatic.iter().enumerate() {
        let smaller: Vec<usize> = aromatic
            .iter()
            .enumerate()
            .filter_map(|(other, candidate)| {
                (other != index && candidate.len() < ring.len()).then_some(other)
            })
            .collect();
        is_envelope[index] = (2..=4).any(|terms| {
            ring_xor_combination_matches(&bond_sets[index], bond_sets, &smaller, terms)
        });
    }
}

fn ring_xor_combination_matches(
    target: &[BondIdx],
    bond_sets: &[Vec<BondIdx>],
    candidates: &[usize],
    terms: usize,
) -> bool {
    fn search(
        target: &[BondIdx],
        bond_sets: &[Vec<BondIdx>],
        candidates: &[usize],
        terms: usize,
        start: usize,
        accumulated: &[BondIdx],
    ) -> bool {
        if terms == 0 {
            return accumulated == target;
        }
        (start..=candidates.len().saturating_sub(terms)).any(|position| {
            let bonds = &bond_sets[candidates[position]];
            let next = if accumulated.is_empty() {
                bonds.clone()
            } else {
                bond_sym_diff(accumulated, bonds)
            };
            search(
                target,
                bond_sets,
                candidates,
                terms - 1,
                position + 1,
                &next,
            )
        })
    }

    candidates.len() >= terms
        && search(target, bond_sets, candidates, terms, 0, &[])
}

pub fn count_aromatic_rings(mol: &Molecule) -> usize {
    let perceived = aromatic_ring_input(mol);
    let Some(mol) = perceived.as_ref().or(Some(mol)) else {
        unreachable!();
    };
    if perceived.is_none()
        && !mol.atoms().any(|(_, atom)| atom.aromatic)
        && crate::sssr::sssr_ring_count(mol) == 0
    {
        return 0;
    }
    let mut components = RingComponentSummary::new(mol);
    let (simple_count, mut keep) = components.count_simple_components(mol);
    if !keep.iter().any(|&value| value) {
        return simple_count;
    }
    components.expand_component_mask(&mut keep);
    simple_count + count_complex_aromatic_rings(mol, &keep)
}

fn aromatic_ring_input(mol: &Molecule) -> Option<Molecule> {
    if mol.atoms().any(|(_, atom)| atom.aromatic) {
        return None;
    }
    if crate::sssr::sssr_ring_count(mol) == 0 {
        return None;
    }
    let perceived = apply_aromaticity(mol);
    if perceived.atom_count() == mol.atom_count() && perceived.bond_count() == mol.bond_count() {
        crate::sssr::seed_ring_data_from(&perceived, mol);
    }
    Some(perceived)
}

struct RingComponentSummary {
    cyclic: std::sync::Arc<Vec<bool>>,
    parent: Vec<u32>,
    ring_atom: Vec<bool>,
    edges: Vec<u32>,
    atoms: Vec<u32>,
    aromatic_atoms: Vec<u32>,
}

impl RingComponentSummary {
    fn new(mol: &Molecule) -> Self {
        let size = mol.atom_count();
        let cyclic = crate::sssr::ring_bond_flags_shared(mol);
        let mut summary = Self {
            cyclic,
            parent: (0..size as u32).collect(),
            ring_atom: vec![false; size],
            edges: vec![0; size],
            atoms: vec![0; size],
            aromatic_atoms: vec![0; size],
        };
        summary.connect(mol);
        summary.collect_counts(mol);
        summary
    }

    fn root(&mut self, mut atom: u32) -> u32 {
        while self.parent[atom as usize] != atom {
            let parent = self.parent[atom as usize];
            self.parent[atom as usize] = self.parent[parent as usize];
            atom = parent;
        }
        atom
    }

    fn connect(&mut self, mol: &Molecule) {
        for (index, bond) in mol.bonds() {
            if !self.cyclic[index.0 as usize] {
                continue;
            }
            self.ring_atom[bond.atom1.0 as usize] = true;
            self.ring_atom[bond.atom2.0 as usize] = true;
            let (left, right) = (self.root(bond.atom1.0), self.root(bond.atom2.0));
            if left != right {
                self.parent[left as usize] = right;
            }
        }
    }

    fn collect_counts(&mut self, mol: &Molecule) {
        for (index, bond) in mol.bonds() {
            if self.cyclic[index.0 as usize] {
                let root = self.root(bond.atom1.0) as usize;
                self.edges[root] += 1;
            }
        }
        for atom in 0..self.ring_atom.len() {
            if self.ring_atom[atom] {
                let root = self.root(atom as u32) as usize;
                self.atoms[root] += 1;
                self.aromatic_atoms[root] +=
                    u32::from(mol.atom(AtomIdx(atom as u32)).aromatic);
            }
        }
    }

    fn count_simple_components(&mut self, mol: &Molecule) -> (usize, Vec<bool>) {
        let size = self.parent.len();
        let root_of: Vec<u32> = (0..size as u32)
            .map(|atom| {
                if self.ring_atom[atom as usize] {
                    self.root(atom)
                } else {
                    u32::MAX
                }
            })
            .collect();
        let mut count = 0;
        let mut keep = vec![false; size];
        for (root, keep_root) in keep.iter_mut().enumerate() {
            if !self.is_populated_root(root) {
                continue;
            }
            if self.edges[root] == self.atoms[root] {
                count += usize::from(self.aromatic_atoms[root] == self.atoms[root]);
            } else if self.edges[root] == self.atoms[root] + 1 {
                count += self.count_rank_two_component(mol, &root_of, root as u32);
            } else {
                *keep_root = true;
            }
        }
        (count, keep)
    }

    fn is_populated_root(&mut self, root: usize) -> bool {
        self.ring_atom[root]
            && self.root(root as u32) as usize == root
            && self.aromatic_atoms[root] > 0
    }

    fn count_rank_two_component(
        &self,
        mol: &Molecule,
        root_of: &[u32],
        root: u32,
    ) -> usize {
        let mut rings = Vec::new();
        if !crate::sssr::unique_basis_of_small_component(
            mol,
            &self.cyclic,
            root_of,
            root,
            &mut rings,
        ) {
            return 0;
        }
        rings
            .iter()
            .filter(|ring| ring.iter().all(|&atom| mol.atom(atom).aromatic))
            .count()
    }

    fn expand_component_mask(&mut self, keep: &mut [bool]) {
        for atom in 0..self.ring_atom.len() {
            if self.ring_atom[atom] {
                let root = self.root(atom as u32) as usize;
                keep[atom] = keep[root];
            }
        }
    }
}

fn count_complex_aromatic_rings(mol: &Molecule, keep: &[bool]) -> usize {
    let rings = crate::sssr::find_sssr_ring_set_in_components(mol, keep);
    let aromatic: Vec<Vec<AtomIdx>> = augmented_ring_set(mol, &rings)
        .into_iter()
        .filter(|ring| ring.iter().all(|&atom| mol.atom(atom).aromatic))
        .collect();
    if aromatic.len() <= 1 {
        return aromatic.len();
    }
    let bond_sets: Vec<Vec<BondIdx>> = aromatic
        .iter()
        .map(|ring| ring_bond_set(mol, ring))
        .collect();
    let mut is_envelope = vec![false; aromatic.len()];
    strip_envelope_rings(&aromatic, &bond_sets, &mut is_envelope);
    is_envelope.iter().filter(|&&envelope| !envelope).count()
}

// ---------------------------------------------------------------------------
// Per-ring pi electron count
// ---------------------------------------------------------------------------

/// Count pi electrons for a ring atom, returning `None` if the atom is
/// incompatible with aromaticity (e.g. sp3 carbon).
///
/// `aromatic_context`: atoms already confirmed aromatic (from Pass 1 or a
/// previous Pass 2 iteration).  Such atoms contribute 1π unconditionally,
/// without requiring an explicit double bond.
///
/// Rules:
/// - **C**: if already in `aromatic_context` → 1π (confirmed sp2).
///   1. No double bond anywhere: carbanion (`charge == -1`) → 2π (lone pair,
///      e.g. cyclopentadienyl anion); otherwise sp3 → None.
///   2. Two endocyclic double bonds: cumulated ring carbon → None; its
///      orthogonal pi orbitals cannot form one aromatic circuit.
///   3. Has a double bond whose far atom is on NO ring at all (a genuine
///      exocyclic substituent, not a ring-fusion bond into a different ring)
///      and is a more electronegative atom (O/N/S) → 0π (its p-orbital
///      electrons are in the exocyclic π bond, e.g. the carbonyl carbon in
///      tropone/pyridone/pyranone). A double bond whose far atom lies in a
///      DIFFERENT ring (e.g. a fusion carbon whose own Kekule double bond
///      happens to point into the other ring of a fused bicyclic, as in
///      quinazoline/quinoxaline) is a ring bond, not a substituent, and
///      falls through to rule 4 instead — see `all_ring_bonds` below.
///   4. Otherwise (has an endocyclic Double/Aromatic bond, or a double bond
///      into another ring) → 1π.
/// - **N**:
///   1. Has H → 2π (pyrrole-type lone pair).
///   2. Has an explicit `Double` bond → 1π (pyridine-type).
///   3. total_degree == 3 AND ring_degree < total_degree AND no explicit
///      double bond → 2π (lone pair in p orbital): covers both a bridgehead
///      N shared by two fused rings (indolizine) and a substituted
///      pyrrole-type N (N-methylpyrrole, N-glycosylated purine); the overall
///      4n+2 sum, not the substituent, decides ring aromaticity.
///   4. Has in-ring `Aromatic` bond → 1π (pyridine-like aromatic N).
///   5. Already in `aromatic_context` → 1π.
///   6. Otherwise → None.
/// - **O/S**: ring_degree must be 2; contributes 2π (lone pair).
/// - **P (15) / Se (34) / Te (52)**: analogous lone-pair donors; only in
///   [`AromaticityAlgorithm::RdkitLike`] mode.
/// - **Other elements**: None (unsupported).
fn ring_pi_electrons(
    mol: &Molecule,
    ring: &[AtomIdx],
    aromatic_context: &FxHashSet<AtomIdx>,
    algo: AromaticityAlgorithm,
    all_ring_bonds: &FxHashSet<BondIdx>,
) -> Option<u32> {
    let ring_atom_set: FxHashSet<AtomIdx> = ring.iter().copied().collect();
    let mut total_pi: u32 = 0;
    for &atom_idx in ring {
        if aromatic_context.contains(&atom_idx) {
            total_pi += 1;
            continue;
        }
        let (electrons, _) = evaluate_atom_pi_contribution_inner(
            mol,
            atom_idx,
            &ring_atom_set,
            algo,
            all_ring_bonds,
        );
        total_pi += u32::from(electrons?);
    }
    Some(total_pi)
}
