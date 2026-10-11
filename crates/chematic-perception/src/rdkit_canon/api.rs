fn rdkit_bond_type(order: BondOrder) -> u32 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => BT_AROMATIC,
        BondOrder::Dative => 17,
        BondOrder::Zero => 21,
        _ => 0,
    }
}

/// Tetrahedral tag relative to adjacency order, in RDKit's encoding.
fn adjacency_chiral_tag(mol: &Molecule, a: AtomIdx) -> u8 {
    let atom = mol.atom(a);
    let ccw = match atom.chirality {
        Chirality::None => return 0,
        Chirality::CounterClockwise => true,
        Chirality::Clockwise => false,
        _ => return 3,
    };
    let adj: Vec<u32> = mol.neighbors(a).map(|(nb, _)| nb.0).collect();
    let declared: Vec<u32> = match mol.stereo_neighbor_order(a) {
        Some(order) => order.to_vec(),
        None if atom.hydrogen_count.is_some_and(|h| h > 0) => return 3,
        None => adj.clone(),
    };
    if declared.len() != adj.len() || declared.iter().any(|x| !adj.contains(x)) {
        return 3;
    }
    let odd = count_swaps(&adj, &declared) % 2 == 1;
    match (ccw, odd) {
        (true, false) | (false, true) => 2,
        _ => 1,
    }
}

fn rings_of(mol: &Molecule) -> Vec<Vec<usize>> {
    match crate::rdkit_sssr_ring_order(mol) {
        Some(r) => r
            .into_iter()
            .map(|ring| ring.into_iter().map(|a| a.0 as usize).collect())
            .collect(),
        None => crate::find_symmetrized_sssr(mol)
            .rings()
            .iter()
            .map(|ring| ring.iter().map(|a| a.0 as usize).collect())
            .collect(),
    }
}

fn build_input(mol: &Molecule, rings: &[Vec<usize>], ez: &[(BondIdx, bool)]) -> Input {
    let n = mol.atom_count();
    let mut nrings = vec![0u32; n];
    for ring in rings {
        for &a in ring {
            nrings[a] += 1;
        }
    }
    let mut inp = Input {
        anum: Vec::with_capacity(n),
        isotope: Vec::with_capacity(n),
        charge: Vec::with_capacity(n),
        map: Vec::with_capacity(n),
        chiral: Vec::with_capacity(n),
        nrings,
        ring_stereo: vec![false; n],
        nbrs: vec![SmallVec::new(); n],
        bonds: Vec::with_capacity(mol.bond_count()),
        atom_symbols: None,
        bond_symbols: None,
        break_ties: true,
        special_chirality: true,
        in_play: None,
    };
    for (idx, atom) in mol.atoms() {
        inp.anum.push(if atom.wildcard {
            0
        } else {
            u32::from(atom.element.atomic_number())
        });
        inp.isotope.push(atom.isotope.map_or(0, u32::from));
        inp.charge.push(i32::from(atom.charge));
        inp.map.push(atom.atom_map.map_or(0, i32::from));
        inp.chiral.push(adjacency_chiral_tag(mol, idx));
        inp.nbrs[idx.0 as usize] = mol.neighbors(idx).map(|(nb, _)| nb.0).collect();
    }
    for (bidx, bond) in mol.bonds() {
        let stereo = ez
            .iter()
            .find(|(b, _)| *b == bidx)
            .map_or(STEREONONE, |&(_, e)| if e { STEREOE } else { STEREOZ });
        inp.bonds.push((
            bond.atom1.0,
            bond.atom2.0,
            rdkit_bond_type(bond.order),
            stereo,
            None,
        ));
    }
    inp
}

/// RDKit's canonical atom ranks (`Chem.CanonicalRankAtoms(mol, breakTies=True)`)
/// for the invariants chematic carries; see the module documentation.
pub fn rdkit_canonical_atom_ranks(mol: &Molecule) -> Vec<u32> {
    let rings = rings_of(mol);
    rank_with_rings(mol, &rings, &[])
}

/// [`rdkit_canonical_atom_ranks`] with double-bond stereo: `ez` lists the
/// bonds RDKit's stereo perception marks `STEREOE` (`true`) or `STEREOZ`
/// (`false`), which RDKit's bond invariants compare (`CanonicalRankAtoms`
/// on a molecule read from SMILES with `/` `\`). chematic's legacy E/Z
/// assignment gives RDKit's labels on 1,454 of 1,458 stereo bonds of the
/// exposed 10k and ChEMBL 5k corpora (four oxime/ylidene CIP flips).
pub fn rdkit_canonical_atom_ranks_with_bond_stereo(
    mol: &Molecule,
    ez: &[(BondIdx, bool)],
) -> Vec<u32> {
    let rings = rings_of(mol);
    rank_with_rings(mol, &rings, ez)
}

pub(crate) fn rank_with_rings(
    mol: &Molecule,
    rings: &[Vec<usize>],
    ez: &[(BondIdx, bool)],
) -> Vec<u32> {
    let inp = build_input(mol, rings, ez);
    let mut ranker = Ranker::new(&inp);
    for i in 0..inp.anum.len() {
        ranker.atoms[i].total_hs = u32::from(implicit_hcount(mol, AtomIdx(i as u32)));
    }
    ranker.rank()
}

/// One atom as RDKit's `Canon::rankMolAtoms` sees it (see [`rdkit_rank_mol_atoms`]).
#[derive(Debug, Clone, Default)]
pub struct RdkitRankAtom {
    /// `getAtomicNum()` (0 for dummies).
    pub atomic_num: u32,
    /// `getIsotope()` (0: none).
    pub isotope: u32,
    /// `getFormalCharge()`.
    pub formal_charge: i32,
    /// `molAtomMapNumber` (0: none).
    pub atom_map: i32,
    /// 0 `CHI_UNSPECIFIED`, 1 `CHI_TETRAHEDRAL_CW`, 2 `CHI_TETRAHEDRAL_CCW`
    /// (relative to the atom's bonds in bond-index order), 3 any other tag.
    pub chiral_tag: u8,
    /// `getTotalNumHs()`.
    pub total_num_hs: u32,
    /// `RingInfo::numAtomRings`.
    pub num_rings: u32,
    /// A `CW`/`CCW` atom carrying `_ringStereoAtoms`.
    pub ring_stereo: bool,
}

/// One bond as RDKit's `Canon::rankMolAtoms` sees it.
#[derive(Debug, Clone, Default)]
pub struct RdkitRankBond {
    /// Begin atom index.
    pub begin: u32,
    /// End atom index.
    pub end: u32,
    /// RDKit `Bond::BondType` value (`AROMATIC` = 12 for aromatic bonds).
    pub bond_type: u32,
    /// RDKit `Bond::BondStereo` value.
    pub stereo: u32,
    /// `getStereoAtoms()` (begin side, end side) for cis/trans stereo.
    pub stereo_atoms: Option<(u32, u32)>,
}

/// RDKit `Canon::rankMolAtoms(mol, ranks, breakTies=true, includeChirality=true,
/// includeIsotopes=true, includeAtomMaps=true, includeChiralPresence=false,
/// includeStereoGroups=true, useNonStereoRanks=false)` for a molecule without
/// stereo groups, described directly in RDKit's terms. Bonds are listed in
/// bond-index order; every atom's neighbour order follows it.
pub fn rdkit_rank_mol_atoms(atoms: &[RdkitRankAtom], bonds: &[RdkitRankBond]) -> Vec<u32> {
    let n = atoms.len();
    let mut nbrs: Vec<SmallVec<[u32; 4]>> = vec![SmallVec::new(); n];
    for b in bonds {
        nbrs[b.begin as usize].push(b.end);
        nbrs[b.end as usize].push(b.begin);
    }
    let inp = Input {
        anum: atoms.iter().map(|a| a.atomic_num).collect(),
        isotope: atoms.iter().map(|a| a.isotope).collect(),
        charge: atoms.iter().map(|a| a.formal_charge).collect(),
        map: atoms.iter().map(|a| a.atom_map).collect(),
        chiral: atoms.iter().map(|a| a.chiral_tag).collect(),
        nrings: atoms.iter().map(|a| a.num_rings).collect(),
        ring_stereo: atoms.iter().map(|a| a.ring_stereo).collect(),
        nbrs,
        bonds: bonds
            .iter()
            .map(|b| (b.begin, b.end, b.bond_type, b.stereo, b.stereo_atoms))
            .collect(),
        atom_symbols: None,
        bond_symbols: None,
        break_ties: true,
        special_chirality: true,
        in_play: None,
    };
    let mut ranker = Ranker::new(&inp);
    for (i, a) in atoms.iter().enumerate() {
        ranker.atoms[i].total_hs = a.total_num_hs;
    }
    ranker.rank()
}

/// `Canon::rankFragmentAtoms(mol, ranks, atomsInPlay, bondsInPlay, nullptr,
/// nullptr, breakTies=true, includeChirality, includeIsotopes,
/// includeAtomMaps=true, includeChiralPresence=false,
/// includeRingStereo=true)`, as `MolFragmentToSmiles` ranks. `atoms` and
/// `bonds` describe the whole molecule (chirality, isotopes and bond stereo
/// already cleared by the caller where they are not included; `num_rings`
/// from the ring perception `rankFragmentAtoms` sees); only the atoms and
/// bonds in play take part in the neighbourhoods, and an atom in play counts
/// its bonds leaving the fragment as hydrogens (RDKit's github #1567).
pub fn rdkit_rank_fragment_atoms(
    atoms: &[RdkitRankAtom],
    bonds: &[RdkitRankBond],
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    include_chirality: bool,
) -> Vec<u32> {
    let n = atoms.len();
    let mut nbrs: Vec<SmallVec<[u32; 4]>> = vec![SmallVec::new(); n];
    for b in bonds {
        nbrs[b.begin as usize].push(b.end);
        nbrs[b.end as usize].push(b.begin);
    }
    let inp = Input {
        anum: atoms.iter().map(|a| a.atomic_num).collect(),
        isotope: atoms.iter().map(|a| a.isotope).collect(),
        charge: atoms.iter().map(|a| a.formal_charge).collect(),
        map: atoms.iter().map(|a| a.atom_map).collect(),
        chiral: atoms.iter().map(|a| a.chiral_tag).collect(),
        nrings: atoms.iter().map(|a| a.num_rings).collect(),
        ring_stereo: atoms.iter().map(|a| a.ring_stereo).collect(),
        nbrs,
        bonds: bonds
            .iter()
            .map(|b| (b.begin, b.end, b.bond_type, b.stereo, b.stereo_atoms))
            .collect(),
        atom_symbols: None,
        bond_symbols: None,
        break_ties: true,
        special_chirality: include_chirality,
        in_play: Some((atoms_in_play.to_vec(), bonds_in_play.to_vec())),
    };
    let mut ranker = Ranker::new(&inp);
    for (i, a) in atoms.iter().enumerate() {
        if atoms_in_play[i] {
            let full_degree = inp.nbrs[i].len() as u32;
            ranker.atoms[i].total_hs = a.total_num_hs + full_degree - ranker.atoms[i].degree;
        }
    }
    ranker.rank()
}

/// `Canon::rankFragmentAtoms(mol, ranks, all atoms, all bonds, &atomSymbols,
/// &bondSymbols, breakTies=false, includeChirality=false,
/// includeIsotopes=false, includeAtomMaps=false, includeChiralPresence=false,
/// includeRingStereo=false)`, the symmetry classes
/// `Chirality::findPotentialStereo` ranks with: atoms compare by degree and
/// `atom_symbols`, bonds by `bond_symbols` and type. Only `begin`, `end`,
/// `bond_type` of the bonds and `num_rings`, `total_num_hs` of the atoms are
/// read.
pub fn rdkit_rank_fragment_atoms_with_symbols(
    atoms: &[RdkitRankAtom],
    bonds: &[RdkitRankBond],
    atom_symbols: &[String],
    bond_symbols: &[String],
) -> Vec<u32> {
    let n = atoms.len();
    let mut nbrs: Vec<SmallVec<[u32; 4]>> = vec![SmallVec::new(); n];
    for b in bonds {
        nbrs[b.begin as usize].push(b.end);
        nbrs[b.end as usize].push(b.begin);
    }
    let inp = Input {
        anum: atoms.iter().map(|a| a.atomic_num).collect(),
        isotope: vec![0; n],
        charge: atoms.iter().map(|a| a.formal_charge).collect(),
        map: vec![0; n],
        chiral: vec![0; n],
        nrings: atoms.iter().map(|a| a.num_rings).collect(),
        ring_stereo: vec![false; n],
        nbrs,
        bonds: bonds
            .iter()
            .map(|b| (b.begin, b.end, b.bond_type, STEREONONE, None))
            .collect(),
        atom_symbols: Some(atom_symbols.to_vec()),
        bond_symbols: Some(bond_symbols.to_vec()),
        break_ties: false,
        special_chirality: false,
        in_play: None,
    };
    let mut ranker = Ranker::new(&inp);
    for (i, a) in atoms.iter().enumerate() {
        ranker.atoms[i].total_hs = a.total_num_hs;
    }
    ranker.rank()
}

/// The Kekulé structure RDKit's `Kekulize(mol, canonical=true)` gives `mol`.
///
/// `valid` is any Kekulé structure of `mol` (from `chematic_core::kekulize`);
/// it decides which aromatic atoms take a double bond. Returns `None` when
/// RDKit's walk does not reach a structure for those candidates within its
/// back-tracking limit, or `mol` has aromatic bonds outside its rings; the
/// caller keeps `valid` then.
pub fn rdkit_canonical_kekule(mol: &Molecule, valid: &KekuleResult) -> Option<KekuleResult> {
    if valid.is_empty() {
        return Some(KekuleResult::default());
    }
    rdkit_canonical_kekule_in_rings(mol, valid, &rings_of(mol))
}

/// [`rdkit_canonical_kekule`] with `mol`'s rings already in hand (RDKit's
/// `RingInfo`: [`crate::rdkit_sssr_ring_order`] where it applies).
pub fn rdkit_canonical_kekule_with_rings(
    mol: &Molecule,
    valid: &KekuleResult,
    rings: &[Vec<AtomIdx>],
) -> Option<KekuleResult> {
    if valid.is_empty() {
        return Some(KekuleResult::default());
    }
    let rings: Vec<Vec<usize>> = rings
        .iter()
        .map(|r| r.iter().map(|a| a.0 as usize).collect())
        .collect();
    rdkit_canonical_kekule_in_rings(mol, valid, &rings)
}

fn rdkit_canonical_kekule_in_rings(
    mol: &Molecule,
    valid: &KekuleResult,
    rings: &[Vec<usize>],
) -> Option<KekuleResult> {
    let context = CanonicalKekuleContext::new(mol, valid, rings)?;
    let mut order: Vec<BondOrder> = mol.bonds().map(|(_, b)| b.order).collect();
    let mut covered = vec![false; order.len()];
    let mut fus_done = vec![false; context.rings.len()];
    let mut curr = 0usize;
    while curr < context.rings.len() {
        let fused = context.take_fused_group(curr, &mut fus_done);
        if !context.kekulize_fused_group(mol, &fused, &mut covered, &mut order) {
            return None;
        }
        match (0..context.rings.len()).find(|&ring| !fus_done[ring]) {
            Some(r) => curr = r,
            None => break,
        }
    }
    let mut out = KekuleResult::default();
    for (bi, &arom) in context.aromatic_bond.iter().enumerate() {
        if arom {
            if !covered[bi] {
                return None;
            }
            out.insert(BondIdx(bi as u32), order[bi]);
        }
    }
    Some(out)
}

struct CanonicalKekuleContext {
    rings: Vec<Vec<usize>>,
    ranks: Vec<u32>,
    aromatic_bond: Vec<bool>,
    aromatic_atom: Vec<bool>,
    candidates: Vec<bool>,
    bond_rings: Vec<Vec<usize>>,
    ring_neighbors: Vec<Vec<usize>>,
}

impl CanonicalKekuleContext {
    fn new(mol: &Molecule, valid: &KekuleResult, rings: &[Vec<usize>]) -> Option<Self> {
        let rings = rings.to_vec();
        let aromatic_bond: Vec<bool> = mol
            .bonds()
            .map(|(_, bond)| bond.order == BondOrder::Aromatic)
            .collect();
        let mut aromatic_atom = vec![false; mol.atom_count()];
        for (index, atom) in mol.atoms() {
            aromatic_atom[index.0 as usize] |= atom.aromatic;
        }
        for (_, bond) in mol.bonds().filter(|(_, bond)| bond.order == BondOrder::Aromatic) {
            aromatic_atom[bond.atom1.0 as usize] = true;
            aromatic_atom[bond.atom2.0 as usize] = true;
        }
        let mut candidates = vec![false; mol.atom_count()];
        for (&index, &order) in valid {
            if order == BondOrder::Double && aromatic_bond[index.0 as usize] {
                let bond = mol.bond(index);
                candidates[bond.atom1.0 as usize] = true;
                candidates[bond.atom2.0 as usize] = true;
            }
        }
        let bond_rings = ring_bonds(mol, &rings)?;
        let ring_neighbors = fused_ring_neighbors(&bond_rings);
        Some(Self {
            ranks: rank_with_rings(mol, &rings, &[]),
            rings,
            aromatic_bond,
            aromatic_atom,
            candidates,
            bond_rings,
            ring_neighbors,
        })
    }

    fn take_fused_group(&self, start: usize, done: &mut [bool]) -> Vec<usize> {
        let mut fused = Vec::new();
        let mut stack = vec![start];
        while let Some(ring) = stack.pop() {
            if done[ring] {
                continue;
            }
            done[ring] = true;
            fused.push(ring);
            for &neighbor in self.ring_neighbors[ring].iter().rev() {
                if !done[neighbor] {
                    stack.push(neighbor);
                }
            }
        }
        fused
    }

    fn kekulize_fused_group(
        &self,
        mol: &Molecule,
        fused: &[usize],
        covered: &mut [bool],
        order: &mut [BondOrder],
    ) -> bool {
        let (atoms, in_group) = self.fused_atoms(fused);
        if !atoms.iter().any(|&atom| self.aromatic_atom[atom]) {
            return true;
        }
        for &ring in fused {
            for &bond in &self.bond_rings[ring] {
                covered[bond] = true;
            }
        }
        let done = atoms
            .iter()
            .copied()
            .filter(|&atom| !self.aromatic_atom[atom])
            .collect();
        for (index, (_, bond)) in mol.bonds().enumerate() {
            if self.aromatic_bond[index]
                && in_group[bond.atom1.0 as usize]
                && in_group[bond.atom2.0 as usize]
            {
                order[index] = BondOrder::Single;
            }
        }
        kekulize_worker(
            mol,
            &atoms,
            &in_group,
            self.candidates.clone(),
            done,
            &self.ranks,
            &self.aromatic_bond,
            order,
        )
    }

    fn fused_atoms(&self, fused: &[usize]) -> (Vec<usize>, Vec<bool>) {
        let mut atoms = Vec::new();
        let mut present = vec![false; self.aromatic_atom.len()];
        for &ring in fused {
            for &atom in &self.rings[ring] {
                if !present[atom] {
                    present[atom] = true;
                    atoms.push(atom);
                }
            }
        }
        (atoms, present)
    }
}

fn ring_bonds(mol: &Molecule, rings: &[Vec<usize>]) -> Option<Vec<Vec<usize>>> {
    rings
        .iter()
        .map(|ring| {
            let mut bonds: Vec<usize> = (0..ring.len())
                .map(|index| {
                    mol.bond_between(
                        AtomIdx(ring[index] as u32),
                        AtomIdx(ring[(index + 1) % ring.len()] as u32),
                    )
                    .map(|(bond, _)| bond.0 as usize)
                })
                .collect::<Option<_>>()?;
            bonds.sort_unstable();
            Some(bonds)
        })
        .collect()
}

fn fused_ring_neighbors(bond_rings: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut neighbors = vec![Vec::new(); bond_rings.len()];
    for left in 0..bond_rings.len() {
        for right in left + 1..bond_rings.len() {
            if bond_rings[left]
                .iter()
                .any(|bond| bond_rings[right].binary_search(bond).is_ok())
            {
                neighbors[left].push(right);
                neighbors[right].push(left);
            }
        }
    }
    neighbors
}

/// RDKit `kekulizeWorker` (no wedged bonds; `maxBackTracks` 100).
#[allow(clippy::too_many_arguments)]
fn kekulize_worker(
    mol: &Molecule,
    all_atms: &[usize],
    in_all: &[bool],
    cands: Vec<bool>,
    done: Vec<usize>,
    ranks: &[u32],
    aromatic_bond: &[bool],
    order: &mut [BondOrder],
) -> bool {
    KekulizeWorker::new(mol, in_all, cands, done, ranks, aromatic_bond, order).run(all_atms)
}

struct KekulizeWorker<'a> {
    mol: &'a Molecule,
    in_all: &'a [bool],
    ranks: &'a [u32],
    aromatic_bond: &'a [bool],
    order: &'a mut [BondOrder],
    candidates: Vec<bool>,
    done: Vec<usize>,
    is_done: Vec<bool>,
    stack: VecDeque<usize>,
    options: std::collections::HashMap<usize, VecDeque<usize>>,
    last_option: Option<usize>,
    local_added: Vec<bool>,
    added: Vec<bool>,
    backtrack_moves: Vec<usize>,
    backtracks: usize,
}

impl<'a> KekulizeWorker<'a> {
    fn new(
        mol: &'a Molecule,
        in_all: &'a [bool],
        candidates: Vec<bool>,
        done: Vec<usize>,
        ranks: &'a [u32],
        aromatic_bond: &'a [bool],
        order: &'a mut [BondOrder],
    ) -> Self {
        let mut is_done = vec![false; in_all.len()];
        for &atom in &done {
            is_done[atom] = true;
        }
        let bond_count = order.len();
        Self {
            mol,
            in_all,
            ranks,
            aromatic_bond,
            order,
            candidates,
            done,
            is_done,
            stack: VecDeque::new(),
            options: Default::default(),
            last_option: None,
            local_added: vec![false; bond_count],
            added: vec![false; bond_count],
            backtrack_moves: Vec::new(),
            backtracks: 0,
        }
    }

    fn run(mut self, atoms: &[usize]) -> bool {
        let mut sorted = atoms.to_vec();
        sorted.sort_by_key(|&atom| (self.ranks[atom], atom));
        while self.done.len() < sorted.len() || !self.stack.is_empty() {
            let Some(current) = self.next_atom(&sorted) else {
                return false;
            };
            self.done.push(current);
            self.is_done[current] = true;
            if !self.candidates[current] {
                continue;
            }
            let mut options = self.options_for(current);
            if let Some(candidate) = options.pop_front() {
                self.select_candidate(current, candidate, options);
            } else if self.requires_match(current) && !self.backtrack() {
                self.reset_local_bonds();
                return false;
            }
        }
        true
    }

    fn next_atom(&mut self, sorted: &[usize]) -> Option<usize> {
        self.stack
            .pop_front()
            .or_else(|| sorted.iter().copied().find(|&atom| !self.is_done[atom]))
    }

    fn options_for(&mut self, current: usize) -> VecDeque<usize> {
        if let Some(options) = self.options.get(&current) {
            return options.clone();
        }
        let mut neighbors: Vec<usize> = self
            .mol
            .neighbors(AtomIdx(current as u32))
            .map(|(atom, _)| atom.0 as usize)
            .filter(|&atom| self.in_all[atom] && !self.is_done[atom])
            .collect();
        neighbors.sort_by_key(|&atom| (self.ranks[atom], atom));
        let mut options = VecDeque::new();
        for neighbor in neighbors {
            if !self.stack.contains(&neighbor) {
                self.stack.push_back(neighbor);
            }
            let bond = self.bond_index(current, neighbor);
            if self.candidates[neighbor] && self.aromatic_bond[bond] {
                options.push_back(neighbor);
            }
        }
        options
    }

    fn select_candidate(
        &mut self,
        current: usize,
        candidate: usize,
        remaining: VecDeque<usize>,
    ) {
        let revisiting = self.options.contains_key(&current);
        let bond = self.bond_index(current, candidate);
        self.order[bond] = BondOrder::Double;
        self.candidates[current] = false;
        self.candidates[candidate] = false;
        self.added[bond] = true;
        self.local_added[bond] = true;
        if remaining.is_empty() {
            self.options.remove(&current);
            if revisiting {
                self.backtrack_moves.pop();
                self.last_option = self.backtrack_moves.last().copied();
            }
        } else {
            if !revisiting {
                self.backtrack_moves.push(current);
            }
            self.last_option = Some(current);
            self.options.insert(current, remaining);
        }
    }

    fn backtrack(&mut self) -> bool {
        const MAX_BACKTRACKS: usize = 100;
        let Some(last_option) = self.last_option.filter(|_| self.backtracks < MAX_BACKTRACKS)
        else {
            return false;
        };
        let first = self
            .done
            .iter()
            .position(|&atom| atom == last_option)
            .expect("last option was visited");
        let last = self
            .done
            .iter()
            .rposition(|&atom| atom == last_option)
            .expect("last option was visited");
        for &atom in self.done[last..].iter().rev() {
            self.stack.push_front(atom);
        }
        let retained = self.done[..first].to_vec();
        self.rollback_bonds(&retained);
        self.done = retained;
        self.is_done.fill(false);
        for &atom in &self.done {
            self.is_done[atom] = true;
        }
        self.backtracks += 1;
        true
    }

    fn rollback_bonds(&mut self, retained: &[usize]) {
        let mut retained_atom = vec![false; self.in_all.len()];
        for &atom in retained {
            retained_atom[atom] = true;
        }
        for bond_index in 0..self.added.len() {
            if !self.added[bond_index] {
                continue;
            }
            let bond = self.mol.bond(BondIdx(bond_index as u32));
            let (left, right) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
            if !retained_atom[left] && !retained_atom[right] {
                self.added[bond_index] = false;
                self.order[bond_index] = BondOrder::Single;
                self.candidates[left] = true;
                self.candidates[right] = true;
            }
        }
    }

    fn reset_local_bonds(&mut self) {
        for (index, &added) in self.local_added.iter().enumerate() {
            if added {
                self.order[index] = BondOrder::Single;
            }
        }
    }

    fn requires_match(&self, atom: usize) -> bool {
        let atom = self.mol.atom(AtomIdx(atom as u32));
        atom.element.atomic_number() != 0 && !atom.wildcard
    }

    fn bond_index(&self, left: usize, right: usize) -> usize {
        self.mol
            .bond_between(AtomIdx(left as u32), AtomIdx(right as u32))
            .map(|(bond, _)| bond.0 as usize)
            .expect("ring neighbours are bonded")
    }
}
