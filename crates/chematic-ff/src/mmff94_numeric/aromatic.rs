/// Partial, behaviorally-calibrated port of RDKit's MMFF aromaticity model.
///
/// The returned molecule carries MMFF-specific aromatic flags and Kekule
/// bond orders. It is intentionally separate from chematic's general Huckel
/// model. Failure to Kekulize is returned as a typed error rather than
/// falling back to the caller's aromatic view.
pub fn compute_mmff94_aromatic_view(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
) -> Result<Molecule, NumericTypeError> {
    if rings.is_empty() {
        return Ok(mol.clone());
    }
    // RDKit's MMFF pass consumes RingInfo in size/insertion order. The
    // symmetrized SSSR canonicalizes its public ring list, so restore the
    // equivalent compatibility order here: smaller rings first and, within
    // one size, the earliest atom encountered by the input graph. This is
    // deliberately local to the RDKit-compatibility layer; general ring
    // perception remains atom-order independent.
    //
    // RDKit's pass stops once every ring atom has been looked at, so a ring
    // deferred until its neighbours are typed is decided by where it sits in
    // that order. Rings RDKit lists are taken in RDKit's own order
    // (`rdkit_sssr_ring_order`, a port of its `findSSSR`/`symmetrizeSSSR`);
    // the size/first-atom order above is the fallback.
    // The caller passes RDKit's own ring list where the port has it, so the
    // input order is RDKit's; the sort below only orders a fallback list.
    let rdkit_position: std::collections::HashMap<Vec<u32>, usize> = rings
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, ring)| {
            let mut key: Vec<u32> = ring.iter().map(|a| a.0).collect();
            key.sort_unstable();
            (key, i)
        })
        .collect();
    let mut rings = rings.to_vec();
    rings.sort_by_key(|ring| {
        let mut key: Vec<u32> = ring.iter().map(|a| a.0).collect();
        key.sort_unstable();
        (
            rdkit_position.get(&key).copied().unwrap_or(usize::MAX),
            ring.len(),
            ring.iter().map(|atom| atom.0).min().unwrap_or(u32::MAX),
        )
    });
    // RDKit kekulizes with canonical atom ranks before this pass; where a
    // system has several Kekule structures (fullerene cages) the pass reads
    // a different set of hexagons as aromatic depending on which one it
    // gets, so take RDKit's (`rdkit_canonical_kekule`), falling back to
    // chematic's own valid structure.
    //
    // `MolFromSmiles` perceives aromaticity from scratch whatever the input
    // spelling, so a Kekule-written ring system gets the same structure as
    // its aromatic spelling: kekulize RDKit's aromatic view, not the input's
    // own double bonds.
    let rdkit_view = chematic_perception::apply_aromaticity_rdkit_parity_experimental(mol).ok();
    let base: &Molecule = rdkit_view.as_ref().unwrap_or(mol);
    let kmol = match chematic_core::kekulize(base) {
        Ok(kek) if kek.is_empty() => base.clone(),
        Ok(kek) => {
            // Where every aromatic system is one ring that no other ring
            // shares a bond with, the Kekule structures only move that ring's
            // double bonds around: once the pass accepts the ring its bonds
            // are aromatic whichever one was used, so chematic's own structure
            // gives RDKit's result without the canonical ranking. A ring the
            // pass rejects keeps its double bonds, so it falls through.
            if let Some(rings_of_systems) = isolated_aromatic_rings(base, &rings) {
                let kmol = chematic_core::apply_kekule(base, &kek);
                let (view, accepted) = mmff_aromatic_pass(&kmol, &rings);
                if rings_of_systems.iter().all(|&r| accepted[r]) {
                    return Ok(view);
                }
            }
            let kek = chematic_perception::rdkit_canonical_kekule_with_rings(base, &kek, &rings)
                .unwrap_or(kek);
            chematic_core::apply_kekule(base, &kek)
        }
        Err(e) => {
            return Err(NumericTypeError(format!(
                "MMFF94 aromaticity re-perception failed at the Kekulization stage: {} \
                 -- refusing to type-assign this molecule rather than silently reusing its \
                 pre-Kekulization (possibly-aromatic-perceived) form as the MMFF view",
                e.detail
            )));
        }
    };

    Ok(mmff_aromatic_pass(&kmol, &rings).0)
}

/// For each aromatic system of `base` (connected aromatic bonds), the index
/// in `rings` of the one ring it is, when every system is a single ring
/// whose bonds lie on no other ring; `None` otherwise.
fn isolated_aromatic_rings(base: &Molecule, rings: &[Vec<AtomIdx>]) -> Option<Vec<usize>> {
    let key = |a: AtomIdx, b: AtomIdx| (a.0.min(b.0), a.0.max(b.0));
    let aromatic: std::collections::HashSet<(u32, u32)> = base
        .bonds()
        .filter(|(_, b)| b.order == BondOrder::Aromatic)
        .map(|(_, b)| key(b.atom1, b.atom2))
        .collect();
    if aromatic.is_empty() {
        return None;
    }
    let mut covered: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
    let mut out = Vec::new();
    for (ri, ring) in rings.iter().enumerate() {
        let len = ring.len();
        let bonds: Vec<(u32, u32)> = (0..len)
            .map(|j| key(ring[j], ring[(j + 1) % len]))
            .collect();
        let n_arom = bonds.iter().filter(|b| aromatic.contains(b)).count();
        if n_arom == 0 {
            continue;
        }
        // A ring partly aromatic, or sharing an aromatic bond with a ring
        // already taken: not an isolated system.
        if n_arom != len || bonds.iter().any(|b| !covered.insert(*b)) {
            return None;
        }
        out.push(ri);
    }
    // Every aromatic bond on some ring, and no aromatic atom joining two
    // rings (a spiro or ring-assembly system through aromatic bonds).
    if covered.len() != aromatic.len() {
        return None;
    }
    let mut seen = std::collections::HashSet::new();
    for &ri in &out {
        for &a in &rings[ri] {
            if !seen.insert(a) {
                return None;
            }
        }
    }
    Some(out)
}

/// RDKit's `setMMFFAromaticity` on a kekulized molecule: the MMFF view and,
/// per ring, whether the ring itself passed.
fn mmff_aromatic_pass(kmol: &Molecule, rings: &[Vec<AtomIdx>]) -> (Molecule, Vec<bool>) {
    let mut pass = MmffAromaticityPass::new(kmol.atom_count(), rings.len());
    pass.resolve(kmol, rings);
    let view = build_mmff_aromatic_view(kmol, rings, &pass.aromatic, &pass.ring_accepted);
    (view, pass.ring_accepted)
}

struct MmffAromaticityPass {
    resolved: Vec<bool>,
    aromatic: Vec<bool>,
    ring_accepted: Vec<bool>,
}

impl MmffAromaticityPass {
    fn new(atom_count: usize, ring_count: usize) -> Self {
        Self {
            resolved: vec![false; atom_count],
            aromatic: vec![false; atom_count],
            ring_accepted: vec![false; ring_count],
        }
    }

    fn resolve(&mut self, mol: &Molecule, rings: &[Vec<AtomIdx>]) {
        let mut previous = None;
        for _ in 0..rings.len() + 2 {
            let resolved_before = self.resolved_count(rings);
            if previous.is_some_and(|count| resolved_before <= count) {
                break;
            }
            previous = Some(resolved_before);
            for (ring_index, ring) in rings.iter().enumerate() {
                let Some(accepted) = self.evaluate_ring(mol, rings, ring) else {
                    continue;
                };
                for atom in ring {
                    self.resolved[atom.0 as usize] = true;
                    if accepted {
                        self.aromatic[atom.0 as usize] = true;
                    }
                }
                self.ring_accepted[ring_index] |= accepted;
            }
            if rings
                .iter()
                .all(|ring| ring.iter().all(|atom| self.resolved[atom.0 as usize]))
            {
                break;
            }
        }
    }

    fn evaluate_ring(
        &self,
        mol: &Molecule,
        rings: &[Vec<AtomIdx>],
        ring: &[AtomIdx],
    ) -> Option<bool> {
        let (mut electrons, contains_nos, exocyclic_double) =
            self.ring_pi_electrons(mol, rings, ring)?;
        if !ring.iter().all(|&atom| {
            !matches!(mol.atom(atom).element, Element::C | Element::N)
                || total_degree(mol, atom) <= 3
        }) {
            return Some(false);
        }
        if contains_nos && !exocyclic_double && ring.len() % 2 == 1 {
            electrons += 2;
        }
        Some(electrons > 2 && (electrons - 2) % 4 == 0)
    }

    fn ring_pi_electrons(
        &self,
        mol: &Molecule,
        rings: &[Vec<AtomIdx>],
        ring: &[AtomIdx],
    ) -> Option<(i32, bool, bool)> {
        let mut electrons = 0;
        let mut contains_nos = false;
        let mut exocyclic_double = false;
        for (position, &atom_index) in ring.iter().enumerate() {
            let atom = mol.atom(atom_index);
            contains_nos |= matches!(atom.element, Element::N | Element::O)
                || (atom.element == Element::S && total_degree(mol, atom_index) == 2);
            let next = ring[(position + 1) % ring.len()];
            let ring_order = mol
                .bond_between(atom_index, next)
                .map(|(_, bond)| bond.order)
                .unwrap_or(BondOrder::Single);
            if ring_order == BondOrder::Double {
                electrons += 2;
                continue;
            }
            let candidate = atom.element == Element::C
                || (atom.element == Element::N && total_valence(mol, atom_index) == 4);
            if candidate {
                self.add_exocyclic_contribution(
                    mol,
                    rings,
                    ring,
                    atom_index,
                    &mut electrons,
                    &mut exocyclic_double,
                )?;
            }
        }
        Some((electrons, contains_nos, exocyclic_double))
    }

    fn add_exocyclic_contribution(
        &self,
        mol: &Molecule,
        rings: &[Vec<AtomIdx>],
        ring: &[AtomIdx],
        atom: AtomIdx,
        electrons: &mut i32,
        exocyclic_double: &mut bool,
    ) -> Option<()> {
        for bond in bonds_of(mol, atom) {
            if ring.contains(&bond.neighbor)
                || matches!(bond.order, BondOrder::Single | BondOrder::Up | BondOrder::Down)
            {
                continue;
            }
            let neighbor_in_ring = rings.iter().any(|candidate| candidate.contains(&bond.neighbor));
            if neighbor_in_ring && !self.resolved[bond.neighbor.0 as usize] {
                return None;
            }
            if bond.order == BondOrder::Double {
                if self.aromatic[bond.neighbor.0 as usize] {
                    *electrons += 1;
                } else {
                    *exocyclic_double = true;
                }
            }
        }
        Some(())
    }

    fn resolved_count(&self, rings: &[Vec<AtomIdx>]) -> i64 {
        rings
            .iter()
            .map(|ring| {
                ring.iter()
                    .filter(|atom| self.resolved[atom.0 as usize])
                    .count() as i64
            })
            .sum()
    }
}

fn build_mmff_aromatic_view(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    aromatic_atoms: &[bool],
    accepted_rings: &[bool],
) -> Molecule {
    if !mol.adjacency_in_bond_order() {
        return rebuild_mmff_aromatic_view(mol, rings, aromatic_atoms, accepted_rings);
    }
    let mut view = mol.clone();
    for (index, _) in mol.atoms() {
        if mol.atom(index).aromatic != aromatic_atoms[index.0 as usize] {
            view.set_atom_aromatic(index, aromatic_atoms[index.0 as usize]);
        }
    }
    for (ring_index, ring) in rings.iter().enumerate() {
        if !accepted_rings[ring_index] {
            continue;
        }
        for (position, &atom) in ring.iter().enumerate() {
            let next = ring[(position + 1) % ring.len()];
            let bond = mol
                .bond_between(atom, next)
                .expect("consecutive ring atoms are bonded")
                .0;
            if view.bond(bond).order != BondOrder::Aromatic {
                view.set_bond_order(bond, BondOrder::Aromatic);
            }
        }
    }
    view
}

fn rebuild_mmff_aromatic_view(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    aromatic_atoms: &[bool],
    accepted_rings: &[bool],
) -> Molecule {
    let mut aromatic_bonds = vec![false; mol.bond_count()];
    for (ring_index, ring) in rings.iter().enumerate() {
        if !accepted_rings[ring_index] {
            continue;
        }
        for (position, &atom) in ring.iter().enumerate() {
            if let Some((bond, _)) = mol.bond_between(atom, ring[(position + 1) % ring.len()]) {
                aromatic_bonds[bond.0 as usize] = true;
            }
        }
    }
    let mut builder =
        chematic_core::MoleculeBuilder::with_capacity(mol.atom_count(), mol.bond_count());
    for (index, atom) in mol.atoms() {
        let mut atom = atom.clone();
        atom.aromatic = aromatic_atoms[index.0 as usize];
        builder.add_atom(atom);
    }
    for (index, bond) in mol.bonds() {
        let order = if aromatic_bonds[index.0 as usize] {
            BondOrder::Aromatic
        } else {
            bond.order
        };
        builder
            .add_bond(bond.atom1, bond.atom2, order)
            .expect("duplicate bond during MMFF94 aromaticity re-perception");
    }
    builder.build()
}

/// RDKit's `isAtomNOxide`: a >=3-connected nitrogen with a terminal
/// (degree-1) oxygen neighbor.
fn is_atom_n_oxide(mol: &Molecule, idx: AtomIdx) -> bool {
    if mol.atom(idx).element != Element::N || total_degree(mol, idx) < 3 {
        return false;
    }
    bonds_of(mol, idx)
        .iter()
        .any(|b| mol.atom(b.neighbor).element == Element::O && total_degree(mol, b.neighbor) == 1)
}

/// The alpha/beta heteroatom classification RDKit's `setMMFFHeavyAtomType`
/// computes once per aromatic 5-ring C/N atom and reuses across several
/// branches. `alpha`/`beta` neighbors are ring-adjacent O/S/(non-N-oxide,
/// 3-connected) N in the *same* 5-membered aromatic ring as `atom`.
struct AlphaBetaHeteroatoms {
    alpha: Vec<AtomIdx>,
    beta: Vec<AtomIdx>,
    is_alpha_os: bool,
    is_beta_os: bool,
    alpha_or_beta_in_same_ring: bool,
}

fn is_alpha_beta_heteroatom_candidate(mol: &Molecule, idx: AtomIdx) -> bool {
    let atom = mol.atom(idx);
    matches!(atom.element, Element::O | Element::S)
        || (atom.element == Element::N && total_degree(mol, idx) == 3 && !is_atom_n_oxide(mol, idx))
}

fn find_alpha_beta_heteroatoms(
    mol: &Molecule,
    rings: &[Vec<AtomIdx>],
    atom: AtomIdx,
) -> AlphaBetaHeteroatoms {
    let mut alpha = Vec::new();
    let mut beta = Vec::new();

    if matches!(mol.atom(atom).element, Element::C | Element::N) {
        for nb in bonds_of(mol, atom) {
            let nbr = nb.neighbor;
            if !atom_in_aromatic_ring_of_size(mol, rings, nbr, 5) {
                continue;
            }
            if atoms_share_aromatic_ring_of_size(mol, rings, atom, nbr, 5)
                && is_alpha_beta_heteroatom_candidate(mol, nbr)
            {
                alpha.push(nbr);
            }
            for nb2 in bonds_of(mol, nbr) {
                let nbr2 = nb2.neighbor;
                if nbr2 == atom {
                    continue;
                }
                if !atom_in_aromatic_ring_of_size(mol, rings, nbr2, 5) {
                    continue;
                }
                if atoms_share_aromatic_ring_of_size(mol, rings, atom, nbr2, 5)
                    && is_alpha_beta_heteroatom_candidate(mol, nbr2)
                {
                    beta.push(nbr2);
                }
            }
        }
    }

    let is_alpha_os = alpha
        .iter()
        .any(|&a| matches!(mol.atom(a).element, Element::O | Element::S));
    let is_beta_os = beta
        .iter()
        .any(|&b| matches!(mol.atom(b).element, Element::O | Element::S));
    let alpha_or_beta_in_same_ring = !alpha.is_empty()
        && !beta.is_empty()
        && alpha.iter().any(|&a| {
            beta.iter()
                .any(|&b| atoms_share_aromatic_ring_of_size(mol, rings, a, b, 5))
        });

    AlphaBetaHeteroatoms {
        alpha,
        beta,
        is_alpha_os,
        is_beta_os,
        alpha_or_beta_in_same_ring,
    }
}
