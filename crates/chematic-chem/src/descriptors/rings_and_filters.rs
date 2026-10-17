// ---------------------------------------------------------------------------
// 13. Basic count descriptors
// ---------------------------------------------------------------------------

/// Number of heteroatoms (non-C, non-H heavy atoms).
pub fn num_heteroatoms(mol: &Molecule) -> usize {
    // RDKit's `CalcNumHeteroatoms`: atomic number other than 1 and 6, so
    // dummy atoms (atomic number 0) count.
    mol.atoms()
        .filter(|(_, a)| {
            let an = if a.wildcard {
                0
            } else {
                a.element.atomic_number()
            };
            an != 1 && an != 6
        })
        .count()
}

/// RDKit's `CalcNumRings`: the number of rings in RDKit's ring information,
/// its symmetrized SSSR. It exceeds [`ring_count`] (the SSSR size, i.e. the
/// cycle rank) for cages and bridged systems such as cubane or adamantane.
#[inline(always)]
pub fn rdkit_num_rings(mol: &Molecule) -> usize {
    mol.derived_rdkit_num_rings(|| rdkit_num_rings_uncached(mol))
}

fn rdkit_num_rings_uncached(mol: &Molecule) -> usize {
    // The SMILES parser seeds an all-false RingBondFlags value when it has
    // already proved that the graph is a forest.  Reuse that proof instead
    // of allocating union-find state and rediscovering the same fact.  Every
    // graph mutation clears the derived cache, so edited molecules cannot
    // observe a stale zero.
    if mol
        .derived_if_computed::<Vec<bool>>(chematic_core::DerivedSlot::RingBondFlags)
        .is_some_and(|flags| !flags.iter().any(|&is_ring| is_ring))
    {
        return 0;
    }
    if let Some(count) = single_smiles_closure_ring_count(mol) {
        return count;
    }
    match simple_ring_system_count(mol) {
        Some(RingCountDecision::Count(count)) => count,
        decision => {
            // RDKit 2026.09 switched its default symmetrized-ring backend to
            // RingDecomposerLib. Its ring families are the relevant cycles;
            // this matters for large macrocycles with several equivalent
            // paths, where the legacy SSSR under-counts. Keep the legacy port
            // as a bounded fallback and for organometallic graphs whose
            // sanitize-time dative rewrite is not represented in the input.
            const MAX_RELEVANT_RINGS: usize = 1_000_000;
            // The RDKit 2026.09 rebaseline changes exposed by the pinned 10k
            // corpus all have legacy macrocycles of 28-32 atoms. Smaller
            // rings are handled by the equivalent, much faster legacy port.
            const MACROCYCLE_MIN_SIZE: usize = 24;
            let (legacy_count, may_need_ring_families) = match decision {
                Some(RingCountDecision::Relevant { legacy_count }) => (legacy_count, true),
                Some(RingCountDecision::Count(_)) => unreachable!(),
                None => {
                    let legacy = rdkit_ring_list(mol);
                    let may_need = legacy.iter().any(|ring| ring.len() >= MACROCYCLE_MIN_SIZE);
                    (legacy.len(), may_need)
                }
            };
            if may_need_ring_families
                && !chematic_perception::rdkit_organometallic_dative_bonds(mol)
                    .iter()
                    .any(|&dative| dative)
                && let Ok(rings) = chematic_perception::relevant_cycles(mol, MAX_RELEVANT_RINGS)
            {
                // RingDecomposerLib keeps the independent relevant families,
                // not every combination of equal shortest paths.  The graph
                // backend intentionally expands those combinations for its
                // public `relevant_cycles` API, which can over-count highly
                // repetitive cyclic oligosaccharides (for example 70 versus
                // RDKit's legacy-compatible 22 or 16).  Genuine symmetric
                // macrocycle changes in the pinned RDKit corpus stay below
                // this bounded expansion (10 -> 22).
                let family_count = rings.len();
                if family_count <= legacy_count.saturating_mul(3) {
                    return family_count;
                }
            }
            legacy_count
        }
    }
}

/// Fast path for an unmodified SMILES graph with one ring-closure bond.
///
/// A ring digit only establishes a graph cycle when the closure's endpoints
/// are already connected through RDKit's ring graph. It may instead join two
/// dot-separated components, or the only alternate path may contain a metal
/// coordination bond that RDKit excludes from ring finding.
fn single_smiles_closure_ring_count(mol: &Molecule) -> Option<usize> {
    if mol.smiles_ring_closure_bond_count() != 1 {
        return None;
    }
    let (closure_idx, closure) = mol
        .bonds()
        .find(|(bond_idx, _)| mol.is_smiles_ring_closure(*bond_idx))?;
    if mol.bonds().any(|(_, bond)| {
        matches!(
            bond.order,
            BondOrder::QueryAny
                | BondOrder::QuerySingleOrDouble
                | BondOrder::QuerySingleOrAromatic
                | BondOrder::QueryDoubleOrAromatic
        )
    }) {
        return None;
    }

    let excluded = chematic_perception::rdkit_organometallic_dative_bonds(mol);
    if excluded[closure_idx.0 as usize] || closure.order == BondOrder::Dative {
        return Some(0);
    }
    // Most parsed molecules fit in two visited words and a short pending
    // stack. SmallVec retains the same traversal with heap fallback for
    // large rings or many branches.
    let mut seen = smallvec::SmallVec::<[u64; 2]>::from_elem(
        0,
        mol.atom_count().div_ceil(u64::BITS as usize),
    );
    let mut stack = smallvec::SmallVec::<[AtomIdx; 16]>::new();
    stack.push(closure.atom1);
    let start = closure.atom1.0 as usize;
    seen[start / 64] |= 1u64 << (start % 64);
    while let Some(atom) = stack.pop() {
        for (neighbor, bond_idx) in mol.neighbors(atom) {
            if bond_idx == closure_idx
                || excluded[bond_idx.0 as usize]
                || mol.bond(bond_idx).order == BondOrder::Dative
            {
                continue;
            }
            if neighbor == closure.atom2 {
                return Some(1);
            }
            let index = neighbor.0 as usize;
            let visited = &mut seen[index / 64];
            let mask = 1u64 << (index % 64);
            if *visited & mask == 0 {
                *visited |= mask;
                stack.push(neighbor);
            }
        }
    }
    Some(0)
}

/// The number of ring systems when every ring system is a single simple
/// cycle (as many ring bonds as ring atoms): then RDKit's symmetrized SSSR
/// is exactly those cycles. `None` for fused, spiro or bridged systems, and
/// where RDKit's ring graph differs from chematic's (organometallic dative
/// bonds, query bond orders).
enum RingCountDecision {
    Count(usize),
    Relevant { legacy_count: usize },
}

fn simple_ring_system_count(mol: &Molecule) -> Option<RingCountDecision> {
    if mol.bonds().any(|(_, b)| {
        matches!(
            b.order,
            BondOrder::QueryAny
                | BondOrder::QuerySingleOrDouble
                | BondOrder::QuerySingleOrAromatic
                | BondOrder::QueryDoubleOrAromatic
        )
    }) || chematic_perception::rdkit_organometallic_dative_bonds(mol)
        .iter()
        .any(|&d| d)
    {
        return None;
    }
    let flags = chematic_perception::ring_bond_flags_shared(mol);
    // The bridge pass that produced `flags` also labels the cyclic
    // components. Reuse those labels instead of rebuilding a union-find.
    let components = chematic_perception::sssr::ring_components_shared(mol);
    let labels = &components.label;
    let mut balance = vec![0i32; components.count as usize];
    let mut on_ring = vec![false; mol.atom_count()];
    let mut roots = Vec::new();
    for (bi, bond) in mol.bonds() {
        if !flags[bi.0 as usize] {
            continue;
        }
        let r = labels[bond.atom1.0 as usize] as usize;
        // Before counting atoms, a zero balance means this is the first
        // cyclic bond in the component. Non-ring singleton labels are never
        // added, so they cannot be mistaken for a simple cycle below.
        if balance[r] == 0 {
            roots.push(r);
        }
        balance[r] -= 1;
        on_ring[bond.atom1.0 as usize] = true;
        on_ring[bond.atom2.0 as usize] = true;
    }
    for (a, &is_on_ring) in on_ring.iter().enumerate() {
        if is_on_ring {
            balance[labels[a] as usize] += 1;
        }
    }
    let mut rings = 0usize;
    let mut needs_relevant_cycles = false;
    for &r in &roots {
        match balance[r] {
            0 => rings += 1,
            -1 => rings += two_cycle_system_rings(mol, &flags, labels, r)?,
            _ => {
                let (count, has_macrocycle) = system_rings_rdkit(mol, &flags, labels, r)?;
                // RingDecomposerLib differs from the legacy symmetrized SSSR
                // in the exposed corpus only for large, symmetric ring
                // families. Keep the fast legacy result for compact cages
                // and fused systems; route macrocycles to the relevant-cycle
                // backend in `rdkit_num_rings`.
                if has_macrocycle {
                    needs_relevant_cycles = true;
                }
                rings += count;
            }
        }
    }
    Some(if needs_relevant_cycles {
        RingCountDecision::Relevant {
            legacy_count: rings,
        }
    } else {
        RingCountDecision::Count(rings)
    })
}

/// Legacy symmetrized-SSSR count for one ring system, plus whether that
/// result contains a macrocycle that may require RDKit 2026.09's relevant
/// cycle families.
fn system_rings_rdkit(
    mol: &Molecule,
    flags: &[bool],
    labels: &[u32],
    root: usize,
) -> Option<(usize, bool)> {
    let mut local = vec![u32::MAX; mol.atom_count()];
    let mut k = 0usize;
    let mut bonds = Vec::new();
    for bi in mol.rdkit_bond_order() {
        if !flags[bi.0 as usize] {
            continue;
        }
        let bond = mol.bond(bi);
        if labels[bond.atom1.0 as usize] as usize != root {
            continue;
        }
        let mut id = |a: AtomIdx| {
            let slot = &mut local[a.0 as usize];
            if *slot == u32::MAX {
                *slot = k as u32;
                k += 1;
            }
            *slot as usize
        };
        let (a, b) = (id(bond.atom1), id(bond.atom2));
        bonds.push((a, b, true));
    }
    let rings = chematic_perception::rdkit_symmetrized_sssr(k, &bonds)?;
    const MACROCYCLE_MIN_SIZE: usize = 24;
    let has_macrocycle = rings.iter().any(|ring| ring.len() >= MACROCYCLE_MIN_SIZE);
    Some((rings.len(), has_macrocycle))
}

/// RDKit's ring count for a ring system of cycle rank 2 rooted at `root`:
/// two cycles sharing one atom give 2; two branch atoms joined by three
/// paths of `a <= b <= c` bonds give the rings `a+b` and `a+c`, plus `b+c`
/// when `a == b` (it then has the size of `a+c` and can stand in for it).
fn two_cycle_system_rings(
    mol: &Molecule,
    flags: &[bool],
    labels: &[u32],
    root: usize,
) -> Option<usize> {
    let ring_neighbors = |a: AtomIdx| -> smallvec::SmallVec<[AtomIdx; 4]> {
        mol.neighbors(a)
            .filter(|(_, bi)| flags[bi.0 as usize])
            .map(|(nb, _)| nb)
            .collect()
    };
    let mut branch = smallvec::SmallVec::<[AtomIdx; 2]>::new();
    for (idx, _) in mol.atoms() {
        if labels[idx.0 as usize] as usize != root {
            continue;
        }
        match ring_neighbors(idx).len() {
            2 => {}
            3 => branch.push(idx),
            4 => return Some(2),
            _ => return None,
        }
    }
    let &[u, v] = branch.as_slice() else {
        return None;
    };
    let mut lengths = smallvec::SmallVec::<[usize; 3]>::new();
    for start in ring_neighbors(u) {
        let (mut prev, mut cur, mut len) = (u, start, 1usize);
        while cur != v {
            let next = ring_neighbors(cur).into_iter().find(|&x| x != prev)?;
            prev = cur;
            cur = next;
            len += 1;
            if len > mol.atom_count() {
                return None;
            }
        }
        lengths.push(len);
    }
    if lengths.len() != 3 {
        return None;
    }
    lengths.sort_unstable();
    Some(if lengths[0] == lengths[1] { 3 } else { 2 })
}

/// RDKit's ring information (`SanitizeMol`'s symmetrized SSSR), falling back
/// to chematic's symmetrized SSSR where RDKit itself would use its
/// approximate ring finder.
pub(crate) fn rdkit_ring_list(mol: &Molecule) -> Vec<Vec<AtomIdx>> {
    chematic_perception::rdkit_sssr_ring_order(mol).unwrap_or_else(|| {
        chematic_perception::find_symmetrized_sssr(mol)
            .rings()
            .to_vec()
    })
}

/// Total number of rings (SSSR count).
pub fn ring_count(mol: &Molecule) -> usize {
    // The SSSR size is the cycle rank; no need to build the rings.
    chematic_perception::sssr_ring_count(mol)
}

/// Number of distinct connected ring systems.
///
/// Two SSSR rings belong to the same ring system when they share at least
/// one atom (fused, bridged, or spiro connections).  This differs from
/// [`ring_count`] which counts every SSSR ring individually.
///
/// Examples:
/// - benzene: 1 ring system (simple 6-membered ring)
/// - naphthalene: 1 ring system (two fused rings → one system)
/// - biphenyl: 2 ring systems (two independent benzene rings connected by a bond)
///
/// Useful for scaffold complexity assessment (ScaffoldGraph, Chemical-Descriptors).
pub fn ring_system_count(mol: &Molecule) -> usize {
    use chematic_perception::find_ring_families;
    let sssr = find_sssr(mol);
    find_ring_families(mol, &sssr).len()
}

/// Lipinski (1997) HBA count — total number of N and O heavy atoms.
///
/// The original Rule-of-Five definition simply counts all nitrogen and oxygen
/// atoms regardless of hybridisation or substitution.  For the more
/// chemically accurate Ertl (2000) definition see [`hba_count`].
///
/// Used by Chemical-Descriptors, DOPtools, MolVS, and classic implementations.
pub fn hba_count_lipinski(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(_, a)| {
            let an = a.element.atomic_number();
            an == 7 || an == 8
        })
        .count()
}

/// Fraction of heavy atoms that participate in a rotatable bond.
///
/// `fraction_rotatable_bonds = rotatable_bond_count / heavy_atom_count`.
/// Returns `0.0` for molecules with zero heavy atoms.
/// Normalises molecular flexibility for cross-scaffold comparison (DOPtools).
pub fn fraction_rotatable_bonds(mol: &Molecule) -> f64 {
    let hac = heavy_atom_count(mol);
    if hac == 0 {
        return 0.0;
    }
    rotatable_bond_count(mol) as f64 / hac as f64
}

/// RDKit's ring classification (Lipinski.cpp) over its own ring
/// information: each ring of [`rdkit_ring_list`] as `(bonds aromatic?,
/// bonds single non-aromatic?, any non-carbon atom on a ring bond?)`, on the
/// RDKit-perceived aromaticity of `mol`.
fn rdkit_ring_classes(mol: &Molecule) -> Vec<(bool, bool, bool)> {
    let view = descriptor_aromaticity(mol);
    let view: &Molecule = &view;
    rdkit_ring_list(view)
        .iter()
        .map(|ring| {
            let n = ring.len();
            let (mut all_aromatic, mut all_single, mut hetero) = (true, true, false);
            for i in 0..n {
                let (a, b) = (ring[i], ring[(i + 1) % n]);
                let Some((bidx, _)) = view.bond_between(a, b) else {
                    continue;
                };
                let order = view.bond(bidx).order;
                let aromatic = order == BondOrder::Aromatic;
                all_aromatic &= aromatic;
                all_single &= !aromatic
                    && matches!(order, BondOrder::Single | BondOrder::Up | BondOrder::Down);
                // Dummy atoms (atomic number 0 in RDKit) count as non-carbon.
                let non_carbon = |x: AtomIdx| {
                    let atom = view.atom(x);
                    atom.wildcard || atom.element.atomic_number() != 6
                };
                hetero |= non_carbon(a) || non_carbon(b);
            }
            (all_aromatic, all_single, hetero)
        })
        .collect()
}

/// Number of aliphatic rings, RDKit's `CalcNumAliphaticRings`: rings of
/// RDKit's ring information with at least one non-aromatic bond.
pub fn num_aliphatic_rings(mol: &Molecule) -> usize {
    rdkit_ring_classes(mol).iter().filter(|c| !c.0).count()
}

/// Number of saturated rings, RDKit's `CalcNumSaturatedRings`: every ring
/// bond single and non-aromatic (exocyclic double bonds do not count).
pub fn num_saturated_rings(mol: &Molecule) -> usize {
    rdkit_ring_classes(mol).iter().filter(|c| c.1).count()
}

/// Returns true if all bonds *within* the ring are single bonds (non-aromatic, non-double).
fn is_ring_saturated(mol: &Molecule, ring: &[AtomIdx]) -> bool {
    let n = ring.len();
    (0..n).all(|i| {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        mol.bond_between(a, b)
            .map(|(bidx, _)| {
                !matches!(
                    mol.bond(bidx).order,
                    BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
                )
            })
            .unwrap_or(true)
    })
}

/// Number of aromatic heterocycles, RDKit's `CalcNumAromaticHeterocycles`:
/// rings whose bonds are all aromatic, with a non-carbon atom.
pub fn num_aromatic_heterocycles(mol: &Molecule) -> usize {
    rdkit_ring_classes(mol)
        .iter()
        .filter(|c| c.0 && c.2)
        .count()
}

/// Number of aliphatic heterocycles, RDKit's `CalcNumAliphaticHeterocycles`:
/// rings with a non-aromatic bond and a non-carbon atom.
pub fn num_aliphatic_heterocycles(mol: &Molecule) -> usize {
    rdkit_ring_classes(mol)
        .iter()
        .filter(|c| !c.0 && c.2)
        .count()
}

/// Number of saturated heterocycles, RDKit's `CalcNumSaturatedHeterocycles`.
pub fn num_saturated_heterocycles(mol: &Molecule) -> usize {
    rdkit_ring_classes(mol)
        .iter()
        .filter(|c| c.1 && c.2)
        .count()
}

/// Number of spiro atoms.
///
/// A spiro atom is in ≥ 2 rings and is the sole shared atom between at least one pair
/// of those rings.  The "exactly 2 rings" shortcut is wrong for molecules where the
/// augmented ring set contains an XOR ring that also passes through the spiro centre
/// (e.g. tropane bridge, peroxide cages).
/// Example: spiro[4.5]decane (`C1CCCCC11CCCC1`) has 1 spiro atom.
fn num_spiro_atoms_from(mol: &Molecule, rings: &[Vec<AtomIdx>]) -> usize {
    mol.atoms()
        .filter(|(idx, _)| {
            let member: Vec<_> = rings.iter().filter(|r| r.contains(idx)).collect();
            if member.len() < 2 {
                return false;
            }
            // Spiro: exists any pair of rings that shares ONLY this one atom.
            (0..member.len()).any(|i| {
                (i + 1..member.len())
                    .any(|j| member[i].iter().filter(|a| member[j].contains(a)).count() == 1)
            })
        })
        .count()
}

pub fn num_spiro_atoms(mol: &Molecule) -> usize {
    num_spiro_atoms_from(mol, &all_ring_list(mol))
}

/// Number of bridgehead atoms.
///
/// An atom is a bridgehead if it lies at the junction of a bridged ring system:
/// for some pair of rings that share ≥ 2 bonds, the atom is incident to exactly
/// one of those shared bonds.
///
/// Example: norbornane (`C1CC2CCC1C2`) has 2 bridgehead atoms.
/// Naphthalene has 0 (junction atoms share exactly 1 bond — fused, not bridged).
/// Spiro[4.5]decane has 0 (no shared bonds between the two rings).
fn num_bridgehead_atoms_from(mol: &Molecule, rings: &[Vec<AtomIdx>]) -> usize {
    // Pre-compute sorted bond-index list for each ring (sorted so intersection is a merge).
    let bond_lists: Vec<Vec<BondIdx>> = rings
        .iter()
        .map(|r| {
            let mut bonds: Vec<BondIdx> = (0..r.len())
                .filter_map(|i| {
                    let a = r[i];
                    let b = r[(i + 1) % r.len()];
                    mol.bond_between(a, b).map(|(bidx, _)| bidx)
                })
                .collect();
            bonds.sort_unstable_by_key(|b| b.0);
            bonds
        })
        .collect();

    let mut bridgehead_set: FxHashSet<AtomIdx> = FxHashSet::default();

    for i in 0..rings.len() {
        for j in (i + 1)..rings.len() {
            // Intersect sorted bond lists.
            let mut shared_bonds: Vec<BondIdx> = Vec::new();
            let (mut ai, mut bi) = (0usize, 0usize);
            while ai < bond_lists[i].len() && bi < bond_lists[j].len() {
                match bond_lists[i][ai].0.cmp(&bond_lists[j][bi].0) {
                    std::cmp::Ordering::Equal => {
                        shared_bonds.push(bond_lists[i][ai]);
                        ai += 1;
                        bi += 1;
                    }
                    std::cmp::Ordering::Less => ai += 1,
                    std::cmp::Ordering::Greater => bi += 1,
                }
            }
            if shared_bonds.len() < 2 {
                continue; // fused (1 shared bond) or disjoint (0)
            }

            let ring_j_set: FxHashSet<AtomIdx> = rings[j].iter().copied().collect();
            for &atom in rings[i].iter().filter(|a| ring_j_set.contains(a)) {
                let incident = shared_bonds
                    .iter()
                    .filter(|&&b| {
                        let bond = mol.bond(b);
                        bond.atom1 == atom || bond.atom2 == atom
                    })
                    .count();
                if incident == 1 {
                    bridgehead_set.insert(atom);
                }
            }
        }
    }

    bridgehead_set.len()
}

pub fn num_bridgehead_atoms(mol: &Molecule) -> usize {
    num_bridgehead_atoms_from(mol, &all_ring_list(mol))
}

// ---------------------------------------------------------------------------
// Ring descriptor bundle — find_sssr called exactly once
// ---------------------------------------------------------------------------

/// All ring-derived descriptor values for a molecule, computed with a single `find_sssr` call.
///
/// Use [`ring_bundle`] to obtain this struct. Individual descriptor functions remain
/// available for single-value use but each re-invoke `find_sssr` independently.
#[derive(Debug, Clone)]
pub struct RingBundle {
    pub ring_count: usize,
    pub ring_system_count: usize,
    pub aromatic_ring_count: usize,
    pub num_aliphatic_rings: usize,
    pub num_saturated_rings: usize,
    pub num_aromatic_heterocycles: usize,
    pub num_aliphatic_heterocycles: usize,
    pub num_saturated_heterocycles: usize,
    pub num_spiro_atoms: usize,
    pub num_bridgehead_atoms: usize,
    pub rotatable_bond_count: usize,
    pub hba_count: usize,
    pub fraction_rotatable_bonds: f64,
}

/// Compute all ring-derived descriptors with a single `find_sssr` call.
///
/// When `mol.descriptors()` or `molecule_report()` needs multiple ring values,
/// call this once and read fields — saves ~10-15 redundant SSSR computations per molecule.
pub fn ring_bundle(mol: &Molecule) -> RingBundle {
    let sssr = find_sssr(mol);
    let rings = sssr.rings();
    let all_rings = all_ring_list(mol);

    // Aromatic ring count: use aromatic_ring_list (augmented + envelope-stripped).
    let aromatic_ring_count = aromatic_ring_list(mol).len();

    let rotatable_bond_count = rotatable_bond_count(mol);
    let hba_count = hba_count(mol);
    let hac = heavy_atom_count(mol);
    let fraction_rotatable_bonds = if hac == 0 {
        0.0
    } else {
        rotatable_bond_count as f64 / hac as f64
    };

    RingBundle {
        ring_count: rings.len(),
        ring_system_count: find_ring_families(mol, &sssr).len(),
        aromatic_ring_count,
        num_aliphatic_rings: all_rings
            .iter()
            .filter(|r| r.iter().any(|&i| !mol.atom(i).aromatic))
            .count(),
        num_saturated_rings: all_rings
            .iter()
            .filter(|r| is_ring_saturated(mol, r))
            .count(),
        num_aromatic_heterocycles: all_rings
            .iter()
            .filter(|r| {
                r.iter().all(|&i| mol.atom(i).aromatic)
                    && r.iter().any(|&i| {
                        let an = mol.atom(i).element.atomic_number();
                        an != 6 && an != 1
                    })
            })
            .count(),
        num_aliphatic_heterocycles: all_rings
            .iter()
            .filter(|r| {
                r.iter().any(|&i| !mol.atom(i).aromatic)
                    && r.iter().any(|&i| {
                        let an = mol.atom(i).element.atomic_number();
                        an != 6 && an != 1
                    })
            })
            .count(),
        num_saturated_heterocycles: all_rings
            .iter()
            .filter(|r| {
                is_ring_saturated(mol, r)
                    && r.iter().any(|&i| {
                        let an = mol.atom(i).element.atomic_number();
                        an != 6 && an != 1
                    })
            })
            .count(),
        num_spiro_atoms: num_spiro_atoms_from(mol, &all_rings),
        num_bridgehead_atoms: num_bridgehead_atoms_from(mol, &all_rings),
        rotatable_bond_count,
        hba_count,
        fraction_rotatable_bonds,
    }
}

/// Number of potential tetrahedral stereocenters (specified + unspecified).
///
/// Matches RDKit `CalcNumAtomStereoCenters`: counts all atoms whose CIP priorities
/// for their four substituents are all distinct, regardless of whether @/@@ is
/// specified in the input SMILES.
pub fn num_stereocenters(mol: &Molecule) -> usize {
    // RDKit's count (`CalcNumAtomStereoCenters`) comes from its legacy
    // stereo perception, which e.g. leaves symmetric bridgehead amines
    // unflagged; use the port of it, falling back to the CIP-based
    // potential centers if the molecule is outside the port's scope.
    match chematic_smiles::rdkit_atom_stereocenter_counts(mol) {
        Ok((total, _)) => total,
        Err(_) => potential_stereocenter_indices(mol).len(),
    }
}

/// Return the atom indices of potential tetrahedral stereocenters.
///
/// This is the atom-level counterpart of [`num_stereocenters`].  Keeping the
/// selected atoms available prevents a count-only comparison from hiding a
/// false positive and false negative that cancel each other out.
pub fn potential_stereocenter_indices(mol: &Molecule) -> Vec<AtomIdx> {
    use chematic_core::CipCode;
    use std::collections::HashMap;

    // Pass 1: graph-only provisional R/S assignments.
    let provisional: HashMap<_, CipCode> = crate::cip::assign_cip(mol)
        .assignments
        .into_iter()
        .filter(|(_, c)| matches!(c, CipCode::R | CipCode::S))
        .collect();

    // Pass 2: count with Rule 5 tie-breaking for graph-tied atoms.
    mol.atoms()
        .filter(|(idx, _)| crate::cip::is_potential_stereocenter_rule5(mol, *idx, &provisional))
        .map(|(idx, _)| idx)
        .collect()
}

/// Number of unspecified (undefined) stereocenters.
///
/// Counts potential tetrahedral stereocenters whose chirality is not specified
/// (no @/@@ in the input SMILES). Reuse the same Rule-5-aware potential-center
/// calculation as [`num_stereocenters`] so ordinary CH2/CH3 atoms with
/// repeated substituents are not counted as stereocenters.
pub fn num_unspecified_stereocenters(mol: &Molecule) -> usize {
    use chematic_core::Chirality;
    // As [`num_stereocenters`]: RDKit's legacy perception when possible.
    if let Ok((_, unspecified)) = chematic_smiles::rdkit_atom_stereocenter_counts(mol) {
        return unspecified;
    }
    potential_stereocenter_indices(mol)
        .into_iter()
        .filter(|idx| mol.atom(*idx).chirality == Chirality::None)
        .count()
}

// ---------------------------------------------------------------------------
// 14. Drug-likeness filters
// ---------------------------------------------------------------------------

/// Veber (2002) oral bioavailability filter.
///
/// Passes when TPSA ≤ 140 Å² **and** rotatable bonds ≤ 10.
pub fn veber_passes(mol: &Molecule) -> bool {
    tpsa(mol) <= 140.0 && rotatable_bond_count(mol) <= 10
}

/// Med-Chem Friendly (MCF) composite filter.
///
/// Returns `true` when the compound passes all common quality gates used in
/// medicinal chemistry triage, mirroring the "MCF" concept in the
/// `medchemfilters` Python library:
///   - No PAINS structural alerts
///   - No Brenk reactive-group alerts
///   - Lipinski Ro5 (MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10)
///   - Veber oral bioavailability (TPSA ≤ 140 Å², RotBonds ≤ 10)
pub fn mcf_passes(mol: &Molecule) -> bool {
    crate::pains_passes(mol)
        && crate::brenk_passes(mol)
        && lipinski_passes(mol)
        && veber_passes(mol)
}

/// Egan (2000) absorption/permeability filter ("Egg model").
///
/// Passes when TPSA ≤ 131.6 Å² **and** Crippen LogP ≤ 5.88.
pub fn egan_passes(mol: &Molecule) -> bool {
    tpsa(mol) <= 131.6 && logp_crippen(mol) <= 5.88
}

/// REOS (Rapid Elimination Of Swill) filter.
///
/// Six criteria for hit-identification library quality:
/// MW 200–500, LogP −5 to 5, HBD 0–5, HBA 0–10, formal charge −2 to 2,
/// rotatable bonds 0–8, heavy atoms 15–50.
pub fn reos_passes(mol: &Molecule) -> bool {
    let mw = molecular_weight(mol);
    let lp = logp_crippen(mol);
    let hbd = hbd_count(mol) as i32;
    let hba = hba_count(mol) as i32;
    let fc = formal_charge_sum(mol);
    let rotb = rotatable_bond_count(mol) as i32;
    let hac = heavy_atom_count(mol) as i32;

    (200.0..=500.0).contains(&mw)
        && (-5.0..=5.0).contains(&lp)
        && (0..=5).contains(&hbd)
        && (0..=10).contains(&hba)
        && (-2..=2).contains(&fc)
        && (0..=8).contains(&rotb)
        && (15..=50).contains(&hac)
}

/// Ghose (1999) drug-likeness filter.
///
/// Four criteria: MW 160–480, LogP −0.4 to 5.6, heavy atoms 20–70,
/// Molar Refractivity 40–130.
pub fn ghose_passes(mol: &Molecule) -> bool {
    let mw = molecular_weight(mol);
    let lp = logp_crippen(mol);
    let hac = heavy_atom_count(mol) as f64;
    let mr = molar_refractivity(mol);

    (160.0..=480.0).contains(&mw)
        && (-0.4..=5.6).contains(&lp)
        && (20.0..=70.0).contains(&hac)
        && (40.0..=130.0).contains(&mr)
}

// ---------------------------------------------------------------------------
// Fragment / Lead-like filters
// ---------------------------------------------------------------------------

/// Rule of Three (Ro3) — fragment-based drug discovery filter (Congreve 2003).
///
/// Passes when MW ≤ 300, LogP ≤ 3, HBD ≤ 3, HBA ≤ 3, and RotBonds ≤ 3.
/// Used to screen fragment libraries for FBDD campaigns.
pub fn ro3_passes(mol: &Molecule) -> bool {
    molecular_weight(mol) <= 300.0
        && logp_crippen(mol) <= 3.0
        && hbd_count(mol) <= 3
        && hba_count(mol) <= 3
        && rotatable_bond_count(mol) <= 3
}

/// Lead-like filter (Oprea 2001).
///
/// Passes when MW ≤ 450, LogP −3.5 to 4.5, RotBonds ≤ 10, and RingCount 1–4.
/// Lead-like compounds have lower MW and LogP than typical drug candidates,
/// leaving room for optimisation-related property increases.
pub fn lead_like_passes(mol: &Molecule) -> bool {
    let mw = molecular_weight(mol);
    let lp = logp_crippen(mol);
    let rotb = rotatable_bond_count(mol);
    let rings = ring_count(mol);

    mw <= 450.0 && (-3.5..=4.5).contains(&lp) && rotb <= 10 && (1..=4).contains(&rings)
}

/// Pfizer 3/75 filter (Leeson & Springthorpe 2007).
///
/// Returns `true` when the compound is NOT in the high-metabolic-liability zone,
/// i.e. the combination `LogP > 3 AND TPSA < 75` is **absent**.
/// Compounds in the "3/75" zone are more likely to be CYP3A4 substrates with
/// high metabolic clearance.
pub fn pfizer_3_75_passes(mol: &Molecule) -> bool {
    !(logp_crippen(mol) > 3.0 && tpsa(mol) < 75.0)
}

// ---------------------------------------------------------------------------
// CNS MPO Score
// ---------------------------------------------------------------------------

/// Central Nervous System Multi-Parameter Optimisation (CNS MPO) score (Wager 2010).
///
/// Combines six desirability functions, each returning 0–1, for a total of 0–6.
/// Higher scores indicate better predicted CNS drug-like properties.
///
/// Component properties and desirability thresholds:
/// - cLogP: 1.0 if ≤ 3; 0.0 if ≥ 5; linear between 3–5
/// - cLogD (pH 7.4): 1.0 if ≤ 2; 0.0 if ≥ 4; linear between 2–4
/// - MW: 1.0 if ≤ 360; 0.0 if ≥ 500; linear between 360–500
/// - TPSA: 1.0 if 40–90; 0.0 if <0 or >120; linear in 0–40 and 90–120 shoulders
/// - HBD: 1.0 if 0; 0.0 if ≥ 2; linear between 0–2
/// - pKa (most basic site): 1.0 if ≤ 8; 0.0 if ≥ 10; linear between 8–10
///
/// CNS MPO score accepting pre-computed descriptor values.
///
/// Equivalent to [`cns_mpo_score`] but reuses `logp`, `tpsa`, `mw`, `hbd`, and
/// `pka_b` that the caller has already computed, saving a full Crippen SMARTS pass
/// and a pKa scan.  `mol` is still needed for LogD (ionisation class).
///
/// - `logp`: Crippen LogP (from `logp_crippen` or `logp_and_mr`)
/// - `tpsa`: topological polar surface area
/// - `mw`: molecular weight
/// - `hbd`: H-bond donor count
/// - `pka_b`: most basic pKa (`pka_base(mol).unwrap_or(0.0)`)
pub fn cns_mpo_from_parts(
    mol: &Molecule,
    logp: f64,
    tpsa: f64,
    mw: f64,
    hbd: usize,
    pka_b: f64,
) -> f64 {
    #[inline]
    fn ld(val: f64, lo: f64, hi: f64) -> f64 {
        if val <= lo {
            1.0
        } else if val >= hi {
            0.0
        } else {
            (hi - val) / (hi - lo)
        }
    }
    let d_tpsa = if !(0.0..=120.0).contains(&tpsa) {
        0.0
    } else if tpsa <= 40.0 {
        tpsa / 40.0
    } else if tpsa <= 90.0 {
        1.0
    } else {
        (120.0 - tpsa) / 30.0
    };
    ld(logp, 3.0, 5.0)
        + ld(crate::logd::logd_from_logp(logp, mol, 7.4), 2.0, 4.0)
        + ld(mw, 360.0, 500.0)
        + d_tpsa
        + (1.0 - hbd as f64 / 2.0).clamp(0.0, 1.0)
        + ld(pka_b, 8.0, 10.0)
}

/// Reference: Wager T.T. et al., ACS Chem. Neurosci. 2010, 1, 435–449.
pub fn cns_mpo_score(mol: &Molecule) -> f64 {
    #[inline]
    fn linear_desirability(val: f64, lo: f64, hi: f64) -> f64 {
        if val <= lo {
            1.0
        } else if val >= hi {
            0.0
        } else {
            (hi - val) / (hi - lo)
        }
    }

    // Compute LogP once; reuse it for LogD to avoid double Crippen matching.
    let logp = logp_crippen(mol);
    let d_logp = linear_desirability(logp, 3.0, 5.0);

    // cLogD pH 7.4: optimal ≤ 2, zero ≥ 4
    let d_logd = linear_desirability(crate::logd::logd_from_logp(logp, mol, 7.4), 2.0, 4.0);

    // MW: optimal ≤ 360, zero ≥ 500
    let d_mw = linear_desirability(molecular_weight(mol), 360.0, 500.0);

    // TPSA: plateau 40–90, shoulders 0–40 and 90–120
    let psa = tpsa(mol);
    let d_tpsa = if !(0.0..=120.0).contains(&psa) {
        0.0
    } else if psa <= 40.0 {
        psa / 40.0
    } else if psa <= 90.0 {
        1.0
    } else {
        (120.0 - psa) / 30.0
    };

    // HBD: optimal 0, zero ≥ 2
    let hbd = hbd_count(mol) as f64;
    let d_hbd = (1.0 - hbd / 2.0).clamp(0.0, 1.0);

    // pKa (most basic): optimal ≤ 8, zero ≥ 10; None (no basic site) → pKa ≈ 0 → score 1.0
    let pka_b = crate::pka::pka_base(mol).unwrap_or(0.0);
    let d_pka = linear_desirability(pka_b, 8.0, 10.0);

    d_logp + d_logd + d_mw + d_tpsa + d_hbd + d_pka
}

// ---------------------------------------------------------------------------
// Per-atom property vectors
// ---------------------------------------------------------------------------

/// Per-atom hybridization state, in RDKit's model.
///
/// Returns a `Vec<u8>` of length `mol.atom_count()` with:
///   `1` = sp, `2` = sp2, `3` = sp3, `0` = other (s, sp3d, sp3d2, metals
///   without a defined state, wildcard atoms).
///
/// Values follow RDKit 2026.03's `ConjugHybrid` (the same model as the SMARTS
/// `^n` primitive, [`chematic_smarts::rdkit_hybridization`], evaluated on
/// RDKit's aromaticity view of `mol`): a lone pair
/// counts as an orbital, so an amide or enamine N, an aryl ether O and a
/// carboxylate O⁻ are sp2; five- and six-orbital centres (PF₅, SF₆) are
/// sp3d/sp3d2 and report `0`. Before 1.0.35 this used bond orders only (any
/// atom without a multiple bond was sp3).
pub fn hybridization_per_atom(mol: &Molecule) -> Vec<u8> {
    // Where chematic's perception may disagree with RDKit's sanitization
    // (radicals, unusual valences, ...), take RDKit's hybridization.
    if chematic_perception::rdkit_model_may_disagree(mol)
        && let Some(h) = chematic_smiles::rdkit_hybridizations(mol)
    {
        return h;
    }
    // On RDKit's aromaticity view, as RDKit's sanitized molecule: a
    // Kekulé-written thiophene S is aromatic, hence sp2.
    let view = descriptor_aromaticity(mol);
    let view: &Molecule = &view;
    view.atoms()
        .map(|(idx, _)| {
            chematic_smarts::rdkit_hybridization(view, idx)
                .filter(|h| (1..=3).contains(h))
                .unwrap_or(0)
        })
        .collect()
}

/// Per-atom formal charge.
///
/// Returns a `Vec<i8>` of length `mol.atom_count()` where index `i` is the
/// formal charge of `AtomIdx(i as u32)`.  Neutral atoms have value `0`.
pub fn formal_charge_per_atom(mol: &Molecule) -> Vec<i8> {
    mol.atoms().map(|(_, a)| a.charge).collect()
}

/// Per-atom implicit hydrogen count.
///
/// Returns a `Vec<u8>` of length `mol.atom_count()` where index `i` is the
/// number of implicit H atoms on `AtomIdx(i as u32)`.
/// Reuses [`chematic_core::implicit_hcount`] per atom.
pub fn implicit_hcount_per_atom(mol: &Molecule) -> Vec<u8> {
    mol.atoms()
        .map(|(idx, _)| implicit_hcount(mol, idx))
        .collect()
}

// ---------------------------------------------------------------------------
// TPSA per-atom
// ---------------------------------------------------------------------------

/// Per-atom TPSA contributions (Ertl 2000).
///
/// Returns a `Vec<f64>` of length `mol.atom_count()` where index `i` holds the
/// TPSA contribution of `AtomIdx(i as u32)`.  Atoms that contribute nothing
/// (C, halogens, metals, …) have value `0.0`.  The sum equals `tpsa(mol)`.
///
/// Mirrors the pattern of [`logp_crippen_per_atom`].
pub fn tpsa_per_atom(mol: &Molecule) -> Vec<f64> {
    tpsa_contributions(mol)
}

#[cfg(test)]
#[path = "ring_component_count_tests.rs"]
mod ring_component_count_tests;
