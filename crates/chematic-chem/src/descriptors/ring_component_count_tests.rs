//! Differential oracle for cyclic-component reuse.
use super::*;

fn reference_count(mol: &Molecule) -> usize {
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
    let mut seen = vec![false; mol.atom_count()];
    let mut stack = vec![closure.atom1];
    seen[closure.atom1.0 as usize] = true;
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
            let visited = &mut seen[neighbor.0 as usize];
            if !*visited {
                *visited = true;
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
    let n = mol.atom_count();
    let mut parent: Vec<u32> = (0..n as u32).collect();
    fn find(parent: &mut [u32], mut x: u32) -> u32 {
        while parent[x as usize] != x {
            let p = parent[parent[x as usize] as usize];
            parent[x as usize] = p;
            x = p;
        }
        x
    }
    let mut ring_bonds = 0usize;
    for (bi, bond) in mol.bonds() {
        if !flags[bi.0 as usize] {
            continue;
        }
        ring_bonds += 1;
        let (a, b) = (
            find(&mut parent, bond.atom1.0),
            find(&mut parent, bond.atom2.0),
        );
        if a != b {
            parent[a as usize] = b;
        }
    }
    if ring_bonds == 0 {
        return Some(RingCountDecision::Count(0));
    }
    // Per system: ring atoms minus ring bonds (0 for a simple cycle).
    let mut balance = vec![0i32; n];
    let mut on_ring = vec![false; n];
    for (bi, bond) in mol.bonds() {
        if flags[bi.0 as usize] {
            on_ring[bond.atom1.0 as usize] = true;
            on_ring[bond.atom2.0 as usize] = true;
            let r = find(&mut parent, bond.atom1.0);
            balance[r as usize] -= 1;
        }
    }
    let mut roots = Vec::new();
    for (a, &is_on_ring) in on_ring.iter().enumerate() {
        if is_on_ring {
            let r = find(&mut parent, a as u32) as usize;
            balance[r] += 1;
            if r == a {
                roots.push(a);
            }
        }
    }
    let mut rings = 0usize;
    let mut needs_relevant_cycles = false;
    for &r in &roots {
        match balance[r] {
            0 => rings += 1,
            -1 => rings += two_cycle_system_rings(mol, &flags, &mut parent, r)?,
            _ => {
                let (count, has_macrocycle) = system_rings_rdkit(mol, &flags, &mut parent, r)?;
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
    parent: &mut [u32],
    root: usize,
) -> Option<(usize, bool)> {
    fn find(parent: &mut [u32], mut x: u32) -> u32 {
        while parent[x as usize] != x {
            let p = parent[parent[x as usize] as usize];
            parent[x as usize] = p;
            x = p;
        }
        x
    }
    let mut local = vec![u32::MAX; mol.atom_count()];
    let mut k = 0usize;
    let mut bonds = Vec::new();
    for bi in mol.rdkit_bond_order() {
        if !flags[bi.0 as usize] {
            continue;
        }
        let bond = mol.bond(bi);
        if find(parent, bond.atom1.0) as usize != root {
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
    parent: &mut [u32],
    root: usize,
) -> Option<usize> {
    fn find(parent: &mut [u32], mut x: u32) -> u32 {
        while parent[x as usize] != x {
            let p = parent[parent[x as usize] as usize];
            parent[x as usize] = p;
            x = p;
        }
        x
    }
    let ring_neighbors = |a: AtomIdx| -> smallvec::SmallVec<[AtomIdx; 4]> {
        mol.neighbors(a)
            .filter(|(_, bi)| flags[bi.0 as usize])
            .map(|(nb, _)| nb)
            .collect()
    };
    let mut branch = smallvec::SmallVec::<[AtomIdx; 2]>::new();
    for (idx, _) in mol.atoms() {
        if find(parent, idx.0) as usize != root {
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

fn compare(input: &Molecule) {
    // Fresh clones drop derived caches, then exercise both cold and warm
    // calls. Compare again after graph edits below to cover invalidation.
    let expected = reference_count(&input.clone());
    assert_eq!(rdkit_num_rings(input), expected);
    assert_eq!(rdkit_num_rings(input), expected);
}

#[test]
fn fixed_ten_thousand_ring_counts_match_reference() {
    let corpus = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
    ));
    let mut rows = 0;
    for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
        let mol = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
        compare(&mol);
        rows += 1;
    }
    assert_eq!(rows, 10000);
}

#[test]
fn edited_graphs_keep_ring_counts_and_invalidate_components() {
    for smiles in [
        "CCC",
        "C1CCCCC1",
        "C1CCCCC1.C1CC1",
        "C1CCC2(CC1)CCCC2",
        "C1CC2CCC1C2",
        "C12C3C4C1C5C2C3C45",
        "C1=CC2=CC=CC=C2C=C1",
        "N1CCN->[Cu]1",
        "[Fe]1(Cl)(Cl)Cl.Cl1",
    ] {
        let mut mol = chematic_smiles::parse_template(smiles).unwrap();
        compare(&mol);
        let before = chematic_perception::sssr::ring_components_shared(&mol);
        let branch = mol.add_atom(chematic_core::Atom::new(Element::C));
        mol.add_bond(AtomIdx(0), branch, BondOrder::Single).unwrap();
        compare(&mol);
        let after = chematic_perception::sssr::ring_components_shared(&mol);
        assert!(!std::sync::Arc::ptr_eq(&before, &after));
        if mol.atom_count() > 3 {
            mol.add_bond(branch, AtomIdx(2), BondOrder::Single).unwrap();
            compare(&mol);
        }
        let bond = mol.bond(BondIdx(0)).clone();
        mol.remove_bond(BondIdx(0));
        compare(&mol);
        mol.add_bond(bond.atom1, bond.atom2, BondOrder::QueryAny)
            .unwrap();
        compare(&mol);
    }
}

#[test]
fn seeded_bond_flags_without_components_keep_the_partition() {
    let mol = chematic_smiles::parse_template("C1CCC2(CC1)CCC2.CC").unwrap();
    let flags = chematic_perception::ring_bond_flags_shared(&mol);
    let copy = mol.clone();
    copy.seed_derived(chematic_core::DerivedSlot::RingBondFlags, flags);
    compare(&copy);
}

#[test]
fn every_five_atom_simple_graph_matches_reference() {
    let edges: Vec<_> = (0..5)
        .flat_map(|a| (a + 1..5).map(move |b| (a, b)))
        .collect();
    for mask in 0u16..(1 << edges.len()) {
        let mut builder = chematic_core::MoleculeBuilder::new();
        for _ in 0..5 {
            builder.add_atom(chematic_core::Atom::new(Element::C));
        }
        for (i, &(a, b)) in edges.iter().enumerate() {
            if mask & (1 << i) != 0 {
                builder
                    .add_bond(AtomIdx(a), AtomIdx(b), BondOrder::Single)
                    .unwrap();
            }
        }
        let mol = builder.build();
        compare(&mol);
    }
}
