//! Diagnostic-only probe for issue #337's symmetrized-SSSR boundary.
//!
//! It compares the existing RDKit-D2-root candidate population with the
//! candidate population obtained by probing every atom as a root. This does
//! not change production SSSR behavior; it tests whether the current D2-root
//! reduction is the reason the six MMFF94 residual macrocycles are missing.
//!
//! Run with:
//! `cargo run --release -p chematic-perception --example issue337_all_root_cycle_probe`

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use chematic_core::{AtomIdx, BondIdx, BondOrder};
use chematic_perception::{find_smallest_rings_bfs, find_sssr, select_rdkit_d2_roots};

const FIXTURES: &[(&str, &str)] = &[
    (
        "chembl_tier_b_0009",
        "c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2",
    ),
    (
        "chembl_tier_b_0023",
        "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)-c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
    ),
    (
        "chembl_tier_b_0028",
        "C1=C\\c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCNc2cc[n+](c3ccccc23)Cc2ccc/1cc2",
    ),
    (
        "chembl_tier_b_0029",
        "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)CCc2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
    ),
    (
        "chembl_tier_b_0030",
        "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)Cc2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
    ),
    (
        "chembl_tier_b_0034",
        "c1ccc2c(c1)c1cc[n+]2Cc2ccc3c(c2)Cc2cc(ccc2-3)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
    ),
];

fn canonical_atom_set(ring: &[AtomIdx]) -> Vec<u32> {
    let mut atoms: Vec<u32> = ring.iter().map(|atom| atom.0).collect();
    atoms.sort_unstable();
    atoms
}

fn all_root_candidates(mol: &chematic_core::Molecule) -> BTreeSet<Vec<u32>> {
    (0..mol.atom_count())
        .flat_map(|raw| find_smallest_rings_bfs(mol, AtomIdx(raw as u32)))
        .map(|ring| canonical_atom_set(&ring))
        .collect()
}

fn d2_root_candidates(mol: &chematic_core::Molecule) -> BTreeSet<Vec<u32>> {
    select_rdkit_d2_roots(mol)
        .into_iter()
        .flat_map(|root| find_smallest_rings_bfs(mol, root))
        .map(|ring| canonical_atom_set(&ring))
        .collect()
}

fn shortest_path_without_bond(
    mol: &chematic_core::Molecule,
    start: AtomIdx,
    goal: AtomIdx,
    blocked: BondIdx,
) -> Option<Vec<AtomIdx>> {
    let mut parent = vec![None; mol.atom_count()];
    let mut seen = vec![false; mol.atom_count()];
    let mut queue = VecDeque::new();
    seen[start.0 as usize] = true;
    queue.push_back(start);
    while let Some(current) = queue.pop_front() {
        if current == goal {
            break;
        }
        for (next, bond) in mol.neighbors(current) {
            if bond == blocked || seen[next.0 as usize] {
                continue;
            }
            seen[next.0 as usize] = true;
            parent[next.0 as usize] = Some(current);
            queue.push_back(next);
        }
    }
    if !seen[goal.0 as usize] {
        return None;
    }
    let mut path = vec![goal];
    let mut current = goal;
    while current != start {
        current = parent[current.0 as usize]?;
        path.push(current);
    }
    path.reverse();
    Some(path)
}

fn edge_exchange_same_size_cycles(
    mol: &chematic_core::Molecule,
    base_ring: &[AtomIdx],
) -> BTreeSet<Vec<u32>> {
    let mut alternatives = BTreeSet::new();
    for i in 0..base_ring.len() {
        let left = base_ring[i];
        let right = base_ring[(i + 1) % base_ring.len()];
        let Some((blocked, _)) = mol.bond_between(left, right) else {
            continue;
        };
        let Some(path) = shortest_path_without_bond(mol, left, right, blocked) else {
            continue;
        };
        if path.len() == base_ring.len() {
            alternatives.insert(canonical_atom_set(&path));
        }
    }
    alternatives
}

fn ring_bonds(mol: &chematic_core::Molecule, ring: &[AtomIdx]) -> BTreeSet<BondIdx> {
    (0..ring.len())
        .filter_map(|i| {
            mol.bond_between(ring[i], ring[(i + 1) % ring.len()])
                .map(|(b, _)| b)
        })
        .collect()
}

fn accepted_by_existing_symmetrized_contract(
    mol: &chematic_core::Molecule,
    base: &[Vec<AtomIdx>],
    candidate: &[u32],
) -> bool {
    let candidate_atoms: Vec<_> = candidate.iter().copied().map(AtomIdx).collect();
    let candidate_bonds = ring_bonds(mol, &candidate_atoms);
    let base_bonds: Vec<_> = base.iter().map(|ring| ring_bonds(mol, ring)).collect();
    let mut bond_ring_count = std::collections::BTreeMap::<BondIdx, usize>::new();
    for ring in &base_bonds {
        for &bond in ring {
            *bond_ring_count.entry(bond).or_default() += 1;
        }
    }
    base_bonds.iter().any(|basis| {
        basis.iter().any(|bond| candidate_bonds.contains(bond))
            && basis.iter().all(|bond| {
                bond_ring_count.get(bond).copied().unwrap_or(0) != 1
                    || candidate_bonds.contains(bond)
            })
    }) && base.iter().any(|ring| ring.len() == candidate.len())
}

fn gf2_rank(rows: &[BTreeSet<BondIdx>]) -> usize {
    let mut pivots = BTreeMap::<BondIdx, BTreeSet<BondIdx>>::new();
    let mut rank = 0;
    for row in rows {
        let mut current = row.clone();
        while let Some(&pivot) = current.iter().next_back() {
            if let Some(existing) = pivots.get(&pivot) {
                let mut reduced = current;
                for bond in existing {
                    if !reduced.insert(*bond) {
                        reduced.remove(bond);
                    }
                }
                current = reduced;
            } else {
                pivots.insert(pivot, current);
                rank += 1;
                break;
            }
        }
    }
    rank
}

fn has_basis_exchange(
    mol: &chematic_core::Molecule,
    base: &[Vec<AtomIdx>],
    candidate: &[u32],
) -> bool {
    let base_bonds: Vec<_> = base.iter().map(|ring| ring_bonds(mol, ring)).collect();
    let candidate_atoms: Vec<_> = candidate.iter().copied().map(AtomIdx).collect();
    let candidate_bonds = ring_bonds(mol, &candidate_atoms);
    let base_rank = gf2_rank(&base_bonds);
    base.iter().enumerate().any(|(index, ring)| {
        if ring.len() != candidate.len() {
            return false;
        }
        let mut replaced = base_bonds.clone();
        replaced[index] = candidate_bonds.clone();
        gf2_rank(&replaced) == base_rank
    })
}

fn enumerate_exact_cycles(
    mol: &chematic_core::Molecule,
    target_len: usize,
    max_cycles: usize,
) -> (BTreeSet<Vec<u32>>, bool) {
    #[allow(clippy::too_many_arguments)]
    fn dfs(
        mol: &chematic_core::Molecule,
        start: AtomIdx,
        current: AtomIdx,
        target_len: usize,
        path: &mut Vec<AtomIdx>,
        seen: &mut Vec<bool>,
        found: &mut BTreeSet<Vec<u32>>,
        max_cycles: usize,
    ) -> bool {
        if found.len() >= max_cycles {
            return true;
        }
        if path.len() == target_len {
            let closes = mol.neighbors(current).any(|(next, bond)| {
                next == start
                    && !matches!(mol.bond(bond).order, BondOrder::Zero | BondOrder::Dative)
            });
            if closes && path.iter().all(|atom| atom.0 >= start.0) {
                found.insert(canonical_atom_set(path));
            }
            return false;
        }
        for (next, bond) in mol.neighbors(current) {
            if matches!(mol.bond(bond).order, BondOrder::Zero | BondOrder::Dative)
                || seen[next.0 as usize]
            {
                continue;
            }
            seen[next.0 as usize] = true;
            path.push(next);
            if dfs(mol, start, next, target_len, path, seen, found, max_cycles) {
                return true;
            }
            path.pop();
            seen[next.0 as usize] = false;
        }
        false
    }

    let mut found = BTreeSet::new();
    let mut capped = false;
    for raw in 0..mol.atom_count() {
        let start = AtomIdx(raw as u32);
        let mut path = vec![start];
        let mut seen = vec![false; mol.atom_count()];
        seen[raw] = true;
        if dfs(
            mol, start, start, target_len, &mut path, &mut seen, &mut found, max_cycles,
        ) {
            capped = true;
            break;
        }
    }
    (found, capped)
}

fn main() {
    println!("issue #337 all-root cycle probe (diagnostic only)");
    for &(name, smiles) in FIXTURES {
        let mol = chematic_smiles::parse(smiles).unwrap_or_else(|err| panic!("{name}: {err}"));
        let base = find_sssr(&mol);
        let all = all_root_candidates(&mol);
        let d2 = d2_root_candidates(&mol);
        let missing_from_d2: BTreeSet<_> = all.difference(&d2).cloned().collect();
        let base_sets: BTreeSet<_> = base
            .rings()
            .iter()
            .map(|ring| canonical_atom_set(ring))
            .collect();
        let base_macrocycles: Vec<_> = base.rings().iter().filter(|ring| ring.len() >= 9).collect();
        let exchanged: BTreeSet<_> = base_macrocycles
            .iter()
            .flat_map(|ring| edge_exchange_same_size_cycles(&mol, ring))
            .filter(|candidate| !base_sets.contains(candidate))
            .collect();
        let accepted_exchanged: Vec<_> = exchanged
            .iter()
            .filter(|candidate| {
                accepted_by_existing_symmetrized_contract(&mol, base.rings(), candidate)
            })
            .map(|candidate| candidate.len())
            .collect();
        let mut exact_cycle_counts = Vec::new();
        let mut exact_cycle_capped = false;
        for ring in &base_macrocycles {
            let (cycles, capped) = enumerate_exact_cycles(&mol, ring.len(), 2_000_000);
            let alternatives: Vec<_> = cycles
                .iter()
                .filter(|candidate| !base_sets.contains(*candidate))
                .collect();
            let accepted = alternatives
                .iter()
                .filter(|candidate| {
                    accepted_by_existing_symmetrized_contract(&mol, base.rings(), candidate)
                })
                .count();
            let basis_exchange = alternatives
                .iter()
                .filter(|candidate| has_basis_exchange(&mol, base.rings(), candidate))
                .count();
            exact_cycle_counts.push((
                ring.len(),
                cycles.len(),
                alternatives.len(),
                accepted,
                basis_exchange,
            ));
            exact_cycle_capped |= capped;
        }
        let new_same_size: Vec<_> = missing_from_d2
            .iter()
            .filter(|candidate| base_sets.iter().any(|base| base.len() == candidate.len()))
            .collect();
        println!(
            "{name}: base_sizes={:?} d2_candidates={} all_root_candidates={} missing_from_d2={} missing_same_size={:?} edge_exchange_same_size={:?} accepted_by_existing_contract={:?} exact_cycle_counts={:?} exact_cycle_capped={exact_cycle_capped}",
            base.rings().iter().map(Vec::len).collect::<Vec<_>>(),
            d2.len(),
            all.len(),
            missing_from_d2.len(),
            new_same_size
                .iter()
                .map(|ring| ring.len())
                .collect::<Vec<_>>(),
            exchanged.iter().map(|ring| ring.len()).collect::<Vec<_>>(),
            accepted_exchanged,
            exact_cycle_counts,
        );
    }
}
