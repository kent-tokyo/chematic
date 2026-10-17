//! Ordered-emission preservation against c84c68f297f95562d1d29c75234f2190cd5f9d75.
#[cfg(test)]
mod tests {
    use super::super::*;
    use chematic_core::{Atom, Element, MoleculeBuilder};

    fn compare(mol: &Molecule) {
        let bonds: Vec<_> = mol
            .bonds()
            .map(|(_, b)| checked_bond_invariant(b.order).unwrap())
            .collect();
        let chirality = MorganChirality {
            tagged: vec![true; mol.atom_count()],
            code: (0..mol.atom_count()).map(|i| 1 + (i % 3) as u32).collect(),
        };
        for rings in [
            vec![false; mol.atom_count()],
            vec![true; mol.atom_count()],
            ring_atom_flags(mol),
        ] {
            for radius in 0..=3 {
                for stereo in [None, Some(&chirality)] {
                    let mut actual = Vec::new();
                    let mut expected = Vec::new();
                    let result = expand_emissions_small_with_chirality(
                        mol,
                        &rings,
                        &bonds,
                        radius,
                        stereo,
                        |a, r, h| actual.push((a, r, h)),
                    );
                    let previous =
                        frozen_binary_search(mol, &rings, &bonds, radius, stereo, |a, r, h| {
                            expected.push((a, r, h))
                        });
                    assert_eq!(result, previous);
                    assert_eq!(
                        actual,
                        expected,
                        "radius={radius}, bonds={}",
                        mol.bond_count()
                    );
                }
            }
        }
    }

    #[test]
    fn ordered_emissions_match_all_five_vertex_graphs() {
        // Every undirected simple graph, including disconnected and symmetric graphs.
        for mask in 0u16..1024 {
            let mut builder = MoleculeBuilder::new();
            let atoms: Vec<_> = (0..5)
                .map(|_| builder.add_atom(Atom::new(Element::C)))
                .collect();
            let mut bit = 0;
            for a in 0..5 {
                for b in a + 1..5 {
                    if mask & (1 << bit) != 0 {
                        builder
                            .add_bond(atoms[a], atoms[b], BondOrder::Single)
                            .unwrap();
                    }
                    bit += 1;
                }
            }
            compare(&builder.build());
        }
    }

    #[test]
    fn ordered_emissions_preserve_spill_stereo_and_fallback_boundaries() {
        compare(&MoleculeBuilder::new().build());
        for smiles in [
            "[He]",
            "C.C.C",
            "CC",
            "CC(=O)O",
            "C1CCCCC1",
            "C1CC2CCC1C2",
            "C[C@H](F)[C@@H](Cl)Br",
            "C[C@H](O)/C=C/[C@@H](N)C(=O)O",
            "[Na+].CC(=O)[O-]",
            "c1ccc2ccccc2c1",
        ] {
            compare(&chematic_smiles::parse(smiles).unwrap());
        }
        // Atom buffers spill above 64; 128 bonds still use the compact path.
        for n in [64, 65, 128, 129, 130] {
            compare(&chematic_smiles::parse(&"C".repeat(n)).unwrap());
        }
    }

    // Complete previous implementation, only the function name/visibility changed.
    fn frozen_binary_search<F>(
        mol: &Molecule,
        ring_atoms: &[bool],
        bond_invariants: &[u32],
        max_radius: u32,
        chirality: Option<&MorganChirality>,
        mut emit: F,
    ) -> Option<()>
    where
        F: FnMut(u32, u32, u32),
    {
        let n = mol.atom_count();
        if mol.bond_count() > 128 {
            return None;
        }
        if n == 0 {
            return Some(());
        }
        let mut current: SmallVec<[u32; 64]> = SmallVec::with_capacity(n);
        for i in 0..n {
            let id = connectivity_invariant(mol, AtomIdx(i as u32), ring_atoms);
            current.push(id);
            emit(i as u32, 0, id);
        }
        // Most drug-like molecules fit in these inline buffers. Larger inputs
        // transparently spill to the heap with unchanged indexing and ordering.
        let mut next: SmallVec<[u32; 64]> = SmallVec::from_elem(0, n);
        let mut dead: SmallVec<[bool; 64]> = SmallVec::from_elem(false, n);
        let mut chiral_atoms: SmallVec<[bool; 64]> = SmallVec::from_elem(false, n);
        let mut env: SmallVec<[u128; 64]> = SmallVec::from_elem(0, n);
        let mut next_env: SmallVec<[u128; 64]> = SmallVec::from_elem(0, n);
        let mut seen: SmallVec<[u128; 128]> = SmallVec::new();
        let mut round: SmallVec<[(u128, u32, u32); 64]> = SmallVec::with_capacity(n);
        let mut pairs: SmallVec<[(u32, u32); 6]> = SmallVec::new();

        for layer in 0..max_radius {
            round.clear();
            for i in 0..n {
                // Atoms that are (or become) dead keep their environment and
                // carry a zero invariant into the next round, as in the general
                // path (whose next-invariant vector starts zeroed every round).
                if dead[i] {
                    next_env[i] = env[i];
                    next[i] = 0;
                    continue;
                }
                let nbrs = mol.neighbor_slice(AtomIdx(i as u32));
                if nbrs.is_empty() {
                    dead[i] = true;
                    next_env[i] = env[i];
                    next[i] = 0;
                    continue;
                }
                let mut bond_env = 0u128;
                pairs.clear();
                for &(nb, bond) in nbrs {
                    bond_env |= 1u128 << bond.0;
                    bond_env |= env[nb.0 as usize];
                    pairs.push((bond_invariants[bond.0 as usize], current[nb.0 as usize]));
                }
                pairs.sort_unstable();
                let mut invar = hash_combine(layer, current[i]);
                let mut looks_chiral = chirality.is_some_and(|c| c.tagged[i]);
                let mut previous_neighbor = None;
                for &(bond_inv, nb_inv) in &pairs {
                    invar = hash_combine(invar, hash_pair(bond_inv, nb_inv));
                    if looks_chiral
                        && !chiral_atoms[i]
                        && (bond_inv != 1 || previous_neighbor == Some(nb_inv))
                    {
                        looks_chiral = false;
                    }
                    previous_neighbor = Some(nb_inv);
                }
                if looks_chiral && let Some(c) = chirality {
                    chiral_atoms[i] = true;
                    invar = hash_combine(invar, c.code[i]);
                }
                next[i] = invar;
                next_env[i] = bond_env;
                round.push((bond_env, invar, i as u32));
            }
            round.sort_unstable();
            let seen_before = seen.len();
            let mut g = 0;
            while g < round.len() {
                let bond_env = round[g].0;
                let mut end = g + 1;
                while end < round.len() && round[end].0 == bond_env {
                    end += 1;
                }
                // Sorted by (env, invariant, atom): round[g] is the group winner.
                if seen[..seen_before].binary_search(&bond_env).is_ok() {
                    for &(_, _, atom) in &round[g..end] {
                        dead[atom as usize] = true;
                    }
                } else {
                    let (_, winner_invariant, winner) = round[g];
                    emit(winner, layer + 1, winner_invariant);
                    seen.push(bond_env);
                    for &(_, _, atom) in &round[g + 1..end] {
                        dead[atom as usize] = true;
                    }
                }
                g = end;
            }
            seen.sort_unstable();
            std::mem::swap(&mut current, &mut next);
            std::mem::swap(&mut env, &mut next_env);
        }
        Some(())
    }
}
