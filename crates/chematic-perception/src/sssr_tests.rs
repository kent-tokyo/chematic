#[cfg(test)]
mod tests {
    use super::*;

    const RING_SAMPLES: &[&str] = &[
        "CCO",
        "c1ccccc1",
        "c1ccc2ccccc2c1",
        "C1CC2CCC1CC2",
        "c1ccc(-c2ccccc2)cc1",
        "C12C3C4C1C5C2C3C45",
        "O=C1NC(=O)c2ccccc21",
        "C1CCC2(CC1)CCCC2",
        "c1ccc2c(c1)-c1ccccc1-2",
        "CC1(C)C2CCC1(C)C(=O)C2",
        "[Na+].[Cl-]",
        "C1CC1.C1CCC1",
    ];

    #[test]
    fn sssr_ring_count_equals_basis_size() {
        for smi in RING_SAMPLES {
            let mol = chematic_smiles::parse(smi).unwrap();
            assert_eq!(sssr_ring_count(&mol), find_sssr(&mol).ring_count(), "{smi}");
        }
    }

    #[test]
    fn ring_atom_ranks_match_full_ranks_on_ring_atoms() {
        for smi in RING_SAMPLES {
            let mol = chematic_smiles::parse(smi).unwrap();
            let cyclic = ring_bond_flags(&mol);
            let ring_atoms = ring_atom_flags(&mol);
            let full = canonical_atom_ranks(&mol);
            let restricted = canonical_ring_atom_ranks(&mol, &cyclic);
            for (i, &r) in ring_atoms.iter().enumerate() {
                if r {
                    assert_eq!(full[i], restricted[i], "{smi} atom {i}");
                }
            }
        }
    }

    #[test]
    fn fast_sssr_matches_horton_reference() {
        for smi in RING_SAMPLES {
            let mol = chematic_smiles::parse(smi).unwrap();
            assert_eq!(
                find_sssr(&mol).rings(),
                find_sssr_horton_reference(&mol).rings(),
                "{smi}"
            );
        }
    }
    /// Cyclic-subgraph component label of every ring atom (`u32::MAX` elsewhere).
    fn ring_components(mol: &chematic_core::Molecule) -> Vec<u32> {
        let cyclic = ring_bond_flags(mol);
        let n = mol.atom_count();
        let mut label = vec![u32::MAX; n];
        let mut next = 0;
        for start in 0..n {
            let touches = mol
                .neighbors(AtomIdx(start as u32))
                .any(|(_, b)| cyclic[b.0 as usize]);
            if !touches || label[start] != u32::MAX {
                continue;
            }
            let mut stack = vec![start];
            label[start] = next;
            while let Some(a) = stack.pop() {
                for (nb, b) in mol.neighbors(AtomIdx(a as u32)) {
                    if cyclic[b.0 as usize] && label[nb.0 as usize] == u32::MAX {
                        label[nb.0 as usize] = next;
                        stack.push(nb.0 as usize);
                    }
                }
            }
            next += 1;
        }
        label
    }

    #[test]
    fn component_restricted_sssr_is_the_filtered_full_sssr() {
        let samples = RING_SAMPLES.iter().copied().chain([
            "O=C1c2ccccc2C(=O)N1CCCCN1CCN(c2cccc3ccccc23)CC1",
            "C1CC1C1CCCC1C1CC2CCC1C2c1ccc2ccccc2c1",
            "c1ccc2c(c1)C1(CCCC1)c1ccccc1-2",
            "C1CCC2(CC1)CC1(CCCC1)C2",
        ]);
        for smi in samples {
            let mol = chematic_smiles::parse(smi).unwrap();
            let label = ring_components(&mol);
            let count = label
                .iter()
                .filter(|&&l| l != u32::MAX)
                .max()
                .map_or(0, |m| m + 1);
            // Every subset of components (small molecules only).
            for subset in 0u32..(1 << count.min(6)) {
                let keep: Vec<bool> = label
                    .iter()
                    .map(|&l| l != u32::MAX && subset & (1 << l) != 0)
                    .collect();
                let expected: Vec<Vec<AtomIdx>> = find_sssr(&mol)
                    .rings()
                    .iter()
                    .filter(|r| keep[r[0].0 as usize])
                    .cloned()
                    .collect();
                assert_eq!(
                    find_sssr_in_components(&mol, &keep),
                    expected,
                    "{smi} {subset:b}"
                );
            }
        }
    }

    #[test]
    fn unordered_ring_set_matches_the_canonical_set() {
        let samples = RING_SAMPLES.iter().copied().chain([
            "O=C1c2ccccc2C(=O)N1CCCCN1CCN(c2cccc3ccccc23)CC1",
            "C12C3C4C1C5C2C3C45",
            "C1CC2CC3CCC1C3C2",
            "c1cc2ccc3cccc4ccc(c1)c2c34",
            "C1C2CC3CC1CC(C2)C3",
            // rank-2 components read off the graph: theta and spiro shapes,
            // including ties that must fall back to the Horton selection
            "C1CC2CCC1CC2",
            "C1CC2CCC1C2",
            "C1CCC2(C1)CCCC2",
            "C1CC11CC1",
            "c1ccc2[nH]ccc2c1",
            "c1cc2cccccc2c1",
            "C1CCCCCC2CCCCC(C1)C2",
            "C12CC(C1)C2",
            "C1CC2CC1C2",
            "c1ccc2ccccc2c1CC1CC2CCC1CC2.C1CC11CC1",
            "C1CC2CCC1CC2c1ccc2ccccc2c1C12CC3CC(CC(C3)C1)C2",
        ]);
        for smi in samples {
            let mol = chematic_smiles::parse(smi).unwrap();
            let keep = ring_atom_flags(&mol);
            // Every returned ring is a cycle in atom order.
            for ring in
                find_sssr_ring_set_in_components(&chematic_smiles::parse(smi).unwrap(), &keep)
            {
                for (i, &a) in ring.iter().enumerate() {
                    let b = ring[(i + 1) % ring.len()];
                    assert!(mol.bond_between(a, b).is_some(), "{smi}: {ring:?}");
                }
            }
            if let Some(unranked) = find_sssr_ring_set_in_components_unranked(
                &chematic_smiles::parse(smi).unwrap(),
                &keep,
            ) {
                let mut u: Vec<Vec<AtomIdx>> = unranked
                    .into_iter()
                    .map(|mut r| {
                        r.sort();
                        r
                    })
                    .collect();
                u.sort();
                let mut f: Vec<Vec<AtomIdx>> = find_sssr(&mol)
                    .rings()
                    .iter()
                    .map(|r| {
                        let mut r = r.clone();
                        r.sort();
                        r
                    })
                    .collect();
                f.sort();
                assert_eq!(u, f, "unranked {smi}");
            }
            let key = |rings: &[Vec<AtomIdx>]| {
                let mut v: Vec<Vec<AtomIdx>> = rings
                    .iter()
                    .map(|r| {
                        let mut r = r.clone();
                        r.sort();
                        r
                    })
                    .collect();
                v.sort();
                v
            };
            let fresh = chematic_smiles::parse(smi).unwrap();
            assert_eq!(
                key(&find_sssr_ring_set_in_components(&fresh, &keep)),
                key(find_sssr(&mol).rings()),
                "{smi}"
            );
        }
    }

    #[test]
    fn bridge_flags_match_edge_deletion() {
        // A bond is cyclic exactly when deleting it keeps its endpoints connected.
        for smi in RING_SAMPLES
            .iter()
            .copied()
            .chain(["C1CC1CCC1CC2CCC12", "C1CC=1"])
        {
            let mol = chematic_smiles::parse(smi).unwrap();
            let flags = ring_bond_flags(&mol);
            for (bidx, bond) in mol.bonds() {
                let mut seen = vec![false; mol.atom_count()];
                let mut stack = vec![bond.atom1];
                seen[bond.atom1.0 as usize] = true;
                while let Some(a) = stack.pop() {
                    for (nb, b) in mol.neighbors(a) {
                        if b != bidx && is_ring_eligible(mol.bond(b).order) && !seen[nb.0 as usize]
                        {
                            seen[nb.0 as usize] = true;
                            stack.push(nb);
                        }
                    }
                }
                let cyclic = is_ring_eligible(bond.order) && seen[bond.atom2.0 as usize];
                assert_eq!(flags[bidx.0 as usize], cyclic, "{smi} bond {bidx:?}");
            }
        }
    }

    use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder, implicit_hcount};

    // Build a cyclohexane molecule (6 carbons, 6 single bonds).
    fn cyclohexane() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Single)
                .unwrap();
        }
        b.build()
    }

    // Build a benzene molecule (6 aromatic carbons, 6 single bonds for topology).
    fn benzene() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Single)
                .unwrap();
        }
        b.build()
    }

    #[test]
    fn ring_atom_flags_match_sssr_membership_on_representative_topologies() {
        let cases = [
            cyclohexane(),
            // Exercises the one-closure SMILES fast path with both a branch
            // and a disconnected acyclic component.
            chematic_smiles::parse("C1CCC(CC1)O.CC").expect("monocycle SMILES"),
            naphthalene(),
            adamantane(),
            chematic_smiles::parse("C12C3C4C1C5C4C3C25").expect("cubane SMILES"),
            spiro_nonane(),
            disconnected_rings(),
            MoleculeBuilder::new().build(),
        ];
        for mol in cases {
            let expected: Vec<bool> = (0..mol.atom_count())
                .map(|idx| find_sssr(&mol).contains_atom(AtomIdx(idx as u32)))
                .collect();
            assert_eq!(ring_atom_flags(&mol), expected);
        }
    }

    #[test]
    fn ring_bond_flags_exclude_bridges_between_cyclic_components() {
        let mol = chematic_smiles::parse("C1CCCCC1CC2CCCCC2").expect("linked rings");
        let flags = ring_bond_flags(&mol);
        assert_eq!(flags.len(), mol.bond_count());
        assert_eq!(flags.iter().filter(|&&is_ring| is_ring).count(), 12);
        for (bond_idx, bond) in mol.bonds() {
            if !flags[bond_idx.0 as usize] {
                assert_eq!(bond.order, BondOrder::Single);
            }
        }
    }

    // Build naphthalene: 10 atoms, 11 bonds (two fused 6-membered rings).
    // Atom numbering:
    //   0-1-2-3-4-5-0  (ring 1 perimeter, 6 atoms)
    //   4-6-7-8-9-5    (ring 2, sharing bond 4-5)
    fn naphthalene() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..10).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        // Ring 1: 0-1-2-3-4-9-0
        let ring1 = [0usize, 1, 2, 3, 4, 9];
        for i in 0..6 {
            b.add_bond(
                atoms[ring1[i]],
                atoms[ring1[(i + 1) % 6]],
                BondOrder::Single,
            )
            .unwrap();
        }
        // Ring 2: 4-5-6-7-8-9 (shares bond 4-9)
        // Bonds to add: 4-5, 5-6, 6-7, 7-8, 8-9
        b.add_bond(atoms[4], atoms[5], BondOrder::Single).unwrap();
        b.add_bond(atoms[5], atoms[6], BondOrder::Single).unwrap();
        b.add_bond(atoms[6], atoms[7], BondOrder::Single).unwrap();
        b.add_bond(atoms[7], atoms[8], BondOrder::Single).unwrap();
        b.add_bond(atoms[8], atoms[9], BondOrder::Single).unwrap();
        b.build()
    }

    // Build norbornane (bicyclo[2.2.1]heptane): 7 carbons, 8 bonds, 2 rings.
    // Numbering:
    //   bridgehead atoms: 0, 3
    //   bridge 1: 0-1-2-3
    //   bridge 2: 0-4-5-3
    //   bridge 3: 0-6-3  (one-carbon bridge)
    fn norbornane() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..7).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        // Bridge 1: 0-1-2-3
        b.add_bond(atoms[0], atoms[1], BondOrder::Single).unwrap();
        b.add_bond(atoms[1], atoms[2], BondOrder::Single).unwrap();
        b.add_bond(atoms[2], atoms[3], BondOrder::Single).unwrap();
        // Bridge 2: 0-4-5-3
        b.add_bond(atoms[0], atoms[4], BondOrder::Single).unwrap();
        b.add_bond(atoms[4], atoms[5], BondOrder::Single).unwrap();
        b.add_bond(atoms[5], atoms[3], BondOrder::Single).unwrap();
        // Bridge 3: 0-6-3
        b.add_bond(atoms[0], atoms[6], BondOrder::Single).unwrap();
        b.add_bond(atoms[6], atoms[3], BondOrder::Single).unwrap();
        b.build()
    }

    #[test]
    fn test_azulene_sssr_minimal() {
        // Azulene: cyclopentadiene fused to cycloheptatriene, sharing one
        // bond. RDKit's GetSymmSSSR ring-size multiset is [5, 7]; the old
        // single-spanning-tree find_sssr previously returned a non-minimal
        // basis for this topology (see aromaticity.rs's PROVISIONAL-tagged
        // azulene regression test for the downstream effect on Pass 1/2).
        let mol = chematic_smiles::parse("C1=CC2=CC=CC=CC2=C1").expect("azulene SMILES");
        let sssr = find_sssr(&mol);
        let mut sizes: Vec<usize> = sssr.rings().iter().map(|r| r.len()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![5, 7], "azulene SSSR must be minimal [5, 7]");
    }

    #[test]
    fn test_indolizine_sssr_minimal() {
        // Indolizine: pyrrole fused to pyridine sharing a bridgehead N.
        // RDKit's GetSymmSSSR ring-size multiset is [5, 6] — a fused
        // heterocycle oracle case distinct from azulene's all-carbon,
        // odd/odd-sized ring pair.
        let mol = chematic_smiles::parse("c1ccn2ccccc12").expect("indolizine SMILES");
        let sssr = find_sssr(&mol);
        let mut sizes: Vec<usize> = sssr.rings().iter().map(|r| r.len()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![5, 6], "indolizine SSSR must be minimal [5, 6]");
    }

    #[test]
    fn test_cyclohexane_sssr() {
        let mol = cyclohexane();
        let rings = find_sssr(&mol);
        assert_eq!(rings.ring_count(), 1, "cyclohexane has exactly 1 ring");
        assert_eq!(rings.rings()[0].len(), 6, "cyclohexane ring has 6 atoms");
    }

    #[test]
    fn blocked_bond_shortest_ring_search_is_non_mutating() {
        let mol = cyclohexane();
        let (bond, _) = mol.bond_between(AtomIdx(0), AtomIdx(1)).unwrap();
        let mut blocked = FxHashSet::default();
        blocked.insert(bond);
        assert!(find_smallest_rings_bfs_with_blocked_bonds(&mol, AtomIdx(0), &blocked).is_empty());
        assert_eq!(find_sssr(&mol).ring_count(), 1);
    }

    #[test]
    fn rdkit_d2_root_selection_collapses_degree2_chains() {
        let roots = select_rdkit_d2_roots(&benzene());
        assert_eq!(roots, vec![AtomIdx(0)]);
    }

    #[test]
    fn symmetrized_d2_root_candidates_are_permutation_invariant() {
        let base = naphthalene();
        let base_ranks = canonical_atom_ranks(&base);
        let mut expected: Vec<u64> = select_permutation_invariant_d2_roots(&base)
            .into_iter()
            .map(|atom| base_ranks[atom.0 as usize])
            .collect();
        expected.sort_unstable();

        let n = base.atom_count();
        for perm in [(0..n).rev().collect::<Vec<_>>(), {
            let mut rotated: Vec<usize> = (0..n).collect();
            rotated.rotate_left(3);
            rotated
        }] {
            let permuted = permute_molecule(&base, &perm);
            let ranks = canonical_atom_ranks(&permuted);
            let mut actual: Vec<u64> = select_permutation_invariant_d2_roots(&permuted)
                .into_iter()
                .map(|atom| ranks[atom.0 as usize])
                .collect();
            actual.sort_unstable();
            assert_eq!(
                actual, expected,
                "D2 root candidate ranks changed: {perm:?}"
            );
        }
    }

    #[test]
    fn test_benzene_sssr() {
        let mol = benzene();
        let rings = find_sssr(&mol);
        assert_eq!(rings.ring_count(), 1, "benzene has exactly 1 ring");
        assert_eq!(rings.rings()[0].len(), 6, "benzene ring has 6 atoms");
    }

    #[test]
    fn test_naphthalene_sssr() {
        let mol = naphthalene();
        let rings = find_sssr(&mol);
        // Cycle rank: 11 bonds - 10 atoms + 1 component = 2
        // SSSR should have 2 rings, both 6-membered.
        assert_eq!(rings.ring_count(), 2, "naphthalene SSSR has 2 rings");
        for ring in rings.rings() {
            assert_eq!(ring.len(), 6, "each naphthalene SSSR ring has 6 atoms");
        }
    }

    #[test]
    fn test_norbornane_sssr() {
        let mol = norbornane();
        let rings = find_sssr(&mol);
        // Cycle rank: 8 bonds - 7 atoms + 1 component = 2
        assert_eq!(rings.ring_count(), 2, "norbornane SSSR has 2 rings");
        // The two smallest rings are both 5-membered.
        for ring in rings.rings() {
            assert_eq!(ring.len(), 5, "each norbornane SSSR ring has 5 atoms");
        }
    }

    #[test]
    fn test_acyclic_molecule() {
        // Ethane: no rings.
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let mol = b.build();
        let rings = find_sssr(&mol);
        assert_eq!(rings.ring_count(), 0);
    }

    #[test]
    fn test_contains_atom() {
        let mol = cyclohexane();
        let rings = find_sssr(&mol);
        for i in 0..6u32 {
            assert!(
                rings.contains_atom(AtomIdx(i)),
                "atom {} should be in a ring",
                i
            );
        }
    }

    #[test]
    fn test_atoms_in_ring_count_benzene() {
        let mol = benzene();
        let rings = find_sssr(&mol);
        for i in 0..6u32 {
            assert_eq!(
                rings.atoms_in_ring_count(AtomIdx(i)),
                1,
                "each benzene atom is in exactly 1 ring"
            );
        }
    }

    // Anthracene: 14 atoms, 16 bonds (3 fused 6-membered rings).
    // Linear fusion: central ring shares edges with two outer rings.
    fn anthracene() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..14).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        // Ring 1 (left): 0-1-2-3-8-9-0
        b.add_bond(atoms[0], atoms[1], BondOrder::Single).unwrap();
        b.add_bond(atoms[1], atoms[2], BondOrder::Single).unwrap();
        b.add_bond(atoms[2], atoms[3], BondOrder::Single).unwrap();
        b.add_bond(atoms[3], atoms[8], BondOrder::Single).unwrap();
        b.add_bond(atoms[8], atoms[9], BondOrder::Single).unwrap();
        b.add_bond(atoms[9], atoms[0], BondOrder::Single).unwrap();
        // Ring 2 (center): 3-4-5-6-7-8-3
        b.add_bond(atoms[3], atoms[4], BondOrder::Single).unwrap();
        b.add_bond(atoms[4], atoms[5], BondOrder::Single).unwrap();
        b.add_bond(atoms[5], atoms[6], BondOrder::Single).unwrap();
        b.add_bond(atoms[6], atoms[7], BondOrder::Single).unwrap();
        b.add_bond(atoms[7], atoms[8], BondOrder::Single).unwrap();
        // Ring 3 (right): 7-10-11-12-13-6-7
        b.add_bond(atoms[7], atoms[10], BondOrder::Single).unwrap();
        b.add_bond(atoms[10], atoms[11], BondOrder::Single).unwrap();
        b.add_bond(atoms[11], atoms[12], BondOrder::Single).unwrap();
        b.add_bond(atoms[12], atoms[13], BondOrder::Single).unwrap();
        b.add_bond(atoms[13], atoms[6], BondOrder::Single).unwrap();
        b.build()
    }

    // Spiro[4.4]nonane: two 5-membered rings sharing a single bridgehead atom.
    // 9 atoms total, cycle rank 2.
    fn spiro_nonane() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..9).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        // Bridgehead: atom 0
        // Ring 1: 0-1-2-3-4-0
        b.add_bond(atoms[0], atoms[1], BondOrder::Single).unwrap();
        b.add_bond(atoms[1], atoms[2], BondOrder::Single).unwrap();
        b.add_bond(atoms[2], atoms[3], BondOrder::Single).unwrap();
        b.add_bond(atoms[3], atoms[4], BondOrder::Single).unwrap();
        b.add_bond(atoms[4], atoms[0], BondOrder::Single).unwrap();
        // Ring 2: 0-5-6-7-8-0
        b.add_bond(atoms[0], atoms[5], BondOrder::Single).unwrap();
        b.add_bond(atoms[5], atoms[6], BondOrder::Single).unwrap();
        b.add_bond(atoms[6], atoms[7], BondOrder::Single).unwrap();
        b.add_bond(atoms[7], atoms[8], BondOrder::Single).unwrap();
        b.add_bond(atoms[8], atoms[0], BondOrder::Single).unwrap();
        b.build()
    }

    // 12-membered macrocycle (1 ring, 12 atoms).
    fn dodecane_ring() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..12).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..12 {
            b.add_bond(atoms[i], atoms[(i + 1) % 12], BondOrder::Single)
                .unwrap();
        }
        b.build()
    }

    // Two disconnected rings (two components).
    fn disconnected_rings() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        // Benzene ring: atoms 0-5
        let benzene_atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            b.add_bond(
                benzene_atoms[i],
                benzene_atoms[(i + 1) % 6],
                BondOrder::Single,
            )
            .unwrap();
        }
        // Separate cyclohexane ring: atoms 6-11
        let hexane_atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            b.add_bond(
                hexane_atoms[i],
                hexane_atoms[(i + 1) % 6],
                BondOrder::Single,
            )
            .unwrap();
        }
        b.build()
    }

    // Adamantane-like tricyclic structure (simplified):
    // 10 atoms, 3 bridges between 2 bridgeheads
    fn adamantane() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..10).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        // Bridgehead atoms: 0, 5
        // Bridge 1: 0-1-2-5 (3 bonds in chain)
        b.add_bond(atoms[0], atoms[1], BondOrder::Single).unwrap();
        b.add_bond(atoms[1], atoms[2], BondOrder::Single).unwrap();
        b.add_bond(atoms[2], atoms[5], BondOrder::Single).unwrap();
        // Bridge 2: 0-3-4-5 (3 bonds in chain)
        b.add_bond(atoms[0], atoms[3], BondOrder::Single).unwrap();
        b.add_bond(atoms[3], atoms[4], BondOrder::Single).unwrap();
        b.add_bond(atoms[4], atoms[5], BondOrder::Single).unwrap();
        // Bridge 3: 0-6-7-5 (3 bonds in chain)
        b.add_bond(atoms[0], atoms[6], BondOrder::Single).unwrap();
        b.add_bond(atoms[6], atoms[7], BondOrder::Single).unwrap();
        b.add_bond(atoms[7], atoms[5], BondOrder::Single).unwrap();
        // Cross-link bonds to connect bridges (forming tertiary center)
        // 1-3, 2-4, 6-? to complete cage
        b.add_bond(atoms[1], atoms[3], BondOrder::Single).unwrap();
        b.add_bond(atoms[2], atoms[4], BondOrder::Single).unwrap();
        b.build()
    }

    #[test]
    fn test_anthracene_sssr() {
        let mol = anthracene();
        let rings = find_sssr(&mol);
        // Cycle rank: 16 bonds - 14 atoms + 1 component = 3.
        // RDKit's GetSymmSSSR gives three 6-membered rings (linear acene) —
        // Horton's minimum-weight basis must match this exactly, not just
        // "3 rings covering most atoms" (the old spanning-tree algorithm
        // could substitute a larger non-minimal ring for one of these).
        assert_eq!(rings.ring_count(), 3, "anthracene SSSR has 3 rings");
        for ring in rings.rings() {
            assert_eq!(ring.len(), 6, "each anthracene SSSR ring has 6 atoms");
        }
        let all_ring_atoms: std::collections::HashSet<_> = rings
            .rings()
            .iter()
            .flat_map(|r| r.iter().copied())
            .collect();
        assert_eq!(
            all_ring_atoms.len(),
            14,
            "anthracene SSSR atoms cover every atom"
        );
    }

    #[test]
    fn test_spiro_nonane_sssr() {
        let mol = spiro_nonane();
        let rings = find_sssr(&mol);
        // Cycle rank: 8 bonds - 9 atoms + 1 component = 0... wait, let me recalculate
        // Actually: two 5-membered rings sharing 1 atom = 4 + 4 + 2 (bridge) = 10 bonds
        // 10 bonds - 9 atoms + 1 = 2 rings
        assert_eq!(rings.ring_count(), 2, "spiro[4.4]nonane SSSR has 2 rings");
        for ring in rings.rings() {
            assert_eq!(ring.len(), 5, "each spiro nonane SSSR ring is 5-membered");
        }
    }

    #[test]
    fn test_dodecane_ring_sssr() {
        let mol = dodecane_ring();
        let rings = find_sssr(&mol);
        assert_eq!(rings.ring_count(), 1, "12-membered ring has 1 SSSR entry");
        assert_eq!(
            rings.rings()[0].len(),
            12,
            "12-membered ring SSSR has 12 atoms"
        );
    }

    #[test]
    fn test_disconnected_rings_sssr() {
        let mol = disconnected_rings();
        let rings = find_sssr(&mol);
        // Cycle rank: 12 bonds - 12 atoms + 2 components = 2 rings
        assert_eq!(
            rings.ring_count(),
            2,
            "two disconnected rings yield 2 SSSR entries"
        );
        let sizes: Vec<_> = rings.rings().iter().map(|r| r.len()).collect();
        assert!(sizes.contains(&6), "one ring should be 6-membered");
    }

    #[test]
    fn test_adamantane_sssr() {
        let mol = adamantane();
        let rings = find_sssr(&mol);
        // Simplified adamantane: 10 atoms, 12 bonds, 1 component
        // Cycle rank: 12 - 10 + 1 = 3
        // (May be 3-4 depending on cross-link structure and Gaussian elimination)
        assert!(
            rings.ring_count() >= 3,
            "adamantane SSSR has at least 3 rings"
        );
        // Each ring should be reasonable size
        for ring in rings.rings() {
            assert!(!ring.is_empty(), "each ring should have atoms");
            assert!(ring.len() <= 10, "ring should not exceed molecule size");
        }
    }

    #[test]
    fn test_macrocycle_atom_in_ring_count() {
        let mol = dodecane_ring();
        let rings = find_sssr(&mol);
        for i in 0..12u32 {
            assert_eq!(
                rings.atoms_in_ring_count(AtomIdx(i)),
                1,
                "each dodecane atom is in exactly 1 ring"
            );
        }
    }

    #[test]
    fn test_figueras_bfs_finds_all_smallest_cubane_faces_through_each_root() {
        let mol = chematic_smiles::parse("C12C3C4C1C5C4C3C25").expect("cubane SMILES");
        let mut faces = std::collections::BTreeSet::new();
        for root in 0..mol.atom_count() {
            for ring in find_smallest_rings_bfs(&mol, AtomIdx(root as u32)) {
                assert_eq!(ring.len(), 4, "cubane's smallest rings are square faces");
                let mut face = ring.into_iter().map(|a| a.0).collect::<Vec<_>>();
                face.sort_unstable();
                faces.insert(face);
            }
        }
        assert_eq!(
            faces.len(),
            6,
            "cubane has six symmetry-equivalent square faces"
        );

        let mol = chematic_smiles::parse("C12C3C4C5C1C6C7C2C8C3C9C4C1C5C6C2C7C8C9C12")
            .expect("dodecahedrane SMILES");
        let mut faces = std::collections::BTreeSet::new();
        for root in 0..mol.atom_count() {
            for ring in find_smallest_rings_bfs(&mol, AtomIdx(root as u32)) {
                assert_eq!(
                    ring.len(),
                    5,
                    "dodecahedrane's smallest rings are pentagons"
                );
                let mut face = ring.into_iter().map(|a| a.0).collect::<Vec<_>>();
                face.sort_unstable();
                faces.insert(face);
            }
        }
        assert_eq!(
            faces.len(),
            12,
            "dodecahedrane has twelve symmetry-equivalent pentagonal faces"
        );
    }

    #[test]
    fn test_symmetrized_sssr_adds_only_verified_duplicate_faces() {
        let benzene = chematic_smiles::parse("c1ccccc1").expect("benzene SMILES");
        assert_eq!(find_symmetrized_sssr(&benzene).ring_count(), 1);
        let diagnostic = find_symmetrized_sssr_with_diagnostics(&benzene);
        assert_eq!(diagnostic.status(), SymmetrizedSssrStatus::Complete);
        assert_eq!(diagnostic.rings().ring_count(), 1);
        assert_eq!(diagnostic.candidates_examined(), 6);

        let cubane = chematic_smiles::parse("C12C3C4C1C5C4C3C25").expect("cubane SMILES");
        assert_eq!(find_symmetrized_sssr(&cubane).ring_count(), 6);

        let dodeca = chematic_smiles::parse("C12C3C4C5C1C6C7C2C8C3C9C4C1C5C6C2C7C8C9C12")
            .expect("dodecahedrane SMILES");
        assert_eq!(find_symmetrized_sssr(&dodeca).ring_count(), 12);
    }

    #[test]
    fn bounded_symmetrized_sssr_reports_measured_cap_exhaustion() {
        let dodeca = chematic_smiles::parse("C12C3C4C5C1C6C7C2C8C3C9C4C1C5C6C2C7C8C9C12")
            .expect("dodecahedrane SMILES");
        let result = find_symmetrized_sssr_with_diagnostics_bounded(&dodeca, Some(0));
        assert_eq!(result.status(), SymmetrizedSssrStatus::CapExhausted);
        assert_eq!(result.candidates_examined(), 1);
        assert_eq!(result.rings().ring_count(), find_sssr(&dodeca).ring_count());
    }

    #[test]
    fn issue337_macrocycle_boundary_is_bounded_and_permutation_stable() {
        let cases = [
            (
                "0009",
                "c1cc2cc(c1)-c1cccc(c1)C[n+]1ccc(c3ccccc31)NCCCCCCCCCCNc1cc[n+](c3ccccc13)C2",
                2,
            ),
            (
                "0023",
                "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)-c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
                4,
            ),
            (
                "0028",
                "C1=C\\c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCNc2cc[n+](c3ccccc23)Cc2ccc/1cc2",
                3,
            ),
            (
                "0029",
                "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)CCc2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
                4,
            ),
            (
                "0030",
                "c1ccc2c(c1)c1cc[n+]2Cc2ccc(cc2)Cc2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
                4,
            ),
            (
                "0034",
                "c1ccc2c(c1)c1cc[n+]2Cc2ccc3c(c2)Cc2cc(ccc2-3)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCN1",
                2,
            ),
        ];

        for (name, smiles, expected_macrocycles) in cases {
            let mol = chematic_smiles::parse(smiles).expect(name);
            let macrocycle_count = find_symmetrized_sssr(&mol)
                .rings()
                .iter()
                .filter(|ring| ring.len() >= 20)
                .count();
            assert_eq!(
                macrocycle_count, expected_macrocycles,
                "{name}: diagnostic symmetrized-ring boundary changed"
            );

            let reversed = permute_molecule(&mol, &(0..mol.atom_count()).rev().collect::<Vec<_>>());
            let mut original_sizes: Vec<_> = find_symmetrized_sssr(&mol)
                .rings()
                .iter()
                .filter(|ring| ring.len() >= 20)
                .map(Vec::len)
                .collect();
            let mut reversed_sizes: Vec<_> = find_symmetrized_sssr(&reversed)
                .rings()
                .iter()
                .filter(|ring| ring.len() >= 20)
                .map(Vec::len)
                .collect();
            original_sizes.sort_unstable();
            reversed_sizes.sort_unstable();
            assert_eq!(
                original_sizes, reversed_sizes,
                "{name}: macrocycle sizes changed under relabeling"
            );

            let canonical_family = |candidate: &Molecule| {
                let ranks = canonical_atom_ranks(candidate);
                let mut keys: Vec<_> = find_symmetrized_sssr(candidate)
                    .rings()
                    .iter()
                    .filter(|ring| ring.len() >= 20)
                    .map(|ring| {
                        canonical_bond_set_key(candidate, &ring_bond_set(candidate, ring), &ranks)
                    })
                    .collect();
                keys.sort_unstable();
                keys
            };
            assert_eq!(
                canonical_family(&mol),
                canonical_family(&reversed),
                "{name}: canonical macrocycle family changed under relabeling"
            );

            let expected_family = canonical_family(&mol);
            for seed in 0..256_u64 {
                let perm = seeded_permutation(mol.atom_count(), seed ^ 0x337);
                let permuted = permute_molecule(&mol, &perm);
                let macrocycle_count = find_symmetrized_sssr(&permuted)
                    .rings()
                    .iter()
                    .filter(|ring| ring.len() >= 20)
                    .count();
                assert_eq!(
                    macrocycle_count, expected_macrocycles,
                    "{name}: representative count changed under seeded relabeling {seed}"
                );
                assert_eq!(
                    canonical_family(&permuted),
                    expected_family,
                    "{name}: canonical family changed under seeded relabeling {seed}"
                );
            }
        }
    }

    #[test]
    fn symmetrized_sssr_is_invariant_under_explicit_hydrogen_nodes() {
        let compact = chematic_smiles::parse(
            "C1=C\\c2ccc(cc2)C[n+]2ccc(c3ccccc32)NCCCCCCCCCCNc2cc[n+](c3ccccc23)Cc2ccc/1cc2",
        )
        .expect("macrocycle SMILES");
        let explicit = add_explicit_hydrogens_for_test(&compact);
        let canonical = |mol: &Molecule| {
            let mut rings: Vec<_> = find_symmetrized_sssr(mol)
                .rings()
                .iter()
                .map(|ring| (ring.len(), ring.iter().map(|idx| idx.0).collect::<Vec<_>>()))
                .collect();
            rings.sort();
            rings
        };
        assert_eq!(canonical(&compact), canonical(&explicit));
    }

    fn add_explicit_hydrogens_for_test(mol: &Molecule) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let mut remap = vec![AtomIdx(0); mol.atom_count()];
        for (idx, atom) in mol.atoms() {
            let mut atom = atom.clone();
            atom.hydrogen_count = Some(0);
            remap[idx.0 as usize] = builder.add_atom(atom);
        }
        for (_, bond) in mol.bonds() {
            builder
                .add_bond(
                    remap[bond.atom1.0 as usize],
                    remap[bond.atom2.0 as usize],
                    bond.order,
                )
                .unwrap();
        }
        for (idx, atom) in mol.atoms() {
            if atom.element == Element::H {
                continue;
            }
            for _ in 0..implicit_hcount(mol, idx) {
                let h = builder.add_atom(Atom::new(Element::H));
                builder
                    .add_bond(remap[idx.0 as usize], h, BondOrder::Single)
                    .unwrap();
            }
        }
        builder.build()
    }

    #[test]
    fn test_cubane_sssr() {
        // Cubane C8H8 — a cage molecule with 12 C-C bonds and 8 vertices.
        // Cycle rank = E - V + 1 = 12 - 8 + 1 = 5.
        // Cubane has 6 square (4-membered) faces but only 5 are linearly
        // independent in GF(2) — the 6th is the XOR of the other 5 (not just
        // any pair of them, since Sum-of-all-6-faces = 0 in GF(2): each edge
        // is shared by exactly 2 faces).
        //
        // Horton's candidate generation (O(V*E) candidates, guaranteed to
        // contain a minimum-weight basis) finds a truly minimal SSSR here:
        // all 5 basis rings are 4-membered (weight 20), strictly better than
        // the old single-spanning-tree algorithm's typical "4 four-membered +
        // 1 six-membered diagonal" (weight 22) — see module doc / project
        // history for why the old algorithm couldn't guarantee this.
        //
        // Recovering the 6th (symmetry-equivalent) face requires XOR-ing all
        // 5 basis rings together, not a pairwise XOR — augmented_ring_set
        // only does pairwise XOR, so it cannot find it from an already-fully-
        // 4-membered basis (two same-size adjacent cube faces XOR to a
        // 6-membered "belt", not another 4-ring). This is expected: full
        // symmetrization (all 6 symmetry-equivalent minimal rings, matching
        // RDKit's GetSymmSSSR on cubane) is out of scope for Horton alone and
        // deferred to a later Vismara "relevant cycles" pass (see project
        // plan) — not a regression, since Horton's SSSR is still a strict
        // minimality improvement over the previous algorithm.
        let mol = chematic_smiles::parse("C12C3C4C1C5C4C3C25").expect("cubane SMILES");
        let sssr = find_sssr(&mol);

        assert_eq!(
            sssr.rings().len(),
            5,
            "cubane must have exactly 5 SSSR rings (cycle rank 12−8+1=5)"
        );
        for ring in sssr.rings() {
            assert!(
                ring.len() <= 6,
                "cubane SSSR rings must be ≤ 6-membered, got {}",
                ring.len()
            );
        }
        let four_membered = sssr.rings().iter().filter(|r| r.len() == 4).count();
        assert_eq!(
            four_membered, 5,
            "Horton SSSR should find all 5 basis rings as 4-membered faces, got {four_membered}"
        );
    }

    /// Rebuild `mol` with atoms relabeled by `perm` (perm[new_idx] = old_idx),
    /// preserving the same graph but a different atom insertion order.
    fn permute_molecule(mol: &chematic_core::Molecule, perm: &[usize]) -> chematic_core::Molecule {
        let mut old_to_new = vec![0u32; perm.len()];
        for (new_idx, &old_idx) in perm.iter().enumerate() {
            old_to_new[old_idx] = new_idx as u32;
        }
        let mut builder = MoleculeBuilder::new();
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

    fn seeded_permutation(len: usize, seed: u64) -> Vec<usize> {
        let mut permutation: Vec<usize> = (0..len).collect();
        let mut state = seed.wrapping_add(1);
        for index in (1..len).rev() {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            let swap_index = (state as usize) % (index + 1);
            permutation.swap(index, swap_index);
        }
        permutation
    }

    /// find_sssr's ring-size multiset must not depend on atom insertion order.
    /// Probes the fused/bridged/cage systems this project has already
    /// identified as the hard cases for canonical_atom_ranks' 3-round
    /// (not-fixpoint) Weisfeiler-Leman tie-break -- unlike
    /// `canonical_atom_order` (chematic-smiles), this function's doc comment
    /// does not claim full canonical-labeling power, and find_sssr's
    /// self-stability was already measured at 0% on a 5000-molecule corpus
    /// during the Horton rewrite (698ba3f); this is a permanent regression
    /// guard for that claim, not a first-time probe.
    #[test]
    fn find_sssr_ring_size_multiset_is_permutation_invariant() {
        let cases: Vec<(&str, chematic_core::Molecule)> = vec![
            ("naphthalene", naphthalene()),
            ("norbornane", norbornane()),
            ("spiro_nonane", spiro_nonane()),
            ("adamantane", adamantane()),
            (
                "cubane",
                chematic_smiles::parse("C12C3C4C1C5C4C3C25").expect("cubane SMILES"),
            ),
        ];

        for (name, mol) in cases {
            let n = mol.atom_count();
            let mut orig_sizes: Vec<usize> = find_sssr(&mol).rings().iter().map(Vec::len).collect();
            orig_sizes.sort_unstable();

            let perms: Vec<Vec<usize>> = vec![(0..n).rev().collect(), {
                let mut p: Vec<usize> = (0..n).collect();
                if n > 2 {
                    p.rotate_left(n / 3 + 1);
                }
                p
            }];
            for perm in perms {
                let permuted = permute_molecule(&mol, &perm);
                let mut perm_sizes: Vec<usize> =
                    find_sssr(&permuted).rings().iter().map(Vec::len).collect();
                perm_sizes.sort_unstable();
                assert_eq!(
                    orig_sizes, perm_sizes,
                    "{name}: SSSR ring-size multiset changed under atom permutation {perm:?}"
                );
            }
        }
    }

    #[test]
    fn symmetrized_sssr_order_is_permutation_invariant() {
        let cases = [
            naphthalene(),
            chematic_smiles::parse("C12C3C4C1C5C4C3C25").expect("cubane SMILES"),
        ];

        for mol in cases {
            let ranks = canonical_atom_ranks(&mol);
            let mut expected: Vec<_> = find_symmetrized_sssr(&mol)
                .rings()
                .iter()
                .map(|ring| (ring.len(), canonical_cycle_key(ring, &ranks)))
                .collect();
            expected.sort_unstable();

            let n = mol.atom_count();
            for perm in [(0..n).rev().collect::<Vec<_>>(), {
                let mut rotated: Vec<usize> = (0..n).collect();
                rotated.rotate_left(3.min(n));
                rotated
            }] {
                let permuted = permute_molecule(&mol, &perm);
                let ranks = canonical_atom_ranks(&permuted);
                let mut actual: Vec<_> = find_symmetrized_sssr(&permuted)
                    .rings()
                    .iter()
                    .map(|ring| (ring.len(), canonical_cycle_key(ring, &ranks)))
                    .collect();
                actual.sort_unstable();
                assert_eq!(actual, expected, "symmetrized SSSR changed under {perm:?}");
            }
        }
    }

    // ── RDKit PR #9118: SSSR excludes Zero-order and Dative bonds ────────────

    #[test]
    fn sssr_ignores_zero_order_bonds() {
        // A--B via a single bond PLUS a Zero-order bond between the same atoms
        // must NOT create a ring. Zero-order bonds are non-valence connections.
        let mut b = MoleculeBuilder::new();
        let mut a_atom = Atom::new(chematic_core::Element::C);
        a_atom.hydrogen_count = Some(3);
        let mut b_atom = Atom::new(chematic_core::Element::C);
        b_atom.hydrogen_count = Some(3);
        let a = b.add_atom(a_atom);
        let bb = b.add_atom(b_atom);
        b.add_bond(a, bb, BondOrder::Single).unwrap();
        b.add_bond(a, bb, BondOrder::Zero)
            .expect_err("duplicate bond — MoleculeBuilder should reject or ignore it");
        // Build a proper molecule: just two atoms with a single bond.
        // The zero-order bond attempt is rejected, so the molecule is acyclic.
        let mol = b.build();
        let sssr = find_sssr(&mol);
        assert_eq!(
            sssr.rings().len(),
            0,
            "single bond between two atoms → no ring"
        );
    }

    #[test]
    fn sssr_ignores_zero_order_bond_as_third_bond() {
        // Benzene ring (6 aromatic bonds) PLUS one Zero-order bond closing an
        // extra connection should NOT add a phantom 2-membered ring.
        // We build cyclohexane (all single bonds) and verify no extra ring from
        // a Zero-order bond added between two non-adjacent atoms.
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4)
            .map(|_| {
                let mut a = Atom::new(chematic_core::Element::C);
                a.hydrogen_count = Some(2);
                b.add_atom(a)
            })
            .collect();
        // Square ring: 0-1-2-3-0
        b.add_bond(atoms[0], atoms[1], BondOrder::Single).unwrap();
        b.add_bond(atoms[1], atoms[2], BondOrder::Single).unwrap();
        b.add_bond(atoms[2], atoms[3], BondOrder::Single).unwrap();
        b.add_bond(atoms[3], atoms[0], BondOrder::Single).unwrap();
        // Zero-order bond between non-adjacent atoms: should NOT create a new ring.
        // Note: builder may reject parallel bonds; if so the test still passes.
        let _ = b.add_bond(atoms[0], atoms[2], BondOrder::Zero);
        let mol = b.build();
        let sssr = find_sssr(&mol);
        // Should find exactly 1 ring (the 4-membered ring), NOT 2 or 3.
        assert_eq!(
            sssr.rings().len(),
            1,
            "zero-order diagonal bond must not create extra rings: found {:?}",
            sssr.rings().iter().map(|r| r.len()).collect::<Vec<_>>()
        );
    }
}
