#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Atom;
    use crate::element::Element;
    use crate::molecule::{MoleculeBuilder, STEREO_H_SENTINEL};
    use crate::stereo_group::StereoGroup;

    /// Build benzene as a fully aromatic SMILES-style molecule.
    fn benzene() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// Build pyridine (5 aromatic C + 1 aromatic N, ring).
    fn pyridine() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::aromatic(Element::C));
        let c2 = b.add_atom(Atom::aromatic(Element::C));
        let c3 = b.add_atom(Atom::aromatic(Element::C));
        let n = b.add_atom(Atom::aromatic(Element::N));
        let c4 = b.add_atom(Atom::aromatic(Element::C));
        let c5 = b.add_atom(Atom::aromatic(Element::C));
        let atoms = [c1, c2, c3, n, c4, c5];
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// Build furan (4 aromatic C + 1 aromatic O, ring).
    fn furan() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::aromatic(Element::O));
        let c1 = b.add_atom(Atom::aromatic(Element::C));
        let c2 = b.add_atom(Atom::aromatic(Element::C));
        let c3 = b.add_atom(Atom::aromatic(Element::C));
        let c4 = b.add_atom(Atom::aromatic(Element::C));
        let atoms = [o, c1, c2, c3, c4];
        for i in 0..5 {
            b.add_bond(atoms[i], atoms[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// Build pyrrole ([nH] + 4 aromatic C, ring).
    fn pyrrole() -> Molecule {
        let mut b = MoleculeBuilder::new();
        // N with 1 H — bracket-style (hydrogen_count = Some(1))
        let mut n_atom = Atom::aromatic(Element::N);
        n_atom.hydrogen_count = Some(1);
        let n = b.add_atom(n_atom);
        let c1 = b.add_atom(Atom::aromatic(Element::C));
        let c2 = b.add_atom(Atom::aromatic(Element::C));
        let c3 = b.add_atom(Atom::aromatic(Element::C));
        let c4 = b.add_atom(Atom::aromatic(Element::C));
        let atoms = [n, c1, c2, c3, c4];
        for i in 0..5 {
            b.add_bond(atoms[i], atoms[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    #[test]
    fn test_kekulize_benzene() {
        let mol = benzene();
        let result = kekulize(&mol).expect("benzene kekulization failed");
        assert_eq!(result.len(), 6); // all 6 aromatic bonds assigned

        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        let singles = result.values().filter(|&&o| o == BondOrder::Single).count();
        assert_eq!(doubles, 3, "benzene must have 3 double bonds");
        assert_eq!(singles, 3, "benzene must have 3 single bonds");
    }

    #[test]
    fn test_kekulize_pyridine() {
        let mol = pyridine();
        let result = kekulize(&mol).expect("pyridine kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 3, "pyridine must have 3 double bonds");
    }

    #[test]
    fn test_kekulize_furan() {
        let mol = furan();
        let result = kekulize(&mol).expect("furan kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 2, "furan must have 2 double bonds");
    }

    #[test]
    fn test_kekulize_pyrrole() {
        let mol = pyrrole();
        let result = kekulize(&mol).expect("pyrrole kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 2, "pyrrole must have 2 double bonds");
    }

    #[test]
    fn test_kekulize_naphthalene() {
        // 10 aromatic C, 11 aromatic bonds (fused bicyclic)
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..10)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        // Ring 1: 0-1-2-3-4-9
        let ring1 = [0, 1, 2, 3, 4, 9];
        for i in 0..ring1.len() {
            b.add_bond(
                atoms[ring1[i]],
                atoms[ring1[(i + 1) % ring1.len()]],
                BondOrder::Aromatic,
            )
            .unwrap();
        }
        // Ring 2: 4-5-6-7-8-9 (shares bond 4-9)
        let ring2 = [4, 5, 6, 7, 8, 9];
        for i in 0..ring2.len() {
            let a = atoms[ring2[i]];
            let bb = atoms[ring2[(i + 1) % ring2.len()]];
            // Skip already-added bond (4-9)
            if mol_has_no_bond_yet(&b, a, bb) {
                b.add_bond(a, bb, BondOrder::Aromatic).unwrap();
            }
        }
        let mol = b.build();
        let result = kekulize(&mol).expect("naphthalene kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 5, "naphthalene must have 5 double bonds");
    }

    #[test]
    fn test_apply_kekule() {
        let mol = benzene();
        let kekule = kekulize(&mol).unwrap();
        let kekule_mol = apply_kekule(&mol, &kekule);

        // After applying, no aromatic bonds should remain.
        for (_, bond) in kekule_mol.bonds() {
            assert_ne!(
                bond.order,
                BondOrder::Aromatic,
                "apply_kekule should remove all aromatic bonds"
            );
        }
    }

    #[test]
    fn test_no_aromatic_bonds_noop() {
        // A molecule with no aromatic bonds should return empty result.
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let mol = b.build();
        let result = kekulize(&mol).unwrap();
        assert!(result.is_empty());
    }

    // Helper: check that the builder does not yet have a bond between a and b.
    fn mol_has_no_bond_yet(b: &MoleculeBuilder, a: AtomIdx, bb: AtomIdx) -> bool {
        for (_, partner) in b.atom_neighbors(a) {
            if partner == bb {
                return false;
            }
        }
        true
    }

    // B6 coverage: exotic ring systems that require correct matching on non-bipartite graphs.

    /// Azulene: fused 5+7 bicyclic, 10 aromatic C, 11 bonds.
    /// Valid Kekulé form exists → must yield 5 double bonds.
    #[test]
    fn test_kekulize_azulene() {
        let mut b = MoleculeBuilder::new();
        let a: Vec<_> = (0..10)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        // 5-ring: 0-1-2-3-4-0
        for i in 0..5 {
            b.add_bond(a[i], a[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        // 7-ring extras (shares bond 0-4): 4-5-6-7-8-9-0
        for (x, y) in [(4usize, 5usize), (5, 6), (6, 7), (7, 8), (8, 9), (9, 0)] {
            if mol_has_no_bond_yet(&b, a[x], a[y]) {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
        }
        let mol = b.build();
        let result = kekulize(&mol).expect("azulene kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 5, "azulene needs 5 double bonds");
    }

    /// Acenaphthylene: naphthalene + fused cyclopentadiene bridge, 12 aromatic C.
    #[test]
    fn test_kekulize_acenaphthylene() {
        // Connectivity: naphthalene core (atoms 0-9) + bridge atoms 10,11
        // Naphthalene: ring1=0-1-2-3-4-9, ring2=4-5-6-7-8-9
        // Bridge: 0-11-10-1 (the 5-ring across the 1,8 positions of naphthalene)
        let mut b = MoleculeBuilder::new();
        let a: Vec<_> = (0..12)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        // Naphthalene ring 1: 0-1-2-3-4-9-0
        for (x, y) in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 9), (9, 0)] {
            b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
        }
        // Naphthalene ring 2: 4-5-6-7-8-9 (4-9 shared)
        for (x, y) in [(4, 5), (5, 6), (6, 7), (7, 8), (8, 9)] {
            b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
        }
        // Bridge: 0-11-10-1
        for (x, y) in [(0, 11), (11, 10), (10, 1)] {
            b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
        }
        let mol = b.build();
        let result = kekulize(&mol).expect("acenaphthylene kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 6, "acenaphthylene needs 6 double bonds");
    }

    // ---- Edge-case tests added v0.1.100 ----
    //
    // Build complex PAH manually to test edge cases in the matching algorithm.

    /// Biphenylene: two benzene rings connected by a cyclobutadiene bridge.
    ///
    /// Topology: Ring A (0–5), Ring B (6–11), bridge bonds (0–11) and (5–6).
    /// The 4-membered ring 0-5-6-11-0 is the challenging part.
    fn biphenylene() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let a: Vec<_> = (0..12)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        // Ring A (6): 0-1-2-3-4-5-0
        for i in 0..6 {
            b.add_bond(a[i], a[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        // Ring B (6): 6-7-8-9-10-11-6
        for i in 0..6 {
            b.add_bond(a[6 + i], a[6 + (i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        // 4-membered bridge: closes ring 0-5-6-11-0
        b.add_bond(a[5], a[6], BondOrder::Aromatic).unwrap();
        b.add_bond(a[0], a[11], BondOrder::Aromatic).unwrap();
        b.build()
    }

    /// Naphtho-4-ring: three fused 6-membered rings sharing a common bond sequence.
    /// Linear anthracene topology (0-1-2-3-4-5, 5-6-7-8-9-4, 6-10-11-12-13-7).
    fn anthracene() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let a: Vec<_> = (0..14)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        // Ring A: 0-1-2-3-4-5-0
        for i in 0..6 {
            b.add_bond(a[i], a[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        // Ring B: 5-4-9-8-7-6-5 (shares bond 4-5 with Ring A)
        for (x, y) in [(4, 9), (9, 8), (8, 7), (7, 6), (6, 5)] {
            b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
        }
        // Ring C: 6-7-13-12-11-10-6 (shares bond 6-7 with Ring B)
        for (x, y) in [(7, 13), (13, 12), (12, 11), (11, 10), (10, 6)] {
            b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
        }
        b.build()
    }

    #[test]
    fn test_kekulize_biphenylene() {
        // Biphenylene: 4-membered cyclobutadiene bridge between two benzenes.
        // 12 aromatic C → 6 double bonds.
        let mol = biphenylene();
        let result = kekulize(&mol).expect("biphenylene kekulization should succeed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 6, "biphenylene needs 6 double bonds");
    }

    #[test]
    fn test_kekulize_anthracene() {
        // Anthracene (C14H10): 3 linearly fused 6-membered rings.  7 double bonds.
        let mol = anthracene();
        let result = kekulize(&mol).expect("anthracene kekulization should succeed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(doubles, 7, "anthracene needs 7 double bonds");
    }

    #[test]
    fn test_kekulize_biphenylene_double_bond_count() {
        // Cross-check: all 14 bonds should be assigned single or double.
        let mol = biphenylene();
        let result = kekulize(&mol).expect("biphenylene kekulization");
        let singles = result.values().filter(|&&o| o == BondOrder::Single).count();
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(singles + doubles, 14, "biphenylene has 14 aromatic bonds");
    }

    #[test]
    fn test_kekulize_large_fused_6rings() {
        // 4 fused 6-membered rings (pyrene-like) in a simple linear topology.
        // Simulates a large all-even-ring PAH that the BFS should handle easily.
        let mol = {
            let mut b = MoleculeBuilder::new();
            let a: Vec<_> = (0..16)
                .map(|_| b.add_atom(Atom::aromatic(Element::C)))
                .collect();
            // Ring 0-1-2-3-4-5
            for i in 0..6 {
                b.add_bond(a[i], a[(i + 1) % 6], BondOrder::Aromatic)
                    .unwrap();
            }
            // Ring 5-4-9-8-7-6
            for (x, y) in [(4, 9), (9, 8), (8, 7), (7, 6), (6, 5)] {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
            // Ring 1-2-11-10-13-12
            for (x, y) in [(2, 11), (11, 10), (10, 13), (13, 12), (12, 1)] {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
            // Ring 6-7-15-14-11-2 (closed)
            for (x, y) in [(7, 15), (15, 14), (14, 11)] {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
            b.build()
        };
        let result = kekulize(&mol).expect("4-ring PAH kekulization should succeed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert!(
            doubles >= 6,
            "4-ring PAH needs at least 6 double bonds, got {doubles}"
        );
    }

    #[test]
    fn test_kekulize_deterministic() {
        // kekulize() is deterministic: rebuilding the same molecule gives the same count.
        let mol1 = biphenylene();
        let mol2 = biphenylene();
        let r1 = kekulize(&mol1).expect("pass1");
        let r2 = kekulize(&mol2).expect("pass2");
        assert_eq!(
            r1.values().filter(|&&o| o == BondOrder::Double).count(),
            r2.values().filter(|&&o| o == BondOrder::Double).count(),
            "kekulization must be deterministic"
        );
    }

    /// Fluoranthene: 16 aromatic C in a fused 5+6+6+6 system.
    #[test]
    fn test_kekulize_fluoranthene() {
        // Simplified: 3 fused 6-rings + 1 fused 5-ring sharing atoms
        // Build as corannulene-like: two naphthalene units bridged by a 5-ring
        // Atoms 0-15 (16 aromatic C), total bonds = 19 (for fluoranthene)
        // Use a simplified topology that gives the right ring count
        // 6-ring A: 0-1-2-3-4-5
        // 6-ring B: 0-5-6-7-8-9
        // 6-ring C: 2-3-10-11-12-13
        // 5-ring D: 0-9-14-15-1 (bridge)
        let mut b = MoleculeBuilder::new();
        let a: Vec<_> = (0..16)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        // Ring A (6-ring): 0-1-2-3-4-5-0
        for i in 0..6 {
            b.add_bond(a[i], a[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        // Ring B (6-ring): 0-5-6-7-8-9-0
        for (x, y) in [(5, 6), (6, 7), (7, 8), (8, 9), (9, 0)] {
            if mol_has_no_bond_yet(&b, a[x], a[y]) {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
        }
        // Ring C (6-ring): 1-2-10-11-12-13-1
        for (x, y) in [(2, 10), (10, 11), (11, 12), (12, 13), (13, 1)] {
            if mol_has_no_bond_yet(&b, a[x], a[y]) {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
        }
        // Ring D (5-ring): 9-8-14-15-13-9
        for (x, y) in [(8, 14), (14, 15), (15, 13), (13, 9)] {
            if mol_has_no_bond_yet(&b, a[x], a[y]) {
                b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
            }
        }
        let mol = b.build();
        let result = kekulize(&mol).expect("fluoranthene-like kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(
            doubles, 8,
            "fluoranthene-like structure needs 8 double bonds"
        );
    }

    /// Indolizine (`c1ccn2cccc2c1`) — bridgehead N at C9a (aromatic degree 3).
    /// 9 atoms: N contributes lone pair, 8 C atoms form 4 double bonds.
    /// Requires Pass 3 (bridgehead-N exclusion).
    #[test]
    fn kekulize_indolizine() {
        // Bonds: 6-ring (0-1-2-3-7-8-0) ∪ 5-ring (3-4-5-6-7-3), fused at edge 3-7.
        let mut b = MoleculeBuilder::new();
        let c: Vec<_> = (0..9)
            .map(|i| {
                if i == 3 {
                    b.add_atom(Atom::aromatic(Element::N))
                } else {
                    b.add_atom(Atom::aromatic(Element::C))
                }
            })
            .collect();
        for (x, y) in [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 3),
            (7, 8),
            (8, 0),
        ] {
            b.add_bond(c[x], c[y], BondOrder::Aromatic).unwrap();
        }
        let mol = b.build();
        let result = kekulize(&mol);
        assert!(
            result.is_ok(),
            "indolizine kekulization failed: {:?}",
            result.err()
        );
        let doubles = result
            .unwrap()
            .values()
            .filter(|&&o| o == BondOrder::Double)
            .count();
        assert_eq!(doubles, 4, "indolizine: 4 double bonds (N lone-pair donor)");
    }

    /// Quinolizine (`c1ccn2ccccc2c1`) — bridgehead N in two 6-membered rings.
    /// 10 atoms (even), bipartite graph → passes 1/2 handle it; Pass 3 must not break it.
    #[test]
    fn kekulize_quinolizine() {
        // 6-ring A: 0-1-2-3-8-9-0  6-ring B: 3-4-5-6-7-8-3  fused at edge 3-8.
        let mut b = MoleculeBuilder::new();
        let c: Vec<_> = (0..10)
            .map(|i| {
                if i == 3 {
                    b.add_atom(Atom::aromatic(Element::N))
                } else {
                    b.add_atom(Atom::aromatic(Element::C))
                }
            })
            .collect();
        for (x, y) in [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 3),
            (8, 9),
            (9, 0),
        ] {
            b.add_bond(c[x], c[y], BondOrder::Aromatic).unwrap();
        }
        let mol = b.build();
        let result = kekulize(&mol);
        assert!(
            result.is_ok(),
            "quinolizine kekulization failed: {:?}",
            result.err()
        );
        let doubles = result
            .unwrap()
            .values()
            .filter(|&&o| o == BondOrder::Double)
            .count();
        assert_eq!(doubles, 5, "quinolizine: 5 double bonds");
    }

    /// Corannulene (C₂₀H₁₀) — bowl-shaped PAH with five 5-membered rings fused to
    /// five 6-membered rings.  The C-only aromatic subgraph is non-bipartite (five odd
    /// cycles), so passes 1–3 fail; requires Pass 4 (Edmonds' blossom).
    #[test]
    fn kekulize_corannulene() {
        // Inner 5-ring hub: 0-1-2-3-4-0
        // Spokes to outer ring: 0-5, 1-7, 2-9, 3-11, 4-13
        // Outer 10-ring: 5-6-7-8-9-10-11-12-13-14-5
        // Outer 5 "cap" pairs: (5,15),(6,15),(7,16),(8,16),(9,17),(10,17),(11,18),(12,18),(13,19),(14,19)
        // 20 vertices, 25 edges, 10 double bonds expected.
        let mut b = MoleculeBuilder::new();
        let a: Vec<_> = (0..20)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        let edges: &[(usize, usize)] = &[
            // inner 5-ring
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 0),
            // spokes
            (0, 5),
            (1, 7),
            (2, 9),
            (3, 11),
            (4, 13),
            // outer 10-ring
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 11),
            (11, 12),
            (12, 13),
            (13, 14),
            (14, 5),
            // outer "cap" bonds
            (5, 15),
            (6, 15),
            (7, 16),
            (8, 16),
            (9, 17),
            (10, 17),
            (11, 18),
            (12, 18),
            (13, 19),
            (14, 19),
        ];
        for &(x, y) in edges {
            b.add_bond(a[x], a[y], BondOrder::Aromatic).unwrap();
        }
        let mol = b.build();
        let result = kekulize(&mol);
        assert!(
            result.is_ok(),
            "corannulene kekulization failed: {:?}",
            result.err()
        );
        let doubles = result
            .unwrap()
            .values()
            .filter(|&&o| o == BondOrder::Double)
            .count();
        assert_eq!(doubles, 10, "corannulene: 10 double bonds");
    }

    /// 1-borazarobenzene (`b1ccccn1`) — aromatic B has an empty p orbital
    /// (electron acceptor), not a lone pair. B must be in a double bond.
    /// This was the last remaining kekulization failure in the 5000-molecule corpus.
    #[test]
    fn kekulize_boron_azine() {
        // 6-ring: B(0)-C(1)-C(2)-C(3)-C(4)-N(5)-B(0)
        // Valid Kekulé: B=C, C=C, C=N  (3 double bonds)
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6)
            .map(|i| {
                if i == 0 {
                    b.add_atom(Atom::aromatic(Element::B))
                } else if i == 5 {
                    b.add_atom(Atom::aromatic(Element::N))
                } else {
                    b.add_atom(Atom::aromatic(Element::C))
                }
            })
            .collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        let mol = b.build();
        let result = kekulize(&mol);
        assert!(
            result.is_ok(),
            "b1ccccn1 kekulization failed: {:?}",
            result.err()
        );
        let doubles = result
            .unwrap()
            .values()
            .filter(|&&o| o == BondOrder::Double)
            .count();
        assert_eq!(doubles, 3, "b1ccccn1: 3 double bonds");
    }

    // -----------------------------------------------------------------
    // Kekule-S0: apply_kekule must preserve stereo side channels.
    //
    // `@`/`@@` (`Atom::chirality`) is a *relative* descriptor over the
    // SMILES neighbor-encounter order recorded in `stereo_neighbor_order`
    // -- copying `chirality` alone without that order is not enough to
    // reproduce the same tetrahedral configuration downstream.
    // -----------------------------------------------------------------

    /// A chiral aromatic ring (5-membered, one sp3 substituent) with a
    /// stereo group, a stashed bond direction, and a `stereo_neighbor_order`
    /// entry set on the sp3 atom -- exercises all three side channels
    /// `apply_kekule` must preserve.
    fn chiral_aromatic_with_stereo_metadata() -> (Molecule, AtomIdx, BondIdx) {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::aromatic(Element::C));
        let c2 = b.add_atom(Atom::aromatic(Element::C));
        let c3 = b.add_atom(Atom::aromatic(Element::C));
        let c4 = b.add_atom(Atom::aromatic(Element::C));
        let o = b.add_atom(Atom::aromatic(Element::O));
        for i in 0..4 {
            let ring = [c1, c2, c3, c4, o];
            b.add_bond(ring[i], ring[i + 1], BondOrder::Aromatic)
                .unwrap();
        }
        b.add_bond(o, c1, BondOrder::Aromatic).unwrap();

        // sp3 chiral substituent on c1: -CH(F)(Cl), a real stereocenter.
        let mut chiral = Atom::new(Element::C);
        chiral.chirality = crate::atom::Chirality::CounterClockwise;
        let ch = b.add_atom(chiral);
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let exocyclic_bond = b.add_bond(c1, ch, BondOrder::Single).unwrap();
        b.add_bond(ch, f, BondOrder::Single).unwrap();
        b.add_bond(ch, cl, BondOrder::Single).unwrap();

        let mut mol = b.build();
        // SMILES-text-order neighbor sequence: c1, implicit-H sentinel, F, Cl.
        mol.set_stereo_neighbor_order(ch, vec![c1.0, STEREO_H_SENTINEL, f.0, cl.0]);
        mol.add_stereo_group(StereoGroup::new(
            crate::stereo_group::StereoGroupKind::Absolute,
            vec![ch],
        ));
        mol.set_bond_direction(exocyclic_bond, BondOrder::Up);

        (mol, ch, exocyclic_bond)
    }

    #[test]
    fn apply_kekule_preserves_stereo_neighbor_order() {
        let (mol, ch, _) = chiral_aromatic_with_stereo_metadata();
        let before = mol.stereo_neighbor_order(ch).map(|s| s.to_vec());
        assert!(before.is_some(), "test setup sanity");

        let result = kekulize(&mol).expect("kekulizable");
        assert!(
            !result.is_empty(),
            "test setup sanity: ring must need kekulization"
        );
        let kekulized = apply_kekule(&mol, &result);

        assert_eq!(
            kekulized.stereo_neighbor_order(ch).map(|s| s.to_vec()),
            before,
            "stereo_neighbor_order must survive apply_kekule verbatim"
        );
    }

    #[test]
    fn apply_kekule_preserves_stereo_groups() {
        let (mol, _, _) = chiral_aromatic_with_stereo_metadata();
        let before = mol.stereo_groups().to_vec();
        assert!(!before.is_empty(), "test setup sanity");

        let result = kekulize(&mol).expect("kekulizable");
        let kekulized = apply_kekule(&mol, &result);

        assert_eq!(
            kekulized.stereo_groups(),
            before.as_slice(),
            "stereo_groups must survive apply_kekule verbatim"
        );
    }

    #[test]
    fn apply_kekule_preserves_bond_directions() {
        let (mol, _, exocyclic_bond) = chiral_aromatic_with_stereo_metadata();
        let before = mol.bond_direction(exocyclic_bond);
        assert!(before.is_some(), "test setup sanity");

        let result = kekulize(&mol).expect("kekulizable");
        let kekulized = apply_kekule(&mol, &result);

        assert_eq!(
            kekulized.bond_direction(exocyclic_bond),
            before,
            "bond_directions must survive apply_kekule verbatim"
        );
    }

    #[test]
    fn apply_kekule_preserves_atom_and_bond_index_mapping() {
        let (mol, _, _) = chiral_aromatic_with_stereo_metadata();
        let result = kekulize(&mol).expect("kekulizable");
        let kekulized = apply_kekule(&mol, &result);

        assert_eq!(kekulized.atom_count(), mol.atom_count());
        assert_eq!(kekulized.bond_count(), mol.bond_count());
        for (idx, atom) in mol.atoms() {
            let after = kekulized.atom(idx);
            assert_eq!(atom.element, after.element, "atom {idx:?} element moved");
            assert_eq!(
                atom.chirality, after.chirality,
                "atom {idx:?} chirality moved"
            );
        }
        for (bidx, bond) in mol.bonds() {
            let after = kekulized.bond(bidx);
            assert_eq!(bond.atom1, after.atom1, "bond {bidx:?} atom1 moved");
            assert_eq!(bond.atom2, after.atom2, "bond {bidx:?} atom2 moved");
        }
    }

    #[test]
    fn apply_kekule_empty_result_is_full_clone() {
        // No aromatic bonds at all -- kekulize() returns an empty map, and
        // apply_kekule must take the early-return clone path, which trivially
        // preserves every side channel (it's the same molecule).
        let mut b = MoleculeBuilder::new();
        let mut chiral = Atom::new(Element::C);
        chiral.chirality = crate::atom::Chirality::Clockwise;
        let c = b.add_atom(chiral);
        let f = b.add_atom(Atom::new(Element::F));
        let cl = b.add_atom(Atom::new(Element::CL));
        let br = b.add_atom(Atom::new(Element::BR));
        let h_bond = b.add_atom(Atom::new(Element::I));
        b.add_bond(c, f, BondOrder::Single).unwrap();
        b.add_bond(c, cl, BondOrder::Single).unwrap();
        b.add_bond(c, br, BondOrder::Single).unwrap();
        b.add_bond(c, h_bond, BondOrder::Single).unwrap();

        let mut mol = b.build();
        mol.set_stereo_neighbor_order(c, vec![f.0, cl.0, br.0, h_bond.0]);

        let result = kekulize(&mol).expect("no aromatic bonds -- trivially kekulizable");
        assert!(result.is_empty(), "test setup sanity: no aromatic bonds");

        let kekulized = apply_kekule(&mol, &result);
        assert_eq!(
            kekulized.stereo_neighbor_order(c).map(|s| s.to_vec()),
            mol.stereo_neighbor_order(c).map(|s| s.to_vec())
        );
        assert_eq!(kekulized.atom_count(), mol.atom_count());
        assert_eq!(kekulized.bond_count(), mol.bond_count());
    }

    // -----------------------------------------------------------------
    // K1: charge-aware `atom_must_be_matched` -- fixtures from
    // docs/rfcs/aromaticity_rdkit_parity_rfc.md section 1a, previously hard
    // failures. Builders mirror `furan`/`pyrrole` above: an explicit-H
    // bracket atom is built via `Atom::aromatic` + a field override, matching
    // the SMILES bracket notation named in each doc comment.
    // -----------------------------------------------------------------

    /// Tropylium cation (`c1ccc[cH+]cc1`): 7-ring, 6 neutral aromatic CH +
    /// 1 cationic aromatic CH+. Returns the cation's `AtomIdx` alongside the
    /// molecule so tests can check its bonds specifically.
    fn tropylium_cation() -> (Molecule, AtomIdx) {
        let mut b = MoleculeBuilder::new();
        let mut atoms: Vec<AtomIdx> = (0..4)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        let mut cation = Atom::aromatic(Element::C);
        cation.charge = 1;
        cation.hydrogen_count = Some(1);
        let cation_idx = b.add_atom(cation); // ring position matches c1ccc[cH+]cc1
        atoms.push(cation_idx);
        atoms.extend((0..2).map(|_| b.add_atom(Atom::aromatic(Element::C))));
        for i in 0..7 {
            b.add_bond(atoms[i], atoms[(i + 1) % 7], BondOrder::Aromatic)
                .unwrap();
        }
        (b.build(), cation_idx)
    }

    /// Imidazolium (`c1c[nH+]c[nH]1`): 5-ring, 3 aromatic C + 1 protonated
    /// `[nH+]` + 1 neutral `[nH]`.
    fn imidazolium() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let c0 = b.add_atom(Atom::aromatic(Element::C));
        let c1 = b.add_atom(Atom::aromatic(Element::C));
        let mut n_plus = Atom::aromatic(Element::N);
        n_plus.charge = 1;
        n_plus.hydrogen_count = Some(1);
        let n2 = b.add_atom(n_plus);
        let c3 = b.add_atom(Atom::aromatic(Element::C));
        let mut n_neutral = Atom::aromatic(Element::N);
        n_neutral.hydrogen_count = Some(1);
        let n4 = b.add_atom(n_neutral);
        let atoms = [c0, c1, n2, c3, n4];
        for i in 0..5 {
            b.add_bond(atoms[i], atoms[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// Pyridinium (`c1cc[nH+]cc1`): 6-ring, 5 aromatic C + 1 protonated `[nH+]`.
    fn pyridinium() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let mut atoms: Vec<AtomIdx> = (0..5)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        let mut n_plus = Atom::aromatic(Element::N);
        n_plus.charge = 1;
        n_plus.hydrogen_count = Some(1);
        atoms.insert(3, b.add_atom(n_plus)); // ring position matches c1cc[nH+]cc1
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// Pyrylium (`c1cc[o+]cc1`): 6-ring, 5 aromatic C + 1 cationic aromatic O+.
    fn pyrylium() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let mut atoms: Vec<AtomIdx> = (0..5)
            .map(|_| b.add_atom(Atom::aromatic(Element::C)))
            .collect();
        let mut o_plus = Atom::aromatic(Element::O);
        o_plus.charge = 1;
        atoms.insert(3, b.add_atom(o_plus)); // ring position matches c1cc[o+]cc1
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// Tellurophene (`c1cc[te]c1`): 5-ring, 4 aromatic C + 1 neutral aromatic Te.
    fn tellurophene() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let te = b.add_atom(Atom::aromatic(Element::TE));
        let c1 = b.add_atom(Atom::aromatic(Element::C));
        let c2 = b.add_atom(Atom::aromatic(Element::C));
        let c3 = b.add_atom(Atom::aromatic(Element::C));
        let c4 = b.add_atom(Atom::aromatic(Element::C));
        let ring = [te, c1, c2, c3, c4];
        for i in 0..5 {
            b.add_bond(ring[i], ring[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// Phosphole (`c1cc[pH]c1`): 5-ring, 4 aromatic C + 1 neutral `[pH]`.
    fn phosphole() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let mut p_atom = Atom::aromatic(Element::P);
        p_atom.hydrogen_count = Some(1);
        let p = b.add_atom(p_atom);
        let c1 = b.add_atom(Atom::aromatic(Element::C));
        let c2 = b.add_atom(Atom::aromatic(Element::C));
        let c3 = b.add_atom(Atom::aromatic(Element::C));
        let c4 = b.add_atom(Atom::aromatic(Element::C));
        let ring = [p, c1, c2, c3, c4];
        for i in 0..5 {
            b.add_bond(ring[i], ring[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// Thiophene (4 aromatic C + 1 aromatic S, ring) -- regression coverage:
    /// neutral chalcogen donor rule must still exempt S after the K1 charge-aware fix.
    fn thiophene() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let s = b.add_atom(Atom::aromatic(Element::S));
        let c1 = b.add_atom(Atom::aromatic(Element::C));
        let c2 = b.add_atom(Atom::aromatic(Element::C));
        let c3 = b.add_atom(Atom::aromatic(Element::C));
        let c4 = b.add_atom(Atom::aromatic(Element::C));
        let ring = [s, c1, c2, c3, c4];
        for i in 0..5 {
            b.add_bond(ring[i], ring[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// Selenophene (4 aromatic C + 1 aromatic Se, ring) -- regression coverage,
    /// same rationale as `thiophene` above.
    fn selenophene() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let se = b.add_atom(Atom::aromatic(Element::SE));
        let c1 = b.add_atom(Atom::aromatic(Element::C));
        let c2 = b.add_atom(Atom::aromatic(Element::C));
        let c3 = b.add_atom(Atom::aromatic(Element::C));
        let c4 = b.add_atom(Atom::aromatic(Element::C));
        let ring = [se, c1, c2, c3, c4];
        for i in 0..5 {
            b.add_bond(ring[i], ring[(i + 1) % 5], BondOrder::Aromatic)
                .unwrap();
        }
        b.build()
    }

    /// An aromatic N cation that already reaches valence 4 (`C[n+]1(C)cccc1`)
    /// needs no ring double bond; pyridinium-type ones still do. Matches
    /// RDKit's SMILES parser on each case.
    #[test]
    fn test_kekulize_aromatic_nitrogen_cation_follows_rdkit_valence() {
        // (ring size, substituent count on N, N's H count, kekulizes)
        for (ring, subs, h, ok) in [
            (5usize, 2usize, 0u8, true), // C[n+]1(C)cccc1
            (5, 0, 2, true),             // [nH2+]1cccc1
            (6, 0, 1, true),             // c1cc[nH+]cc1
            (6, 1, 0, true),             // C[n+]1ccccc1
            (5, 1, 0, false),            // C[n+]1cccc1
        ] {
            let mut b = MoleculeBuilder::new();
            let mut n = Atom::aromatic(Element::N);
            n.charge = 1;
            n.hydrogen_count = Some(h);
            let mut atoms = vec![b.add_atom(n)];
            atoms.extend((1..ring).map(|_| b.add_atom(Atom::aromatic(Element::C))));
            for i in 0..ring {
                b.add_bond(atoms[i], atoms[(i + 1) % ring], BondOrder::Aromatic)
                    .unwrap();
            }
            for _ in 0..subs {
                let c = b.add_atom(Atom::new(Element::C));
                b.add_bond(atoms[0], c, BondOrder::Single).unwrap();
            }
            assert_eq!(
                kekulize(&b.build()).is_ok(),
                ok,
                "ring {ring}, {subs} substituents, H {h}"
            );
        }
    }

    /// `[c+]1ccccc1` and `[c-]1ccccc1` (no H) kekulize as in RDKit;
    /// `[cH+]1ccccc1` and `[cH-]1ccccc1` do not (RDKit rejects them too).
    #[test]
    fn test_kekulize_six_ring_charged_carbon_follows_h_count() {
        for (charge, h, ok) in [
            (1i8, 0u8, true),
            (1, 1, false),
            (-1, 0, true),
            (-1, 1, false),
        ] {
            let mut b = MoleculeBuilder::new();
            let mut cation = Atom::aromatic(Element::C);
            cation.charge = charge;
            cation.hydrogen_count = Some(h);
            let mut atoms = vec![b.add_atom(cation)];
            atoms.extend((0..5).map(|_| b.add_atom(Atom::aromatic(Element::C))));
            for i in 0..6 {
                b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Aromatic)
                    .unwrap();
            }
            assert_eq!(
                kekulize(&b.build()).is_ok(),
                ok,
                "charge {charge} H count {h}"
            );
        }
    }

    #[test]
    fn test_kekulize_tropylium_cation() {
        let (mol, cation_idx) = tropylium_cation();
        let result = kekulize(&mol).expect("tropylium cation kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(
            doubles, 3,
            "tropylium: 6 neutral C alternate into 3 double bonds"
        );
        // The cationic carbon must end up with only single bonds.
        for (bidx, order) in &result {
            let bond = mol.bond(*bidx);
            if bond.atom1 == cation_idx || bond.atom2 == cation_idx {
                assert_eq!(
                    *order,
                    BondOrder::Single,
                    "cationic C must not get a double bond"
                );
            }
        }
    }

    #[test]
    fn test_kekulize_imidazolium() {
        let mol = imidazolium();
        let result = kekulize(&mol).expect("imidazolium kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(
            doubles, 2,
            "imidazolium: 2 double bonds ([nH+] matched, [nH] not)"
        );
    }

    #[test]
    fn test_kekulize_pyridinium() {
        let mol = pyridinium();
        let result = kekulize(&mol).expect("pyridinium kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(
            doubles, 3,
            "pyridinium: 3 double bonds, same as neutral pyridine"
        );
    }

    #[test]
    fn test_kekulize_pyrylium() {
        let mol = pyrylium();
        let result = kekulize(&mol).expect("pyrylium kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(
            doubles, 3,
            "pyrylium: 3 double bonds, O+ matched like pyridine's N"
        );
    }

    #[test]
    fn test_kekulize_tellurophene() {
        let mol = tellurophene();
        let result = kekulize(&mol).expect("tellurophene kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(
            doubles, 2,
            "tellurophene: 2 double bonds, Te donates its lone pair"
        );
    }

    #[test]
    fn test_kekulize_phosphole() {
        let mol = phosphole();
        let result = kekulize(&mol).expect("phosphole kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(
            doubles, 2,
            "phosphole: 2 double bonds, [pH] donates its lone pair like [nH]"
        );
    }

    #[test]
    fn test_kekulize_thiophene_regression() {
        let mol = thiophene();
        let result = kekulize(&mol).expect("thiophene kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(
            doubles, 2,
            "thiophene must still kekulize after the K1 charge-aware fix"
        );
    }

    #[test]
    fn test_kekulize_selenophene_regression() {
        let mol = selenophene();
        let result = kekulize(&mol).expect("selenophene kekulization failed");
        let doubles = result.values().filter(|&&o| o == BondOrder::Double).count();
        assert_eq!(
            doubles, 2,
            "selenophene must still kekulize after the K1 charge-aware fix"
        );
    }
}
