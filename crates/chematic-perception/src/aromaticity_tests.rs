#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};

    // =========================================================================
    // Molecule builder helpers (kekulized, manually constructed)
    // =========================================================================

    fn benzene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(atoms[i], atoms[(i + 1) % 6], order).unwrap();
        }
        b.build()
    }

    fn cyclohexane() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..6).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..6 {
            b.add_bond(atoms[i], atoms[(i + 1) % 6], BondOrder::Single)
                .unwrap();
        }
        b.build()
    }

    fn pyridine_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let n = b.add_atom(Atom::new(Element::N));
        let atoms_c: Vec<_> = (0..5).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        let ring = [
            n, atoms_c[0], atoms_c[1], atoms_c[2], atoms_c[3], atoms_c[4],
        ];
        for i in 0..6 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(ring[i], ring[(i + 1) % 6], order).unwrap();
        }
        b.build()
    }

    fn furan_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::new(Element::O));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        let ring = [o, c1, c2, c3, c4];
        b.add_bond(ring[0], ring[1], BondOrder::Single).unwrap();
        b.add_bond(ring[1], ring[2], BondOrder::Double).unwrap();
        b.add_bond(ring[2], ring[3], BondOrder::Single).unwrap();
        b.add_bond(ring[3], ring[4], BondOrder::Double).unwrap();
        b.add_bond(ring[4], ring[0], BondOrder::Single).unwrap();
        b.build()
    }

    fn pyrrole_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let mut n_atom = Atom::new(Element::N);
        n_atom.hydrogen_count = Some(1);
        let n = b.add_atom(n_atom);
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        let ring = [n, c1, c2, c3, c4];
        b.add_bond(ring[0], ring[1], BondOrder::Single).unwrap();
        b.add_bond(ring[1], ring[2], BondOrder::Double).unwrap();
        b.add_bond(ring[2], ring[3], BondOrder::Single).unwrap();
        b.add_bond(ring[3], ring[4], BondOrder::Double).unwrap();
        b.add_bond(ring[4], ring[0], BondOrder::Single).unwrap();
        b.build()
    }

    /// Same ring as `pyrrole_kekule()`, but the N has NO explicit
    /// `hydrogen_count` — matching how the SMILES parser actually builds a
    /// bare, non-bracket `N` (e.g. from `Chem.Kekulize` + non-canonical
    /// `MolToSmiles(kekuleSmiles=True)` round-tripping an `[nH]`-written
    /// pyrrole/imidazole/purine nitrogen). `pyrrole_kekule()` above sidesteps
    /// the bug this reproduces by setting `hydrogen_count` manually.
    fn pyrrole_kekule_implicit_h() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let n = b.add_atom(Atom::new(Element::N));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        let ring = [n, c1, c2, c3, c4];
        b.add_bond(ring[0], ring[1], BondOrder::Single).unwrap();
        b.add_bond(ring[1], ring[2], BondOrder::Double).unwrap();
        b.add_bond(ring[2], ring[3], BondOrder::Single).unwrap();
        b.add_bond(ring[3], ring[4], BondOrder::Double).unwrap();
        b.add_bond(ring[4], ring[0], BondOrder::Single).unwrap();
        b.build()
    }

    fn naphthalene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..10).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        let ring1 = [0usize, 1, 2, 3, 4, 9];
        let orders1 = [
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ];
        for i in 0..6 {
            b.add_bond(atoms[ring1[i]], atoms[ring1[(i + 1) % 6]], orders1[i])
                .unwrap();
        }
        let ring2_extra = [(4, 5), (5, 6), (6, 7), (7, 8), (8, 9)];
        let orders2 = [
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ];
        for (i, &(a, bb)) in ring2_extra.iter().enumerate() {
            b.add_bond(atoms[a], atoms[bb], orders2[i]).unwrap();
        }
        b.build()
    }

    fn cyclobutadiene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..4 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(atoms[i], atoms[(i + 1) % 4], order).unwrap();
        }
        b.build()
    }

    fn cyclooctatetraene_kekule() -> chematic_core::Molecule {
        let mut b = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..8).map(|_| b.add_atom(Atom::new(Element::C))).collect();
        for i in 0..8 {
            let order = if i % 2 == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            b.add_bond(atoms[i], atoms[(i + 1) % 8], order).unwrap();
        }
        b.build()
    }

    /// Helper: parse an aromatic SMILES and return the molecule with aromatic bonds
    /// (no kekulization).  Use for compounds where kekulization is unsupported.
    #[cfg(test)]
    fn mol_aromatic(smiles: &str) -> chematic_core::Molecule {
        chematic_smiles::parse(smiles).expect("valid SMILES")
    }

    /// Helper: parse SMILES and kekulize.  Panics if kekulization fails.
    #[cfg(test)]
    fn mol_kekulized(smiles: &str) -> chematic_core::Molecule {
        let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
        let k = chematic_core::kekulize(&mol).expect("kekulizable");
        chematic_core::apply_kekule(&mol, &k)
    }

    // =========================================================================
    // Regression: kekulized single-ring aromatics (Pass 1 only, no context)
    // =========================================================================

    #[test]
    fn test_benzene_is_aromatic() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "all 6 benzene atoms aromatic"
        );
        for i in 0..6u32 {
            assert!(model.is_atom_aromatic(AtomIdx(i)));
        }
    }

    #[test]
    fn test_cyclohexane_not_aromatic() {
        let mol = cyclohexane();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 0, "cyclohexane not aromatic");
    }

    #[test]
    fn test_pyridine_is_aromatic() {
        let mol = pyridine_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6);
    }

    #[test]
    fn test_furan_is_aromatic() {
        let mol = furan_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5);
    }

    #[test]
    fn test_pyrrole_is_aromatic() {
        let mol = pyrrole_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5);
    }

    #[test]
    fn test_apply_aromaticity_preserves_pyrrole_nh_implicit_hydrogen() {
        // Regression test: apply_aromaticity_ex() normalizes all aromatic-
        // model ring bonds to BondOrder::Aromatic, which discards the
        // Kekule Single/Double pattern that distinguishes a pyrrole-type N
        // (needs 1 implicit H) from a pyridine-type N (needs 0) once both
        // have exactly 2 aromatic-order ring bonds and no explicit bracket
        // H count. Without preserving the pre-normalization value,
        // implicit_hcount() on the perceived molecule silently returns 0
        // instead of 1 for the unsubstituted pyrrole N -- wrong molecular
        // formula/weight, and a representation-dependent divergence from
        // the same molecule parsed directly from aromatic-written SMILES
        // (where `[nH]`'s bracket H count is correct by construction).
        let mol = pyrrole_kekule_implicit_h();
        let n_idx = AtomIdx(0);
        assert_eq!(
            implicit_hcount(&mol, n_idx),
            1,
            "pre-normalization: bare N with 2 single ring bonds must show 1 implicit H"
        );

        let perceived = apply_aromaticity(&mol);
        assert!(perceived.atom(n_idx).aromatic, "ring N must be aromatic");
        assert_eq!(
            implicit_hcount(&perceived, n_idx),
            1,
            "post-apply_aromaticity: pyrrole N must still show 1 implicit H, not 0"
        );
    }

    #[test]
    fn test_apply_aromaticity_does_not_add_h_to_pyridine_type_n() {
        // Sibling check to the pyrrole regression above: a pyridine-type
        // ring N (1 ring single + 1 ring double pre-normalization, no H)
        // must NOT gain a spurious implicit H from the preservation logic --
        // its pre- and post-normalization implicit_hcount already agree
        // (both 0), so it must be left untouched.
        let mol = pyridine_kekule();
        let n_idx = AtomIdx(0);
        assert_eq!(implicit_hcount(&mol, n_idx), 0);

        let perceived = apply_aromaticity(&mol);
        assert!(perceived.atom(n_idx).aromatic);
        assert_eq!(implicit_hcount(&perceived, n_idx), 0);
        assert_eq!(
            perceived.atom(n_idx).hydrogen_count,
            None,
            "pyridine N must not gain an explicit hydrogen_count -- would force spurious bracket notation"
        );
    }

    #[test]
    fn test_naphthalene_both_rings_aromatic() {
        let mol = naphthalene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "all 10 naphthalene atoms aromatic"
        );
    }

    #[test]
    fn test_bond_aromaticity_benzene() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        let count = mol
            .bonds()
            .filter(|(b, _)| model.is_bond_aromatic(*b))
            .count();
        assert_eq!(count, 6);
    }

    #[test]
    fn test_apply_aromaticity_benzene() {
        let mol = benzene_kekule();
        let aromatic = apply_aromaticity(&mol);
        for (_, atom) in aromatic.atoms() {
            assert!(atom.aromatic, "every benzene carbon should be aromatic");
        }
        let aromatic_bond_count = aromatic
            .bonds()
            .filter(|(_, b)| b.order == BondOrder::Aromatic)
            .count();
        assert_eq!(aromatic_bond_count, 6);
    }

    #[test]
    fn test_apply_aromaticity_cyclohexane_unchanged() {
        let mol = cyclohexane();
        let result = apply_aromaticity(&mol);
        for (_, atom) in result.atoms() {
            assert!(!atom.aromatic);
        }
        for (_, bond) in result.bonds() {
            assert_ne!(bond.order, BondOrder::Aromatic);
        }
    }

    // =========================================================================
    // Antiaromaticity
    // =========================================================================

    #[test]
    fn test_cyclobutadiene_antiaromatic() {
        let mol = cyclobutadiene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            0,
            "cyclobutadiene not aromatic"
        );
        assert!(model.has_antiaromaticity(), "cyclobutadiene antiaromatic");
        assert_eq!(model.antiaromatic_rings().len(), 1);
        let classifications = model.ring_classifications();
        assert_eq!(classifications.len(), 1);
        assert_eq!(classifications[0].1, RingAromaticity::Antiaromatic);
        assert_eq!(classifications[0].2, 4);
    }

    #[test]
    fn test_cyclooctatetraene_antiaromatic() {
        let mol = cyclooctatetraene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 0, "COT not aromatic");
        assert!(model.has_antiaromaticity(), "COT antiaromatic");
        assert_eq!(model.antiaromatic_rings().len(), 1);
        let cls = &model.ring_classifications()[0];
        assert_eq!(cls.1, RingAromaticity::Antiaromatic);
        assert_eq!(cls.2, 8);
    }

    // =========================================================================
    // Ring classifications
    // =========================================================================

    #[test]
    fn test_ring_classifications_benzene() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        let classifications = model.ring_classifications();
        assert_eq!(classifications.len(), 1);
        assert_eq!(classifications[0].1, RingAromaticity::Aromatic);
        assert_eq!(classifications[0].2, 6);
    }

    #[test]
    fn test_ring_classifications_naphthalene() {
        let mol = naphthalene_kekule();
        let model = assign_aromaticity(&mol);
        let classifications = model.ring_classifications();
        assert_eq!(classifications.len(), 2, "naphthalene has two rings");
        for (_, classification, count) in classifications {
            assert_eq!(*classification, RingAromaticity::Aromatic);
            assert_eq!(*count, 6);
        }
    }

    #[test]
    fn test_non_aromatic_cyclohexane() {
        let mol = cyclohexane();
        let model = assign_aromaticity(&mol);
        for (_, classification, _) in model.ring_classifications() {
            assert_ne!(*classification, RingAromaticity::Aromatic);
            assert_ne!(*classification, RingAromaticity::Antiaromatic);
        }
    }

    // =========================================================================
    // Electron distribution
    // =========================================================================

    #[test]
    fn test_thiophene_aromatic() {
        let mut b = MoleculeBuilder::new();
        let s = b.add_atom(Atom::new(Element::S));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        let ring = [s, c1, c2, c3, c4];
        b.add_bond(ring[0], ring[1], BondOrder::Single).unwrap();
        b.add_bond(ring[1], ring[2], BondOrder::Double).unwrap();
        b.add_bond(ring[2], ring[3], BondOrder::Single).unwrap();
        b.add_bond(ring[3], ring[4], BondOrder::Double).unwrap();
        b.add_bond(ring[4], ring[0], BondOrder::Single).unwrap();
        let mol = b.build();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5);
        assert_eq!(model.ring_classifications()[0].2, 6);
    }

    #[test]
    fn issue_767_default_api_rejects_neutral_hypercoordinate_chalcogen_donors() {
        for smi in [
            "C1=CS(O)(O)C=C1",
            "C1=CS(O)C=C1",
            "c1s(O)(O)ccc1",
            "O=s1cccc1",
            "O=s1(=O)cccc1",
        ] {
            let mol = mol_aromatic(smi);
            let perceived = apply_aromaticity(&mol);
            assert_eq!(
                perceived.atoms().filter(|(_, atom)| atom.aromatic).count(),
                0,
                "{smi}: stale aromatic atom flags must not survive"
            );
            assert_eq!(
                perceived
                    .bonds()
                    .filter(|(_, bond)| bond.order == BondOrder::Aromatic)
                    .count(),
                0,
                "{smi}: stale aromatic bond orders must not survive"
            );
        }
    }

    #[test]
    fn test_electron_distribution_tracking() {
        let mol = benzene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.ring_classifications()[0].2, 6, "benzene: 6 × 1π = 6");

        let mol = pyrrole_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.ring_classifications()[0].2,
            6,
            "pyrrole: N(2π) + 4C(1π) = 6"
        );

        let mol = furan_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.ring_classifications()[0].2,
            6,
            "furan: O(2π) + 4C(1π) = 6"
        );
    }

    // =========================================================================
    // Aromatic-SMILES input (BondOrder::Aromatic, no kekulization)
    // Verifies that assign_aromaticity works on pre-kekulization molecules.
    // =========================================================================

    #[test]
    fn test_benzene_aromatic_smiles() {
        // c1ccccc1 — parsed with BondOrder::Aromatic bonds
        let mol = mol_aromatic("c1ccccc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "benzene from aromatic SMILES"
        );
    }

    #[test]
    fn test_naphthalene_aromatic_smiles() {
        let mol = mol_aromatic("c1ccc2ccccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "naphthalene from aromatic SMILES"
        );
    }

    #[test]
    fn test_pyridine_aromatic_smiles() {
        let mol = mol_aromatic("c1ccncc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "pyridine from aromatic SMILES"
        );
    }

    #[test]
    fn test_furan_aromatic_smiles() {
        let mol = mol_aromatic("c1ccoc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5, "furan from aromatic SMILES");
    }

    #[test]
    fn test_pyrrole_aromatic_smiles() {
        // [nH] bracket atom: hydrogen_count = Some(1)
        let mol = mol_aromatic("c1cc[nH]c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            5,
            "pyrrole from aromatic SMILES"
        );
    }

    #[test]
    fn test_thiophene_aromatic_smiles() {
        let mol = mol_aromatic("c1ccsc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            5,
            "thiophene from aromatic SMILES"
        );
    }

    // =========================================================================
    // Fused-ring kekulized systems (Pass 2 propagation)
    // =========================================================================

    #[test]
    fn test_indole_aromatic() {
        // c1ccc2[nH]ccc2c1 — indole (9 atoms, 5-ring + 6-ring fused)
        let mol = mol_kekulized("c1ccc2[nH]ccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 indole atoms aromatic"
        );
    }

    #[test]
    fn test_benzimidazole_aromatic() {
        // Two N atoms in fused 5+6 ring system
        let mol = mol_kekulized("c1ccc2[nH]cnc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 9, "all 9 benzimidazole atoms");
    }

    #[test]
    fn test_quinoline_aromatic() {
        let mol = mol_kekulized("c1ccc2ncccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 10, "all 10 quinoline atoms");
    }

    #[test]
    fn test_acridine_aromatic() {
        // 3 fused 6-membered rings, central N: 13 atoms
        let mol = mol_kekulized("c1ccc2nc3ccccc3cc2c1");
        let model = assign_aromaticity(&mol);
        // acridine is C13H9N → 14 heavy atoms (13 C + 1 N), all aromatic
        assert_eq!(model.aromatic_atom_count(), 14, "all 14 acridine atoms");
    }

    // =========================================================================
    // Fused-ring aromatic-SMILES input (BondOrder::Aromatic, kekulize fails)
    // =========================================================================

    #[test]
    fn test_indolizine_aromatic() {
        // c1ccn2cccc2c1 — indolizine: bridgehead N, kekulization unsupported.
        // The SSSR finds a 6-ring and a 9-ring; the 5-ring is recovered via
        // augmentation (XOR of 6- and 9-ring).
        // Pass 1: 5-ring (augmented) detected via bridgehead-N rule → 6π.
        // Pass 2: 6-ring detected using N already aromatic from 5-ring → 6π.
        // The 9-ring (SSSR artifact) is NonAromatic (9π ≠ 4n+2), but all
        // 9 atoms are correctly flagged aromatic via the 5- and 6-ring.
        let mol = mol_aromatic("c1ccn2cccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 indolizine atoms aromatic"
        );
        // At least the 6-ring should be classified as Aromatic in the SSSR set.
        let has_aromatic_ring = model
            .ring_classifications()
            .iter()
            .any(|(_, cls, _)| *cls == RingAromaticity::Aromatic);
        assert!(has_aromatic_ring, "at least one SSSR ring aromatic");
    }

    #[test]
    #[ignore = "PROVISIONAL: regressed by the Horton SSSR fix, see comment below"]
    fn test_purine_aromatic() {
        // c1cnc2[nH]cnc2n1 — purine: 9 atoms, kekulizable
        //
        // Regressed by the Horton SSSR rewrite (confirmed passing on the old
        // single-spanning-tree find_sssr, failing only after Horton; see
        // debug dump captured during diagnosis). Root cause, empirically
        // confirmed: the 6-membered ring (pyrimidine-type) passes Pass 1
        // alone (6π) and marks its atoms aromatic. The 5-membered ring
        // (imidazole-type) evaluates to 4π in isolation — its two fusion
        // carbons each have their only double bond exocyclic to a ring N,
        // which the exocyclic-to-heteroatom rule scores as 0π — and 4π trips
        // `classify_ring_aromaticity`'s "4n → Antiaromatic" branch. Pass 1
        // treats Antiaromatic as definitive and never retries it in Pass 2,
        // even though the fusion carbons would each contribute 1π (not 0π)
        // once `aromatic_context` recognizes them as already-aromatic — that
        // recount gives 6π (aromatic). The old, non-minimal SSSR never hit
        // this path because it fed a different (structurally wrong) ring set
        // into Pass 1 in the first place.
        //
        // Fix belongs in the aromatic_context-removal PR (see
        // greedy-hopping-crescent.md step 5), not here: retrying
        // Antiaromatic rings in Pass 2 is a real fix, but must not be
        // bundled into the SSSR PR per the "measure free recoveries with
        // zero aromaticity.rs changes" staging requirement.
        let mol = mol_kekulized("c1cnc2[nH]cnc2n1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 purine atoms aromatic"
        );
    }

    #[test]
    fn test_purine_aromatic_from_aromatic_smiles() {
        let mol = mol_aromatic("c1cnc2[nH]cnc2n1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "purine from aromatic SMILES"
        );
    }

    #[test]
    fn test_2_pyridinone_aromatic() {
        // O=c1ccncc1 — 2-pyridinone (aromatic SMILES, N without H, exo C=O).
        // Kekulization fails; tested on the aromatic-bond form directly.
        // The exo C=O gives the C atom has_double_any=true → 1π.
        // N has Aromatic bonds in ring → 1π (pyridine-like).
        // Total: 6 × 1π = 6π → aromatic.
        let mol = mol_aromatic("O=c1ccncc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "all 6 ring atoms of 2-pyridinone aromatic"
        );
    }

    #[test]
    fn test_quinolone_aromatic() {
        // O=c1ccc2ncccc2c1 — quinolone: fused 6+6 with exo C=O, kekulize fails
        let mol = mol_aromatic("O=c1ccc2ncccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "all 10 quinolone ring atoms aromatic"
        );
        assert_eq!(
            model.ring_classifications().len(),
            2,
            "two rings classified"
        );
    }

    #[test]
    fn test_indole_aromatic_smiles() {
        let mol = mol_aromatic("c1ccc2[nH]ccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "indole from aromatic SMILES"
        );
    }

    // =========================================================================
    // Bridgehead N rule: specifically test that the rule fires correctly
    // =========================================================================

    #[test]
    fn test_bridgehead_n_contributes_lone_pair() {
        // Indolizine: the bridgehead N (degree 3, no H, no explicit double bond)
        // must be detected as a 2π contributor for the 5-membered ring.
        // We verify by checking the 5-ring classification (if accessible).
        let mol = mol_aromatic("c1ccn2cccc2c1");
        let model = assign_aromaticity(&mol);
        // All 9 atoms aromatic: both rings must be aromatic.
        assert_eq!(model.aromatic_atom_count(), 9);
        // The bridgehead N itself must be in the aromatic set.
        // In the SMILES c1ccn2cccc2c1, n is atom index 3.
        assert!(
            model.is_atom_aromatic(AtomIdx(3)),
            "bridgehead N must be aromatic"
        );
    }

    #[test]
    fn test_non_bridgehead_n_no_false_positive() {
        // Pyrimidine: two N atoms in a 6-membered ring, no bridgehead.
        // Both N have ring_degree == total_degree == 2.
        // Should be detected as aromatic via has_aromatic_in_ring (Aromatic bonds).
        let mol = mol_aromatic("c1ccncn1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6, "pyrimidine is aromatic");
    }

    #[test]
    fn test_imidazole_aromatic() {
        // c1cn[nH]c1 / c1c[nH]cn1 — imidazole: one pyridine-type N, one pyrrole-type N
        let mol = mol_aromatic("c1cn[nH]c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 5, "imidazole is aromatic");
    }

    // =========================================================================
    // Pass 2 specifically: rings that need fused-ring context
    // =========================================================================

    #[test]
    fn test_pass2_needed_for_indolizine_6ring() {
        // The augmented 5-ring (XOR of SSSR 6-ring and 9-ring) is detected aromatic in Pass 1.
        // The SSSR 6-ring is then detected aromatic in Pass 2 (N already aromatic → 1π).
        // The SSSR 9-ring (9π) remains NonAromatic per Hückel.
        // Key assertion: all 9 atoms are aromatic (correct overall perception).
        let mol = mol_aromatic("c1ccn2cccc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 indolizine atoms aromatic"
        );
        // The bridgehead N must be aromatic.
        assert!(
            model.is_atom_aromatic(AtomIdx(3)),
            "bridgehead N is aromatic"
        );
        // The 6-ring (SSSR ring, improved by Pass 2) should be classified Aromatic.
        let aromatic_count = model
            .ring_classifications()
            .iter()
            .filter(|(_, cls, _)| *cls == RingAromaticity::Aromatic)
            .count();
        assert!(aromatic_count >= 1, "at least one SSSR ring is aromatic");
    }

    #[test]
    fn test_no_pass2_needed_for_naphthalene() {
        // Naphthalene: both rings pass independently in Pass 1.
        // Verifies Pass 2 doesn't break things that already work.
        let mol = naphthalene_kekule();
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 10);
        let classes = model.ring_classifications();
        assert_eq!(classes.len(), 2);
        for (_, cls, _) in classes {
            assert_eq!(*cls, RingAromaticity::Aromatic);
        }
    }

    #[test]
    fn test_anthracene_aromatic() {
        // c1ccc2cc3ccccc3cc2c1 — anthracene: 3 linearly fused 6-rings, 14 atoms
        let mol = mol_kekulized("c1ccc2cc3ccccc3cc2c1");
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 14, "all 14 anthracene atoms");
    }

    // =========================================================================
    // Regression: aromatic-bond path must not perturb kekulized correctness
    // =========================================================================

    #[test]
    fn test_kekulized_path_unaffected_by_aromatic_bond_changes() {
        // Kekulized benzene: bonds are Double/Single, not Aromatic.
        // The new Aromatic-bond branches must stay dormant.
        let mol = benzene_kekule();
        // Verify no aromatic bonds in input.
        for (_, bond) in mol.bonds() {
            assert_ne!(bond.order, BondOrder::Aromatic, "input must be kekulized");
        }
        let model = assign_aromaticity(&mol);
        assert_eq!(model.aromatic_atom_count(), 6);
        // All 6 bonds in benzene ring should be aromatic.
        let aromatic_bonds = mol
            .bonds()
            .filter(|(b, _)| model.is_bond_aromatic(*b))
            .count();
        assert_eq!(aromatic_bonds, 6);
    }

    #[test]
    fn test_keto_pyridinone_aromatic() {
        // O=C1NC=CC=C1 — 2-pyridinone keto form with N-H.
        // π count: C(=O)(0π, exocyclic-only double bond to O) + N-H(2π) +
        // 4×C in 2 ring C=C (1π each) = 6π → aromatic. Matches RDKit, which
        // marks all 6 ring atoms aromatic (exocyclic O stays non-aromatic).
        let mol = mol_kekulized("O=C1NC=CC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "keto pyridinone ring is Hückel aromatic (6π = 4n+2)"
        );
    }

    #[test]
    fn test_tropone_aromatic() {
        // O=C1C=CC=CC=C1 — tropone (cycloheptatrienone), Kekulized input.
        // Carbonyl C contributes 0π (exocyclic-only double bond to O); the
        // other 6 ring carbons contribute 1π each from 3 endocyclic C=C.
        // Total 6π → aromatic, matching RDKit (all 7 ring atoms aromatic).
        let mol = mol_kekulized("O=C1C=CC=CC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            7,
            "all 7 tropone ring atoms aromatic"
        );
    }

    #[test]
    fn test_4_pyridone_aromatic() {
        // O=C1C=CNC=C1 — 4-pyridone, Kekulized input. Same 6π accounting as
        // 2-pyridone, just with N para to the carbonyl. Matches RDKit.
        let mol = mol_kekulized("O=C1C=CNC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "all 6 4-pyridone ring atoms aromatic"
        );
    }

    #[test]
    fn test_pyranone_aromatic() {
        // O=C1C=COC=C1 — 4H-pyran-4-one, Kekulized input. Ring O contributes
        // 2π (lone pair), carbonyl C contributes 0π, remaining 4 ring carbons
        // contribute 1π each from 2 endocyclic C=C. Total 6π. Matches RDKit.
        let mol = mol_kekulized("O=C1C=COC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "all 6 pyranone ring atoms aromatic"
        );
    }

    #[test]
    fn test_cyclopentadienyl_anion_aromatic() {
        // [CH-]1C=CC=C1 — cyclopentadienyl anion. The carbanion carbon has no
        // double bond but contributes 2π (lone pair); the other 4 carbons
        // contribute 1π each from 2 endocyclic C=C. Total 6π. Matches RDKit
        // (all 5 atoms aromatic).
        let mol = mol_kekulized("[CH-]1C=CC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            5,
            "all 5 cyclopentadienyl anion atoms aromatic"
        );
    }

    // ── K2a: charge-aware ring_pi_electrons -- tropylium/imidazolium/
    // pyridinium/pyrylium now genuinely confirmed aromatic by the raw
    // Huckel model itself (not just a stale parser flag surviving), under
    // BOTH documented calling conventions (`apply_aromaticity`'s own doc
    // comment: "may be kekulized... or may retain Aromatic bond orders from
    // the SMILES parser"). RDKit-verified: all four are aromatic cations,
    // all-atom/all-bond, per rdkit==2026.03.3 (see
    // docs/rfcs/aromaticity_rdkit_parity_rfc.md and the K2a PR description for
    // the full 40-fixture oracle re-run against a live RDKit).
    //
    // K1 (fix/kekulize-charge-aware-k1, already merged) made
    // chematic_core::kekulize() succeed for all four; this fix is the
    // separate, independent charge-blindness bug in the Huckel
    // pi-electron-counting layer (`ring_pi_electrons`) that K1 explicitly
    // did not touch. Deliberately does NOT touch `build_molecule_from_model`
    // (that promote-only-vs-demote question is tracked separately as K2b) --
    // these four fixtures need no demotion at all: their atom flags were
    // already `true` from the aromatic-notation parse, and once the model
    // itself confirms the ring, the EXISTING promote-only bond loop already
    // correctly promotes their bonds to `Aromatic` for the first time. That
    // is what actually fixes the pre-existing atom/bond flag inconsistency
    // for these four -- no demotion capability required.
    fn assert_fully_aromatic(mol: &Molecule, n: usize, label: &str) {
        let applied = apply_aromaticity(mol);
        for (idx, atom) in applied.atoms() {
            assert!(atom.aromatic, "{label}: atom {idx:?} should be aromatic");
        }
        assert_eq!(applied.atom_count(), n, "{label}: unexpected atom count");
        for (_, bond) in applied.bonds() {
            assert_eq!(
                bond.order,
                BondOrder::Aromatic,
                "{label}: every ring bond should end up Aromatic order"
            );
        }
    }

    #[test]
    fn test_tropylium_cation_aromatic_raw_and_kekulized() {
        let raw = chematic_smiles::parse("c1ccc[cH+]cc1").expect("valid SMILES");
        assert_fully_aromatic(&raw, 7, "tropylium (raw)");
        let kek = mol_kekulized("c1ccc[cH+]cc1");
        assert_fully_aromatic(&kek, 7, "tropylium (kekulized)");
        assert_eq!(
            assign_aromaticity(&raw).aromatic_atom_count(),
            7,
            "tropylium: raw model itself must confirm all 7 atoms, not rely on a stale flag"
        );
        assert_eq!(
            assign_aromaticity(&kek).aromatic_atom_count(),
            7,
            "tropylium: kekulized model itself must confirm all 7 atoms"
        );
    }

    #[test]
    fn test_imidazolium_aromatic_raw_and_kekulized() {
        let raw = chematic_smiles::parse("c1c[nH+]c[nH]1").expect("valid SMILES");
        assert_fully_aromatic(&raw, 5, "imidazolium (raw)");
        let kek = mol_kekulized("c1c[nH+]c[nH]1");
        assert_fully_aromatic(&kek, 5, "imidazolium (kekulized)");
        assert_eq!(assign_aromaticity(&raw).aromatic_atom_count(), 5);
        assert_eq!(assign_aromaticity(&kek).aromatic_atom_count(), 5);
    }

    #[test]
    fn test_pyridinium_aromatic_raw_and_kekulized() {
        let raw = chematic_smiles::parse("c1cc[nH+]cc1").expect("valid SMILES");
        assert_fully_aromatic(&raw, 6, "pyridinium (raw)");
        let kek = mol_kekulized("c1cc[nH+]cc1");
        assert_fully_aromatic(&kek, 6, "pyridinium (kekulized)");
        assert_eq!(assign_aromaticity(&raw).aromatic_atom_count(), 6);
        assert_eq!(assign_aromaticity(&kek).aromatic_atom_count(), 6);
    }

    #[test]
    fn test_pyrylium_aromatic_raw_and_kekulized() {
        let raw = chematic_smiles::parse("c1cc[o+]cc1").expect("valid SMILES");
        assert_fully_aromatic(&raw, 6, "pyrylium (raw)");
        let kek = mol_kekulized("c1cc[o+]cc1");
        assert_fully_aromatic(&kek, 6, "pyrylium (kekulized)");
        assert_eq!(assign_aromaticity(&raw).aromatic_atom_count(), 6);
        assert_eq!(assign_aromaticity(&kek).aromatic_atom_count(), 6);
    }

    // ── K2a scope guard: tellurophene/phosphole are explicitly NOT fixed by
    // the charge-aware change above (they need real Se/Te/P electron-donor
    // support in the default Huckel engine, out of scope -- see the K2a/K2b
    // PR descriptions). Pin the current (still-gap) count so a future
    // change to this area doesn't silently start claiming these are fixed
    // without an explicit, source-grounded review.
    #[test]
    fn test_tellurophene_and_phosphole_still_unsupported_under_default_huckel() {
        let te = mol_kekulized("c1cc[te]c1");
        assert_eq!(
            assign_aromaticity(&te).aromatic_atom_count(),
            0,
            "tellurophene: still unsupported under default Huckel (K2a does not add Te support)"
        );
        let p = mol_kekulized("c1cc[pH]c1");
        assert_eq!(
            assign_aromaticity(&p).aromatic_atom_count(),
            0,
            "phosphole: still unsupported under default Huckel (K2a does not add P support)"
        );
    }

    // ── K2b fused-diazine fix (fix/aromaticity-flag-demotion-k2b follow-up) ─
    //
    // Opt-in only, via `assign_aromaticity_authoritative_experimental` --
    // per coordinator decision, `apply_aromaticity`/`apply_aromaticity_ex`
    // (and the plain `assign_aromaticity`/`assign_aromaticity_ex` they call)
    // stay byte-identical to their pre-K2b behavior. The
    // `test_known_gap_fused_diazine_exocyclic_misfire_antiaromatic` pin that
    // used to live here (asserting the DEFAULT engine's wrong 6/10 count) is
    // superseded by `test_default_engine_unaffected_by_fused_diazine_fix`
    // below (same assertion, renamed for clarity: this is now a permanent
    // "default stays reverted" guard, not a "known gap" pin -- the gap is
    // only closed for the opt-in engine, not fixed in the default at all).
    // The azulene pin further below is untouched either way (separate,
    // still-open, out-of-scope mechanism, never affected by this fix in
    // ANY engine).

    #[test]
    fn test_authoritative_experimental_fixes_fused_diazine_ring_fusion() {
        // c1cnc2ccccc2n1 -- a bare, unsubstituted naphthyridine isomer (15
        // chars, no substituents). RDKit: fully aromatic, all 10 atoms/bonds
        // (verified live against rdkit==2026.03.3). Under the DEFAULT engine
        // (`assign_aromaticity`), chematic confirms only 6/10 (the
        // pyridine-type ring) -- see
        // `test_default_engine_unaffected_by_fused_diazine_fix` below: the
        // benzo ring's Pass 1 evaluation wrongly zeroes out BOTH its fusion
        // carbons via `CarbonExocyclicHeteroatomDouble` (each fusion
        // carbon's own Kekule double bond points into the OTHER (pyridine)
        // ring, toward a nitrogen there -- from the benzo ring's own,
        // single-ring-only perspective using only `ring_atom_set`, that bond
        // looks exactly like a genuine exocyclic C=O/C=N substituent
        // (tropone's shape), which is what that rule is actually meant to
        // catch). Landing on EXACTLY pi=4 classifies the ring `Antiaromatic`,
        // which Pass 2 never retries ("definitive, do not retry").
        //
        // Fixed under the OPT-IN `assign_aromaticity_authoritative_experimental`
        // engine by making the rule bond-level (`all_ring_bonds`, built once
        // from every SSSR/augmented ring): a double bond whose far atom sits
        // on a DIFFERENT ring is a ring-fusion bond, not a substituent, so it
        // no longer zeroes the atom -- both fusion carbons now fall through
        // to the ordinary sp2 default (1π each), the benzo ring lands on
        // pi=6 (Aromatic) directly in Pass 1, and Pass 2 promotes the
        // pyridine-type ring via `AlreadyAromaticContext` as before.
        // Confirmed Kekule-choice-dependent, not shape-dependent: plain
        // quinoxaline and quinazoline (`c1ccc2nccnc2c1`, `c1ccc2ncncc2c1`)
        // never reproduced this in the first place (only ONE fusion carbon
        // was affected for those, landing on the retryable odd pi=5
        // NonAromatic case). This molecule was constructed as a minimal
        // repro for the dominant pattern seen in 33/84 corpus regressions
        // K2b's demotion fix surfaced (fused quinazoline/quinoxaline/
        // purine-shaped bicyclics with an N-substituent elsewhere in the
        // molecule); it is not itself one of the 84 (it is unsubstituted).
        let mol = mol_kekulized("c1cnc2ccccc2n1");
        let model = assign_aromaticity_authoritative_experimental(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "all 10 atoms should be aromatic under the opt-in engine, matching RDKit"
        );
        assert!(
            mol.atoms().all(|(idx, _)| model.is_atom_aromatic(idx)),
            "every atom should be aromatic"
        );
    }

    #[test]
    fn test_default_engine_unaffected_by_fused_diazine_fix() {
        // Same molecule as above, through the DEFAULT engine
        // (`assign_aromaticity`) -- must stay exactly as it was before the
        // K2b fused-diazine follow-up fix existed (6/10, still wrong vs
        // RDKit), confirming `apply_aromaticity`/`apply_aromaticity_ex`
        // remain byte-identical to pre-K2b behavior per the coordinator
        // decision to ship this as opt-in only.
        let mol = mol_kekulized("c1cnc2ccccc2n1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "default engine must stay unaffected: only the pyridine-type ring \
             (6/10 atoms) confirmed, matching pre-K2b behavior"
        );
    }

    /// A handful of the 33-molecule `fused_diazine_quinazoline_quinoxaline_purine`
    /// corpus cluster (K2b's own diagnosis; see the PR description), pinned
    /// as permanent regression tests against the OPT-IN engine now that this
    /// fix resolves them there. Not exhaustive -- the fixed corpus-vs-RDKit
    /// comparison (`scripts/aromaticity_atom_parity.py` equivalent run
    /// against `scripts/descriptor_census_corpus.smi`) is the authoritative
    /// check; these are a stable, minimal sample.
    #[test]
    fn test_authoritative_experimental_fused_diazine_cluster_sample_matches_rdkit() {
        // (smiles, expected RDKit-aromatic atom count, all-aromatic?)
        let cases: &[(&str, usize)] = &[
            ("COc1cccc2nc(N3CCNCC3)cnc12", 10),
            ("Fc1cccc2nc(N3CCNCC3)cnc12", 10),
            ("Clc1cccc2nc(N3CCNCC3)cnc12", 10),
            ("CN1CCN(c2cnc3cc(Cl)ccc3n2)CC1", 10),
            ("Clc1cc2ncc(N3CCNCC3)nc2cc1Cl", 10),
            ("O=C(O)C1CN(c2cnc3ccccc3n2)CCN1", 10),
        ];
        for (smi, expected) in cases {
            let mol = mol_kekulized(smi);
            let model = assign_aromaticity_authoritative_experimental(&mol);
            assert_eq!(
                model.aromatic_atom_count(),
                *expected,
                "{smi}: expected {expected} aromatic atoms (the fused \
                 quinoxaline/naphthyridine core) under the opt-in engine, matching RDKit"
            );
        }
    }

    #[test]
    fn test_known_gap_azulene_nonalternant_odd_odd_split() {
        // c1ccc2cccc-2cc1 -- azulene itself (already the canonical example
        // in this codebase and in docs/rfcs/aromaticity_a1_rfc.md). RDKit: fully
        // aromatic, all 10 atoms, 9/10 bonds (the explicit fusion bond the
        // SMILES itself writes non-aromatic, `-2`, stays a formal single
        // bond even in RDKit's own answer). chematic: 0/10 -- both the
        // 5-ring and 7-ring independently get an ODD pi count (5 and 7) in
        // Pass 1, so neither is Aromatic nor Antiaromatic (both
        // `NonAromatic`), and Pass 2 never seeds because seeding requires
        // an ALREADY-aromatic adjacent ring, which neither ring is able to
        // become on its own. The default model now applies a deliberately
        // narrow all-carbon odd/odd fused-envelope fallback for this case.
        let mol = mol_kekulized("c1ccc2cccc-2cc1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "default Hückel's bounded fused-envelope fallback recognizes azulene"
        );
    }

    // ── N-substituted pyrrole-type N: bridgehead-branch guard removal ────────
    //
    // The bridgehead-N branch used to require the exocyclic substituent to be
    // sp2, to defensively block imide N (phthalimide). That guard also
    // blocked the much more common case of a plain alkyl/aryl/sugar
    // substituent on an otherwise-aromatic pyrrole-type N. It was removed;
    // these tests cover both the newly-fixed cases and the phthalimide
    // regression it was guarding against (which stays correct via the
    // overall 4n+2 sum, not the substituent).

    #[test]
    fn test_n_methylpyrrole_aromatic() {
        let mol = mol_kekulized("CN1C=CC=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            5,
            "all 5 N-methylpyrrole ring atoms aromatic"
        );
    }

    #[test]
    fn test_n_methylimidazole_aromatic() {
        let mol = mol_kekulized("CN1C=CN=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            5,
            "all 5 N-methylimidazole ring atoms aromatic"
        );
    }

    #[test]
    fn test_n_methylindole_aromatic() {
        let mol = mol_kekulized("CN1C=CC2=CC=CC=C21");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 N-methylindole ring atoms aromatic"
        );
    }

    #[test]
    fn test_9_methylpurine_aromatic() {
        let mol = mol_kekulized("CN1C=NC2=NC=NC=C21");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            9,
            "all 9 9-methylpurine ring atoms aromatic"
        );
    }

    #[test]
    fn test_phthalimide_5ring_not_aromatic() {
        // O=C1NC(=O)c2ccccc21 — only the fused benzo ring is aromatic (6
        // atoms); the imide 5-ring (2 carbonyl C + N) is not: carbonyl
        // carbons contribute 0π each (exocyclic C=O rule), N contributes 2π,
        // the two ring-fusion carbons contribute 1π each — 4π total, not
        // 4n+2. Regression guard for the bridgehead-N guard removal above.
        let mol = mol_kekulized("O=C1NC(=O)c2ccccc21");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "only the 6 benzo atoms of phthalimide are aromatic"
        );
    }

    #[test]
    fn test_n_methylphthalimide_5ring_not_aromatic() {
        // O=C1N(C)C(=O)c2ccccc21 — same as phthalimide but N-methylated;
        // same accounting applies (N still contributes 2π regardless of
        // substituent), 5-ring still non-aromatic.
        let mol = mol_kekulized("O=C1N(C)C(=O)c2ccccc21");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            6,
            "only the 6 benzo atoms of N-methylphthalimide are aromatic"
        );
    }

    #[test]
    fn test_azulene_kekulized_aromatic() {
        // C1=CC2=CC=CC=CC2=C1 — non-alternant fused bicyclic, all 10 atoms
        // aromatic per RDKit. Regression coverage: this was previously
        // (incorrectly) believed to need a ring-system rewrite, based on a
        // test that never called apply_aromaticity() on Kekulized input.
        //
        // Regressed by the Horton SSSR rewrite (confirmed passing on the old
        // single-spanning-tree find_sssr, failing only after Horton). Root
        // cause, empirically confirmed via debug dump: Horton's correct,
        // minimal SSSR is exactly the 5-ring + 7-ring (matches RDKit). Each
        // evaluated standalone has an ODD pi-electron count (5-ring: 5pi,
        // 7-ring: 7pi — every ring atom contributes 1pi via a double bond,
        // whether the double bond is endo- or exocyclic-to-a-carbon), so
        // neither passes Pass 1 and neither can seed Pass 2's
        // aromatic_context bootstrap. Azulene's aromaticity is a genuinely
        // non-alternant, whole-perimeter (10-atom, 10pi) delocalized system
        // — it needs the full-ring-system envelope as a Hückel candidate,
        // which `augmented_ring_set` deliberately excludes (its docstring
        // names naphthalene's spurious 10-ring as the exact case to avoid).
        // The old, non-minimal SSSR happened to hand a large fundamental
        // cycle straight to Pass 1 that included the whole perimeter,
        // papering over this gap by coincidence.
        //
        // The bounded all-carbon odd/odd fused-envelope fallback now handles
        // this case without changing the broader default model.
        let mol = mol_kekulized("C1=CC2=CC=CC=CC2=C1");
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            10,
            "all 10 azulene atoms aromatic"
        );
    }

    // ── RDKit #9271: charged / zwitterionic aromatic systems ─────────────────

    #[test]
    fn test_fluorescein_dianion_aromatic() {
        // Fluorescein dianion: RDKit #9271 incorrectly marked xanthene bonds as
        // single instead of aromatic. Verify chematic parses and identifies
        // aromatic atoms correctly (two benzene rings + xanthene O-bridge ring).
        // Kekulé-form SMILES: all atoms uppercase.
        let smi = "C1=CC=C(C(=C1)C2=C3C=CC(=O)C=C3OC4=C2C=CC(=C4)[O-])C(=O)[O-]";
        let mol = chematic_smiles::parse(smi).expect("fluorescein dianion should parse");
        // The molecule should parse without panic. Verify aromatic ring count:
        // fluorescein has 3 aromatic rings (2 benzene + xanthene core).
        let arc = count_aromatic_rings(&mol);
        assert!(
            arc >= 2,
            "fluorescein dianion: expected ≥2 aromatic rings, got {arc} \
             (RDKit #9271: charged aromatics may be misclassified)"
        );
    }

    #[test]
    fn test_rhodamine_zwitterion_parses() {
        // Rhodamine-type zwitterion with N+ and bridging O (RDKit #9271).
        // Must parse cleanly and produce a valid aromatic ring count.
        let smi = "CCN(CC)c1ccc2c(-c3ccccc3C(=O)O)c3ccc(=[N+](CC)CC)cc-3oc2c1";
        let mol = chematic_smiles::parse(smi).expect("rhodamine zwitterion should parse");
        let arc = count_aromatic_rings(&mol);
        assert!(arc >= 3, "rhodamine: expected ≥3 aromatic rings, got {arc}");
    }

    #[test]
    fn test_cyclopentadienyl_not_aromatic_kekulized() {
        // C1=CC=CC1 — cyclopentadiene (4 C with doubles + 1 sp3 CH2): not aromatic.
        let mut b = MoleculeBuilder::new();
        let c0 = b.add_atom(Atom::new(Element::C)); // sp3
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let c3 = b.add_atom(Atom::new(Element::C));
        let c4 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c0, c1, BondOrder::Single).unwrap();
        b.add_bond(c1, c2, BondOrder::Double).unwrap();
        b.add_bond(c2, c3, BondOrder::Single).unwrap();
        b.add_bond(c3, c4, BondOrder::Double).unwrap();
        b.add_bond(c4, c0, BondOrder::Single).unwrap();
        let mol = b.build();
        let model = assign_aromaticity(&mol);
        assert_eq!(
            model.aromatic_atom_count(),
            0,
            "cyclopentadiene not aromatic"
        );
    }

    // =========================================================================
    // RdkitLike mode: P/Se/Te heteroaromatics
    // =========================================================================

    #[test]
    fn test_phosphole_rdkit_aromatic() {
        // c1cc[pH]c1 — P donates its lone pair in the RDKit-compatible mode.
        let mol = mol_aromatic("c1cc[pH]c1");
        let m = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            m.aromatic_atom_count(),
            5,
            "phosphole: all 5 atoms aromatic in RdkitLike"
        );
    }

    #[test]
    fn test_azulene_rdkit_like_uses_whole_perimeter() {
        // The strict per-ring Hückel pass sees azulene as an odd/odd fused
        // split. RDKit evaluates the connected 10π perimeter instead.
        let mol = mol_kekulized("C1=CC2=CC=CC=CC2=C1");
        let m = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            m.aromatic_atom_count(),
            10,
            "azulene: whole perimeter must be aromatic in RdkitLike"
        );
    }

    #[test]
    fn test_selenophene_huckel_not_aromatic() {
        // c1cc[se]c1 — in strict Hückel mode, Se is unsupported → 0 aromatic atoms
        // (assign_aromaticity_ex re-derives from scratch, ignoring parser's aromatic flags)
        let mol = mol_aromatic("c1cc[se]c1");
        let m = assign_aromaticity(&mol); // default Hückel
        assert_eq!(
            m.aromatic_atom_count(),
            0,
            "selenophene: Se not aromatic in Hückel mode"
        );
    }

    #[test]
    fn test_selenophene_rdkit_aromatic() {
        // c1cc[se]c1 — in RdkitLike mode, Se donates 2π → 6π total → aromatic
        let mol = mol_aromatic("c1cc[se]c1");
        let m = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            m.aromatic_atom_count(),
            5,
            "selenophene: all 5 atoms aromatic in RdkitLike"
        );
    }

    #[test]
    fn test_tellurophene_rdkit_aromatic() {
        // c1cc[te]c1 — Te analogous to Se (2π donor)
        let mol = mol_aromatic("c1cc[te]c1");
        let m = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            m.aromatic_atom_count(),
            5,
            "tellurophene: all 5 atoms aromatic in RdkitLike"
        );
    }

    #[test]
    fn test_benzoselenophene_rdkit() {
        // Fused benzene + selenophene
        let mol = mol_aromatic("c1ccc2[se]ccc2c1");
        let m = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            m.aromatic_atom_count(),
            9,
            "benzoselenophene: 9 atoms aromatic"
        );
    }

    #[test]
    fn test_rdkit_mode_does_not_break_benzene() {
        // Benzene must give same result in both modes
        let mol = mol_aromatic("c1ccccc1");
        let m_h = assign_aromaticity(&mol);
        let m_r = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(m_h.aromatic_atom_count(), m_r.aromatic_atom_count());
    }

    #[test]
    fn test_rdkit_mode_does_not_break_thiophene() {
        let mol = mol_aromatic("c1ccsc1");
        let m_h = assign_aromaticity(&mol);
        let m_r = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            m_h.aromatic_atom_count(),
            m_r.aromatic_atom_count(),
            "thiophene same in both modes"
        );
    }

    // ── Known regressions from fix #2 (bridgehead-N guard removal) ──────────
    //
    // Re-measured after the Horton SSSR rewrite landed (find_sssr is now
    // minimal and deterministic, 0% self-instability on the 5000-molecule
    // corpus): all 32 counts below are UNCHANGED under the DEFAULT engine.
    // Zero free recoveries there.
    //
    // These 32 molecules share one root cause: a "fake bridgehead" N (same
    // local shape as a genuine bridgehead or N-substituted azole) feeds a
    // central ring that only closes via the `aromatic_context` bypass reusing
    // an unrelated ring's atoms. Fixing this requires removing the bypass in
    // favor of proper ring-system candidate enumeration (see project plan/
    // issue tracker).
    //
    // RESOLVED, but only under the OPT-IN `assign_aromaticity_authoritative_experimental`
    // engine (K2b fused-diazine follow-up fix; see
    // `test_authoritative_experimental_fixes_bridgehead_n_false_positives`
    // below): all 32 of these benzo-fused bridgehead-N tricyclics
    // (`...C3=NCCCN23`-shaped) ALSO have a fusion carbon whose own Kekule
    // double bond points into the adjacent ring at a heteroatom -- the exact
    // same misclassification the fused-diazine fix targets, just in a
    // three-ring rather than two-ring shape. Spot-checked live against
    // rdkit==2026.03.3 for 4 of the 32 (the shortest, a 15/12 case, and two
    // of the 28/24 cases), all matching. The originally-suspected root cause
    // above (the `aromatic_context`/`AlreadyAromaticContext` bypass) was
    // evidently either wrong or not the operative mechanism for this
    // specific molecule class -- not re-investigated further, since the fix
    // that resolved it was general (scoped to the fused-diazine cluster) and
    // not bridgehead-N-specific. This is opt-in only: the DEFAULT engine
    // (`assign_aromaticity`) is unaffected and still shows the original
    // `expected_wrong` counts below (see the coordinator decision requiring
    // `apply_aromaticity`/`apply_aromaticity_ex` to stay byte-identical to
    // pre-K2b behavior).
    // (kekulized SMILES, current chematic aromatic_atom_count() under the
    // default engine, RDKit's correct count).
    // Named at module level (not a local in the test below) so
    // Aromaticity-A1-0's corpus tests, further down this module, can reuse
    // the identical pinned data instead of re-deriving a copy that could
    // silently drift out of sync with it.
    const KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES: &[(&str, usize, usize)] = &[
        ("C[Si](C)(C)C1=CC=C(C2=CC3=CC=CC=C3C3=NCCCN23)C=C1", 16, 12),
        (
            "C1=C(C2=CC=C(CCC3=CC=CC=C3)C=C2)N2CCCN=C2C2=CC=CC=C12",
            22,
            18,
        ),
        ("ClC1=CC=C(OCC2=CC3=CC=CC=C3C3=NCCCN23)C=C1", 16, 12),
        ("N[C@@H](CC1=CC=CC=C1)C1=CC2=CC=CC=C2C2=NCCCN12", 16, 12),
        (
            "CC(C)(C)C1=CC=C(C2=C(CC3=CC=CC=C3)C3=CC=CC=C3C3=NCCCN32)C=C1",
            22,
            18,
        ),
        (
            "C[Si](C)(C)C1=CC=C(C2=C(CC3=CC=CC=C3)C3=CC=CC=C3C3=NCCCN32)C=C1",
            22,
            18,
        ),
        (
            "C1=C(C2=CC=C(C3=CC=CC=C3)C=C2)N2CCCN=C2C2=CC=CC=C12",
            22,
            18,
        ),
        (
            "C1=C(C2=CC=C(OCC3=CC=CC=C3)C=C2)N2CCCN=C2C2=CC=CC=C12",
            22,
            18,
        ),
        ("COC1=C(OC)C(OC)=CC(C2=CC3=CC=CC=C3C3=NCCCN23)=C1", 16, 12),
        ("CC1=CC2=CC=CC=C2C2=NCCCN12", 10, 6),
        (
            "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(NC(=O)NC4CCCCC4)C=C3)C3=NCCCN23)C=C1",
            16,
            12,
        ),
        (
            "C1=CC=C(CCC2=CC=C(C3=C(CC4=CC=CC=C4)C4=CC=CC=C4C4=NCCCN43)C=C2)C=C1",
            28,
            24,
        ),
        (
            "CCCCC1=C(C2=CC=C(CCC3=CC=CC=C3)C=C2)N2CCCN=C2C2=CC=CC=C12",
            22,
            18,
        ),
        (
            "CCCCC1=C(C2=CC=C(C(C)(C)C)C=C2)N2CCCN=C2C2=CC=CC=C12",
            16,
            12,
        ),
        ("CCCCCCC1=CC2=CC=CC=C2C2=NCCCN12", 10, 6),
        (
            "CCOC1=CC=C(CC2=C(CCCC3=CC=CC4=CC=CC=C34)N3CCCN=C3C3=CC=CC=C23)C=C1",
            26,
            22,
        ),
        (
            "CCOC1=CC=C(CC2=C(C3=CC=C(CCC4=CC=CC=C4)C=C3)N3CCCN=C3C3=CC=CC=C23)C=C1",
            28,
            24,
        ),
        (
            "CN(C)CCC1=C(C2=CC=C(C(C)(C)C)C=C2)N2CCCN=C2C2=CC=CC=C12",
            16,
            12,
        ),
        (
            "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(N/C(S)=N/C4CCCCC4)C=C3)C3=NCCCN23)C=C1",
            16,
            12,
        ),
        ("C1=C(/C=C/C2=CC=CC=C2)N2CCCN=C2C2=CC=CC=C12", 16, 12),
        ("CC(C)(C)C1=CC=C(C2=CC3=CC=CC=C3C3=NCCCN23)C=C1", 16, 12),
        (
            "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(NC(=O)CC4=CC=CC=N4)C=C3)C3=NCCCN23)C=C1",
            22,
            18,
        ),
        (
            "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(NC(=O)NC4=C(Cl)C=C(Cl)C=C4)C=C3)C3=NCCCN23)C=C1",
            22,
            18,
        ),
        ("C1=C(CC2=CC=CC=C2)C2=CC=CC=C2C2=NCCCN12", 16, 12),
        ("ClC1=CC=C(C2=CC3=CC=CC=C3C3=NCCCN23)C=C1", 16, 12),
        ("C1=C(C2=CC=CC=C2)N2CCCN=C2C2=CC=CC=C12", 16, 12),
        (
            "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(N(CC4=CC=CC=C4)CC4=CC=CC=C4)C=C3)C3=NCCCN23)C=C1",
            28,
            24,
        ),
        (
            "CC(C)(C)C1=CC=C(C2=CC3=C(C=C(N)C=C3)C3=NCCCN23)C=C1",
            16,
            12,
        ),
        ("CC1=C2C(=NC=C1)N(C1CC1)C1=NC=CC=C1C(=O)N2C", 15, 12),
        ("CC(=O)N1C2=NC=CC=C2C(=O)N(C)C2=CC=CN=C21", 15, 12),
        ("CN1C(=O)C2=CC=CN=C2N(C(C)(C)C)C2=NC=CC=C21", 15, 12),
        ("CCCN1C2=NC=CC=C2C(=O)N(C)C2=CC=CN=C21", 15, 12),
    ];

    #[test]
    fn test_known_regressions_from_bridgehead_n_fix() {
        for (smi, expected_wrong, rdkit_correct) in KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES {
            let mol = mol_kekulized(smi);
            let model = assign_aromaticity(&mol);
            assert_eq!(
                model.aromatic_atom_count(),
                *expected_wrong,
                "{smi}: expected current (wrong) count {expected_wrong} under the default \
                 engine (RDKit correct: {rdkit_correct})"
            );
        }
    }

    #[test]
    fn test_authoritative_experimental_fixes_bridgehead_n_false_positives() {
        // Beneficial, unattempted side effect of the K2b fused-diazine
        // follow-up fix, now reachable only via the opt-in engine -- see
        // this const's preceding doc comment.
        for (smi, _expected_wrong, rdkit_correct) in KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES {
            let mol = mol_kekulized(smi);
            let model = assign_aromaticity_authoritative_experimental(&mol);
            assert_eq!(
                model.aromatic_atom_count(),
                *rdkit_correct,
                "{smi}: expected {rdkit_correct} aromatic atoms under the opt-in \
                 authoritative-experimental engine, matching RDKit"
            );
        }
    }

    // ── Known order-dependence: same molecule, different Kekulized traversal ─
    //
    // Originally found because these 3 molecules passed with RDKit's
    // canonical Kekulized SMILES but failed with at least one other valid
    // Kekulized ordering of the identical structure -- confirmed via
    // atom-map-number alignment (no substructure matching). Root cause was
    // NOT Pass 1/Pass 2 (verified order-invariant by construction) -- it was
    // `find_sssr` itself, non-deterministic and non-minimal.
    //
    // Re-measured after the Horton SSSR rewrite (find_sssr is now
    // deterministic and minimal, 0% self-instability on the 5000-molecule
    // corpus): the 3 pinned failing-traversal counts below are UNCHANGED.
    // The original order-dependence *mechanism* (find_sssr picking a
    // different non-minimal ring depending on traversal) is resolved -- but
    // these 3 specific SMILES still disagree with RDKit's count, so at least
    // one more bug (likely `aromatic_context`, same as the 32-molecule
    // corpus above) also affects this molecule class. Not re-diagnosed here;
    // a fresh worst-of-N run against the full corpus would confirm whether
    // order-dependence itself (canonical vs. this pinned variant disagreeing
    // with each other) is now fully gone, separate from RDKit agreement.
    //
    // The K2b fused-diazine follow-up fix (`assign_aromaticity_authoritative_experimental`)
    // does shift 2 of these 3 counts (16->12) when run through the OPT-IN
    // engine -- confirmed unrelated to and not fixing this bucket (still
    // wrong, by a different amount, a separate multi-causal bug). This test
    // asserts the DEFAULT (`assign_aromaticity`) engine only, which is
    // unaffected by that opt-in fix, so the pinned values below stay as
    // originally measured.
    // Named at module level for the same reason as
    // `KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES` above -- Aromaticity-A1-0's corpus
    // tests reuse this exact pinned data instead of a second copy.
    const KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES: &[(&str, usize, usize)] = &[
        (
            "N1=C2C(N(CC(O)=O)C(=O)N=C2N(C2C=C(C(F)(F)F)C=C(C=2)C(F)(F)F)C2C1=CC=CC=2)=O",
            16,
            20,
        ),
        (
            "[C@H]12N(C([C@H](NC(=O)[C@H]([C@H](OC(=O)[C@@H](N(C)C(CN(C)C1=O)=O)C(C)C)C)NC(=O)C1C=C(OC)C(C)=C3OC4=C(C)C(=O)C(=C(C4=NC=13)C(=O)N[C@H]1C(=O)N[C@@H](C(C)C)C(N3[C@H](C(=O)N(CC(N([C@H](C(C)C)C(O[C@H]1C)=O)C)=O)C)CCC3)=O)N)C(C)C)=O)CCC2",
            6,
            14,
        ),
        ("C12N(C3C=CC=CC=3)C3=NC(=O)N(C)C(C3=NC1=CC=CC=2)=O", 16, 20),
    ];

    #[test]
    fn test_known_order_dependent_regressions() {
        for (smi, expected_wrong, rdkit_correct) in KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES {
            let mol = mol_kekulized(smi);
            let model = assign_aromaticity(&mol);
            assert_eq!(
                model.aromatic_atom_count(),
                *expected_wrong,
                "{smi}: expected current (wrong) count {expected_wrong} (RDKit correct: {rdkit_correct})"
            );
        }
    }

    // ── Aromaticity-A1-0: anti-drift guard for `trace_ring_pi_electrons` ────
    //
    // `trace_ring_pi_electrons` is a deliberately separate implementation
    // from `ring_pi_electrons` (see the doc comment above it) so it can
    // report *why* each atom scored what it did. That separateness is a
    // drift risk: nothing stops the two from silently diverging as either
    // one is edited. This test is the guard -- for every ring in every
    // molecule of the known false-positive/false-negative/negative-control
    // corpus (the same molecules `docs/rfcs/aromaticity_a1_rfc.md`'s diagnostic
    // corpus uses), both functions must agree exactly, in both an empty
    // context (Pass-1-equivalent) and the model's final converged context
    // (an upper-bound Pass-2-equivalent). This does not assert anything
    // about correctness vs RDKit -- only that the trace and the real engine
    // never disagree with each other.
    #[test]
    fn trace_matches_ring_pi_electrons_on_corpus() {
        let smiles: Vec<&str> = KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES
            .iter()
            .map(|(smi, _, _)| *smi)
            .chain(
                KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES
                    .iter()
                    .map(|(smi, _, _)| *smi),
            )
            .chain([
                "C1=CC2=CC=CC=CC2=C1",    // azulene (Kekulized) -- known false negative
                "c1cnc2[nH]cnc2n1",       // purine -- known false negative
                "C1=Cc2ccccc2C2=NCCCN12", // PR #86 minimal false-positive reproducer
                "C1=Cc2ccccc2C2=CCCC12",  // negative control: no bridgehead N
                "C1=Cc2ccccc2C2=CCNC12",  // negative control: N not at bridgehead
                "C1Cc2ccccc2C2=NCCCN12",  // negative control: bridgehead N, no exocyclic C=C
                "c1ccc2[nH]ccc2c1",       // indole -- must stay correct
                "c1ccc2ncccc2c1",         // quinoline -- must stay correct
                "c1ccc2ccccc2c1",         // naphthalene -- must stay correct
            ])
            .collect();

        for algo in [
            AromaticityAlgorithm::Huckel,
            AromaticityAlgorithm::RdkitLike,
        ] {
            for smi in &smiles {
                let mol = mol_kekulized(smi);
                let model = assign_aromaticity_ex(&mol, algo);
                let final_context: FxHashSet<AtomIdx> = mol
                    .atoms()
                    .map(|(idx, _)| idx)
                    .filter(|&idx| model.is_atom_aromatic(idx))
                    .collect();

                let sssr = find_sssr(&mol);
                let rings = augmented_ring_set(&mol, sssr.rings());
                let empty_context: FxHashSet<AtomIdx> = FxHashSet::default();
                let all_ring_bonds: FxHashSet<BondIdx> =
                    rings.iter().flat_map(|r| ring_bond_set(&mol, r)).collect();

                for ring in &rings {
                    for ctx in [&empty_context, &final_context] {
                        let expected = ring_pi_electrons(&mol, ring, ctx, algo, &all_ring_bonds);
                        let traced =
                            trace_ring_pi_electrons(&mol, ring, ctx, algo, &all_ring_bonds);
                        assert_eq!(
                            traced.total,
                            expected,
                            "{smi} (algo={algo:?}, ring={ring:?}, ctx_len={}): \
                             trace_ring_pi_electrons diverged from ring_pi_electrons",
                            ctx.len()
                        );
                        // Cross-check the per-atom eligibility bookkeeping too.
                        for a in &traced.atoms {
                            assert_eq!(
                                a.contribution.is_some(),
                                a.reason.is_eligible(),
                                "{smi}: atom {:?} contribution/reason eligibility mismatch",
                                a.atom_idx
                            );
                        }
                    }
                }
            }
        }
    }

    // ── Aromaticity-A1-0: false-positive/false-negative polarity sanity ────
    //
    // These are cheap, structural sanity checks that the corpus buckets are
    // labeled the direction they claim -- not a re-measurement of the full
    // corpus (that's `aromaticity_a1_0_report` + the Python RDKit join, see
    // `docs/rfcs/aromaticity_a1_rfc.md`). Catches an accidental swap or a stale
    // pinned count silently going the other way.
    #[test]
    fn false_positive_corpus_over_counts_vs_rdkit() {
        for (smi, expected_wrong, rdkit_correct) in KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES {
            assert!(
                expected_wrong > rdkit_correct,
                "{smi}: false-positive bucket entry should over-count \
                 (chematic={expected_wrong} should be > rdkit={rdkit_correct})"
            );
        }
    }

    #[test]
    fn false_negative_corpus_under_counts_vs_rdkit() {
        for (smi, expected_wrong, rdkit_correct) in KNOWN_ORDER_DEPENDENT_FALSE_NEGATIVES {
            assert!(
                expected_wrong < rdkit_correct,
                "{smi}: false-negative bucket entry should under-count \
                 (chematic={expected_wrong} should be < rdkit={rdkit_correct})"
            );
        }
    }

    // ── Aromaticity-A1-1a: exhaustive_aromaticity_oracle pinned cases ──────
    //
    // The oracle is a discovery tool, not a correct-answer generator: its
    // candidates are built from the SAME per-atom local rules
    // (`evaluate_atom_pi_contribution`) that are wrong for the false-positive
    // family, so it can't independently arbitrate that family. This test
    // pins what the oracle DOES get right (RDKit-atom-index-verified, not
    // guessed) after two real fixes made during this milestone:
    //
    // 1. Connectivity: `build_conjugated_components`'s conjugation graph
    //    originally only bridged single bonds via a `LonePairDonor` endpoint,
    //    leaving azulene's all-carbon alternating perimeter as 5 disconnected
    //    2-atom pairs (oracle returned an empty set). Fixed: any bond between
    //    two independently-eligible atoms connects (ordinary carbon-carbon
    //    single-bond conjugation, ordinary organic chemistry).
    // 2. Home-ring evaluation: evaluating a multi-ring candidate's electron
    //    sum against its own *flattened* atom set broke the N
    //    bridgehead/substituted-azole rule for any TRUE bridgehead (every
    //    bond looks "in-family" once the family itself is the context) --
    //    indolizine's own bridgehead N came out `Ineligible`, an oracle bug,
    //    not a chematic bug. Fixed via `evaluate_atom_via_home_ring`.
    //
    // Both fixes were originally confirmed correct AND confirmed NOT to
    // silently "fix" the false-positive family by accident.
    //
    // UPDATE (K2b fused-diazine follow-up fix): both the false-positive
    // reproducer AND purine are now RDKit-exact too, as a side effect of the
    // same general `CarbonExocyclicHeteroatomDouble` ring-fusion fix
    // (`evaluate_atom_pi_contribution_inner` mirrors `ring_pi_electrons`'s
    // rule exactly -- see its doc comment). The false-positive reproducer's
    // own fusion carbon (whose double bond points into the bridgehead-N
    // ring's own nitrogen) no longer gets wrongly zeroed, so the oracle
    // stops over-aromatizing into the bridgehead ring and correctly confirms
    // only the plain benzo ring. Purine's 5-ring fusion carbons no longer
    // get wrongly zeroed by the same rule either, so the oracle now confirms
    // all 9 atoms, matching RDKit -- resolving the open finding below.
    // Verified live against rdkit==2026.03.3 for both (not assumed from the
    // fix's general mechanism alone).
    #[test]
    fn exhaustive_oracle_pinned_cases() {
        let algo = AromaticityAlgorithm::RdkitLike;

        // (name, smiles, expected oracle-aromatic atom indices, sorted)
        let matches_rdkit: &[(&str, &str, &[u32])] = &[
            (
                "azulene",
                "C1=CC2=CC=CC=CC2=C1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
            (
                "naphthalene",
                "c1ccc2ccccc2c1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
            (
                "anthracene",
                "c1ccc2cc3ccccc3cc2c1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
            ),
            ("indole", "c1ccc2[nH]ccc2c1", &[0, 1, 2, 3, 4, 5, 6, 7, 8]),
            (
                "quinoline",
                "c1ccc2ncccc2c1",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            ),
            (
                "indolizine (bridgehead N, both rings valid)",
                "c1ccn2ccccc12",
                &[0, 1, 2, 3, 4, 5, 6, 7, 8],
            ),
            ("tropone", "O=c1cccccc1", &[1, 2, 3, 4, 5, 6, 7]),
            ("2-pyridone", "O=c1cccc[nH]1", &[1, 2, 3, 4, 5, 6]),
        ];
        for (name, smi, expected) in matches_rdkit {
            let mol = mol_kekulized(smi);
            let (atoms, _bonds) = exhaustive_aromaticity_oracle(&mol, algo);
            let mut got: Vec<u32> = atoms.iter().map(|a| a.0).collect();
            got.sort();
            assert_eq!(&got, expected, "{name} ({smi}): oracle should match RDKit");
        }

        // Now RDKit-exact -- see this test's doc comment (K2b fused-diazine
        // follow-up fix). RDKit: only the plain benzo ring (6 atoms) is
        // aromatic; the bridgehead-N ring is not (verified live).
        let (fp_atoms, _) =
            exhaustive_aromaticity_oracle(&mol_kekulized("C1=Cc2ccccc2C2=NCCCN12"), algo);
        let mut fp_got: Vec<u32> = fp_atoms.iter().map(|a| a.0).collect();
        fp_got.sort();
        assert_eq!(
            fp_got,
            vec![2, 3, 4, 5, 6, 7],
            "false-positive reproducer: oracle now matches RDKit exactly"
        );

        // Now RDKit-exact -- see this test's doc comment (K2b fused-diazine
        // follow-up fix). RDKit: all 9 atoms aromatic (verified live).
        let (purine_atoms, _) =
            exhaustive_aromaticity_oracle(&mol_kekulized("c1cnc2[nH]cnc2n1"), algo);
        let mut purine_got: Vec<u32> = purine_atoms.iter().map(|a| a.0).collect();
        purine_got.sort();
        assert_eq!(
            purine_got,
            vec![0, 1, 2, 3, 4, 5, 6, 7, 8],
            "purine: oracle now matches RDKit exactly (all 9 atoms)"
        );
    }
}

#[cfg(test)]
mod electron_trace_boundary_contract_tests {
    use super::*;

    #[test]
    fn charged_and_heavier_heterocycles_explain_their_pi_contribution() {
        use ContributionReason::*;
        let cases = [
            ("[CH+]1C=CC=C1", 0, Some(0), CarbonCationVacant),
            ("[CH-]1C=CC=C1", 0, Some(2), CarbonCarbanionLonePair),
            ("[N+]1(C)(C)CCCC1", 0, None, NitrogenIneligible),
            ("[OH+]1CCCC1", 0, None, ChalcogenIneligible),
            ("[O+]1=CC=CC=C1", 0, Some(1), ChalcogenCationPyridineType),
            ("[S+]1=CC=CC=C1", 0, Some(1), ChalcogenCationPyridineType),
            ("O=S1C=CC=C1", 1, None, ChalcogenIneligible),
            ("P1C=CC=C1", 0, Some(2), PnictogenOrChalcogenLonePair),
            ("[Se]1C=CC=C1", 0, Some(2), PnictogenOrChalcogenLonePair),
            ("[Te]1C=CC=C1", 0, Some(2), PnictogenOrChalcogenLonePair),
            ("O=[Se]1C=CC=C1", 1, None, ChalcogenIneligible),
            ("O=[Te]1C=CC=C1", 1, None, ChalcogenIneligible),
            ("[SiH2]1C=CC=C1", 0, None, UnsupportedElement),
        ];
        for (text, atom, electrons, reason) in cases {
            let mol = chematic_smiles::parse(text).unwrap_or_else(|e| panic!("{text}: {e}"));
            let rings = find_sssr(&mol);
            assert_eq!(rings.rings().len(), 1, "{text}");
            let ring = &rings.rings()[0];
            let ring_bonds = ring
                .iter()
                .zip(ring.iter().cycle().skip(1))
                .map(|(&a, &b)| mol.bond_between(a, b).unwrap().0)
                .collect();
            let trace = trace_ring_pi_electrons(
                &mol,
                ring,
                &FxHashSet::default(),
                AromaticityAlgorithm::RdkitLike,
                &ring_bonds,
            );
            let selected = trace
                .atoms
                .iter()
                .find(|a| a.atom_idx == AtomIdx(atom))
                .unwrap();
            assert_eq!(selected.contribution, electrons, "{text}");
            assert_eq!(selected.reason, reason, "{text}");
            assert_eq!(selected.reason.is_eligible(), electrons.is_some(), "{text}");
            assert_eq!(trace.atoms.len(), ring.len());
            let expected_total: Option<u32> = trace
                .atoms
                .iter()
                .map(|a| a.contribution.map(u32::from))
                .sum();
            assert_eq!(trace.total, expected_total, "{text}");
            if mol.atom(AtomIdx(atom)).element.atomic_number() == 15 {
                let strict = trace_ring_pi_electrons(
                    &mol,
                    ring,
                    &FxHashSet::default(),
                    AromaticityAlgorithm::Huckel,
                    &ring_bonds,
                );
                assert_eq!(strict.total, None, "the strict model excludes phosphorus");
            }
        }
    }
}
