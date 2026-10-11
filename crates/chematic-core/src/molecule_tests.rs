#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn molecule_keeps_auto_traits() {
        fn assert_traits<
            T: Send + Sync + std::panic::RefUnwindSafe + std::panic::UnwindSafe + Clone,
        >() {
        }
        assert_traits::<Molecule>();
    }
    use crate::atom::Atom;
    use crate::element::Element;

    #[test]
    fn fragments_report_source_atoms() {
        // O(0) . C(1)-C(2) . N(3) bonded to C(1) makes {1,2,3}; {0} alone.
        let mut b = MoleculeBuilder::new();
        let o = b.add_atom(Atom::new(Element::O));
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        let n = b.add_atom(Atom::new(Element::N));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_bond(n, c1, BondOrder::Single).unwrap();
        let mol = b.build();
        let frags = mol.fragments_with_source_atoms();
        assert_eq!(frags.len(), 2);
        assert_eq!(frags[0].1, vec![o]);
        assert_eq!(frags[1].1, vec![c1, c2, n]);
        for (frag, source) in &frags {
            for (i, &src) in source.iter().enumerate() {
                assert_eq!(frag.atom(AtomIdx(i as u32)).element, mol.atom(src).element);
            }
        }
        let plain = mol.fragments();
        assert_eq!(plain.len(), 2);
        assert_eq!(plain[1].atom_count(), 3);
        assert_eq!(plain[1].bond_count(), 2);
    }

    fn ethane() -> Molecule {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.build()
    }

    #[test]
    fn test_basic_molecule() {
        let mol = ethane();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
    }

    #[test]
    fn rdkit_bond_order_visitor_keeps_chain_then_sorted_closures() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4)
            .map(|_| builder.add_atom(Atom::new(Element::C)))
            .collect();
        let early_closure = builder
            .add_bond(atoms[0], atoms[1], BondOrder::Single)
            .unwrap();
        let chain = builder
            .add_bond(atoms[1], atoms[2], BondOrder::Single)
            .unwrap();
        let late_closure = builder
            .add_bond(atoms[2], atoms[3], BondOrder::Single)
            .unwrap();
        builder.set_smiles_ring_closure(early_closure, 2, 0);
        builder.set_smiles_ring_closure(late_closure, 9, 0);
        let molecule = builder.build();

        let mut visited = Vec::new();
        molecule.for_each_bond_in_rdkit_order(|bond| visited.push(bond));
        assert_eq!(visited, vec![chain, early_closure, late_closure]);
        assert_eq!(molecule.rdkit_bond_order(), visited);
    }

    #[test]
    fn test_adjacency() {
        let mol = ethane();
        let neighbors: Vec<_> = mol.neighbors(AtomIdx(0)).collect();
        assert_eq!(neighbors.len(), 1);
        assert_eq!(neighbors[0].0, AtomIdx(1));
    }

    #[test]
    fn test_bond_between() {
        let mol = ethane();
        assert!(mol.bond_between(AtomIdx(0), AtomIdx(1)).is_some());
        assert!(mol.bond_between(AtomIdx(1), AtomIdx(0)).is_some());
    }

    #[test]
    fn test_duplicate_bond_error() {
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let err = b.add_bond(c1, c2, BondOrder::Double);
        assert!(matches!(err, Err(MolError::DuplicateBond(_, _))));
    }

    #[test]
    fn test_formula() {
        let mut b = MoleculeBuilder::new();
        let c = b.add_atom(Atom::new(Element::C));
        let n = b.add_atom(Atom::new(Element::N));
        b.add_bond(c, n, BondOrder::Single).unwrap();
        let mol = b.build();
        assert_eq!(mol.formula(), "CN");
    }

    #[test]
    fn test_implicit_hydrogen_count() {
        // Isolated C atom (sp3, 4 bonds available): 4 implicit H
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::organic(Element::C));
        let mol = b.build();
        assert_eq!(mol.implicit_hydrogen_count(AtomIdx(0)), 4);
    }

    #[test]
    fn test_total_formula_methane() {
        // Organic C atom with 0 explicit bonds → 4 implicit H → CH4
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::organic(Element::C));
        let mol = b.build();
        assert_eq!(mol.total_formula(), "CH4");
    }

    #[test]
    fn test_total_formula_no_hydrogen() {
        // NaCl — neither Na nor Cl is in the organic subset, no implicit H
        let mut b = MoleculeBuilder::new();
        let na = b.add_atom(Atom::new(Element::NA));
        let cl = b.add_atom(Atom::new(Element::CL));
        b.add_bond(na, cl, BondOrder::Single).unwrap();
        let mol = b.build();
        assert_eq!(mol.total_formula(), "ClNa");
    }

    #[test]
    fn test_with_atom_aromatic() {
        let mol = ethane();
        let updated = mol.with_atom_aromatic(AtomIdx(0), true);
        assert!(updated.atom(AtomIdx(0)).aromatic);
        assert!(!updated.atom(AtomIdx(1)).aromatic);
    }

    #[test]
    fn test_with_bond_order() {
        let mol = ethane();
        let updated = mol.with_bond_order(BondIdx(0), BondOrder::Double);
        assert_eq!(updated.bond(BondIdx(0)).order, BondOrder::Double);
    }

    // --- bond_directions remap correctness (not just presence) ---

    /// 4-atom chain A-B-C-D with a `bond_direction` stash on the LAST bond
    /// (C-D, index 2). Removing the FIRST bond (A-B, index 0) shifts every
    /// surviving bond's index down by one; the stash must follow the C-D
    /// bond to its new index (1), not stay pinned to numeric index 2 (which
    /// would now point at a different physical bond) and not vanish.
    fn chain_with_direction_on_last_bond() -> (Molecule, BondIdx) {
        let mut b = MoleculeBuilder::new();
        let a = b.add_atom(Atom::new(Element::C));
        let bb = b.add_atom(Atom::new(Element::C));
        let c = b.add_atom(Atom::new(Element::C));
        let d = b.add_atom(Atom::new(Element::C));
        b.add_bond(a, bb, BondOrder::Single).unwrap(); // bond 0 (to be removed)
        b.add_bond(bb, c, BondOrder::Single).unwrap(); // bond 1
        let cd = b.add_bond(c, d, BondOrder::Single).unwrap(); // bond 2
        b.set_bond_direction(cd, BondOrder::Up);
        b.set_bond_direction_anchor(cd, c);
        (b.build(), cd)
    }

    #[test]
    fn test_remove_bond_remaps_bond_direction_not_misattributes() {
        let (mut mol, _cd) = chain_with_direction_on_last_bond();
        assert_eq!(mol.bond_count(), 3);
        mol.remove_bond(BondIdx(0)); // remove A-B; C-D shifts from index 2 to 1
        assert_eq!(mol.bond_count(), 2);
        // The stash must have followed C-D to its new index...
        assert_eq!(mol.bond_direction(BondIdx(1)), Some(BondOrder::Up));
        assert_eq!(mol.bond_direction_anchor(BondIdx(1)), Some(AtomIdx(2)));
        // ...and must NOT have leaked onto the bond that shifted into the
        // old numeric slot 2 (which no longer exists) or onto B-C (index 0
        // after the shift), which never had a direction.
        assert_eq!(mol.bond_direction(BondIdx(0)), None);
        assert_eq!(mol.bond_opt(BondIdx(2)), None);
    }

    #[test]
    fn test_remove_bond_drops_direction_for_the_removed_bond_itself() {
        let (mut mol, _cd) = chain_with_direction_on_last_bond();
        mol.remove_bond(BondIdx(2)); // remove C-D itself — its stash must go with it
        assert_eq!(mol.bond_count(), 2);
        assert!(mol.bond_direction(BondIdx(0)).is_none());
        assert!(mol.bond_direction(BondIdx(1)).is_none());
    }

    #[test]
    fn test_with_atom_removed_remaps_bond_direction() {
        let (mol, _cd) = chain_with_direction_on_last_bond();
        // Remove atom A (index 0), unrelated to the C-D bond carrying the
        // stash. Bonds incident to A (A-B) disappear; B-C and C-D survive,
        // renumbered 0 and 1 respectively — the direction must follow C-D.
        let (updated, _atom_remap) = mol.with_atom_removed(AtomIdx(0));
        assert_eq!(updated.bond_count(), 2);
        // Find the surviving C-D bond by scanning for the stash directly,
        // rather than assuming a specific bond index, so this test doesn't
        // depend on internal re-numbering order.
        let has_direction = (0..updated.bond_count())
            .map(|i| BondIdx(i as u32))
            .any(|bidx| updated.bond_direction(bidx) == Some(BondOrder::Up));
        assert!(
            has_direction,
            "bond_direction on C-D must survive atom removal, remapped to its new bond index"
        );
    }

    // --- mutable API ---

    #[test]
    fn test_add_remove_atom() {
        let mut mol = ethane();
        let n_idx = mol.add_atom(Atom::new(Element::N));
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.atom(n_idx).element.atomic_number(), 7);

        let remap = mol.remove_atom(n_idx);
        assert_eq!(mol.atom_count(), 2);
        assert!(remap[n_idx.0 as usize].is_none());
    }

    #[test]
    fn test_remove_atom_preserves_survivor_order() {
        let mut mol = MoleculeBuilder::new().build();
        mol.add_atom(Atom::new(Element::C));
        let removed = mol.add_atom(Atom::new(Element::N));
        mol.add_atom(Atom::new(Element::O));

        let remap = mol.remove_atom(removed);

        assert_eq!(
            mol.atoms()
                .map(|(_, atom)| atom.element)
                .collect::<Vec<_>>(),
            vec![Element::C, Element::O]
        );
        assert_eq!(remap, vec![Some(AtomIdx(0)), None, Some(AtomIdx(1))]);
    }

    #[test]
    fn test_add_remove_bond() {
        let mut mol = ethane();
        let n_idx = mol.add_atom(Atom::new(Element::N));
        let bidx = mol.add_bond(AtomIdx(0), n_idx, BondOrder::Single).unwrap();
        assert_eq!(mol.bond_count(), 2);
        mol.remove_bond(bidx);
        assert_eq!(mol.bond_count(), 1);
    }

    #[test]
    fn test_set_charge_element() {
        let mut mol = ethane();
        mol.set_charge(AtomIdx(0), 1);
        assert_eq!(mol.atom(AtomIdx(0)).charge, 1);
        mol.set_element(AtomIdx(0), Element::N);
        assert_eq!(mol.atom(AtomIdx(0)).element.atomic_number(), 7);
    }

    #[test]
    fn test_is_connected() {
        let mol = ethane();
        assert!(mol.is_connected());

        // Two separate atoms — disconnected
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::new(Element::C));
        b.add_atom(Atom::new(Element::N));
        let disconnected = b.build();
        assert!(!disconnected.is_connected());
    }

    #[test]
    fn test_fragments() {
        // "CC.N" — two components
        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::organic(Element::C));
        let c2 = b.add_atom(Atom::organic(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        b.add_atom(Atom::new(Element::N)); // disconnected N
        let mol = b.build();
        let frags = mol.fragments();
        assert_eq!(frags.len(), 2);
        let sizes: std::collections::HashSet<usize> =
            frags.iter().map(|f| f.atom_count()).collect();
        assert!(sizes.contains(&2));
        assert!(sizes.contains(&1));
    }

    #[test]
    fn test_builder_from_molecule() {
        let mol = ethane();
        let mut b = MoleculeBuilder::from_molecule(&mol);
        b.add_atom(Atom::new(Element::O));
        let mol2 = b.build();
        assert_eq!(mol2.atom_count(), 3);
        assert_eq!(mol2.bond_count(), 1); // original bond preserved
    }

    #[test]
    fn r_group_sidecars_survive_index_preserving_and_remapping_edits() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(Atom::new(Element::C));
        let r1 = builder.add_atom(Atom::wildcard());
        builder.set_r_group(r1, RGroupLabel::numbered(1).unwrap());
        let r2 = builder.add_atom(Atom::wildcard());
        builder.set_r_group(r2, RGroupLabel::numbered(2).unwrap());
        builder.add_bond(carbon, r1, BondOrder::Single).unwrap();
        let mol = builder.build();

        let charged = mol.with_atom_charge(carbon, 1);
        assert_eq!(charged.r_group_label(r1).unwrap().number(), Some(1));
        assert_eq!(charged.r_group_label(r2).unwrap().number(), Some(2));

        let changed = charged.with_atom_element(r1, Element::N);
        assert_eq!(changed.r_group_label(r1), None);
        assert!(!changed.atom(r1).wildcard);
        assert_eq!(changed.r_group_label(r2).unwrap().number(), Some(2));

        let (removed, remap) = mol.with_atom_removed(carbon);
        let remapped_r1 = remap[r1.0 as usize].unwrap();
        let remapped_r2 = remap[r2.0 as usize].unwrap();
        assert_eq!(
            removed.r_group_label(remapped_r1).unwrap().number(),
            Some(1)
        );
        assert_eq!(
            removed.r_group_label(remapped_r2).unwrap().number(),
            Some(2)
        );

        let fragments = mol.fragments();
        assert_eq!(fragments.len(), 2);
        assert!(fragments.iter().any(|fragment| {
            fragment
                .r_group_label(AtomIdx(1))
                .is_some_and(|label| label.number() == Some(1))
        }));
        assert!(fragments.iter().any(|fragment| {
            fragment
                .r_group_label(AtomIdx(0))
                .is_some_and(|label| label.number() == Some(2))
        }));
    }

    // --- safe Option-returning variants ---

    #[test]
    fn test_atom_opt_valid() {
        let mol = ethane();
        assert!(mol.atom_opt(AtomIdx(0)).is_some());
        assert!(mol.atom_opt(AtomIdx(1)).is_some());
        let atom = mol.atom_opt(AtomIdx(0)).unwrap();
        assert_eq!(atom.element.atomic_number(), 6);
    }

    #[test]
    fn test_atom_opt_invalid() {
        let mol = ethane();
        assert!(mol.atom_opt(AtomIdx(2)).is_none());
        assert!(mol.atom_opt(AtomIdx(1000)).is_none());
    }

    #[test]
    fn test_bond_opt_valid() {
        let mol = ethane();
        assert!(mol.bond_opt(BondIdx(0)).is_some());
        let bond = mol.bond_opt(BondIdx(0)).unwrap();
        assert_eq!(bond.order, BondOrder::Single);
    }

    #[test]
    fn test_bond_opt_invalid() {
        let mol = ethane();
        assert!(mol.bond_opt(BondIdx(1)).is_none());
        assert!(mol.bond_opt(BondIdx(1000)).is_none());
    }

    #[test]
    fn test_neighbors_opt_valid() {
        let mol = ethane();
        let neighbors = mol.neighbors_opt(AtomIdx(0)).unwrap();
        assert_eq!(neighbors.len(), 1);
        assert_eq!(neighbors[0].0, AtomIdx(1));
    }

    #[test]
    fn test_neighbors_opt_isolated_atom() {
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::new(Element::C));
        b.add_atom(Atom::new(Element::N));
        let mol = b.build();
        let neighbors = mol.neighbors_opt(AtomIdx(0)).unwrap();
        assert_eq!(neighbors.len(), 0);
    }

    #[test]
    fn test_neighbors_opt_invalid() {
        let mol = ethane();
        assert!(mol.neighbors_opt(AtomIdx(2)).is_none());
        assert!(mol.neighbors_opt(AtomIdx(1000)).is_none());
    }

    #[test]
    fn test_degree_opt_valid() {
        let mol = ethane();
        assert_eq!(mol.degree_opt(AtomIdx(0)), Some(1));
        assert_eq!(mol.degree_opt(AtomIdx(1)), Some(1));
    }

    #[test]
    fn test_degree_opt_isolated_atom() {
        let mut b = MoleculeBuilder::new();
        b.add_atom(Atom::new(Element::C));
        b.add_atom(Atom::new(Element::N));
        let mol = b.build();
        assert_eq!(mol.degree_opt(AtomIdx(0)), Some(0));
        assert_eq!(mol.degree_opt(AtomIdx(1)), Some(0));
    }

    #[test]
    fn test_degree_opt_invalid() {
        let mol = ethane();
        assert!(mol.degree_opt(AtomIdx(2)).is_none());
        assert!(mol.degree_opt(AtomIdx(1000)).is_none());
    }

    #[test]
    fn test_degree_opt_multiple_bonds() {
        // Create a central atom with 3 neighbors
        let mut b = MoleculeBuilder::new();
        let center = b.add_atom(Atom::new(Element::C));
        let n1 = b.add_atom(Atom::new(Element::C));
        let n2 = b.add_atom(Atom::new(Element::N));
        let n3 = b.add_atom(Atom::new(Element::O));
        b.add_bond(center, n1, BondOrder::Single).unwrap();
        b.add_bond(center, n2, BondOrder::Double).unwrap();
        b.add_bond(center, n3, BondOrder::Single).unwrap();
        let mol = b.build();
        assert_eq!(mol.degree_opt(center), Some(3));
        assert_eq!(mol.degree_opt(n1), Some(1));
        assert_eq!(mol.degree_opt(n2), Some(1));
        assert_eq!(mol.degree_opt(n3), Some(1));
    }

    #[test]
    fn set_tag_uses_nonzero_u16_niche() {
        assert_eq!(
            core::mem::size_of::<Option<core::num::NonZeroU16>>(),
            2,
            "tag Option must stay two bytes"
        );
        let mut mol = ethane();
        mol.set_tag(AtomIdx(0), Some(1));
        assert_eq!(mol.atom_tag(AtomIdx(0)).map(|t| t.get()), Some(1));
        mol.set_tag(AtomIdx(0), Some(0)); // zero is absence
        assert!(mol.atom_tag(AtomIdx(0)).is_none());
        mol.set_tag(AtomIdx(1), Some(u16::MAX));
        assert_eq!(mol.atom_tag(AtomIdx(1)).map(|t| t.get()), Some(u16::MAX));
        mol.set_tag(AtomIdx(1), None);
        assert!(mol.atom_tag(AtomIdx(1)).is_none());
    }

    #[test]
    fn atom_tags_survive_core_edits_without_changing_atom_api() {
        // An exhaustive 1.x struct literal must remain source-compatible.
        let atom = Atom {
            element: Element::C,
            isotope: None,
            charge: 0,
            hydrogen_count: None,
            aromatic: false,
            chirality: crate::Chirality::None,
            wildcard: false,
            atom_map: None,
            cip_code: None,
        };
        let mut mol = ethane();
        assert!(mol.atom_tags.is_empty());
        mol.set_tag(AtomIdx(0), None);
        assert!(mol.atom_tags.is_empty());
        mol.set_tag(AtomIdx(0), Some(17));
        mol.set_tag(AtomIdx(1), Some(23));
        assert_eq!(mol.atom(AtomIdx(0)), &atom);
        let expected = vec![Some(17), Some(23)];
        for edited in [
            mol.clone(),
            MoleculeBuilder::from_molecule(&mol).build(),
            mol.with_atom_charge(AtomIdx(0), 1),
            mol.with_atom_element(AtomIdx(0), Element::N),
            mol.with_atom_aromatic(AtomIdx(0), true),
            mol.with_bond_order(BondIdx(0), BondOrder::Double),
            mol.with_bond_removed(BondIdx(0)),
        ] {
            let actual: Vec<_> = edited
                .atoms()
                .map(|(i, _)| edited.atom_tag(i).map(core::num::NonZeroU16::get))
                .collect();
            assert_eq!(actual, expected);
        }
        let (mut extended, born) = mol.with_atom_added(atom);
        assert!(extended.atom_tag(born).is_none());
        extended.set_tag(born, Some(29));
        let connected = extended
            .with_bond_added(AtomIdx(1), born, BondOrder::Single)
            .unwrap()
            .0;
        assert_eq!(connected.atom_tag(born).unwrap().get(), 29);
        let (removed, remap) = extended.with_atom_removed(AtomIdx(0));
        assert_eq!(removed.atom_tag(remap[1].unwrap()).unwrap().get(), 23);
        assert_eq!(removed.atom_tag(remap[2].unwrap()).unwrap().get(), 29);
        extended.remove_atom(AtomIdx(0));
        assert_eq!(extended.atom_tag(AtomIdx(0)).unwrap().get(), 23);
        assert_eq!(extended.atom_tag(AtomIdx(1)).unwrap().get(), 29);
        let new = extended.add_atom(Atom::new(Element::O));
        assert!(extended.atom_tag(new).is_none());
        assert!(extended.atom_tag(AtomIdx(u32::MAX)).is_none());
    }

    #[test]
    fn builder_tags_keep_new_atoms_untagged() {
        let mut mol = ethane();
        mol.set_tag(AtomIdx(0), Some(31));
        let mut builder = MoleculeBuilder::from_molecule(&mol);
        let born = builder.add_atom(Atom::new(Element::O));
        builder.set_tag(born, Some(41));
        builder.set_tag(born, Some(0));
        let built = builder.build();
        assert_eq!(built.atom_tag(AtomIdx(0)), mol.atom_tag(AtomIdx(0)));
        assert!(built.atom_tag(born).is_none());
        let mut larger = MoleculeBuilder::from_molecule(&built);
        larger.copy_atom_tags_from(&mol);
        assert!(larger.build().atom_tag(born).is_none());
    }
}
