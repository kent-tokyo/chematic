use super::*;

#[test]
fn bridge_fragment_counters_and_attachment_order() {
    let source = chematic_smiles::parse("CCOC").unwrap();
    let result = fragment_single_bonds_checked(&source, &[BondIdx(2), BondIdx(0)], None).unwrap();
    assert_eq!(result.cuts_per_atom, [1, 1, 1, 1]);
    assert_eq!(result.molecule.atom_count(), 8);
    assert_eq!(result.molecule.bond_count(), 5);
    assert_eq!(
        result
            .attachments
            .iter()
            .map(|a| (a.anchor.0, a.dummy.0, a.isotope))
            .collect::<Vec<_>>(),
        [(3, 4, 2), (2, 5, 3), (1, 6, 0), (0, 7, 1)]
    );
    assert_eq!(source.atom_count(), 4);
    assert_eq!(source.bond_count(), 3);
}

#[test]
fn bridge_fragment_keeps_original_atoms_and_hydrogens() {
    for input in [
        "CCOC",
        "[13CH3:7]C(=O)O",
        "CC[NH3+]",
        "c1ccccc1CC#N",
        "[2H]OC",
    ] {
        let mol = chematic_smiles::parse(input).unwrap();
        for (index, bond) in mol.bonds() {
            if bond.order != BondOrder::Single || is_ring_bond(&mol, index) {
                continue;
            }
            let res = fragment_single_bonds_checked(&mol, &[index], Some(&[(7, 9)])).unwrap();
            for (a, atom) in mol.atoms() {
                let after = res.molecule.atom(a);
                assert_eq!(atom.element, after.element);
                assert_eq!(atom.charge, after.charge);
                assert_eq!(atom.isotope, after.isotope);
                assert_eq!(atom.atom_map, after.atom_map);
                assert_eq!(
                    mol.implicit_hydrogen_count(a),
                    res.molecule.implicit_hydrogen_count(a)
                );
            }
            assert_eq!(res.attachments[0].isotope, 7);
            assert_eq!(res.attachments[1].isotope, 9);
        }
    }
}

#[test]
fn bridge_fragment_refuses_unsupported_edits() {
    for (input, cuts) in [
        ("CCO", vec![BondIdx(2)]),
        ("CCO", vec![BondIdx(0), BondIdx(0)]),
        ("C=C", vec![BondIdx(0)]),
        ("C1CC1", vec![BondIdx(0)]),
        ("N[C@H](C)O", vec![BondIdx(0)]),
        ("C/C=C/C", vec![BondIdx(0)]),
        ("*CC", vec![BondIdx(1)]),
    ] {
        let mol = chematic_smiles::parse(input).unwrap();
        assert!(
            fragment_single_bonds_checked(&mol, &cuts, None).is_err(),
            "{input}"
        );
    }
    let mol = chematic_smiles::parse("CCO").unwrap();
    assert!(fragment_single_bonds_checked(&mol, &[BondIdx(0)], Some(&[])).is_err());
}

#[test]
fn bridge_fragment_empty_is_a_true_graph_noop() {
    let mol = chematic_smiles::parse("N[C@H](C)O").unwrap();
    let result = fragment_single_bonds_checked(&mol, &[], None).unwrap();
    assert_eq!(result.cuts_per_atom, [0, 0, 0, 0]);
    assert!(result.attachments.is_empty());
    assert_eq!(
        chematic_smiles::canonical_smiles(&mol),
        chematic_smiles::canonical_smiles(&result.molecule)
    );
}
