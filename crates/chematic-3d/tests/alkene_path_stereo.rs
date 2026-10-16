use chematic_3d::rdkit_embed::{AlkenePathStereo as Side, rdkit_alkene_path_stereo};
use chematic_core::AtomIdx;

#[test]
fn all_four_outer_pairs_and_both_path_orders_have_the_declared_relation() {
    for (text, same) in [("F/C(Cl)=C(Br)/I", false), ("F/C(Cl)=C(Br)\\I", true)] {
        let mol = chematic_smiles::parse(text).unwrap();
        for (a, d, cis) in [(0, 5, same), (0, 4, !same), (2, 5, !same), (2, 4, same)] {
            let expected = Some(if cis { Side::Cis } else { Side::Trans });
            assert_eq!(
                rdkit_alkene_path_stereo(&mol, [a, 1, 3, d]).unwrap(),
                expected
            );
            assert_eq!(
                rdkit_alkene_path_stereo(&mol, [d, 3, 1, a]).unwrap(),
                expected
            );
        }
    }
}

#[test]
fn explicit_hydrogen_pairs_have_complementary_relations() {
    let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse("F/C=C/Cl").unwrap());
    for (path, expected) in [
        ([0, 1, 2, 3], Side::Trans),
        ([4, 1, 2, 3], Side::Cis),
        ([0, 1, 2, 5], Side::Cis),
        ([4, 1, 2, 5], Side::Trans),
    ] {
        assert_eq!(
            rdkit_alkene_path_stereo(&mol, path).unwrap(),
            Some(expected)
        );
    }
}

#[test]
fn undeclared_stereo_is_not_guessed() {
    let mol = chematic_smiles::parse("FC=CCl").unwrap();
    assert_eq!(rdkit_alkene_path_stereo(&mol, [0, 1, 2, 3]).unwrap(), None);
    assert_eq!(rdkit_alkene_path_stereo(&mol, [3, 2, 1, 0]).unwrap(), None);
}

#[test]
fn invalid_paths_are_errors_before_any_indexing() {
    let mol = chematic_smiles::parse("F/C(Cl)=C(Br)/I").unwrap();
    for path in [
        [0, 1, 3, usize::MAX],
        [0, 1, 3, 0],
        [0, 2, 3, 5],
        [0, 1, 2, 3],
    ] {
        assert!(rdkit_alkene_path_stereo(&mol, path).is_err(), "{path:?}");
    }
}

#[test]
fn a_single_centered_bond_is_not_an_alkene_path() {
    let mol = chematic_smiles::parse("CCCC").unwrap();
    assert!(rdkit_alkene_path_stereo(&mol, [0, 1, 2, 3]).is_err());
}

#[test]
fn independent_atom_permutations_preserve_the_path_relation() {
    let mol = chematic_smiles::parse("F/C(Cl)=C(Br)/I").unwrap();
    for shift in 0..6 {
        for reversed in [false, true] {
            let mut order: Vec<_> = (0..6).map(|i| AtomIdx((i + shift) % 6)).collect();
            if reversed {
                order.reverse();
            }
            let renamed = mol.renumber_atoms(&order).unwrap();
            let map = |old| order.iter().position(|a| a.0 == old).unwrap();
            assert_eq!(
                rdkit_alkene_path_stereo(&renamed, [map(0), map(1), map(3), map(5)]).unwrap(),
                Some(Side::Trans)
            );
        }
    }
}

#[test]
fn inspection_does_not_mutate_the_original_graph_or_stereo() {
    let mol = chematic_smiles::parse("F/C(Cl)=C(Br)/I").unwrap();
    let before = chematic_smiles::write(&mol);
    for _ in 0..3 {
        rdkit_alkene_path_stereo(&mol, [5, 3, 1, 0]).unwrap();
    }
    assert_eq!(chematic_smiles::write(&mol), before);
}
