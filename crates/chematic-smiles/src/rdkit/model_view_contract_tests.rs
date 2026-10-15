//! Core model views must not depend on unused CIP/potential-center flags.
use super::*;

fn same_core_view(mol: &Molecule) {
    let full =
        rdkit_mol_from_smiles(mol).and_then(|(model, _)| pdb_read::to_chematic(&model, true));
    let optimized = rdkit_model_molecule(mol, true);
    match (full, optimized) {
        (Ok(full), Ok(optimized)) => {
            assert_eq!(
                full.atoms().collect::<Vec<_>>(),
                optimized.atoms().collect::<Vec<_>>()
            );
            assert_eq!(
                full.bonds().collect::<Vec<_>>(),
                optimized.bonds().collect::<Vec<_>>()
            );
            for (idx, _) in full.atoms() {
                assert_eq!(full.neighbor_slice(idx), optimized.neighbor_slice(idx));
                assert_eq!(
                    full.stereo_neighbor_order(idx),
                    optimized.stereo_neighbor_order(idx)
                );
            }
            for (idx, _) in full.bonds() {
                assert_eq!(full.bond_direction(idx), optimized.bond_direction(idx));
            }
        }
        (Err(full), Err(optimized)) => assert_eq!(full.to_string(), optimized.to_string()),
        _ => panic!("model conversion success/error changed"),
    }
}

#[test]
fn core_model_keeps_specified_and_unspecified_stereo() {
    for smiles in [
        "CC(F)C(Cl)Br",
        "N[C@@H](C)C(=O)O",
        "C[C@](C)(F)Cl",
        "F/C=C/F",
        "F/C=C\\F",
        "FC=CF",
        "C1/C=C/CCC1",
        "C1/C=C/CCCCC1",
        "F[Pt@SP1](Cl)(Br)I",
        "N[Co@OH1](F)(Cl)(Br)(I)P",
        "C[C@AL1]=C=CC",
        "[C@H]1(CCCCC1)C",
        "C[C@H]1CCCC[C@@H]1C",
        "[13CH3][C@H](O)C",
    ] {
        same_core_view(&crate::parse(smiles).unwrap());
    }
}

#[test]
fn core_model_keeps_hydrogens_cleanup_and_error_boundaries() {
    for smiles in [
        "[H]O[H]",
        "[2H][C@](F)(Cl)Br",
        "OCl(=O)(=O)=O",
        "[NH4+]",
        "C1C[N+]2=CC3=CC=CC=C3O[Cu]24OC5=CC=CC=C5C=[N+]14",
        "[O]",
        "[CH2]C",
        "[Na+].[Cl-]",
        "c1cccc1",
        "CO(C)C",
    ] {
        same_core_view(&crate::parse_template(smiles).unwrap());
    }
    let mut mol = crate::parse("CO").unwrap();
    let h = mol.add_atom(chematic_core::Atom::new(chematic_core::Element::H));
    mol.add_bond(
        chematic_core::AtomIdx(0),
        h,
        chematic_core::BondOrder::Single,
    )
    .unwrap();
    same_core_view(&mol);
}

#[test]
fn model_callbacks_keep_stereo_after_graph_edits() {
    let mut mol = crate::parse("[H][C@](F)(Cl)Br").unwrap();
    same_core_view(&mol);
    mol.set_element(chematic_core::AtomIdx(2), chematic_core::Element::CL);
    same_core_view(&mol);
    mol.set_hydrogen_count(chematic_core::AtomIdx(1), Some(1));
    same_core_view(&mol);
}

#[test]
fn sanitized_view_cache_is_not_inherited_by_clones_or_edits() {
    let mol = crate::parse("[CH2]C").unwrap();
    let mut view = rdkit_model_molecule(&mol, true).unwrap();
    let slot = chematic_core::DerivedSlot::RdkitParityAromatic;
    type View = Result<Molecule, chematic_perception::AromaticityError>;
    assert!(view.derived_if_computed::<View>(slot).is_some());
    assert!(view.clone().derived_if_computed::<View>(slot).is_none());
    view.set_hydrogen_count(chematic_core::AtomIdx(0), Some(3));
    assert!(view.derived_if_computed::<View>(slot).is_none());
    let plain = rdkit_model_molecule(&mol, false).unwrap();
    assert!(plain.derived_if_computed::<View>(slot).is_none());
}
