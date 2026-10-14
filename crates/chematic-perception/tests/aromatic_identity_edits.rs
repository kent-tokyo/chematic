use chematic_core::{AtomIdx, BondIdx, BondOrder, DerivedSlot};
use chematic_perception::{
    apply_aromaticity_rdkit_parity_experimental, rdkit_model_may_disagree,
    rdkit_parity_view_is_identity,
};

#[test]
fn saturated_identity_does_not_survive_a_conjugating_edit() {
    let mut mol = chematic_smiles::parse("C1CCCCC1").unwrap();
    assert!(rdkit_parity_view_is_identity(&mol));
    for bond in [0, 2, 4] {
        mol.set_bond_order(BondIdx(bond), BondOrder::Double);
    }
    assert!(!rdkit_parity_view_is_identity(&mol));
    let view = apply_aromaticity_rdkit_parity_experimental(&mol).unwrap();
    assert!(view.atoms().all(|(_, atom)| atom.aromatic));
    assert!(
        view.bonds()
            .all(|(_, bond)| bond.order == BondOrder::Aromatic)
    );
    assert!(mol.atoms().all(|(_, atom)| !atom.aromatic));
}

#[test]
fn forest_proof_is_invalidated_when_a_conjugated_ring_is_closed() {
    let mut mol = chematic_smiles::parse("C=CC=CC=C").unwrap();
    assert!(rdkit_parity_view_is_identity(&mol));
    assert!(
        mol.derived_if_computed::<Vec<bool>>(DerivedSlot::RingBondFlags)
            .unwrap()
            .iter()
            .all(|&cyclic| !cyclic)
    );
    mol.add_bond(AtomIdx(0), AtomIdx(5), BondOrder::Single)
        .unwrap();
    assert!(!rdkit_parity_view_is_identity(&mol));
    let view = apply_aromaticity_rdkit_parity_experimental(&mol).unwrap();
    assert!(view.atoms().all(|(_, atom)| atom.aromatic));
}

#[test]
fn nonaromatic_shortcut_keeps_radical_and_cleanup_guards() {
    for smiles in ["[CH2]", "OCl(=O)(=O)=O", "N1CCN->[Cu]1"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        assert!(rdkit_model_may_disagree(&mol), "{smiles}");
    }
}
