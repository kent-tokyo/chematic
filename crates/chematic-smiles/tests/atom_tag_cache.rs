//! Warm-cache regression cases for caller atom labels (PR #668).
use chematic_core::AtomIdx;
use chematic_perception::apply_aromaticity_rdkit_parity_experimental;
use chematic_smiles::parse;

#[test]
fn setting_tag_refreshes_cached_molecule() {
    let mut mol = parse("C1=CC=CC=C1").unwrap();
    let before = apply_aromaticity_rdkit_parity_experimental(&mol).unwrap();
    assert_eq!(before.atom_tag(AtomIdx(0)), None);
    mol.set_tag(AtomIdx(0), Some(42));
    let after = apply_aromaticity_rdkit_parity_experimental(&mol).unwrap();
    assert_eq!(after.atom_tag(AtomIdx(0)), mol.atom_tag(AtomIdx(0)));
}

#[test]
fn clearing_tag_refreshes_cached_molecule() {
    let mut mol = parse("C1=CC=CC=C1").unwrap();
    mol.set_tag(AtomIdx(0), Some(42));
    let before = apply_aromaticity_rdkit_parity_experimental(&mol).unwrap();
    assert_eq!(before.atom_tag(AtomIdx(0)).unwrap().get(), 42);
    mol.set_tag(AtomIdx(0), None);
    let after = apply_aromaticity_rdkit_parity_experimental(&mol).unwrap();
    assert_eq!(after.atom_tag(AtomIdx(0)), None);
}

#[test]
fn tags_survive_perception_and_kekulization() {
    for smiles in ["c1ccccc1", "C1=CC=CC=C1", "c1cc[nH]c1", "CCO"] {
        let mut mol = parse(smiles).unwrap();
        for i in 0..mol.atom_count() {
            mol.set_tag(AtomIdx(i as u32), Some(i as u16 + 1));
        }
        let perceived = chematic_perception::apply_aromaticity(&mol);
        let parity = apply_aromaticity_rdkit_parity_experimental(&mol).unwrap();
        let kekule = chematic_core::apply_kekule(&mol, &chematic_core::kekulize(&mol).unwrap());
        for copy in [mol.clone(), perceived, parity, kekule] {
            for (idx, _) in mol.atoms() {
                assert_eq!(copy.atom_tag(idx), mol.atom_tag(idx), "{smiles}: {idx:?}");
            }
        }
        mol.set_tag(AtomIdx(0), Some(0));
        assert!(
            apply_aromaticity_rdkit_parity_experimental(&mol)
                .unwrap()
                .atom_tag(AtomIdx(0))
                .is_none()
        );
    }
}
