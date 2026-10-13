use chematic_chem::{CipMode, assign_cip_with_mode};
use chematic_core::{AtomIdx, CipCode};

#[test]
fn direct_isotope_priority_precedes_remote_pseudoasymmetry() {
    // At atoms 18 and 23 the two carbon branches are constitutionally tied until
    // the directly labelled 200-C atom is reached. Sequence Rule 2 therefore
    // resolves the branch order before Rule 5 pseudoasymmetry is considered.
    let smiles =
        "C/C(=C\\C1=CC=C(C2=CC=CC=C2)C=C1)C(=O)N[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[200C@H]1O";
    let mol =
        chematic_smiles::rdkit_parsed_molecule(&chematic_smiles::parse(smiles).unwrap()).unwrap();
    let assignment = assign_cip_with_mode(&mol, CipMode::Accurate).unwrap();
    assert_eq!(assignment.get(AtomIdx(18)), Some(CipCode::S));
    assert_eq!(assignment.get(AtomIdx(23)), Some(CipCode::S));
}
