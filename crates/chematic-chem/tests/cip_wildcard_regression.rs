use chematic_chem::{CipMode, CipUnresolvedReason, assign_cip_with_mode};
use chematic_core::AtomIdx;

#[test]
fn wildcard_ligand_is_typed_unresolved_instead_of_guessing_a_cip_label() {
    let smiles = "C12CC[C@](C#CCCCI)(O)[C@@]1(C)CCC1C2CC*2[C@@]1(CCC(N2C)=O)[C]";
    let mol =
        chematic_smiles::rdkit_parsed_molecule(&chematic_smiles::parse(smiles).unwrap()).unwrap();
    let assignment = assign_cip_with_mode(&mol, CipMode::Accurate).unwrap();
    assert_eq!(assignment.get(AtomIdx(20)), None);
    assert!(
        assignment
            .unresolved
            .contains(&(AtomIdx(20), CipUnresolvedReason::OracleUnstable))
    );
}
