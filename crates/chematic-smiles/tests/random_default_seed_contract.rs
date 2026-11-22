use chematic_smiles::{RdkitSmilesParams, parse, rdkit_random_smiles_with};

#[test]
fn first_zero_seed_uses_default_seed_42() {
    let molecule = parse("CC(O)C(=O)N").unwrap();
    let params = RdkitSmilesParams::default();
    let initial = rdkit_random_smiles_with(&molecule, 8, 0, &params).unwrap();
    let explicit = rdkit_random_smiles_with(&molecule, 8, 42, &params).unwrap();
    assert_eq!(initial, explicit);
    assert_eq!(initial.len(), 8);
}
