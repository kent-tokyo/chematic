use chematic_chem::standardize::{cleanup, cleanup_checked, cleanup_unchecked};

#[test]
fn cleanup_alias_and_checked_cleanup_agree_on_ordinary_molecules() {
    for smiles in ["CCO", "CC(=O)O", "[Na+].[Cl-]", "C[N+](=O)[O-]"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let unchecked = cleanup_unchecked(&mol);
        let alias = cleanup(&mol);
        let checked = cleanup_checked(&mol).unwrap();
        let canonical = |mol| chematic_smiles::rdkit_smiles(mol, &Default::default()).unwrap();
        assert_eq!(canonical(&unchecked), canonical(&alias));
        assert_eq!(canonical(&unchecked), canonical(&checked));
    }
}
