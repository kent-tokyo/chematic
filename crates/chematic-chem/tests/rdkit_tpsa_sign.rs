//! `rdMolDescriptors.CalcTPSA` returns +0.0 for a molecule without N or O
//! (RDKit 2026.03.1); the RDKit-default TPSA must keep that sign.

#[test]
fn rdkit_tpsa_without_n_or_o_is_positive_zero() {
    for smiles in ["CCCC", "c1ccccc1", "ClCCl"] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let tpsa = chematic_chem::descriptors::rdkit_tpsa(&mol);
        assert_eq!(tpsa.to_bits(), 0.0f64.to_bits(), "{smiles}: {tpsa}");
    }
}
