use chematic_ff::{assign_mmff94_numeric_types, mmff94_torsion_term_params};

#[test]
fn peroxide_and_sulfenate_torsions_match_rdkit_reference() {
    // RDKit 2025.09.3, MMFF94, GetMMFFTorsionParams for atoms 0-1-2-3.
    for (smiles, expected_v2) in [("COOC", -2.0), ("COSC", -4.0)] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let types = assign_mmff94_numeric_types(&mol).unwrap();
        let torsion = mmff94_torsion_term_params(&mol, 0, 1, 2, 3, &types).unwrap();
        assert_eq!(
            (torsion.v1, torsion.v2, torsion.v3),
            (0.0, expected_v2, 0.0),
            "{smiles}"
        );
    }
}
