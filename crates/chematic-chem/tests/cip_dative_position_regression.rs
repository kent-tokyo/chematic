use chematic_chem::{CipMode, assign_cip_with_mode};
use chematic_core::{AtomIdx, CipCode};

#[test]
fn dative_cleanup_ligand_is_not_a_tetrahedral_position() {
    // RDKit 2026.03.6 cleanup turns the metal bond into a coordinate bond.
    // The original SMILES order consequently has five entries (four physical
    // tetrahedral positions plus the metal), but CIP excludes the dative
    // ligand and labels the carbon from the remaining four positions.
    let fixtures = [
        (
            "C1CC2CN3CCc4cc5c(OCO5)cc4[C@@H]3C[C@H]2([Zn])C[C@@H]1O",
            AtomIdx(18),
            CipCode::S,
        ),
        (
            "O=C(O)[C@@H]1([SnH])C[C@@H](/C=C/CP(=O)(O)O)CCN1",
            AtomIdx(3),
            CipCode::S,
        ),
    ];
    for (smiles, center, expected) in fixtures {
        let parsed = chematic_smiles::parse(smiles).unwrap();
        let mol = chematic_smiles::rdkit_parsed_molecule(&parsed).unwrap();
        let assignment = assign_cip_with_mode(&mol, CipMode::Accurate).unwrap();
        assert_eq!(assignment.get(center), Some(expected), "{smiles}");
    }
}
