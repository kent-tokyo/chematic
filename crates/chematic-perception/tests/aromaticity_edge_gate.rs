//! Public aromaticity regression gate for fused and non-alternant systems.

use chematic_perception::{AromaticityAlgorithm, assign_aromaticity_ex};
use chematic_smiles::parse;

fn kekulized(smiles: &str) -> chematic_core::Molecule {
    let mol = parse(smiles).expect("fixture parses");
    let kekule = chematic_core::kekulize(&mol).expect("fixture kekulizes");
    chematic_core::apply_kekule(&mol, &kekule)
}

#[test]
fn rdkit_like_covers_purine_and_azulene() {
    let cases = [
        ("purine", "c1cnc2[nH]cnc2n1", 9),
        ("azulene", "C1=CC2=CC=CC=CC2=C1", 10),
    ];
    for (name, smiles, expected) in cases {
        let model = assign_aromaticity_ex(&kekulized(smiles), AromaticityAlgorithm::RdkitLike);
        assert_eq!(
            model.aromatic_atom_count(),
            expected,
            "{name}: RDKit-like fused/non-alternant gate"
        );
    }
}
