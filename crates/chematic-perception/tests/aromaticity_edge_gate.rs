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

#[test]
fn rdkit_like_covers_the_fused_nonalternant_holdout() {
    // These are the five historical A1 false-negative fixtures.  Keep the
    // expected atom/bond counts pinned to the RDKit-like model rather than
    // asserting only that the call succeeds; the latter would miss partial
    // aromaticity propagation into a fused component.
    let cases = [
        (
            "fused-macrocycle-1",
            "N1=C2C(N(CC(O)=O)C(=O)N=C2N(C2C=C(C(F)(F)F)C=C(C=2)C(F)(F)F)C2C1=CC=CC=2)=O",
            20,
            21,
        ),
        (
            "fused-macrocycle-2",
            "[C@H]12N(C([C@H](NC(=O)[C@H]([C@H](OC(=O)[C@@H](N(C)C(CN(C)C1=O)=O)C(C)C)C)NC(=O)C1C=C(OC)C(C)=C3OC4=C(C)C(=O)C(=C(C4=NC=13)C(=O)N[C@H]1C(=O)N[C@@H](C(C)C)C(N3[C@H](C(=O)N(CC(N([C@H](C(C)C)C(O[C@H]1C)=O)C)=O)C)CCC3)=O)N)C(C)C)=O)CCC2",
            14,
            15,
        ),
        (
            "fused-macrocycle-3",
            "C12N(C3C=CC=CC=3)C3=NC(=O)N(C)C(C3=NC1=CC=CC=2)=O",
            20,
            21,
        ),
        ("azulene-kekule", "C1=CC2=CC=CC=CC2=C1", 10, 10),
        ("purine", "c1cnc2[nH]cnc2n1", 9, 10),
    ];
    for (name, smiles, atoms, bonds) in cases {
        let mol = kekulized(smiles);
        let model = assign_aromaticity_ex(&mol, AromaticityAlgorithm::RdkitLike);
        let aromatic_bonds = (0..mol.bond_count())
            .filter(|&idx| model.is_bond_aromatic(chematic_core::BondIdx(idx as u32)))
            .count();
        assert_eq!(
            (model.aromatic_atom_count(), aromatic_bonds),
            (atoms, bonds),
            "{name}: fused/non-alternant aromaticity holdout"
        );
    }
}
