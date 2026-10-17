use super::*;

#[test]
fn two_vacant_nontetrahedral_sites_survive_all_smiles_writers() {
    // RDKit 2025.09.3 canonical isomeric SMILES. SP2, TB3, and OH4 have
    // two unoccupied positions, with root/non-root and changed tags.
    for (input, expected) in [
        ("[Pt@SP1](F)Cl", "[F][Pt@SP2][Cl]"),
        ("F[Pt@SP1]Cl", "[F][Pt@SP2][Cl]"),
        ("[Pt@SP2](F)Cl", "[F][Pt@SP3][Cl]"),
        ("Cl[Pt@SP3]F", "[F][Pt@SP3][Cl]"),
        ("[P@TB1](F)(Cl)Br", "F[P@TB20](Cl)Br"),
        ("F[P@TB10](Cl)Br", "F[P@TB10](Cl)Br"),
        ("Br[P@TB19](Cl)F", "F[P@TB16](Cl)Br"),
        ("[S@OH1](F)(Cl)(Br)I", "F[S@OH5](Cl)(Br)I"),
        ("F[S@OH27](Cl)(Br)I", "F[S@OH21](Cl)(Br)I"),
        ("I[S@OH30](Br)(Cl)F", "F[S@OH3](Cl)(Br)I"),
    ] {
        let mol = crate::parse(input).unwrap();
        assert_eq!(
            crate::rdkit_canonical_smiles(&mol).unwrap(),
            expected,
            "{input}"
        );
        assert_eq!(
            mol.stereo_neighbor_order(AtomIdx(if input.starts_with('[') { 0 } else { 1 }))
                .unwrap()
                .iter()
                .filter(|&&id| id == STEREO_VACANT_SENTINEL)
                .count(),
            2,
            "{input}"
        );
        for written in [write(&mol), crate::canonical_smiles(&mol)] {
            let read_back = crate::parse(&written).unwrap();
            assert_eq!(
                crate::rdkit_canonical_smiles(&read_back).unwrap(),
                expected,
                "{input} -> {written}"
            );
        }
        for seed in 1..=16 {
            let written = crate::random_smiles(&mol, seed);
            let read_back = crate::parse(&written).unwrap();
            assert_eq!(
                crate::rdkit_canonical_smiles(&read_back).unwrap(),
                expected,
                "{input} -> {written} (seed {seed})"
            );
        }
    }
}
