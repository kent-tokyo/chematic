use chematic_core::{AtomIdx, Chirality};
use chematic_mol::{parse_mol2, parse_mol2_source_record, write_mol2};

const EXPLICIT_H: &str = include_str!("data/mol2-3d-stereo-0.mol2");
const IMPLICIT_H: &str = include_str!("data/mol2-3d-stereo-noh.mol2");
const AMINO_ACID: &str = include_str!("data/mol2-3d-stereo-1.mol2");
const FOUR_HEAVY: &str = include_str!("data/mol2-3d-stereo-2.mol2");

#[test]
fn mol2_3d_tetrahedral_coordinates_recover_stereo_in_both_readers() {
    for (source, center) in [
        (EXPLICIT_H, 1),
        (IMPLICIT_H, 1),
        (AMINO_ACID, 1),
        (FOUR_HEAVY, 1),
    ] {
        let (mol, coordinates) = parse_mol2(source).unwrap();
        let chirality = mol.atom(AtomIdx(center)).chirality;
        assert!(
            chirality.is_tetrahedral(),
            "3D MOL2 source lacks tetrahedral tag"
        );
        let retained = parse_mol2_source_record(source).unwrap();
        assert_eq!(
            retained.molecule().atom(AtomIdx(center)).chirality,
            chirality
        );

        let mirrored = coordinates
            .iter()
            .map(|&(x, y, z)| (x, y, -z))
            .collect::<Vec<_>>();
        let (opposite, _) = parse_mol2(&write_mol2(&mol, &mirrored)).unwrap();
        assert_eq!(
            opposite.atom(AtomIdx(center)).chirality,
            match chirality {
                Chirality::Clockwise => Chirality::CounterClockwise,
                Chirality::CounterClockwise => Chirality::Clockwise,
                _ => unreachable!(),
            },
            "reflecting genuine 3D coordinates should invert the center"
        );

        let planar = coordinates
            .iter()
            .map(|&(x, y, _)| (x, y, 0.0))
            .collect::<Vec<_>>();
        let (flat, _) = parse_mol2(&write_mol2(&mol, &planar)).unwrap();
        assert_eq!(flat.atom(AtomIdx(center)).chirality, Chirality::None);
    }
}

#[test]
fn mol2_rejects_nonfinite_coordinates_before_stereo_perception() {
    let invalid = EXPLICIT_H.replacen("0.1500", "NaN", 1);
    assert!(parse_mol2(&invalid).is_err());
    assert!(parse_mol2_source_record(&invalid).is_err());
}
