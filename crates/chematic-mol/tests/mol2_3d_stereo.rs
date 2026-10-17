use chematic_core::{AtomIdx, Chirality};
use chematic_mol::{
    parse_mol2, parse_mol2_record, parse_mol2_source_record, write_mol2, write_mol2_record,
};
use chematic_smiles::canonical_smiles;

const EXPLICIT_H: &str = include_str!("data/mol2-3d-stereo-0.mol2");
const IMPLICIT_H: &str = include_str!("data/mol2-3d-stereo-noh.mol2");
const AMINO_ACID: &str = include_str!("data/mol2-3d-stereo-1.mol2");
const FOUR_HEAVY: &str = include_str!("data/mol2-3d-stereo-2.mol2");
const EXPLICIT_H_BRIDGEHEADS: &str = include_str!("data/mol2-3d-bridgehead-explicit-h.mol2");
const EXPLICIT_H_ADAMANTANE: &str = include_str!("data/mol2-3d-adamantane-explicit-h.mol2");
const NH_IMINE_PAIR: [&str; 2] = [
    include_str!("data/mol2-3d-nh-imine-e.mol2"),
    include_str!("data/mol2-3d-nh-imine-z.mol2"),
];
const OPENBABEL_STRESS_STEREO: [(&str, usize); 4] = [
    (include_str!("data/mol2-3d-openbabel-651.mol2"), 2),
    (include_str!("data/mol2-3d-openbabel-658.mol2"), 3),
    (include_str!("data/mol2-3d-openbabel-2601.mol2"), 3),
    (include_str!("data/mol2-3d-openbabel-3599.mol2"), 2),
];
const SINGLE_METHYL_RING: &str = include_str!("data/mol2-3d-methylcyclohexane.mol2");
const SYMMETRIC_PHOSPHORUS: &str = include_str!("data/mol2-3d-symmetric-phosphorus.mol2");
const ALKENE_CASES: [&str; 12] = [
    include_str!("data/mol2-3d-ez-00.mol2"),
    include_str!("data/mol2-3d-ez-01.mol2"),
    include_str!("data/mol2-3d-ez-02.mol2"),
    include_str!("data/mol2-3d-ez-03.mol2"),
    include_str!("data/mol2-3d-ez-04.mol2"),
    include_str!("data/mol2-3d-ez-05.mol2"),
    include_str!("data/mol2-3d-ez-06.mol2"),
    include_str!("data/mol2-3d-ez-07.mol2"),
    include_str!("data/mol2-3d-ez-08.mol2"),
    include_str!("data/mol2-3d-ez-09.mol2"),
    include_str!("data/mol2-3d-ez-10.mol2"),
    include_str!("data/mol2-3d-ez-11.mol2"),
];

#[test]
fn mol2_3d_terminal_nh_imine_geometries_remain_distinct() {
    let mut direct = Vec::new();
    let mut retained = Vec::new();
    for input in NH_IMINE_PAIR {
        let (mol, _) = parse_mol2(input).unwrap();
        direct.push(canonical_smiles(&mol));
        let source = parse_mol2_source_record(input).unwrap();
        retained.push(canonical_smiles(source.molecule()));
    }
    for pair in [&direct, &retained] {
        assert!(
            pair.iter()
                .all(|smiles| smiles.contains('/') || smiles.contains('\\'))
        );
        assert_ne!(pair[0], pair[1]);
    }
}

#[test]
fn mol2_3d_alkenes_recover_ez_in_both_readers_and_after_rotation() {
    let mut smiles = Vec::new();
    let mut source_smiles = Vec::new();
    for source in ALKENE_CASES {
        let (mol, coords) = parse_mol2(source).unwrap();
        let canonical = canonical_smiles(&mol);
        assert!(canonical.contains('/'), "E/Z missing from {canonical}");
        let retained = parse_mol2_source_record(source).unwrap();
        let retained_canonical = canonical_smiles(retained.molecule());
        assert!(retained_canonical.contains('/'));
        source_smiles.push(retained_canonical);

        let mut record = parse_mol2_record(source).unwrap();
        record.coords = coords.iter().map(|&(x, y, z)| (z, x, y)).collect();
        let (rotated_mol, _) = parse_mol2(&write_mol2_record(&record)).unwrap();
        assert_eq!(canonical_smiles(&rotated_mol), canonical);
        smiles.push(canonical);
    }
    for index in (0..smiles.len()).step_by(2) {
        assert_ne!(
            smiles[index],
            smiles[index + 1],
            "opposite geometries were collapsed"
        );
        assert_ne!(
            source_smiles[index],
            source_smiles[index + 1],
            "source reader collapsed opposite geometries"
        );
    }
}

#[test]
fn mol2_alkene_geometry_must_have_a_defined_side() {
    let mut record = parse_mol2_record(ALKENE_CASES[0]).unwrap();
    record.coords = record
        .coords
        .iter()
        .enumerate()
        .map(|(index, _)| (index as f64, 0.0, 0.0))
        .collect();
    let degenerate = write_mol2_record(&record);
    let (mol, _) = parse_mol2(&degenerate).unwrap();
    assert!(!canonical_smiles(&mol).contains('/'));
    let source = parse_mol2_source_record(&degenerate).unwrap();
    assert!(!canonical_smiles(source.molecule()).contains('/'));

    // A planar alkene is still a valid 3D MOL2 conformer for E/Z, even
    // though it cannot specify tetrahedral chirality.
    let mut planar = parse_mol2_record(ALKENE_CASES[0]).unwrap();
    for point in &mut planar.coords {
        point.2 = 0.0;
    }
    let (mol, _) = parse_mol2(&write_mol2_record(&planar)).unwrap();
    assert!(canonical_smiles(&mol).contains('/'));

    // Rotate one alkene end's two substituents by 90 degrees about the
    // double-bond axis. A twisted, orthogonal geometry has no stable E/Z.
    let mut twisted = parse_mol2_record(ALKENE_CASES[0]).unwrap();
    let end = twisted.coords[2];
    let other = twisted.coords[1];
    let axis = (end.0 - other.0, end.1 - other.1, end.2 - other.2);
    let length = (axis.0 * axis.0 + axis.1 * axis.1 + axis.2 * axis.2).sqrt();
    let unit = (axis.0 / length, axis.1 / length, axis.2 / length);
    for index in [3, 5] {
        let p = twisted.coords[index];
        let v = (p.0 - end.0, p.1 - end.1, p.2 - end.2);
        let along = unit.0 * v.0 + unit.1 * v.1 + unit.2 * v.2;
        let cross = (
            unit.1 * v.2 - unit.2 * v.1,
            unit.2 * v.0 - unit.0 * v.2,
            unit.0 * v.1 - unit.1 * v.0,
        );
        twisted.coords[index] = (
            end.0 + unit.0 * along + cross.0,
            end.1 + unit.1 * along + cross.1,
            end.2 + unit.2 * along + cross.2,
        );
    }
    let (mol, _) = parse_mol2(&write_mol2_record(&twisted)).unwrap();
    assert!(!canonical_smiles(&mol).contains('/'));
}

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
fn mol2_3d_bridgeheads_with_explicit_hydrogens_recover_stereo() {
    let (mol, coordinates) = parse_mol2(EXPLICIT_H_BRIDGEHEADS).unwrap();
    let retained = parse_mol2_source_record(EXPLICIT_H_BRIDGEHEADS).unwrap();
    for center in [2, 5] {
        let chirality = mol.atom(AtomIdx(center)).chirality;
        assert!(
            chirality.is_tetrahedral(),
            "bridgehead {center} lost stereo"
        );
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
fn mol2_3d_symmetric_multibridgehead_cage_remains_untagged() {
    let (mol, _) = parse_mol2(EXPLICIT_H_ADAMANTANE).unwrap();
    let source = parse_mol2_source_record(EXPLICIT_H_ADAMANTANE).unwrap();
    for molecule in [&mol, source.molecule()] {
        assert!(
            molecule
                .atoms()
                .all(|(_, atom)| atom.chirality == Chirality::None),
            "symmetric cage acquired an arbitrary tetrahedral tag"
        );
    }
}

#[test]
fn mol2_3d_recovers_paired_ring_and_phosphorus_centers() {
    for (input, expected_centers) in OPENBABEL_STRESS_STEREO {
        let (mol, coordinates) = parse_mol2(input).unwrap();
        let retained = parse_mol2_source_record(input).unwrap();
        let tagged: Vec<_> = mol
            .atoms()
            .filter(|(_, atom)| atom.chirality.is_tetrahedral())
            .map(|(idx, _)| idx)
            .collect();
        assert_eq!(tagged.len(), expected_centers);
        for &center in &tagged {
            assert_eq!(
                retained.molecule().atom(center).chirality,
                mol.atom(center).chirality
            );
        }

        let mirrored = coordinates
            .iter()
            .map(|&(x, y, z)| (x, y, -z))
            .collect::<Vec<_>>();
        let (opposite, _) = parse_mol2(&write_mol2(&mol, &mirrored)).unwrap();
        for &center in &tagged {
            assert_ne!(opposite.atom(center).chirality, mol.atom(center).chirality);
            assert!(opposite.atom(center).chirality.is_tetrahedral());
        }
    }
}

#[test]
fn mol2_3d_single_methyl_ring_remains_untagged() {
    let (mol, _) = parse_mol2(SINGLE_METHYL_RING).unwrap();
    let retained = parse_mol2_source_record(SINGLE_METHYL_RING).unwrap();
    for molecule in [&mol, retained.molecule()] {
        assert!(
            molecule
                .atoms()
                .all(|(_, atom)| atom.chirality == Chirality::None),
            "a ring with one methyl group has no tetrahedral center"
        );
    }
}

#[test]
fn mol2_3d_phosphorus_with_duplicate_ligands_remains_untagged() {
    let (mol, _) = parse_mol2(SYMMETRIC_PHOSPHORUS).unwrap();
    let retained = parse_mol2_source_record(SYMMETRIC_PHOSPHORUS).unwrap();
    for molecule in [&mol, retained.molecule()] {
        assert!(
            molecule
                .atoms()
                .all(|(_, atom)| atom.chirality == Chirality::None),
            "phosphorus with two methyl ligands is not tetrahedrally stereogenic"
        );
    }
}

#[test]
fn mol2_rejects_nonfinite_coordinates_before_stereo_perception() {
    let invalid = EXPLICIT_H.replacen("0.1500", "NaN", 1);
    assert!(parse_mol2(&invalid).is_err());
    assert!(parse_mol2_source_record(&invalid).is_err());
}
