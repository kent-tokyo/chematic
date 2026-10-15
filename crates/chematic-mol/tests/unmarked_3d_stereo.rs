//! RDKit 2025.09.3 embedded 3D MOL blocks with atom parity and bond CFG
//! deliberately removed. The conformer is the only source of stereo.

use chematic_mol::{mol2000::read_mol_with_diagnostics, read_mol_v3000_with_diagnostics};
use chematic_smiles::{canonical_smiles, parse, rdkit_canonical_smiles, write};

const CASES: [(&str, &str, &str); 8] = [
    (
        "V2000 tetrahedral carbon",
        include_str!("fixtures/chematic-unmarked-3d-0-mol.mol"),
        "F[C@](Cl)(Br)I",
    ),
    (
        "V3000 tetrahedral carbon",
        include_str!("fixtures/chematic-unmarked-3d-0-mol_v3000.mol"),
        "F[C@](Cl)(Br)I",
    ),
    (
        "V2000 implicit hydrogen",
        include_str!("fixtures/chematic-unmarked-3d-1-mol.mol"),
        "F[C@H](Cl)Br",
    ),
    (
        "V3000 implicit hydrogen",
        include_str!("fixtures/chematic-unmarked-3d-1-mol_v3000.mol"),
        "F[C@H](Cl)Br",
    ),
    (
        "V2000 alanine",
        include_str!("fixtures/chematic-unmarked-3d-2-mol.mol"),
        "N[C@@H](C)C(=O)O",
    ),
    (
        "V3000 alanine",
        include_str!("fixtures/chematic-unmarked-3d-2-mol_v3000.mol"),
        "N[C@@H](C)C(=O)O",
    ),
    (
        "V2000 sulfoxide",
        include_str!("fixtures/chematic-unmarked-3d-3-mol.mol"),
        "C[S@](=O)c1ccccc1",
    ),
    (
        "V3000 sulfoxide",
        include_str!("fixtures/chematic-unmarked-3d-3-mol_v3000.mol"),
        "C[S@](=O)c1ccccc1",
    ),
];

#[test]
fn unmarked_3d_geometry_recovers_tetrahedral_stereo() {
    for (name, block, expected) in CASES {
        let mol = if name.starts_with("V2000") {
            read_mol_with_diagnostics(block).unwrap().mol
        } else {
            read_mol_v3000_with_diagnostics(block).unwrap().mol
        };
        let expected = rdkit_canonical_smiles(&parse(expected).unwrap()).unwrap();
        for smiles in [write(&mol), canonical_smiles(&mol)] {
            let reparsed = parse(&smiles).unwrap();
            assert_eq!(
                rdkit_canonical_smiles(&reparsed).unwrap(),
                expected,
                "{name}: {smiles}"
            );
        }
    }
}

#[test]
fn equivalent_substituents_do_not_gain_3d_stereo() {
    let block = CASES[0].1.replace(" Cl  ", " F   ");
    let mol = read_mol_with_diagnostics(&block).unwrap().mol;
    assert!(!write(&mol).contains('@'));

    for block in [
        include_str!("fixtures/chematic-unmarked-3d-diethyl-sulfoxide-mol.mol"),
        include_str!("fixtures/chematic-unmarked-3d-diethyl-sulfoxide-mol_v3000.mol"),
    ] {
        let mol = if block.contains("V3000") {
            read_mol_v3000_with_diagnostics(block).unwrap().mol
        } else {
            read_mol_with_diagnostics(block).unwrap().mol
        };
        assert!(!write(&mol).contains('@'));
    }
}

#[test]
fn reflecting_coordinates_reverses_the_inferred_enantiomer() {
    let mut lines: Vec<String> = CASES[0].1.lines().map(str::to_owned).collect();
    for line in &mut lines[4..9] {
        let x: f64 = line[..10].trim().parse().unwrap();
        *line = format!("{:>10.4}{}", -x, &line[10..]);
    }
    let reflected = lines.join("\n") + "\n";
    let mol = read_mol_with_diagnostics(&reflected).unwrap().mol;
    let observed = rdkit_canonical_smiles(&parse(&write(&mol)).unwrap()).unwrap();
    let expected = rdkit_canonical_smiles(&parse("F[C@@](Cl)(Br)I").unwrap()).unwrap();
    assert_eq!(observed, expected);
}

#[test]
fn a_2d_declaration_does_not_infer_stereo_from_nonzero_z() {
    let block = CASES[0].1.replace("RDKit          3D", "RDKit          2D");
    let mol = read_mol_with_diagnostics(&block).unwrap().mol;
    assert!(!write(&mol).contains('@'));
}
