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

// Independent stress-corpus rows, embedded by RDKit 2025.09.3 with all
// CTAB parity/wedge/CFG markers removed. Both readers must infer stereo from
// the coordinates alone.
const RING_CASES: [(&str, &str, &str); 10] = [
    (
        "V2000 paired bridgeheads",
        include_str!("fixtures/chematic-unmarked-3d-bridgehead-v2000.mol"),
        "C1C[C@H]2CC[C@@H]1C2",
    ),
    (
        "V3000 paired bridgeheads",
        include_str!("fixtures/chematic-unmarked-3d-bridgehead-v3000.mol"),
        "C1C[C@H]2CC[C@@H]1C2",
    ),
    (
        "V2000 ring carbons",
        include_str!("fixtures/chematic-unmarked-3d-ring-carbon-v2000.mol"),
        "O=C(O)[C@H]1NCCC[C@H]1C(=O)O",
    ),
    (
        "V3000 ring carbons",
        include_str!("fixtures/chematic-unmarked-3d-ring-carbon-v3000.mol"),
        "O=C(O)[C@H]1NCCC[C@H]1C(=O)O",
    ),
    (
        "V2000 selenium is not tetrahedral",
        include_str!("fixtures/chematic-unmarked-3d-selenium-v2000.mol"),
        "C/C(CN(C)C)=C(/n1ccnc1)[Se]1=C(OCC(C)C)C=CC=C1",
    ),
    (
        "V3000 selenium is not tetrahedral",
        include_str!("fixtures/chematic-unmarked-3d-selenium-v3000.mol"),
        "C/C(CN(C)C)=C(/n1ccnc1)[Se]1=C(OCC(C)C)C=CC=C1",
    ),
    (
        "V2000 nine-membered ring alkene",
        include_str!("fixtures/chematic-unmarked-3d-ring-alkene-9-v2000.mol"),
        "C1=C/CCCCCCC/1",
    ),
    (
        "V3000 nine-membered ring alkene",
        include_str!("fixtures/chematic-unmarked-3d-ring-alkene-9-v3000.mol"),
        "C1=C/CCCCCCC/1",
    ),
    (
        "V2000 ten-membered ring alkene",
        include_str!("fixtures/chematic-unmarked-3d-ring-alkene-10-v2000.mol"),
        "C1=C/CCCCCCCC/1",
    ),
    (
        "V3000 ten-membered ring alkene",
        include_str!("fixtures/chematic-unmarked-3d-ring-alkene-10-v3000.mol"),
        "C1=C/CCCCCCCC/1",
    ),
];

const ADDITIONAL_CASES: [(&str, &str, &str); 9] = [
    (
        "V2000 phosphorus oxide",
        include_str!("fixtures/chematic-unmarked-3d-phosphorus-v2000.mol"),
        "C[P@](=O)(O)c1ccccc1",
    ),
    (
        "V3000 phosphorus oxide",
        include_str!("fixtures/chematic-unmarked-3d-phosphorus-v3000.mol"),
        "C[P@](=O)(O)c1ccccc1",
    ),
    (
        "V2000 paired four-membered ring centers",
        include_str!("fixtures/chematic-unmarked-3d-paired-ring-4-v2000.mol"),
        "C[C@]1(O)C[C@](C)(O)C1",
    ),
    (
        "V3000 paired four-membered ring centers",
        include_str!("fixtures/chematic-unmarked-3d-paired-ring-4-v3000.mol"),
        "C[C@]1(O)C[C@](C)(O)C1",
    ),
    (
        "V3000 dative bond adjacent to center",
        include_str!("fixtures/chematic-unmarked-3d-dative-v3000.mol"),
        "O=C(O)->[C@H](CC#[C][Sn]1=[CH]C=CC=[C]1F)NCP(=O)(O)O",
    ),
    (
        "V2000 paired six-membered ring centers",
        include_str!("fixtures/chematic-unmarked-3d-paired-ring-6-v2000.mol"),
        "C[C@]1(N)CC[C@@](C)(O)CC1",
    ),
    (
        "V3000 paired six-membered ring centers",
        include_str!("fixtures/chematic-unmarked-3d-paired-ring-6-v3000.mol"),
        "C[C@]1(N)CC[C@@](C)(O)CC1",
    ),
    (
        "V2000 paired fluorine/chlorine ring centers",
        include_str!("fixtures/chematic-unmarked-3d-paired-ring-4-fcl-v2000.mol"),
        "C[C@]1(F)C[C@](C)(Cl)C1",
    ),
    (
        "V3000 paired fluorine/chlorine ring centers",
        include_str!("fixtures/chematic-unmarked-3d-paired-ring-4-fcl-v3000.mol"),
        "C[C@]1(F)C[C@](C)(Cl)C1",
    ),
];

#[test]
fn unmarked_3d_geometry_recovers_stereo() {
    for &(name, block, expected) in CASES
        .iter()
        .chain(RING_CASES.iter())
        .chain(ADDITIONAL_CASES.iter())
    {
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
fn symmetric_multibridgehead_cage_does_not_gain_stereo() {
    for block in [
        include_str!("fixtures/chematic-unmarked-3d-adamantane-v2000.mol"),
        include_str!("fixtures/chematic-unmarked-3d-adamantane-v3000.mol"),
    ] {
        let mol = if block.contains("V3000") {
            read_mol_v3000_with_diagnostics(block).unwrap().mol
        } else {
            read_mol_with_diagnostics(block).unwrap().mol
        };
        assert!(!write(&mol).contains('@'));
        assert!(!canonical_smiles(&mol).contains('@'));
    }
}

#[test]
fn one_substituted_ring_center_does_not_gain_stereo() {
    for block in [
        include_str!("fixtures/chematic-unmarked-3d-single-ring-center-v2000.mol"),
        include_str!("fixtures/chematic-unmarked-3d-single-ring-center-v3000.mol"),
    ] {
        let mol = if block.contains("V3000") {
            read_mol_v3000_with_diagnostics(block).unwrap().mol
        } else {
            read_mol_with_diagnostics(block).unwrap().mol
        };
        assert!(!write(&mol).contains('@'));
        assert!(!canonical_smiles(&mol).contains('@'));
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
