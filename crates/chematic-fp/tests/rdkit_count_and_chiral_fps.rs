//! RDKit 2026.03.1 reference values for the Morgan `includeChirality` and
//! `countSimulation` fingerprints and the hashed atom-pair/torsion count
//! fingerprints.

use chematic_fp::{
    RdkitMorganConfig, rdkit_atom_pair_counts, rdkit_morgan_count_simulation,
    rdkit_morgan_fingerprint, rdkit_torsion_counts,
};

fn on_bits(fp: &chematic_fp::BitVecN) -> Vec<usize> {
    (0..fp.bit_width()).filter(|&i| fp.get(i)).collect()
}

fn chiral() -> RdkitMorganConfig {
    RdkitMorganConfig {
        include_chirality: true,
        ..RdkitMorganConfig::default()
    }
}

#[test]
fn morgan_include_chirality_covers_ez_and_tetrahedral() {
    // GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    for (smiles, want) in [
        (
            "F/C=C/C[C@H](N)O",
            vec![
                1, 53, 80, 227, 694, 724, 783, 786, 807, 894, 1002, 1157, 1171, 1589, 1649, 1928,
                1938,
            ],
        ),
        (
            "F/C=C\\C[C@@H](N)O",
            vec![
                1, 50, 80, 227, 532, 694, 699, 786, 797, 807, 938, 1113, 1157, 1171, 1649, 1671,
                1928,
            ],
        ),
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        let fp = rdkit_morgan_fingerprint(&mol, &chiral()).unwrap();
        assert_eq!(on_bits(&fp.fingerprint), want, "{smiles}");
    }
}

#[test]
fn morgan_count_simulation_matches_rdkit() {
    // GetMorganGenerator(radius=2, fpSize=2048, countSimulation=True)
    let mol = chematic_smiles::parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
    let fp = rdkit_morgan_count_simulation(&mol, &RdkitMorganConfig::default()).unwrap();
    assert_eq!(
        on_bits(&fp),
        vec![
            44, 92, 132, 256, 257, 320, 552, 553, 700, 701, 732, 772, 856, 857, 956, 1180, 1181,
            1348, 1349, 1350, 1424, 1425, 1524, 1525, 1544, 1556, 1588, 1692, 1736, 1776, 1820,
            1824, 2020,
        ]
    );
}

#[test]
fn hashed_atom_pair_and_torsion_counts_match_rdkit() {
    let mol = chematic_smiles::parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
    assert_eq!(
        rdkit_atom_pair_counts(&mol, 2048),
        vec![
            (4, 1),
            (71, 2),
            (74, 1),
            (94, 1),
            (95, 1),
            (115, 1),
            (124, 1),
            (134, 1),
            (137, 1),
            (138, 1),
            (158, 1),
            (159, 1),
            (178, 3),
            (191, 3),
            (200, 2),
            (201, 2),
            (241, 3),
            (246, 1),
            (254, 6),
            (255, 1),
            (270, 2),
            (311, 2),
            (312, 5),
            (335, 1),
            (374, 3),
            (377, 2),
            (453, 2),
            (477, 1),
            (704, 2),
            (711, 2),
            (728, 1),
            (776, 1),
            (816, 1),
            (823, 2),
            (1041, 1),
            (1104, 2),
            (1105, 1),
            (1136, 1),
            (1234, 1),
            (1242, 1),
            (1266, 1),
            (1352, 1),
            (1364, 1),
            (1495, 1),
            (1502, 1),
            (1684, 1),
            (1687, 1),
            (1731, 1),
            (1750, 1),
            (1990, 1),
        ]
    );
    assert_eq!(
        rdkit_torsion_counts(&mol, 2048),
        vec![
            (7, 1),
            (94, 3),
            (102, 1),
            (213, 1),
            (409, 1),
            (414, 1),
            (428, 1),
            (429, 1),
            (467, 1),
            (468, 1),
            (491, 1),
            (805, 1),
            (1424, 1),
            (1539, 1),
            (1548, 2),
            (1804, 1),
        ]
    );
}
