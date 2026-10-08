//! Depiction and MOL block outputs pinned to RDKit 2026.03.1
//! (`rdDepictor.Compute2DCoords(Chem.MolFromSmiles(s))`, `Chem.MolToMolBlock`).

use crate::{parse, rdkit_2d_coords, rdkit_mol_block_2d};

fn coords(smiles: &str) -> Vec<[f64; 2]> {
    rdkit_2d_coords(&parse(smiles).unwrap()).unwrap()
}

/// Bit patterns of every coordinate RDKit assigns.
#[test]
fn small_molecules_are_bit_exact() {
    let cases: &[(&str, &[[u64; 2]])] = &[
        (
            "CCO",
            &[
                [0xbff4c8dc2e423980, 0xbfcffffffffffff5],
                [0x0000000000000000, 0x3fe0000000000001],
                [0x3ff4c8dc2e423980, 0xbfd0000000000006],
            ],
        ),
        (
            "c1ccc2ccccc2c1",
            &[
                [0xc004c8dc2e42397f, 0xbfe8000000000004],
                [0xc004c8dc2e42397f, 0x3fe8000000000004],
                [0xbff4c8dc2e42397f, 0x3ff8000000000000],
                [0xbcb0000000000000, 0x3fe8000000000002],
                [0x3ff4c8dc2e423983, 0x3ff7ffffffffffff],
                [0x4004c8dc2e423980, 0x3fe7fffffffffffb],
                [0x4004c8dc2e423980, 0xbfe8000000000006],
                [0x3ff4c8dc2e42397a, 0xbff8000000000000],
                [0x0000000000000000, 0xbfe7ffffffffffff],
                [0xbff4c8dc2e42397f, 0xbff8000000000000],
            ],
        ),
        (
            "C/C=C/C(=O)O",
            &[
                [0xc006ff820c828bed, 0xbfcc1e3fb4f1dc38],
                [0xbff7b8637542fa81, 0x3fd5ab7b669c1b8c],
                [0xbfd3631b00c06309, 0xbfe2d0d9120eea82],
                [0x3ff16dd9e3920497, 0xbf9e716e308cb5c0],
                [0x3ff4d4ddd241402f, 0x3ff7482c7a318bb8],
                [0x400226bb4c52732a, 0xbfee9a2236e15df7],
            ],
        ),
        (
            "C1CC2CCC1C2",
            &[
                [0x3ff8d732e136ae00, 0xbfe8000000000000],
                [0x3ff8d732e136ae00, 0x3fe7fffffffffffe],
                [0x3fc01f42961b720a, 0x3ff36a99b4b1f77e],
                [0xbff4cf623bafd17f, 0x3fe8000000000001],
                [0xbff4cf623bafd17e, 0xbfe8000000000001],
                [0x3fc01f42961b7200, 0xbff36a99b4b1f77e],
                [0xbfe82ee3e1292b0b, 0x3ca8000000000000],
            ],
        ),
        (
            "CC(C)(C)C(C)(C)C",
            &[
                [0xbfe8000000000000, 0x3ff8000000000000],
                [0xbfe8000000000002, 0x3ca0000000000000],
                [0xbfe8000000000002, 0xbff8000000000000],
                [0xc002000000000001, 0x3cb0000000000000],
                [0x3fe8000000000000, 0x0000000000000000],
                [0x3fe8000000000000, 0xbff8000000000000],
                [0x3fe8000000000002, 0x3ff8000000000000],
                [0x4002000000000000, 0x0000000000000000],
            ],
        ),
        (
            "C1CCC2(CC1)CCCC2",
            &[
                [0x4005e7664ee9ed89, 0xbc9999999999999a],
                [0x3fffcecc9dd3db0d, 0xbff4c8dc2e423983],
                [0x3fdf3b32774f6c36, 0xbff4c8dc2e42397f],
                [0xbfd0c4cd88b093c0, 0x3c9b58d8ffa37a7e],
                [0x3fdf3b32774f6c48, 0x3ff4c8dc2e423981],
                [0x3fffcecc9dd3db12, 0x3ff4c8dc2e423980],
                [0xbff24c8da58428b7, 0xbff36a99b4b1f77c],
                [0xc0048fec19fbb43c, 0xbfe7fffffffffffd],
                [0xc0048fec19fbb43a, 0x3fe8000000000001],
                [0xbff24c8da58428b5, 0x3ff36a99b4b1f77e],
            ],
        ),
        (
            "CC.O",
            &[
                [0xbfe8000000000000, 0x3c90000000000000],
                [0x3fe8000000000000, 0xbc90000000000000],
                [0x0000000000000000, 0x3ff0000000000000],
            ],
        ),
    ];
    for (smiles, expected) in cases {
        let got = coords(smiles);
        let got: Vec<[u64; 2]> = got
            .iter()
            .map(|p| [p[0].to_bits(), p[1].to_bits()])
            .collect();
        assert_eq!(&got[..], *expected, "{smiles}");
    }
}

fn fnv(coords: &[[f64; 2]]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for p in coords {
        for v in p {
            h ^= v.to_bits();
            h = h.wrapping_mul(0x0100_0000_01b3);
        }
    }
    h
}

/// Fused, bridged and crowded systems (ring merging, collision removal by
/// bond flips, angle opening and bond shortening): a hash of the bits of
/// every coordinate RDKit assigns.
#[test]
fn larger_molecules_are_bit_exact() {
    let cases: &[(&str, usize, u64)] = &[
        (
            "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C",
            62,
            0xfdb92b7834fa8a27,
        ),
        (
            "c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41",
            60,
            0xb863c165599e8206,
        ),
        (
            "CC(C)(C)C(C(C)(C)C)(C(C)(C)C)C(C)(C)C",
            17,
            0x453b022a5e60cc1e,
        ),
        ("C1CCCCCCCCCCCCCCCCCCCCCCC1", 24, 0xdde03182bade3379),
        (
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            48,
            0x4b8ef7113ea9f4e8,
        ),
        ("C1CCCCC/C=C/CCCC1", 12, 0x92ae30f3501cc442),
    ];
    for &(smiles, n, hash) in cases {
        let got = coords(smiles);
        assert_eq!(got.len(), n, "{smiles}");
        assert_eq!(fnv(&got), hash, "{smiles}");
    }
}

/// Wedges, crossed double bonds, charges, isotopes and a ring closure
/// whose bond symbol is written at the closing digit.
#[test]
fn mol_blocks_match_rdkit() {
    let cases: &[(&str, &str)] = &[
        (
            "C[C@H](O)F",
            "\n     RDKit          2D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2990   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    1.5000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n  2  1  1  1\n  2  3  1  0\n  2  4  1  0\nM  END\n",
        ),
        (
            "C[N+](C)(C)C.[Cl-]",
            "\n     RDKit          2D\n\n  6  4  0  0  0  0  0  0  0  0999 V2000\n   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.2990    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\n  2  4  1  0\n  2  5  1  0\nM  CHG  2   2   1   6  -1\nM  END\n",
        ),
        (
            "CC=CC",
            "\n     RDKit          2D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   -1.9796   -0.1365    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.5994    0.4508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.5994   -0.4508    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.9796    0.1365    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  2  3\n  3  4  1  0\nM  END\n",
        ),
        (
            "c1ccccc1/C=C/C",
            "\n     RDKit          2D\n\n  9  9  0  0  0  0  0  0  0  0999 V2000\n   -0.3293   -1.1192    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.7872   -1.4721    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.8217   -0.3860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.3984    1.0531    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9405    1.4059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0941    0.3198    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5520    0.6727    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5865   -0.4135    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.0444   -0.0606    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  2  0\n  3  4  1  0\n  4  5  2  0\n  5  6  1  0\n  6  7  1  0\n  7  8  2  0\n  8  9  1  0\n  6  1  2  0\nM  END\n",
        ),
        (
            "[13CH3]C(=O)[O-]",
            "\n     RDKit          2D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n   -1.2990   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    1.5000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    1.2990   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  2  0\n  2  4  1  0\nM  CHG  1   4  -1\nM  ISO  1   1  13\nM  END\n",
        ),
        (
            "N#Cc1cc2c(cc1[N+](=O)[O-])=NC(=O)C(=O)N=2",
            "\n     RDKit          2D\n\n 17 18  0  0  0  0  0  0  0  0999 V2000\n    4.3013   -2.4367    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0196   -1.6574    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7380   -0.8780    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4222   -1.5983    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8594   -0.8189    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8253    0.6807    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4905    1.4010    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7721    0.6216    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0879    1.3419    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n    4.3695    0.5625    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1220    2.8415    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1069    1.4601    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.4227    0.7398    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7043    1.5192    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -3.4568   -0.7598    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -4.7726   -1.4801    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -2.1752   -1.5392    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  3  0\n  2  3  1  0\n  3  4  2  0\n  4  5  1  0\n  5  6  1  0\n  6  7  1  0\n  7  8  2  0\n  8  9  1  0\n  9 10  2  0\n  9 11  1  0\n  6 12  2  0\n 12 13  1  0\n 13 14  2  0\n 13 15  1  0\n 15 16  2  0\n 15 17  1  0\n  8  3  1  0\n 17  5  2  0\nM  CHG  2   9   1  11  -1\nM  END\n",
        ),
    ];
    for (smiles, expected) in cases {
        let got = rdkit_mol_block_2d(&parse(smiles).unwrap()).unwrap();
        assert_eq!(&got, expected, "{smiles}");
    }
}

#[test]
fn v3000_and_unsupported_inputs_are_refused() {
    assert!(rdkit_mol_block_2d(&parse("[NH3]->[Pt](Cl)(Cl)<-[NH3]").unwrap()).is_err());
    assert!(rdkit_2d_coords(&parse("[Pt@SP1](Cl)(Cl)(N)N").unwrap()).is_err());
}
