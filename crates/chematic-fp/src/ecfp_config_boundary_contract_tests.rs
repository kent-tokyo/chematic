use crate::*;

#[test]
fn ecfp_bit_origins_follow_size_radius_and_folding_configuration() {
    for smiles in ["", "CCO", "c1ccccc1", "N[C@@H](C)C(=O)O", "[13CH3][NH3+]"] {
        let mol = if smiles.is_empty() {
            chematic_core::MoleculeBuilder::new().build()
        } else {
            chematic_smiles::parse(smiles).unwrap()
        };
        for mode in [EcfpInvariantMode::Chematic, EcfpInvariantMode::RdkitMorgan] {
            for nbits in [1, 64, 2048] {
                for radius in [0, 2, 3, 21] {
                    let mut config = EcfpConfig {
                        radius,
                        nbits,
                        use_chirality: true,
                        use_double_fold: false,
                    };
                    let single = ecfp_with_invariant_mode(&mol, &config, mode);
                    for double_fold in [false, true] {
                        config.use_double_fold = double_fold;
                        let (fp, info) = ecfp_with_bitinfo_and_mode(&mol, &config, mode);
                        assert_eq!(fp, ecfp_with_invariant_mode(&mol, &config, mode));
                        assert_eq!(fp.popcount() as usize, info.len());
                        assert_eq!(single.and(&fp), single);
                        for bit in 0..2048 {
                            assert_eq!(fp.get(bit), info.contains_key(&bit));
                            if bit >= nbits {
                                assert!(!fp.get(bit));
                            }
                        }
                        for origins in info.values() {
                            assert!(!origins.is_empty());
                            for &(atom, r) in origins {
                                assert!((atom as usize) < mol.atom_count());
                                assert!(r <= radius.min(ecfp::MAX_ECFP_RADIUS));
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn ecfp_public_shortcuts_match_explicit_configuration() {
    let mol = chematic_smiles::parse("CC(=O)Nc1ccccc1").unwrap();
    let four = EcfpConfig::default();
    let six = EcfpConfig {
        radius: 3,
        ..four.clone()
    };
    assert_eq!(ecfp4(&mol), ecfp(&mol, &four));
    assert_eq!(ecfp6(&mol), ecfp(&mol, &six));
    assert_eq!(ecfp_with_bitinfo(&mol, &four).0, ecfp4(&mol));
    assert_eq!(
        ecfp4_rdkit_invariants(&mol),
        ecfp_with_invariant_mode(&mol, &four, EcfpInvariantMode::RdkitMorgan)
    );
    assert_eq!(
        ecfp6_rdkit_invariants(&mol),
        ecfp_with_invariant_mode(&mol, &six, EcfpInvariantMode::RdkitMorgan)
    );
    assert_eq!(
        ecfp4_rdkit_environment_experimental(&mol),
        ecfp_with_bitinfo_rdkit_environment_experimental(&mol, &four).0
    );
    assert_eq!(
        ecfp6_rdkit_environment_experimental(&mol),
        ecfp_with_bitinfo_rdkit_environment_experimental(&mol, &six).0
    );
    assert_eq!(tanimoto_ecfp4(&mol, &mol), 1.0);
}
