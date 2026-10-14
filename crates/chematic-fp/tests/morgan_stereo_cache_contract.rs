use chematic_core::{AtomIdx, BondIdx, BondOrder, Chirality};
use chematic_fp::{RdkitMorganConfig, RdkitMorganFpSize, RdkitMorganRadius, rdkit_morgan_bitvec};

fn config() -> RdkitMorganConfig {
    RdkitMorganConfig {
        include_chirality: true,
        ..Default::default()
    }
}

#[test]
fn prepared_stereo_is_independent_of_radius_width_and_plain_calls() {
    let mol = chematic_smiles::parse("C/C=C/[C@H](O)C(=O)O").unwrap();
    for radius in [
        RdkitMorganRadius::R0,
        RdkitMorganRadius::R1,
        RdkitMorganRadius::R2,
        RdkitMorganRadius::R3,
    ] {
        for fp_size in [
            RdkitMorganFpSize::B128,
            RdkitMorganFpSize::B512,
            RdkitMorganFpSize::B2048,
        ] {
            for include_chirality in [true, false, true] {
                let cfg = RdkitMorganConfig {
                    radius,
                    fp_size,
                    include_chirality,
                };
                assert_eq!(
                    rdkit_morgan_bitvec(&mol, &cfg).unwrap(),
                    rdkit_morgan_bitvec(&mol.clone(), &cfg).unwrap()
                );
            }
        }
    }
}

#[test]
fn chirality_and_direction_edits_do_not_reuse_old_preparation() {
    let mut mol = chematic_smiles::parse("C[C@H](N)C(=O)O").unwrap();
    let before = rdkit_morgan_bitvec(&mol, &config()).unwrap();
    mol.set_chirality(AtomIdx(1), Chirality::Clockwise);
    let after = rdkit_morgan_bitvec(&mol, &config()).unwrap();
    assert_ne!(before, after);
    assert_eq!(after, rdkit_morgan_bitvec(&mol.clone(), &config()).unwrap());

    let mut alkene = chematic_smiles::parse("F/C=C/F").unwrap();
    let before = rdkit_morgan_bitvec(&alkene, &config()).unwrap();
    alkene.set_bond_order(BondIdx(2), BondOrder::Down);
    let after = rdkit_morgan_bitvec(&alkene, &config()).unwrap();
    assert_ne!(before, after);
    assert_eq!(
        after,
        rdkit_morgan_bitvec(&alkene.clone(), &config()).unwrap()
    );
}
