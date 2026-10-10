use crate::pharmacophore::*;
use chematic_core::AtomIdx;
#[test]
fn public_features_identify_charges_explicit_thiol_hydrogens_and_aromatic_membership() {
    for (source, atom, kind) in [
        ("N", 0, FeatureType::Donor),
        ("[NH4+]", 0, FeatureType::Positive),
        ("[O-]", 0, FeatureType::Negative),
        ("[H]S[H]", 1, FeatureType::Donor),
        ("S", 0, FeatureType::Acceptor),
        ("CCl", 1, FeatureType::Hydrophobic),
        ("CBr", 1, FeatureType::Hydrophobic),
        ("CI", 1, FeatureType::Hydrophobic),
    ] {
        let mol = chematic_smiles::parse(source).unwrap();
        let features = detect_features(&mol);
        assert!(
            features
                .iter()
                .any(|f| f.atom == AtomIdx(atom) && f.ftype == kind),
            "{source}"
        );
        assert!(
            features
                .iter()
                .all(|f| (f.atom.0 as usize) < mol.atom_count())
        );
    }
    let mol = chematic_smiles::parse("c1ccccc1.O").unwrap();
    for feature in detect_features(&mol)
        .iter()
        .filter(|f| f.ftype == FeatureType::Aromatic)
    {
        assert_eq!(feature.neighbors, (0..6).map(AtomIdx).collect::<Vec<_>>());
    }
}
#[test]
fn feature_bit_encoding_retains_order_and_caps_at_eight_features() {
    let kinds = [
        FeatureType::Donor,
        FeatureType::Acceptor,
        FeatureType::Aromatic,
        FeatureType::Hydrophobic,
        FeatureType::Positive,
        FeatureType::Negative,
    ];
    let mut features = Vec::new();
    assert_eq!(features_to_bitvec(&features), 0);
    for kind in kinds {
        features.push(Feature {
            ftype: kind,
            atom: AtomIdx(0),
            neighbors: vec![],
        });
    }
    assert_eq!(
        features_to_bitvec(&features),
        1 | (2 << 3) | (3 << 6) | (4 << 9) | (5 << 12) | (6 << 15)
    );
    features.extend(features.clone());
    assert_eq!(
        features_to_bitvec(&features),
        features_to_bitvec(&features[..8])
    );
    assert_eq!(features_to_bitvec(&features) >> 24, 0);
}
