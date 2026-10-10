use crate::{
    rdkit_atom_pair::rdkit_atom_pair_sparse_counts,
    rdkit_torsion::{rdkit_legacy_torsion_counts, rdkit_torsion_sparse_counts},
};

#[test]
fn fallible_similarity_search_keeps_identical_ties_stable_and_respects_limits() {
    use crate::search::{FpType, PreparedFingerprintIndex, try_nearest_neighbors};
    let db: Vec<_> = ["CCO", "CCO", "CCN", "c1ccccc1"]
        .iter()
        .map(|s| chematic_smiles::parse(s).unwrap())
        .collect();
    for profile in [
        FpType::Ecfp4,
        FpType::Ecfp6,
        FpType::Ecfp4Chiral,
        FpType::Fcfp4,
        FpType::Maccs,
        FpType::TopoPath,
        FpType::RdkitEcfp4,
    ] {
        let expected = vec![(0, 1.0), (1, 1.0)];
        assert_eq!(
            try_nearest_neighbors(&db[0], &db, 2, profile).unwrap(),
            expected
        );
        assert!(
            try_nearest_neighbors(&db[0], &db, 0, profile)
                .unwrap()
                .is_empty()
        );
        assert!(
            try_nearest_neighbors(&db[0], &[], 2, profile)
                .unwrap()
                .is_empty()
        );
        let prepared = PreparedFingerprintIndex::try_new(&db, profile).unwrap();
        assert_eq!(
            prepared.try_search_threshold(&db[0], 1.0, 10).unwrap(),
            expected
        );
        assert!(
            prepared
                .try_search_threshold(&db[0], 1.0, 0)
                .unwrap()
                .is_empty()
        );
        let all = prepared.try_search_threshold(&db[0], 0.0, 10).unwrap();
        assert_eq!(all.len(), 4);
        assert!(all.windows(2).all(|pair| pair[0].1 >= pair[1].1));
    }
}

#[test]
fn sparse_count_fingerprint_codes_and_multiplicities_match_pinned_rdkit() {
    let fixture: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-sparse-fp-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    for row in fixture["rows"].as_array().unwrap() {
        let text = row["smiles"].as_str().unwrap();
        let mol = if text.is_empty() {
            chematic_core::MoleculeBuilder::new().build()
        } else {
            chematic_smiles::parse(text).unwrap()
        };
        assert_eq!(
            serde_json::to_value(rdkit_atom_pair_sparse_counts(&mol)).unwrap(),
            row["atom_pair"],
            "atom pair: {text}"
        );
        assert_eq!(
            serde_json::to_value(rdkit_torsion_sparse_counts(&mol)).unwrap(),
            row["torsion"],
            "torsion: {text}"
        );
        assert_eq!(
            serde_json::to_value(rdkit_legacy_torsion_counts(&mol)).unwrap(),
            row["legacy_torsion"],
            "legacy torsion: {text}"
        );
    }
}

#[test]
fn fps_errors_explain_headers_hex_lengths_and_bad_digits_before_stopping() {
    use crate::fps::*;
    use std::io::Cursor;
    for (text, context) in [
        ("01\tname\n", "num_bits"),
        ("#num_bits=zero\n", "zero"),
        ("#num_bits=0\n", "num_bits"),
    ] {
        let error = FpsReader::new(Cursor::new(text)).err().unwrap();
        assert!(error.to_string().contains(context));
    }
    for (text, context) in [
        ("missing-tab", "missing-tab"),
        ("0\tname", "expected 2"),
        ("gg\tname", "'g'"),
    ] {
        let mut reader =
            FpsReader::new(Cursor::new(format!("#num_bits=8\n{text}\n01\tgood\n"))).unwrap();
        let error = reader.next().unwrap().unwrap_err();
        assert!(error.to_string().contains(context), "{error}");
        assert!(reader.next().is_none());
    }
}

#[test]
fn prepared_morgan_rejects_query_only_bond_types_with_the_original_index() {
    use crate::rdkit_morgan_ecfp4::*;
    use chematic_core::{BondIdx, BondOrder};
    for order in [
        BondOrder::QueryAny,
        BondOrder::QuerySingleOrDouble,
        BondOrder::QuerySingleOrAromatic,
        BondOrder::QueryDoubleOrAromatic,
    ] {
        let mut mol = chematic_smiles::parse("CCC").unwrap();
        mol.set_bond_order(BondIdx(1), order);
        let e = prepare_rdkit_morgan_ecfp4(&mol).err().unwrap();
        assert!(matches!(
            e,
            RdkitMorganError::UnsupportedBondOrder {
                bond_idx: BondIdx(1),
                ..
            }
        ));
        assert!(
            e.to_string().contains("BondIdx(1)") && e.to_string().contains(&format!("{order:?}"))
        );
    }
}
