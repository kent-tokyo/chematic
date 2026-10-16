use chematic_fp::feature_patterns::FeaturePatternSet;

#[test]
fn masks_use_caller_order_and_all_matched_atoms() {
    let mol = chematic_smiles::parse("CC(=O)O").unwrap();
    let set = FeaturePatternSet::new(["[#8]", "C=O", "[#6]"]).unwrap();
    assert_eq!(set.invariants(&mol, 100_000).unwrap(), [4, 6, 3, 1]);
    assert_eq!(set.patterns(), ["[#8]", "C=O", "[#6]"]);
}

#[test]
fn all_32_bits_are_preserved_and_the_33rd_is_rejected() {
    let mol = chematic_smiles::parse("CO").unwrap();
    let set = FeaturePatternSet::new(["O"; 32]).unwrap();
    assert_eq!(set.invariants(&mol, 100_000).unwrap(), [0, u32::MAX]);
    assert!(FeaturePatternSet::new(["O"; 33]).is_err());
    let mut consumed = 0;
    let source = std::iter::repeat_with(|| {
        consumed += 1;
        "O"
    });
    assert!(FeaturePatternSet::new(source).is_err());
    assert_eq!(consumed, 33);
}

#[test]
fn empty_definitions_mean_zero_features_not_defaults() {
    let set = FeaturePatternSet::new(std::iter::empty::<&str>()).unwrap();
    let mol = chematic_smiles::parse("CCO").unwrap();
    assert_eq!(set.invariants(&mol, 100_000).unwrap(), [0, 0, 0]);
    let empty = chematic_core::MoleculeBuilder::new().build();
    assert!(set.invariants(&empty, 100_000).unwrap().is_empty());
    assert_eq!(
        set.fingerprint(&empty, 2, 2048, 100_000)
            .unwrap()
            .popcount(),
        0
    );
}

#[test]
fn changed_patterns_affect_native_fingerprint_and_upper_bits() {
    let mol = chematic_smiles::parse("CCO").unwrap();
    let oxygen = FeaturePatternSet::new(["O"]).unwrap();
    let carbon = FeaturePatternSet::new(["C"]).unwrap();
    let mut upper = vec!["N"; 31];
    upper.push("O");
    let high = FeaturePatternSet::new(upper).unwrap();
    assert_ne!(
        oxygen.fingerprint(&mol, 2, 2048, 100_000).unwrap(),
        carbon.fingerprint(&mol, 2, 2048, 100_000).unwrap()
    );
    assert_eq!(high.invariants(&mol, 100_000).unwrap(), [0, 0, 1u32 << 31]);
    assert_ne!(
        oxygen.fingerprint(&mol, 2, 2048, 100_000).unwrap(),
        high.fingerprint(&mol, 2, 2048, 100_000).unwrap()
    );
}

#[test]
fn aromatic_and_kekule_forms_share_masks_and_fingerprint() {
    let set = FeaturePatternSet::new(["[a]", "[#7]", "[#6]"]).unwrap();
    let aromatic = chematic_smiles::parse("n1ccccc1").unwrap();
    let kekule = chematic_smiles::parse("N1=CC=CC=C1").unwrap();
    assert_eq!(
        set.invariants(&aromatic, 100_000).unwrap(),
        set.invariants(&kekule, 100_000).unwrap()
    );
    assert_eq!(
        set.fingerprint(&aromatic, 2, 2048, 100_000).unwrap(),
        set.fingerprint(&kekule, 2, 2048, 100_000).unwrap()
    );
}

#[test]
fn budgets_and_invalid_profiles_refuse_partial_results() {
    let mol = chematic_smiles::parse("CCCCCC").unwrap();
    let set = FeaturePatternSet::new(["CC"]).unwrap();
    assert!(set.invariants(&mol, 1).is_err());
    for budget in [0, 1_000_001] {
        assert!(set.invariants(&mol, budget).is_err());
    }
    assert!(set.fingerprint(&mol, 4, 2048, 100_000).is_err());
    assert!(set.fingerprint(&mol, 2, 4096, 100_000).is_err());
    for text in ["", "[invalid", &"C".repeat(33), &"C".repeat(4097)] {
        assert!(FeaturePatternSet::new([text]).is_err());
    }
    let big = chematic_smiles::parse(&"C".repeat(257)).unwrap();
    assert!(set.invariants(&big, 100_000).is_err());
}
