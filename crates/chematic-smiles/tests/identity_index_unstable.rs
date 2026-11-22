use chematic_smiles::{BatchCanonicalization, SmilesBatchCanonicalizer};

#[test]
fn identity_index_rejects_unstable_stereochemistry_without_losing_other_records() {
    let unstable = r"CC1CNC(/C=C\C=C/C=C\C2(C)C(=C(O)C(C2)=O)C(/C=C\C=C/C=C\C=C/C=C1)=O)=O";
    let build = SmilesBatchCanonicalizer::default().build_identity_index(["CCO", unstable, "OCC"]);
    assert_eq!(build.record_count, 3);
    assert_eq!(build.index.record_count(), 2);
    assert_eq!(build.index.unique_key_count(), 1);
    assert_eq!(build.rejected.len(), 1);
    assert_eq!(build.rejected[0].input_index, 1);
    assert_eq!(build.rejected[0].input, unstable);
    assert!(
        matches!(&build.rejected[0].result, BatchCanonicalization::Rejected { error } if error.contains("not stable"))
    );
}
