use chematic_fp::rdkit_legacy_torsion_counts;
use chematic_smiles::parse;

#[test]
fn legacy_torsion_counts_ignore_removable_explicit_hydrogens() {
    // RDKit 2025.09.3 GetTopologicalTorsionFingerprint on hydrogen-suppressed graphs.
    for (plain, explicit, expected) in [
        ("CCCO", "[H]CCCO", vec![(12_893_306_912, 1)]),
        ("CCCC", "[H]CCCC", vec![(4_303_372_320, 1)]),
        ("CCNCC", "[H]CCNCC", vec![(4_445_978_656, 2)]),
    ] {
        assert_eq!(
            rdkit_legacy_torsion_counts(&parse(plain).unwrap()),
            expected,
            "{plain}"
        );
        assert_eq!(
            rdkit_legacy_torsion_counts(&parse(explicit).unwrap()),
            expected,
            "{explicit}"
        );
    }
}
