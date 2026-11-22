use chematic_smiles::{RdkitSmilesError, parse, rdkit_stereoisomer_count};

fn independent_stereocenters(count: usize) -> String {
    std::iter::repeat_n("FC(Cl)(Br)I", count)
        .collect::<Vec<_>>()
        .join(".")
}

#[test]
fn stereoisomer_count_checks_u128_boundary_without_wrapping() {
    // Each component has exactly one independent unspecified tetrahedral center.
    let within = parse(&independent_stereocenters(127)).unwrap();
    assert_eq!(rdkit_stereoisomer_count(&within).unwrap(), 1u128 << 127);
    let beyond = parse(&independent_stereocenters(128)).unwrap();
    let error = rdkit_stereoisomer_count(&beyond).unwrap_err();
    assert!(
        matches!(error, RdkitSmilesError::Unsupported(ref message) if message.contains("2^128"))
    );
}
