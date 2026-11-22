use chematic_crystal::{CrystalError, Lattice, PeriodicStructure};

#[test]
fn huge_finite_cutoff_reports_bounded_search_without_overflow() {
    let structure = PeriodicStructure::new(Lattice::cubic(10.0).unwrap(), vec![]).unwrap();
    let error = structure.neighbors_within(f64::MAX).unwrap_err();
    assert!(matches!(
        error,
        CrystalError::NeighborSearchTooLarge {
            candidate_count: u64::MAX,
            ..
        }
    ));
}
