use chematic_3d::torsion_motif::{TorsionHistogram, fit_von_mises_mixture};

#[test]
fn empty_torsion_fit_returns_zero_density_bins() {
    let histogram = TorsionHistogram::from_angles_deg(&[], 12);
    let fit = fit_von_mises_mixture(&histogram, 2, 10);
    assert!(fit.components.is_empty());
    assert!(fit.probabilities(0).is_empty());
    assert_eq!(fit.probabilities(8), vec![0.0; 8]);
}
