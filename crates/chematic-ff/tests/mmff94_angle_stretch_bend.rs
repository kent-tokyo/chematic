use chematic_ff::{mmff94_angle_energy, mmff94_stbn};

// RDKit 2026.03.6 GetMMFFAngleBendParams on the benzyl-pyridinium
// motif in frozen corpus inputs 89 and 96. Type 58 must substitute
// through its defined equivalence levels, rather than remain 58.
#[test]
fn aromatic_cation_equivalence_selects_the_correct_angle_parameters() {
    for (i, k) in [(37, 58), (58, 37)] {
        let params = mmff94_angle_energy(0, i, 1, k).unwrap();
        assert!((params.ka - 1.16).abs() < 1e-12);
        assert!((params.theta0 - 107.963).abs() < 1e-12);
    }
}

// RDKit 2026.03.6 GetMMFFStretchBendParams at biphenyl's connecting
// single bond: type 1 has no specific 37-37-37 row. The periodic-row
// default is +0.3, not the aromatic-ring type-0 value -0.411.
#[test]
fn biphenyl_stretch_bend_uses_periodic_row_default() {
    assert_eq!(mmff94_stbn(1, 37, 37, 37, 6, 6, 6), Some((0.3, 0.3)));
    assert_eq!(mmff94_stbn(0, 37, 37, 37, 6, 6, 6), Some((-0.411, -0.411)));
}
