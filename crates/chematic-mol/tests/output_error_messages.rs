use chematic_mol::XyzWriteOptions;
use chematic_mol::gro_lammps::{GroLammpsUnits, gro_to_lammps_atomic};

#[test]
fn output_refusals_explain_bad_mass_and_xyz_precision() {
    let error =
        gro_to_lammps_atomic("empty\n0\n2 3 4\n", &[], &[], GroLammpsUnits::Real).unwrap_err();
    assert_eq!(
        error.to_string(),
        "supply 1 to 4096 positive finite masses in grams/mole"
    );
    let options = XyzWriteOptions {
        decimal_places: 18,
        ..Default::default()
    };
    let error = options.validate("geometry").unwrap_err();
    assert!(error.to_string().contains("decimal_places must be 0..=17"));
}

#[test]
fn oversized_grid_error_reports_count_and_limit() {
    let input = "object 1 class gridpositions counts 10000 10000 10000\n";
    let error = chematic_mol::opendx::parse_opendx(input).unwrap_err();
    let message = error.to_string();
    assert!(message.contains("1000000000000 points"), "{message}");
    assert!(message.contains("100000000-point limit"), "{message}");
}

#[test]
fn gaussian_writer_error_reports_invalid_electronic_state() {
    use chematic_core::Element;
    use chematic_mol::{GjfCoordinateUnits, GjfWriteOptions, write_gjf_cartesian};
    let options = GjfWriteOptions::new(0, 0, GjfCoordinateUnits::Angstrom);
    let error = write_gjf_cartesian(&[Element::C], &[[0.; 3]], &[None], &options).unwrap_err();
    assert_eq!(
        error.to_string(),
        "GJF Cartesian output: multiplicity must be positive"
    );
}
