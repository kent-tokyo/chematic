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
