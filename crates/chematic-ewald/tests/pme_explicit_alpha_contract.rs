use chematic_ewald::{BoxVectors, EwaldError, PmeConfig, spme_energy};

#[test]
fn explicit_automatic_alpha_and_zero_charge_particle_agree() {
    let coords = [[1.0, 1.0, 1.0], [3.0, 1.0, 1.0]];
    let charges = [1.0, -1.0];
    let box_vectors = BoxVectors::cubic(10.0);
    let auto = PmeConfig {
        r_cut: 5.0,
        mesh: [8; 3],
        ..Default::default()
    };
    let fixed = PmeConfig {
        alpha: 3.5 / auto.r_cut,
        ..auto.clone()
    };
    let a = spme_energy(&coords, &charges, &box_vectors, &auto).unwrap();
    let b = spme_energy(&coords, &charges, &box_vectors, &fixed).unwrap();
    assert_eq!(a, b);
    let extended = [[1.0, 1.0, 1.0], [3.0, 1.0, 1.0], [7.0, 6.0, 5.0]];
    let with_zero = spme_energy(&extended, &[1.0, -1.0, 0.0], &box_vectors, &fixed).unwrap();
    assert_eq!(b, with_zero);
}

#[test]
fn excessive_grid_is_rejected_before_allocation() {
    let config = PmeConfig {
        mesh: [4097, 4097, 1],
        ..Default::default()
    };
    assert_eq!(
        spme_energy(&[[0.0; 3]], &[1.0], &BoxVectors::cubic(10.0), &config),
        Err(EwaldError::SingularBoxMatrix)
    );
}
