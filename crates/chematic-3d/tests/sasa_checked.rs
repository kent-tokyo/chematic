use chematic_3d::{SasaInputError as Error, sasa_with_radii_checked as sasa};
use chematic_core::{Atom, Element, Molecule, MoleculeBuilder};
use std::f64::consts::PI;

fn atoms(n: usize) -> Molecule {
    let mut b = MoleculeBuilder::new();
    for _ in 0..n {
        b.add_atom(Atom::new(Element::C));
    }
    b.build()
}

#[test]
fn isolated_and_disjoint_spheres_use_supplied_radii_once() {
    let result = sasa(
        &atoms(2),
        &[[0., 0., 0.], [20., 0., 0.]],
        &[1., 2.],
        0.5,
        960,
    )
    .unwrap();
    assert!((result.per_atom[0] - 4. * PI * 1.5_f64.powi(2)).abs() < 1e-12);
    assert!((result.per_atom[1] - 4. * PI * 2.5_f64.powi(2)).abs() < 1e-12);
    assert_eq!(result.total, result.per_atom.iter().sum::<f64>());
}

#[test]
fn overlapping_equal_spheres_match_spherical_caps() {
    let result = sasa(
        &atoms(2),
        &[[0., 0., 0.], [0., 1., 0.]],
        &[1., 1.],
        0.,
        4096,
    )
    .unwrap();
    for area in &result.per_atom {
        assert!((area - 3. * PI).abs() < 1e-12);
    }
}

#[test]
fn multiple_occluders_match_nonoverlapping_spherical_caps() {
    let result = sasa(
        &atoms(3),
        &[[0., 0., 0.], [0., 1., 0.], [0., 2., 0.]],
        &[1.; 3],
        0.,
        4096,
    )
    .unwrap();
    for (area, expected) in result.per_atom.iter().zip([3. * PI, 2. * PI, 3. * PI]) {
        assert!((area - expected).abs() < 1e-12);
    }
    assert!((result.total - 8. * PI).abs() < 1e-12);
}

#[test]
fn unequal_spheres_match_analytic_cap_areas() {
    let result = sasa(
        &atoms(2),
        &[[0., 0., 0.], [2., 0., 0.]],
        &[2., 1.5],
        0.,
        32768,
    )
    .unwrap();
    let expected = [
        4. * PI * 4. - 2. * PI * 2. * 0.5625,
        4. * PI * 2.25 - 2. * PI * 1.5 * 0.9375,
    ];
    for (got, want) in result.per_atom.iter().zip(expected) {
        assert!((got - want).abs() < 0.03);
    }
}

#[test]
fn containment_tangency_and_duplicate_ownership() {
    let mol = atoms(2);
    let nested = sasa(&mol, &[[0.; 3], [0.; 3]], &[1., 2.], 0., 960).unwrap();
    assert_eq!(nested.per_atom[0], 0.);
    assert_eq!(nested.per_atom[1], 16. * PI);
    let tangent = sasa(&mol, &[[0.; 3], [2., 0., 0.]], &[1., 1.], 0., 960).unwrap();
    assert_eq!(tangent.total, 8. * PI);
    assert_eq!(
        sasa(&mol, &[[0.; 3], [0.; 3]], &[1., 1.], 0., 960),
        Err(Error::Coincident {
            first: 0,
            second: 1
        })
    );
}

#[test]
fn strict_counts_and_finite_values() {
    let mol = atoms(1);
    assert!(matches!(
        sasa(&mol, &[], &[1.], 0., 960),
        Err(Error::Count {
            field: "coordinates",
            ..
        })
    ));
    assert!(matches!(
        sasa(&mol, &[[0.; 3]], &[], 0., 960),
        Err(Error::Count { field: "radii", .. })
    ));
    for v in [0., -1., f64::NAN, f64::INFINITY] {
        assert_eq!(
            sasa(&mol, &[[0.; 3]], &[v], 0., 960),
            Err(Error::Radius { index: 0 })
        );
    }
    assert_eq!(
        sasa(&mol, &[[0.; 3]], &[1.], -1., 960),
        Err(Error::ProbeRadius)
    );
    assert_eq!(
        sasa(&mol, &[[f64::NAN, 0., 0.]], &[1.], 0., 960),
        Err(Error::Coordinate { index: 0 })
    );
    assert_eq!(
        sasa(&mol, &[[0.; 3]], &[1.], 0., 0),
        Err(Error::SpherePoints)
    );
}

#[test]
fn extreme_coordinates_and_rescaling_remain_finite() {
    let mol = atoms(2);
    let far = sasa(
        &mol,
        &[[1e308, 0., 0.], [-1e308, 0., 0.]],
        &[1., 1.],
        0.,
        960,
    )
    .unwrap();
    assert_eq!(far.total, 8. * PI);
    for s in [1e-100, 1e100] {
        let result = sasa(&mol, &[[0.; 3], [0., s, 0.]], &[s, s], 0., 4096).unwrap();
        assert!((result.total / (s * s) - 6. * PI).abs() < 1e-12);
    }
    for radius in [1e-300, 1e300] {
        assert_eq!(
            sasa(&atoms(1), &[[0.; 3]], &[radius], 0., 960),
            Err(Error::NumericalRange)
        );
    }
}

#[test]
fn empty_and_resource_limits_are_explicit() {
    assert_eq!(sasa(&atoms(0), &[], &[], 0., 960).unwrap().total, 0.);
    assert_eq!(sasa(&atoms(2049), &[], &[], 0., 960), Err(Error::AtomLimit));
    let xyz = (0..40)
        .map(|i| [i as f64 * 0.01, 0., 0.])
        .collect::<Vec<_>>();
    assert_eq!(
        sasa(&atoms(40), &xyz, &[1.; 40], 0., 65536),
        Err(Error::WorkLimit)
    );
}
