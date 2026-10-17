use chematic_core::Molecule;
use chematic_ff::rdkit_mmff::{MmffPartialChargeError as Error, MmffVariant, RdkitMmffField};

const ELECTRIC: [bool; 7] = [false, false, false, false, false, false, true];

fn mol(smiles: &str) -> Molecule {
    chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).unwrap())
}
fn field(m: &Molecule, xyz: &[[f64; 3]], q: &[f64]) -> RdkitMmffField {
    RdkitMmffField::with_terms_and_partial_charges(
        m,
        xyz,
        MmffVariant::Mmff94,
        100.,
        false,
        ELECTRIC,
        q,
    )
    .unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-12 * b.abs().max(1e-20), "{a} != {b}");
}

#[test]
fn arbitrary_charges_use_buffered_coulomb_and_current_dielectric() {
    let m = mol("[Na+].[Cl-]");
    let xyz = [[0., 0., 0.], [4., 0., 0.]];
    for variant in [MmffVariant::Mmff94, MmffVariant::Mmff94s] {
        let mut f = RdkitMmffField::with_terms_and_partial_charges(
            &m,
            &xyz,
            variant,
            100.,
            false,
            ELECTRIC,
            &[0.3, -0.8],
        )
        .unwrap();
        let pos = xyz.into_iter().flatten().collect::<Vec<_>>();
        for epsilon in [1., 4., 78., 0.5, 1.] {
            f.set_dielectric_constant(epsilon).unwrap();
            near(f.energy(&pos), 332.0716 * (0.3 * -0.8 / epsilon) / 4.05);
            let mut g = [0.; 6];
            f.gradient(&pos, &mut g);
            near(g[0], 332.0716 * (0.3 * -0.8 / epsilon) / 4.05_f64.powi(2));
            near(g[3], -g[0]);
        }
        assert_eq!(m.atom(chematic_core::AtomIdx(0)).charge, 1);
        assert_eq!(m.atom(chematic_core::AtomIdx(1)).charge, -1);
    }
}

#[test]
fn initially_zero_mmff_pairs_and_tiny_supplied_charges_are_not_dropped() {
    let m = mol("C.C");
    let xyz = (0..m.atom_count())
        .map(|i| [i as f64 * 4., 0., 0.])
        .collect::<Vec<_>>();
    let pos = xyz.iter().flatten().copied().collect::<Vec<_>>();
    let default =
        RdkitMmffField::with_terms(&m, &xyz, MmffVariant::Mmff94, 100., false, ELECTRIC).unwrap();
    assert_eq!(default.energy(&pos), 0.);
    for charge in [0.5, 1e-12] {
        let mut charges = vec![0.; m.atom_count()];
        charges[0] = charge;
        charges[1] = -charge;
        near(
            field(&m, &xyz, &charges).energy(&pos),
            -332.0716 * charge * charge / 4.05,
        );
    }
}

#[test]
fn topology_exclusions_scaling_cutoff_and_fragment_policy_are_retained() {
    let m = mol("CCCC");
    let xyz = (0..m.atom_count())
        .map(|i| [i as f64, 0., 0.])
        .collect::<Vec<_>>();
    let pos = xyz.iter().flatten().copied().collect::<Vec<_>>();
    for end in [1, 2, 3] {
        let mut charges = vec![0.; m.atom_count()];
        charges[0] = 1.;
        charges[end] = -1.;
        let expected = if end == 3 {
            -332.0716 * 0.75 / 3.05
        } else {
            0.
        };
        near(field(&m, &xyz, &charges).energy(&pos), expected);
    }
    let ions = mol("[Na+].[Cl-]");
    let xyz = [[0.; 3], [4., 0., 0.]];
    for (cutoff, ignore) in [(3., false), (100., true)] {
        let f = RdkitMmffField::with_terms_and_partial_charges(
            &ions,
            &xyz,
            MmffVariant::Mmff94,
            cutoff,
            ignore,
            ELECTRIC,
            &[0.5, -0.5],
        )
        .unwrap();
        assert_eq!(f.energy(&[0., 0., 0., 4., 0., 0.]), 0.);
    }
}

#[test]
fn validation_rejects_bad_counts_values_geometry_and_coefficients() {
    let m = mol("[Na+].[Cl-]");
    let xyz = [[0.; 3], [4., 0., 0.]];
    let make = |xyz: &[[f64; 3]], q: &[f64], cutoff| {
        RdkitMmffField::with_terms_and_partial_charges(
            &m,
            xyz,
            MmffVariant::Mmff94,
            cutoff,
            false,
            ELECTRIC,
            q,
        )
    };
    assert!(matches!(make(&xyz, &[1.], 100.), Err(Error::Count { .. })));
    for v in [f64::NAN, f64::INFINITY, -f64::INFINITY] {
        assert!(matches!(
            make(&xyz, &[v, 1.], 100.),
            Err(Error::Nonfinite { index: 0 })
        ));
    }
    assert!(matches!(
        make(&[], &[1., 1.], 100.),
        Err(Error::Coordinates)
    ));
    assert!(matches!(
        make(&[[f64::NAN; 3]; 2], &[1., 1.], 100.),
        Err(Error::Coordinates)
    ));
    assert!(matches!(make(&xyz, &[1., 1.], -1.), Err(Error::Threshold)));
    for invalid in [[1e300, 1e300], [1e-300, 1e-300], [1e154, 1e154]] {
        assert!(matches!(
            make(&xyz, &invalid, 100.),
            Err(Error::CoefficientRange)
        ));
    }
}

#[test]
fn supplied_charges_leave_non_electrostatic_terms_identical() {
    let m = mol("CCO");
    let data: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/openbabel_issues/mmff-dielectric/cases.json"
    ))
    .unwrap();
    let xyz: Vec<[f64; 3]> = serde_json::from_value(data["cases"][2]["coords"].clone()).unwrap();
    let pos = xyz.iter().flatten().copied().collect::<Vec<_>>();
    let terms = [true, true, true, true, true, true, false];
    let default =
        RdkitMmffField::with_terms(&m, &xyz, MmffVariant::Mmff94, 100., false, terms).unwrap();
    let custom = RdkitMmffField::with_terms_and_partial_charges(
        &m,
        &xyz,
        MmffVariant::Mmff94,
        100.,
        false,
        terms,
        &vec![0.25; m.atom_count()],
    )
    .unwrap();
    assert_eq!(
        custom.energy(&pos).to_bits(),
        default.energy(&pos).to_bits()
    );
    let mut a = vec![0.; pos.len()];
    let mut b = a.clone();
    custom.gradient(&pos, &mut a);
    default.gradient(&pos, &mut b);
    assert_eq!(a, b);
}
