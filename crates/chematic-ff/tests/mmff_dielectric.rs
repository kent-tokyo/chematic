use chematic_core::Molecule;
use chematic_ff::rdkit_mmff::{MmffDielectricError, MmffVariant, RdkitMmffField};

fn fixture() -> serde_json::Value {
    serde_json::from_str(include_str!(
        "../../../validation/openbabel_issues/mmff-dielectric/cases.json"
    ))
    .unwrap()
}
fn model(
    row: &serde_json::Value,
    variant: MmffVariant,
    terms: [bool; 7],
) -> (RdkitMmffField, Vec<f64>) {
    let mol: Molecule = chematic_chem::add_hydrogens(
        &chematic_smiles::parse(row["smiles"].as_str().unwrap()).unwrap(),
    );
    let coords: Vec<[f64; 3]> = serde_json::from_value(row["coords"].clone()).unwrap();
    let f = RdkitMmffField::with_terms(&mol, &coords, variant, 100., false, terms).unwrap();
    (f, coords.into_iter().flatten().collect())
}
fn gradient(f: &RdkitMmffField, pos: &[f64]) -> Vec<f64> {
    let mut g = vec![0.; pos.len()];
    f.gradient(pos, &mut g);
    g
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-9 + 1e-11 * expected.abs(),
        "{actual} != {expected}"
    );
}
#[test]
fn live_updates_match_pinned_rdkit_energy_and_gradient() {
    let data = fixture();
    assert_eq!(data["rdkit"], "2026.09.1");
    for row in data["cases"].as_array().unwrap() {
        for variant in [MmffVariant::Mmff94, MmffVariant::Mmff94s] {
            let (mut field, pos) = model(row, variant, [true; 7]);
            let name = if variant == MmffVariant::Mmff94 {
                "MMFF94"
            } else {
                "MMFF94s"
            };
            assert_eq!(field.dielectric_constant(), 1.);
            let baseline = field.energy(&pos);
            let baseline_g = gradient(&field, &pos);
            for epsilon in [1., 4., 78., 0.5, 4., 1.] {
                field.set_dielectric_constant(epsilon).unwrap();
                let reference = row["references"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .find(|r| r["variant"] == name && r["epsilon"].as_f64() == Some(epsilon))
                    .unwrap();
                near(field.energy(&pos), reference["energy"].as_f64().unwrap());
                for (a, e) in gradient(&field, &pos)
                    .into_iter()
                    .zip(reference["gradient"].as_array().unwrap())
                {
                    near(a, e.as_f64().unwrap());
                }
            }
            assert_eq!(field.energy(&pos).to_bits(), baseline.to_bits());
            assert_eq!(gradient(&field, &pos), baseline_g);
        }
    }
}
#[test]
fn screening_scales_only_electrostatic_terms() {
    let data = fixture();
    let row = &data["cases"][4];
    let (mut electric, pos) = model(
        row,
        MmffVariant::Mmff94,
        [false, false, false, false, false, false, true],
    );
    let (mut other, _) = model(
        row,
        MmffVariant::Mmff94,
        [true, true, true, true, true, true, false],
    );
    let e = electric.energy(&pos);
    assert!(e.abs() > 0.1);
    let g = gradient(&electric, &pos);
    let old_other = other.energy(&pos);
    let old_g = gradient(&other, &pos);
    for epsilon in [4., 78., 0.5, 1.] {
        electric.set_dielectric_constant(epsilon).unwrap();
        other.set_dielectric_constant(epsilon).unwrap();
        near(electric.energy(&pos), e / epsilon);
        for (a, b) in gradient(&electric, &pos).iter().zip(&g) {
            near(*a, *b / epsilon);
        }
        assert_eq!(other.energy(&pos).to_bits(), old_other.to_bits());
        assert_eq!(gradient(&other, &pos), old_g);
    }
}
#[test]
fn screened_gradient_matches_finite_differences() {
    let data = fixture();
    let (mut field, pos) = model(&data["cases"][4], MmffVariant::Mmff94s, [true; 7]);
    for epsilon in [4., 78.] {
        field.set_dielectric_constant(epsilon).unwrap();
        for (axis, g) in gradient(&field, &pos).into_iter().enumerate() {
            let mut plus = pos.clone();
            let mut minus = pos.clone();
            plus[axis] += 1e-5;
            minus[axis] -= 1e-5;
            let fd = (field.energy(&plus) - field.energy(&minus)) / 2e-5;
            assert!(
                (fd - g).abs() < 2e-5 + g.abs() * 2e-6,
                "axis {axis}: {fd} != {g}"
            );
        }
    }
}
#[test]
fn invalid_updates_are_atomic_and_clones_independent() {
    let data = fixture();
    let (mut field, pos) = model(&data["cases"][0], MmffVariant::Mmff94, [true; 7]);
    field.set_dielectric_constant(4.).unwrap();
    let energy = field.energy(&pos);
    let grad = gradient(&field, &pos);
    let copy = field.clone();
    for value in [0., -0., -1., f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            field.set_dielectric_constant(value),
            Err(MmffDielectricError::InvalidConstant)
        );
        assert_eq!(field.dielectric_constant(), 4.);
        assert_eq!(field.energy(&pos).to_bits(), energy.to_bits());
        assert_eq!(gradient(&field, &pos), grad);
    }
    assert_eq!(
        field.set_dielectric_constant(f64::from_bits(1)),
        Err(MmffDielectricError::CoefficientOverflow)
    );
    assert_eq!(field.dielectric_constant(), 4.);
    assert_eq!(field.energy(&pos).to_bits(), energy.to_bits());
    assert_eq!(gradient(&field, &pos), grad);
    field.set_dielectric_constant(78.).unwrap();
    assert_eq!(copy.dielectric_constant(), 4.);
    assert_eq!(copy.energy(&pos), energy);
}
#[test]
fn minimizer_uses_updated_screening_with_existing_pairs() {
    let data = fixture();
    for variant in [MmffVariant::Mmff94, MmffVariant::Mmff94s] {
        let (mut reused, mut a) = model(&data["cases"][4], variant, [true; 7]);
        let (mut fresh, mut b) = model(&data["cases"][4], variant, [true; 7]);
        let _ = reused.energy(&a);
        let _ = gradient(&reused, &a);
        reused.set_dielectric_constant(4.).unwrap();
        reused.set_dielectric_constant(78.).unwrap();
        fresh.set_dielectric_constant(78.).unwrap();
        let result_a = reused.optimize(&mut a, 10);
        let result_b = fresh.optimize(&mut b, 10);
        assert_eq!(result_a, result_b);
        assert_eq!(a, b);
        assert_eq!(reused.dielectric_constant(), 78.);
    }
}
