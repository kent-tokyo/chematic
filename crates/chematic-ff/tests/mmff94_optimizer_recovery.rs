use chematic_ff::Mmff94EnergyModel;
use chematic_ff::mmff94_energy::mmff94_bond_energy;

#[test]
fn constrained_recovery_on_a_near_equilibrium_bond_preserves_energy_descent() {
    let mol = chematic_smiles::parse("CC").unwrap();
    let model = Mmff94EnergyModel::new(&mol).unwrap();
    let r0 = mmff94_bond_energy(0, 1, 1).unwrap().r0;
    let initial = vec![[0.0, 0.0, 0.0], [r0 + 1e-6, 0.0, 0.0]];
    let gradient = model.bounded_analytic_gradient(&initial);
    let direction: Vec<[f64; 3]> = gradient.iter().map(|g| [-g[0], -g[1], -g[2]]).collect();
    // Exclude every trial on the initial quasi-Newton line search, forcing
    // the public optimizer to use its distinct recovery step. This models
    // a caller geometry boundary, without exposing optimizer internals.
    let mut excluded = Vec::new();
    let mut alpha = 1.0;
    for _ in 0..40 {
        let trial: Vec<[f64; 3]> = initial
            .iter()
            .zip(&direction)
            .map(|(c, p)| {
                [
                    c[0] + alpha * p[0],
                    c[1] + alpha * p[1],
                    c[2] + alpha * p[2],
                ]
            })
            .collect();
        excluded.push(trial);
        alpha *= 0.5;
    }
    let accept =
        |candidate: &[[f64; 3]]| !excluded.iter().any(|point| point.as_slice() == candidate);
    let before = model.energy(&initial);
    let mut coords = initial.clone();
    let result = model
        .minimize_lbfgs_bounded_analytic_with_constraint(&mut coords, 1, accept)
        .unwrap();
    assert_ne!(coords, initial, "an admissible recovery step exists");
    assert!(accept(&coords));
    assert!(
        result.energy < before,
        "recovery must descend: {before} -> {}",
        result.energy
    );
    assert_eq!(result.energy, model.energy(&coords));
}
