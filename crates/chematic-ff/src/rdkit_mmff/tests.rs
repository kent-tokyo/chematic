use super::*;

/// `CC(=O)Nc1ccc(O)cc1` after `AddHs`, RDKit 2026.03.1
/// `EmbedMolecule(randomSeed=42)` coordinates.
pub(crate) const PARACETAMOL: &str = "CC(=O)Nc1ccc(O)cc1";
pub(crate) const COORDS: [[f64; 3]; 20] = [
    [3.762175028044453, 0.30538563337732244, -0.09745815343294277],
    [2.3329369300899403, 0.11610000936619912, -0.5176254421062628],
    [
        2.1306032462731035,
        -0.09593752049575305,
        -1.7310611678300563,
    ],
    [1.3484147775946786, 0.19222362621777822, 0.4794175681721041],
    [
        -0.04651062589597437,
        0.03729534973500948,
        0.2549743150386899,
    ],
    [-0.6601760831775039, -1.210120317536445, 0.3277557282622417],
    [-2.027875532733807, -1.3550238887272026, 0.10667566446336685],
    [-2.822456517029374, -0.2670983319336878, -0.1922446077538659],
    [-4.187714927031642, -0.3854347185261268, -0.416753109003994],
    [-2.215264059602086, 0.990686079372787, -0.2676816459777285],
    [
        -0.8506489653837782,
        1.1179323355827404,
        -0.04458045101770329,
    ],
    [3.8141231841990857, 1.3434050655543417, 0.3024234778452463],
    [4.477057254518884, 0.10425801972676844, -0.9145679992818092],
    [3.909133114355366, -0.42656371299342444, 0.7351586531314641],
    [1.6266288380623064, 0.37726586194308837, 1.4846206751064692],
    [
        -0.03711779601633895,
        -2.0820010342352497,
        0.5648832746587801,
    ],
    [-2.5130630664729114, -2.314409067782831, 0.16030971994515217],
    [-4.834672252546139, -0.3986600801567403, 0.37067475360645646],
    [-2.8406901294426152, 1.8542155191372227, -0.5040194618516627],
    [-0.3648824178056287, 2.096481172374245, -0.10090179197395394],
];

pub(crate) fn paracetamol() -> Molecule {
    chematic_chem::add_hydrogens(&chematic_smiles::parse(PARACETAMOL).unwrap())
}

#[test]
fn interruptible_mmff_preserves_ordinary_results_and_reports_cancellation() {
    let molecule = paracetamol();
    let field = RdkitMmffField::new(&molecule, &COORDS, 100.0, true).unwrap();
    for limit in [0, 1, 200] {
        let mut expected_coords = COORDS.as_flattened().to_vec();
        let mut actual_coords = expected_coords.clone();
        let expected = field.optimize_detailed(&mut expected_coords, limit);
        let actual = field.optimize_interruptible(&mut actual_coords, limit, || false);
        assert_eq!(
            actual,
            RdkitMmffInterruptibleOptimization::Finished(expected)
        );
        assert_eq!(actual_coords, expected_coords);
    }
    let mut coords = COORDS.as_flattened().to_vec();
    let initial = coords.clone();
    let result = field.optimize_interruptible(&mut coords, 200, || true);
    assert_eq!(coords, initial);
    assert_eq!(
        result,
        RdkitMmffInterruptibleOptimization::Interrupted {
            completed_iterations: 0,
            energy_at_position: None,
        }
    );
    let mut coords = COORDS.as_flattened().to_vec();
    let mut polls = 0;
    let result = field.optimize_interruptible(&mut coords, 200, || {
        polls += 1;
        polls == 20
    });
    let RdkitMmffInterruptibleOptimization::Interrupted {
        completed_iterations,
        energy_at_position: Some(energy),
    } = result
    else {
        panic!("expected interruption with a known energy: {result:?}")
    };
    assert!(completed_iterations > 0);
    assert!((energy - field.energy(&coords)).abs() < 1e-10);
}

fn check(variant: MmffVariant, energy: f64, opt_energy: f64, x0: f64) {
    let mol = paracetamol();
    let ff = RdkitMmffField::with_terms(&mol, &COORDS, variant, 100.0, true, [true; 7]).unwrap();
    let mut pos = COORDS.as_flattened().to_vec();
    let initial = ff.energy(&pos);
    assert!(
        (initial - energy).abs() <= 1e-12,
        "initial energy {initial} differs from RDKit {energy}"
    );
    let (status, e) = ff.optimize(&mut pos, 200);
    assert_eq!(status, 0);
    assert!(
        (e - opt_energy).abs() <= 1e-9,
        "optimized energy {e} differs from RDKit {opt_energy} by {}",
        (e - opt_energy).abs()
    );
    assert!(
        (pos[0] - x0).abs() <= 1e-6,
        "first coordinate {} differs from RDKit {x0}",
        pos[0]
    );

    let mut detailed_pos = COORDS.as_flattened().to_vec();
    let detailed = ff.optimize_detailed(&mut detailed_pos, 200);
    assert_eq!(detailed.status, status);
    assert_eq!(detailed.energy, e);
    assert_eq!(detailed_pos, pos);
    assert!(detailed.iterations > 0);
}

/// RDKit 2026.03.1: `MMFFGetMoleculeForceField(...).CalcEnergy()` and
/// `Minimize(maxIts=200)` (status, energy, first coordinate).
#[test]
fn mmff94_matches_rdkit_with_platform_tolerance() {
    check(
        MmffVariant::Mmff94,
        22.565180175565263,
        -12.775974493568183,
        3.6705868960020265,
    );
}

#[test]
fn mmff94s_matches_rdkit_with_platform_tolerance() {
    check(
        MmffVariant::Mmff94s,
        22.959370426002092,
        -11.453744586799015,
        3.66682072353845,
    );
}

#[test]
fn distance_cache_matches_public_energy_and_gradient() {
    let mol = paracetamol();
    let ff = RdkitMmffField::new(&mol, &COORDS, 100.0, true).unwrap();
    let pos = COORDS.as_flattened();
    let cache = RefCell::new(DistanceCache::new(mol.atom_count()));
    assert_eq!(ff.cached_energy(pos, &cache), ff.energy(pos));

    let mut cached_gradient = vec![0.0; pos.len()];
    let mut public_gradient = vec![0.0; pos.len()];
    ff.cached_scaled_gradient(pos, &mut cached_gradient, &cache);
    ff.gradient(pos, &mut public_gradient);
    crate::rdkit_bfgs::scale_gradient(&mut public_gradient);
    assert_eq!(cached_gradient, public_gradient);

    let generation = cache.borrow().generation;
    ff.cached_energy(pos, &cache);
    assert_eq!(cache.borrow().generation, generation);
}
