//! RDKit 2026.09.2 angle and stretch-bend parameters on exposed cage structures.
use chematic_ff::rdkit_mmff::{MmffVariant, RdkitMmffField};

#[test]
fn equal_minimum_cage_cycles_parameterize_all_four_ring_angles() {
    // Synthetic coordinates are deliberately fixed to isolate parameters from
    // embedding. Values and gradients were generated with RDKit 2026.09.2,
    // with one MMFF term class enabled at a time.
    for (smiles, atom, expected) in [
        (
            "COC(=O)[C@]12[C@@H]3[C@@H]4[C@H]1[C@@H]1[C@H]2[C@H]3[C@]41COC(=O)C(Cl)(Cl)Cl",
            4,
            [
                (
                    9972.48369006898,
                    [-34.802686216936735, 253.19253669002177, 11.659983130787126],
                ),
                (
                    -16162.223544065324,
                    [199.7744066260464, 272.38308016630776, -180.9784231905409],
                ),
            ],
        ),
        (
            "COC(=O)[C@]12[C@@H]3[C@@H]4[C@H]1[C@@H]1[C@H]2[C@H]3[C@]41C(=O)OC",
            8,
            [
                (
                    7585.170145722135,
                    [26.59932924758498, -177.6033515272834, -25.15717554228438],
                ),
                (
                    -11785.89001457829,
                    [-474.7768207694628, 565.9695098117993, 361.8951603379655],
                ),
            ],
        ),
    ] {
        let mol = chematic_chem::add_hydrogens(&chematic_smiles::parse(smiles).unwrap());
        let coords: Vec<[f64; 3]> = (0..mol.atom_count())
            .map(|i| {
                [
                    ((i * 37) % 101) as f64 / 17.0,
                    ((i * 53) % 103) as f64 / 19.0,
                    ((i * 67) % 107) as f64 / 23.0,
                ]
            })
            .collect();
        let flat: Vec<f64> = coords.iter().flatten().copied().collect();
        for (term, (expected_energy, expected_gradient)) in [1usize, 2].into_iter().zip(expected) {
            let mut terms = [false; 7];
            terms[term] = true;
            let field =
                RdkitMmffField::with_terms(&mol, &coords, MmffVariant::Mmff94, 100.0, true, terms)
                    .unwrap();
            let actual = field.energy(&flat);
            assert!(
                (actual - expected_energy).abs() < 1e-6,
                "{smiles} term {term}: {actual} vs {expected_energy}"
            );
            let mut gradient = vec![0.0; flat.len()];
            field.gradient(&flat, &mut gradient);
            for axis in 0..3 {
                let actual = gradient[3 * atom + axis];
                assert!(
                    (actual - expected_gradient[axis]).abs() < 1e-6,
                    "{smiles} term {term} axis {axis}: {actual} vs {}",
                    expected_gradient[axis]
                );
            }
        }
    }
}
