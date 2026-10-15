//! Public batch RMSD contracts with analytic two-atom geometry.
use chematic_smiles::{parse, rdkit_all_conformer_best_rms, rdkit_all_conformer_best_rms_to_ref};

fn conformer(length: f64) -> Vec<[f64; 3]> {
    vec![[0.0, 0.0, 0.0], [length, 0.0, 0.0]]
}

fn assert_rms(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert!((actual - expected).abs() < 1e-10, "{actual} != {expected}");
    }
}

#[test]
fn batches_keep_documented_pair_order_and_atom_mapping() {
    // After optimal translation, each atom moves by half the length difference.
    let mol = parse("CO").unwrap();
    let conformers = vec![conformer(1.0), conformer(2.0), conformer(4.0)];
    for symmetrize in [false, true] {
        assert_rms(
            &rdkit_all_conformer_best_rms(&mol, &conformers, 100, symmetrize, None).unwrap(),
            &[0.5, 1.5, 1.0],
        );
        let reference = parse("OC").unwrap();
        let references = vec![conformer(3.0), conformer(6.0)];
        assert_rms(
            &rdkit_all_conformer_best_rms_to_ref(
                &mol,
                &reference,
                &conformers,
                &references,
                100,
                symmetrize,
                None,
            )
            .unwrap(),
            &[1.0, 0.5, 0.5, 2.5, 2.0, 1.0],
        );
    }
}

#[test]
fn empty_and_singleton_batches_have_no_pairs() {
    let mol = parse("CO").unwrap();
    for conformers in [vec![], vec![conformer(1.0)]] {
        assert!(
            rdkit_all_conformer_best_rms(&mol, &conformers, 100, true, None)
                .unwrap()
                .is_empty()
        );
    }
    for (probe, reference) in [
        (vec![], vec![conformer(1.0)]),
        (vec![conformer(1.0)], vec![]),
    ] {
        assert!(
            rdkit_all_conformer_best_rms_to_ref(&mol, &mol, &probe, &reference, 100, true, None)
                .unwrap()
                .is_empty()
        );
    }
}

#[test]
fn batches_reject_bad_coordinates_weights_and_unmatched_graphs() {
    let mol = parse("CO").unwrap();
    let good = vec![conformer(1.0), conformer(2.0)];
    let bad = vec![vec![[0.0; 3]]];
    let error = rdkit_all_conformer_best_rms(&mol, &bad, 100, true, None).unwrap_err();
    assert!(
        error
            .to_string()
            .contains("conformer has 1 positions for 2 atoms")
    );
    for (probe, reference, label) in [(&bad, &good, "probe"), (&good, &bad, "reference")] {
        let error =
            rdkit_all_conformer_best_rms_to_ref(&mol, &mol, probe, reference, 100, true, None)
                .unwrap_err();
        assert!(
            error
                .to_string()
                .contains(&format!("{label} conformer has 1 positions for 2 atoms"))
        );
    }
    let error = rdkit_all_conformer_best_rms(&mol, &good, 100, true, Some(&[1.0])).unwrap_err();
    assert!(error.to_string().contains("Mismatch in number of weights"));
    let different = parse("CN").unwrap();
    let error =
        rdkit_all_conformer_best_rms_to_ref(&mol, &different, &good, &good, 100, true, None)
            .unwrap_err();
    assert!(error.to_string().contains("No sub-structure match"));
}
