use chematic_fp::{
    BitVecN, FoldDensityError, MAX_FOLD_DENSITY_BITS, fold_packed_to_target_density,
    fold_to_target_density,
};

fn fp(width: usize, on: &[usize]) -> BitVecN {
    let mut fp = BitVecN::new(width);
    for &bit in on {
        fp.set(bit);
    }
    fp
}

fn on(fp: &BitVecN) -> Vec<usize> {
    (0..fp.bit_width()).filter(|&bit| fp.get(bit)).collect()
}

#[test]
fn sparse_fingerprints_fold_and_dense_ones_do_not() {
    let sparse = fp(256, &[0]);
    let result = diagnostic(fold_to_target_density(&sparse, 0.3, 64)).unwrap();
    assert_eq!(result.fingerprint.bit_width(), 64);
    assert_eq!(result.folds, 2);
    assert!(!result.target_reached);
    assert_eq!(sparse.bit_width(), 256);
    let dense = fp(256, &(0..256).collect::<Vec<_>>());
    let result = diagnostic(fold_to_target_density(&dense, 0.3, 64)).unwrap();
    assert_eq!(result.fingerprint, dense);
    assert_eq!(result.folds, 0);
    assert!(result.target_reached);
}

#[test]
fn equality_at_density_stops_and_equality_at_minimum_allows_one_fold() {
    let input = fp(128, &(0..16).collect::<Vec<_>>());
    let result = diagnostic(fold_to_target_density(&input, 0.25, 64)).unwrap();
    assert_eq!(result.fingerprint.bit_width(), 64);
    assert_eq!(result.density(), 0.25);
    assert!(result.target_reached);
    let result = diagnostic(fold_to_target_density(&input, 0.125, 64)).unwrap();
    assert_eq!(result.fingerprint, input);
    assert_eq!(result.folds, 0);
}

#[test]
fn collisions_are_or_not_xor() {
    let result = diagnostic(fold_to_target_density(&fp(128, &[0, 64]), 1., 64)).unwrap();
    assert_eq!(on(&result.fingerprint), [0]);
    assert_eq!(result.density(), 1. / 64.);
}

#[test]
fn odd_width_tail_bits_use_modulo_and_can_reduce_density() {
    let result = diagnostic(fold_to_target_density(&fp(5, &[0, 2, 4]), 0.8, 2)).unwrap();
    assert_eq!(result.fingerprint.bit_width(), 2);
    assert_eq!(on(&result.fingerprint), [0]);
    assert_eq!(result.density(), 0.5);
    assert!(!result.target_reached);
    let result = diagnostic(fold_to_target_density(&fp(65, &[0, 32, 64]), 0.1, 32)).unwrap();
    assert_eq!(on(&result.fingerprint), [0]);
    assert_eq!(result.fingerprint.bit_width(), 32);
}

#[test]
fn zero_bits_reach_minimum_but_zero_target_needs_no_fold() {
    let result = diagnostic(fold_to_target_density(&fp(128, &[]), 0.3, 1)).unwrap();
    assert_eq!(result.fingerprint.bit_width(), 1);
    assert_eq!(result.folds, 7);
    assert!(!result.target_reached);
    let result = diagnostic(fold_to_target_density(&fp(128, &[]), 0., 1)).unwrap();
    assert_eq!(result.folds, 0);
    assert!(result.target_reached);
}

#[test]
fn input_shorter_than_minimum_and_single_bit_are_valid() {
    assert_eq!(
        diagnostic(fold_to_target_density(&fp(8, &[]), 1., 64))
            .unwrap()
            .fingerprint
            .bit_width(),
        8
    );
    assert!(
        diagnostic(fold_to_target_density(&fp(1, &[0]), 1., 1))
            .unwrap()
            .target_reached
    );
    assert_eq!(
        diagnostic(fold_to_target_density(&fp(1, &[]), 1., 1))
            .unwrap()
            .folds,
        0
    );
}

#[test]
fn invalid_options_are_typed_errors() {
    for density in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.1, 1.1] {
        assert_eq!(
            diagnostic(fold_to_target_density(&fp(8, &[]), density, 1)),
            Err(FoldDensityError::InvalidDensity)
        );
    }
    assert_eq!(
        diagnostic(fold_to_target_density(&fp(8, &[]), 0.3, 0)),
        Err(FoldDensityError::InvalidMinimumLength)
    );
}

#[test]
fn packed_validation_precedes_decode_and_ignores_no_data() {
    assert_eq!(
        diagnostic(fold_packed_to_target_density(&[], 0, 0.3, 64)),
        Err(FoldDensityError::EmptyFingerprint)
    );
    assert_eq!(
        diagnostic(fold_packed_to_target_density(&[], usize::MAX, 0.3, 64)),
        Err(FoldDensityError::ResourceLimit {
            bits: usize::MAX,
            limit: MAX_FOLD_DENSITY_BITS
        })
    );
    assert_eq!(
        diagnostic(fold_packed_to_target_density(&[0, 0], 8, 0.3, 64)),
        Err(FoldDensityError::PackedLength {
            actual: 2,
            expected: 1
        })
    );
    assert_eq!(
        diagnostic(fold_packed_to_target_density(&[128], 7, 0.3, 64)),
        Err(FoldDensityError::NonzeroPadding)
    );
    assert_eq!(
        diagnostic(fold_packed_to_target_density(&[0], 9, 0.3, 64)),
        Err(FoldDensityError::PackedLength {
            actual: 1,
            expected: 2
        })
    );
}

#[test]
fn maximum_width_accepts_and_next_bit_rejects() {
    let input = fp(MAX_FOLD_DENSITY_BITS, &[0]);
    let result = diagnostic(fold_to_target_density(&input, 1., 1)).unwrap();
    assert_eq!(result.fingerprint.bit_width(), 1);
    assert_eq!(result.folds, 20);
    assert!(result.target_reached);
    assert!(matches!(
        diagnostic(fold_to_target_density(
            &BitVecN::new(MAX_FOLD_DENSITY_BITS + 1),
            1.,
            1
        )),
        Err(FoldDensityError::ResourceLimit { .. })
    ));
}

#[test]
fn packed_and_native_inputs_agree_for_partial_bytes_and_words() {
    for width in [1, 7, 8, 9, 63, 64, 65, 127, 128, 129, 1000, 2048] {
        let input = fp(width, &(0..width).step_by(3).collect::<Vec<_>>());
        for density in [0., 0.3, 0.75, 1.] {
            let direct = diagnostic(fold_to_target_density(&input, density, 1)).unwrap();
            let packed = diagnostic(fold_packed_to_target_density(
                &input.to_le_bytes(),
                width,
                density,
                1,
            ))
            .unwrap();
            assert_eq!(direct, packed);
        }
    }
}

#[path = "support/upstream_diagnostics.rs"]
mod upstream_diagnostics;
use upstream_diagnostics::diagnostic;
