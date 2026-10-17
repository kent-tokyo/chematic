use chematic_chem::descriptor_distances::{
    DescriptorDistanceError as Error, euclidean_distance_matrix_checked as distances,
};

#[test]
fn strict_lower_triangle_order_and_duplicates() {
    assert_eq!(
        distances(&[[0., 0.], [3., 4.], [3., 0.]]).unwrap(),
        vec![5., 3., 4.]
    );
    assert_eq!(distances(&[[1.], [1.], [2.]]).unwrap(), vec![0., 1., 1.]);
}

#[test]
fn empty_singleton_and_zero_feature_inputs() {
    assert_eq!(
        distances(&Vec::<Vec<f64>>::new()).unwrap(),
        Vec::<f64>::new()
    );
    assert!(distances(&[[123.]]).unwrap().is_empty());
    assert_eq!(
        distances(&[Vec::<f64>::new(), vec![], vec![]]).unwrap(),
        vec![0.; 3]
    );
    assert!(matches!(
        distances(&[[f64::NAN]]),
        Err(Error::Nonfinite { row: 0, column: 0 })
    ));
}

#[test]
fn invalid_shapes_and_nonfinite_values_are_refused() {
    assert_eq!(
        distances(&[vec![0., 1.], vec![0.]]),
        Err(Error::Ragged {
            row: 1,
            expected: 2,
            actual: 1
        })
    );
    for value in [f64::NAN, f64::INFINITY, -f64::INFINITY] {
        assert_eq!(
            distances(&[[0., value], [1., 2.]]),
            Err(Error::Nonfinite { row: 0, column: 1 })
        );
    }
}

#[test]
fn scaled_norm_preserves_extreme_representable_distances() {
    for value in [1e-200, 1e200, f64::from_bits(1), f64::MAX] {
        assert_eq!(distances(&[[0.], [value]]).unwrap(), vec![value]);
    }
    for scale in [1e-200, 1., 1e200] {
        let value = distances(&[[0., 0.], [3. * scale, 4. * scale]]).unwrap()[0];
        assert!((value / scale - 5.).abs() < 2e-15);
    }
    let value = distances(&[[0., 0.], [1e308, 1e308]]).unwrap()[0];
    assert!((value / 1e308 - 2_f64.sqrt()).abs() < 1e-15);
}

#[test]
fn unrepresentable_norm_and_subtraction_refuse_partial_output() {
    assert!(matches!(
        distances(&[[0., 0.], [f64::MAX, f64::MAX]]),
        Err(Error::NumericalRange {
            first: 1,
            second: 0
        })
    ));
    assert!(matches!(
        distances(&[[0.], [f64::MAX], [-f64::MAX]]),
        Err(Error::NumericalRange {
            first: 2,
            second: 1
        })
    ));
}

#[test]
fn input_and_work_limits_apply_before_computation() {
    assert_eq!(distances(&vec![vec![]; 1025]), Err(Error::RowLimit));
    assert_eq!(distances(&[vec![0.; 257]]), Err(Error::FeatureLimit));
    assert_eq!(distances(&vec![vec![0.; 256]; 1024]), Err(Error::WorkLimit));
    assert_eq!(distances(&vec![vec![]; 1024]).unwrap().len(), 523776);
    assert_eq!(
        distances(&[vec![0.; 256], vec![1.; 256]]).unwrap(),
        vec![16.]
    );
}

#[test]
fn permutations_and_exact_translations_preserve_correspondence() {
    let rows = [[0., 0.], [3., 4.], [3., 0.]];
    let before = distances(&rows).unwrap();
    let translated = rows.map(|r| r.map(|v| v + 1024.));
    assert_eq!(distances(&translated).unwrap(), before);
    assert_eq!(
        distances(&[rows[2], rows[0], rows[1]]).unwrap(),
        vec![3., 4., 5.]
    );
    assert_eq!(rows, [[0., 0.], [3., 4.], [3., 0.]]);
}
