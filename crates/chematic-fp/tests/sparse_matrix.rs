use chematic_fp::{
    SparseMatrixError as Error, SparseTanimotoMetric as Metric,
    sparse_tanimoto_matrix_checked as matrix,
};

#[test]
fn strict_lower_triangle_and_metric_complement() {
    let rows = vec![vec![], vec![1], vec![1, 2], vec![3]];
    let similarity = matrix(&rows, Metric::Similarity).unwrap();
    assert_eq!(similarity, vec![0., 0., 0.5, 0., 0., 0.]);
    let distance = matrix(&rows, Metric::Distance).unwrap();
    assert_eq!(distance, vec![1., 1., 0.5, 1., 1., 1.]);
    for (s, d) in similarity.iter().zip(distance) {
        assert_eq!(s + d, 1.);
    }
}
#[test]
fn empty_input_and_empty_sets_are_explicit() {
    assert!(matrix(&[], Metric::Distance).unwrap().is_empty());
    assert!(matrix(&[vec![7]], Metric::Similarity).unwrap().is_empty());
    assert_eq!(
        matrix(&[vec![], vec![], vec![7]], Metric::Similarity).unwrap(),
        vec![1., 0., 0.]
    );
    assert_eq!(
        matrix(&[vec![], vec![], vec![7]], Metric::Distance).unwrap(),
        vec![0., 1., 1.]
    );
}
#[test]
fn duplicate_order_and_sources_are_preserved() {
    let rows = vec![vec![9, 1, 9, 5], vec![5, 1], vec![9]];
    let before = rows.clone();
    let s = matrix(&rows, Metric::Similarity).unwrap();
    assert_eq!(s, vec![2. / 3., 1. / 3., 0.]);
    assert_eq!(rows, before);
    let permuted = vec![rows[2].clone(), rows[0].clone(), rows[1].clone()];
    assert_eq!(
        matrix(&permuted, Metric::Similarity).unwrap(),
        vec![s[1], s[2], s[0]]
    );
}
#[test]
fn high_u32_ids_do_not_require_dense_storage() {
    let rows = vec![vec![0, u32::MAX], vec![u32::MAX], vec![u32::MAX - 1]];
    assert_eq!(
        matrix(&rows, Metric::Similarity).unwrap(),
        vec![0.5, 0., 0.]
    );
}
#[test]
fn raw_limits_precede_deduplication_and_exact_caps_work() {
    assert_eq!(
        matrix(&vec![vec![]; 1025], Metric::Distance),
        Err(Error::RowLimit)
    );
    assert_eq!(
        matrix(&[vec![7; 65537]], Metric::Distance),
        Err(Error::BitsPerFingerprint { row: 0 })
    );
    assert!(
        matrix(&[vec![7; 65536]], Metric::Similarity)
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        matrix(&vec![vec![7; 65536]; 16], Metric::Distance),
        Err(Error::TotalIds)
    );
    assert_eq!(
        matrix(&vec![vec![7; 62500]; 16], Metric::Similarity).unwrap(),
        vec![1.; 120]
    );
    assert_eq!(
        matrix(&vec![vec![]; 1024], Metric::Distance).unwrap().len(),
        523776
    );
}
#[test]
fn pair_work_exact_boundary_and_next_bit() {
    let mut rows = vec![(0..200).collect::<Vec<u32>>(); 500];
    rows.push(vec![]);
    let result = matrix(&rows, Metric::Similarity).unwrap();
    assert_eq!(result.len(), 125250);
    assert_eq!(result[0], 1.);
    assert_eq!(*result.last().unwrap(), 0.);
    rows[0].push(300);
    assert_eq!(matrix(&rows, Metric::Similarity), Err(Error::WorkLimit));
}

#[test]
fn agrees_with_existing_native_dense_tanimoto_including_empty_sets() {
    use chematic_fp::BitVec2048;
    let rows = vec![vec![], vec![], vec![1, 3, 2047], vec![3, 5]];
    let dense: Vec<_> = rows
        .iter()
        .map(|ids| {
            let mut fp = BitVec2048::new();
            for &id in ids {
                fp.set(id as usize);
            }
            fp
        })
        .collect();
    let expected: Vec<_> = (1..dense.len())
        .flat_map(|i| {
            (0..i).map({
                let dense = &dense;
                move |j| dense[i].tanimoto(&dense[j])
            })
        })
        .collect();
    assert_eq!(matrix(&rows, Metric::Similarity).unwrap(), expected);
}
