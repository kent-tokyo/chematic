use chematic_3d::rdkit_embed::{BoundsMatrix, triangle_smooth_bounds};

#[test]
fn tolerance_repairs_only_small_relative_bound_inversions() {
    // Upper bounds are above the diagonal; lower bounds are below it.
    let raw = vec![0.0, 2.0, 10.0, 2.0001, 0.0, 10.0, 0.1, 0.1, 0.0];
    let mut tolerated = BoundsMatrix::from_raw(3, raw.clone());
    assert!(triangle_smooth_bounds(&mut tolerated, 0.0001));
    assert_eq!(tolerated.lower(0, 1), 2.0001);
    assert_eq!(tolerated.upper(0, 1), 2.0001);

    for tolerance in [0.0, 0.00001] {
        let mut strict = BoundsMatrix::from_raw(3, raw.clone());
        assert!(!triangle_smooth_bounds(&mut strict, tolerance));
    }
}
