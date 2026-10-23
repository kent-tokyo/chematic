use super::*;

#[test]
fn dense_grid_preserves_cell_membership_and_insertion_order() {
    let empty = DenseGrid::<u32>::new(&[]);
    assert_eq!(empty.get((0, 0)), None);
    let entries = [((-2, 3), 7), ((1, -1), 2), ((-2, 3), 9), ((0, 0), 4)];
    let grid = DenseGrid::new(&entries);
    for x in -3..=2 {
        for y in -2..=4 {
            let expected: Vec<_> = entries
                .iter()
                .filter(|(cell, _)| *cell == (x, y))
                .map(|(_, v)| *v)
                .collect();
            assert_eq!(
                grid.get((x, y)),
                (!expected.is_empty()).then_some(expected.as_slice())
            );
        }
    }
}

#[test]
fn nearest_pair_ties_use_atom_ids_and_ignore_unplaced_atoms() {
    let placed = [
        Some(Point::new(0.0, 0.0)),
        Some(Point::new(1.0, 0.0)),
        Some(Point::new(0.0, 1.0)),
        None,
    ];
    for atoms in [
        [AtomIdx(2), AtomIdx(1), AtomIdx(0), AtomIdx(3)],
        [AtomIdx(3), AtomIdx(0), AtomIdx(1), AtomIdx(2)],
    ] {
        let (a, b) = closest_placed_pair(&atoms, &placed);
        assert_eq!((a.min(b), a.max(b)), (AtomIdx(0), AtomIdx(1)));
    }
}

#[test]
fn stacking_counts_only_nonbonded_pairs_independent_of_order_or_motion() {
    let mol = chematic_smiles::parse("CCC.C").unwrap();
    let points = [
        Point::new(0.0, 0.0),
        Point::new(0.05 * BOND_LEN, 0.0),
        Point::new(0.1 * BOND_LEN, 0.0),
        Point::new(5.0 * BOND_LEN, 0.0),
    ];
    let atoms = [AtomIdx(3), AtomIdx(2), AtomIdx(1), AtomIdx(0)];
    for (dx, dy) in [(0.0, 0.0), (2.0, -3.0)] {
        let placed: Vec<_> = points
            .iter()
            .map(|p| Some(Point::new(p.x + dx, p.y + dy)))
            .collect();
        assert_eq!(stacked_pairs(&mol, &atoms, &placed), 1);
        let mut partial = placed;
        partial[2] = None;
        assert_eq!(stacked_pairs(&mol, &atoms, &partial), 0);
    }
}

#[test]
fn reflection_is_an_involution_and_segment_crossings_exclude_touching() {
    let a = Point::new(1.0, 2.0);
    let b = Point::new(3.0, 4.0);
    let q = Point::new(-2.0, 3.0);
    let mirrored = reflect_point(q, a, b);
    let restored = reflect_point(mirrored, a, b);
    assert!((restored.x - q.x).abs() < 1e-12 && (restored.y - q.y).abs() < 1e-12);
    assert!((mirrored.dist(&a) - q.dist(&a)).abs() < 1e-12);
    let same = reflect_point(q, a, a);
    assert_eq!((same.x, same.y), (q.x, q.y));
    for scale in [1e-3, 1.0, 1e3] {
        let p = |x, y| Point::new(x * scale, y * scale);
        let cases = [
            ([p(0.0, 0.0), p(2.0, 2.0), p(0.0, 2.0), p(2.0, 0.0)], true),
            ([p(0.0, 0.0), p(1.0, 1.0), p(1.0, 1.0), p(2.0, 0.0)], false),
            ([p(0.0, 0.0), p(2.0, 0.0), p(1.0, 0.0), p(3.0, 0.0)], false),
            ([p(0.0, 0.0), p(2.0, 0.0), p(0.0, 1.0), p(2.0, 1.0)], false),
        ];
        for ([a, b, c, d], expected) in cases {
            assert_eq!(segments_intersect(a, b, c, d), expected);
            assert_eq!(segments_intersect(b, a, d, c), expected);
            assert_eq!(segments_intersect(c, d, a, b), expected);
        }
    }
}
