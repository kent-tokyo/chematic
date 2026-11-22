use chematic_fp::{MhfpFingerprint, MhfpLshIndex};

fn fp(hashes: [u64; 4]) -> MhfpFingerprint {
    MhfpFingerprint {
        hashes: hashes.to_vec(),
        num_hashes: 4,
    }
}

#[test]
fn lsh_index_tracks_insertions_and_filters_candidate_similarity() {
    let mut index = MhfpLshIndex::with_bands(2, 2);
    assert!(index.is_empty());
    assert_eq!(index.len(), 0);
    let query = fp([1, 2, 3, 4]);
    assert!(index.query(&query, 0.0).is_empty());
    assert_eq!(index.add(fp([1, 2, 9, 9])), 0);
    assert_eq!(index.add(query.clone()), 1);
    assert_eq!(index.add(fp([5, 6, 7, 8])), 2);
    assert!(!index.is_empty());
    assert_eq!(index.len(), 3);
    assert_eq!(index.query(&query, 0.0), [(1, 1.0), (0, 0.5)]);
    assert_eq!(index.query(&query, 0.6), [(1, 1.0)]);
}
