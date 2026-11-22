use chematic_smarts::cache::SmartsCache;

#[test]
fn clear_discards_patterns_and_allows_reuse() {
    let mol = chematic_smiles::parse("CCO").unwrap();
    let mut cache = SmartsCache::new(2);
    assert!(cache.is_empty());
    assert!(cache.has_match("O", &mol).unwrap());
    assert!(!cache.has_match("N", &mol).unwrap());
    assert_eq!(cache.len(), 2);
    cache.clear();
    assert_eq!(cache.len(), 0);
    assert!(cache.is_empty());
    assert!(cache.has_match("O", &mol).unwrap());
    assert_eq!(cache.len(), 1);
    assert!(cache.has_match("[", &mol).is_err());
    assert_eq!(cache.len(), 1);
}
