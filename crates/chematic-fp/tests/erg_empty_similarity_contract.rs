use chematic_fp::{ERG_VEC_LEN, cosine_erg_vec, tanimoto_erg_vec};

#[test]
fn empty_and_orthogonal_erg_vectors_use_defined_similarity_conventions() {
    let empty = [0.0; ERG_VEC_LEN];
    let mut a = empty;
    let mut b = empty;
    a[0] = 1.0;
    b[1] = 1.0;
    assert_eq!(tanimoto_erg_vec(&empty, &empty), 1.0);
    assert_eq!(cosine_erg_vec(&empty, &empty), 0.0);
    assert_eq!(tanimoto_erg_vec(&empty, &a), 0.0);
    assert_eq!(tanimoto_erg_vec(&a, &b), 0.0);
    assert_eq!(tanimoto_erg_vec(&a, &a), 1.0);
    assert_eq!(cosine_erg_vec(&a, &a), 1.0);
}
