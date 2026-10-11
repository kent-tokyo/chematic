/// Return `true` if the permutation mapping `original` order to `canonical`
/// order has odd parity (i.e. requires an odd number of transpositions).
///
/// Both slices must contain the same multiset of values (no duplicates — a
/// repeated value would collide in the by-value index below). Generic over
/// `T` (any length) so it backs two different callers:
/// `corrected_chirality`'s non-4-element fallback (an allene *end* carbon's
/// 3-element `stereo_neighbor_order` -- see that function's comment; the
/// common 4-element tetrahedral case goes through
/// `chematic_core::remap_tetrahedral_parity` instead, the generalized
/// stereo-geometry module, `docs/rfcs/generalized_stereo_geometry_rfc.md`) and
/// the test-only rank-based tetrahedral fingerprint (`u64` canonical ranks,
/// arbitrary length).
fn permutation_is_odd<T: Eq + std::hash::Hash + Copy>(original: &[T], canonical: &[T]) -> bool {
    let n = original.len();
    let mut pos: HashMap<T, usize> = HashMap::with_capacity(n);
    for (i, &v) in original.iter().enumerate() {
        pos.insert(v, i);
    }
    // perm[i] = position in `original` of the element at `canonical[i]`
    let perm: Vec<usize> = canonical
        .iter()
        .map(|v| *pos.get(v).unwrap_or(&0))
        .collect();

    // Count cycles in the permutation; parity = (n - #cycles) % 2
    let mut visited = vec![false; n];
    let mut num_cycles = 0usize;
    for start in 0..n {
        if !visited[start] {
            num_cycles += 1;
            let mut j = start;
            while !visited[j] {
                visited[j] = true;
                j = perm[j];
            }
        }
    }
    (n - num_cycles) % 2 == 1
}
