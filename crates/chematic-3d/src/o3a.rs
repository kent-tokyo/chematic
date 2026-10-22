//! Open3DAlign (O3A)-style atom correspondence search and alignment for 3D
//! molecules.
//!
//! Unlike simple RMSD superposition ([`crate::align::align_coords`]), which
//! requires atom pairs to already be known, O3A finds the atom correspondence
//! itself ([`correspondence_search`]) and then feeds it into that
//! already-tested Kabsch alignment ([`o3a_align`]).
//!
//! # Algorithm
//! 1. Assign MMFF94 atom types to both molecules; two atoms are considered
//!    "compatible" only when their types match exactly.
//! 2. Pre-align by principal axes (via the inertia-like covariance tensor's
//!    eigenvectors). Principal axes have a sign ambiguity — of the 8 possible
//!    per-axis sign combinations only 4 are proper rotations (det = +1); all
//!    4 are tried as independent seeds, since picking the wrong one can
//!    converge to a mirrored, degenerate correspondence.
//! 3. For each seed: globally pair compatible atoms within a distance cutoff
//!    for up to 64 atoms per molecule, maximizing matched atoms and then
//!    Gaussian overlap. Larger molecules retain the bounded greedy path. Alternate
//!    "refit rotation for the current pairing" (reusing `align_coords`) and
//!    "re-pair against the newly rotated coordinates" until the overlap
//!    score stops improving.
//! 4. Return the correspondence from whichever seed reached the best score.

#[allow(unused_imports)]
use crate::fmath::DetMath;
use chematic_core::Molecule;
use chematic_ff::MMFF94Type;

use crate::align::{align_coords, apply_alignment};
use crate::shape_descriptors::jacobi3;

/// Error returned by O3A correspondence search.
#[derive(Debug)]
pub enum O3AError {
    /// `coords.len()` did not match `mol.atom_count()`.
    CoordinateMismatch,
    /// MMFF94 atom-type assignment failed for one of the molecules.
    TypeAssignment(chematic_ff::AssignError),
}

impl std::fmt::Display for O3AError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            O3AError::CoordinateMismatch => write!(f, "coords length does not match atom count"),
            O3AError::TypeAssignment(e) => write!(f, "MMFF94 type assignment failed: {e}"),
        }
    }
}

impl std::error::Error for O3AError {}

const DIST_CUTOFF: f64 = 4.0; // Å — candidate pairs beyond this are rejected
const MAX_REFINE_ITERS: usize = 20;
const SCORE_SIGMA: f64 = 1.0; // Å, Gaussian overlap score width
const CONVERGENCE_EPS: f64 = 1e-4;
const GLOBAL_ASSIGNMENT_MAX_ATOMS: usize = 64;

fn coordinate_order(a: &[f64; 3], b: &[f64; 3]) -> std::cmp::Ordering {
    a.iter()
        .zip(b)
        .map(|(x, y)| x.total_cmp(y))
        .find(|order| !order.is_eq())
        .unwrap_or(std::cmp::Ordering::Equal)
}

fn canonical_coordinate_order(coords: &[[f64; 3]]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..coords.len()).collect();
    order.sort_by(|&a, &b| coordinate_order(&coords[a], &coords[b]));
    order
}

// ---------------------------------------------------------------------------
// Small linear-algebra helpers (local to this module — not a general
// utility, so kept private rather than widening align.rs's visibility)
// ---------------------------------------------------------------------------

fn centroid(coords: &[[f64; 3]]) -> [f64; 3] {
    let n = coords.len().max(1) as f64;
    let mut c = [0.0; 3];
    for index in canonical_coordinate_order(coords) {
        let p = &coords[index];
        for k in 0..3 {
            c[k] += p[k];
        }
    }
    for v in &mut c {
        *v /= n;
    }
    c
}

/// Unweighted (geometric) covariance tensor of `coords` about `center`.
///
/// Used only to seed a rough pre-alignment via principal axes — not a
/// physically accurate mass-weighted inertia tensor.
fn covariance_tensor(coords: &[[f64; 3]], center: [f64; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0; 3]; 3];
    for index in canonical_coordinate_order(coords) {
        let p = &coords[index];
        let x = p[0] - center[0];
        let y = p[1] - center[1];
        let z = p[2] - center[2];
        t[0][0] += x * x;
        t[1][1] += y * y;
        t[2][2] += z * z;
        t[0][1] += x * y;
        t[0][2] += x * z;
        t[1][2] += y * z;
    }
    t[1][0] = t[0][1];
    t[2][0] = t[0][2];
    t[2][1] = t[1][2];
    t
}

fn det3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Flip the third column's sign if needed so `axes` is a proper rotation
/// (det = +1). `jacobi3`'s eigenvector accumulation is always a proper
/// rotation, but its eigenvalue-ascending sort permutes columns and can
/// introduce an odd (det-flipping) permutation.
fn normalize_proper(mut axes: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    if det3(axes) < 0.0 {
        for row in axes.iter_mut() {
            row[2] *= -1.0;
        }
    }
    axes
}

fn matmul3(a: [[f64; 3]; 3], b_transposed: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for (row, r_row) in r.iter_mut().enumerate() {
        for (col, r_val) in r_row.iter_mut().enumerate() {
            *r_val = (0..3).map(|k| a[row][k] * b_transposed[col][k]).sum();
        }
    }
    r
}

fn apply_rotation_translation(
    coords: &[[f64; 3]],
    center: [f64; 3],
    rot: [[f64; 3]; 3],
    target_center: [f64; 3],
) -> Vec<[f64; 3]> {
    coords
        .iter()
        .map(|p| {
            let cp = [p[0] - center[0], p[1] - center[1], p[2] - center[2]];
            let mut out = [0.0; 3];
            for (row, out_val) in out.iter_mut().enumerate() {
                *out_val = (0..3).map(|k| rot[row][k] * cp[k]).sum::<f64>() + target_center[row];
            }
            out
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Correspondence search
// ---------------------------------------------------------------------------

/// Gaussian overlap score (`Σ exp(-d²/2σ²)`) between two equal-length,
/// already-paired coordinate lists — higher means a tighter, more extensive
/// overlap.
fn overlap_score_paired(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(p, q)| {
            let d2: f64 = (0..3).map(|k| (p[k] - q[k]).powi(2)).sum();
            (-d2 / (2.0 * SCORE_SIGMA * SCORE_SIGMA)).dexp()
        })
        .sum()
}

fn overlap_score(pairs: &[(usize, usize)], coords1: &[[f64; 3]], coords2: &[[f64; 3]]) -> f64 {
    let a: Vec<[f64; 3]> = pairs.iter().map(|&(i, _)| coords1[i]).collect();
    let b: Vec<[f64; 3]> = pairs.iter().map(|&(_, j)| coords2[j]).collect();
    overlap_score_paired(&a, &b)
}

/// Greedily pair each mol1 atom (farthest-from-centroid first) with its
/// nearest not-yet-used, type-compatible mol2 atom within `DIST_CUTOFF`.
///
/// Both tie-breaks below (processing order on an exact centroid-distance
/// tie; candidate choice on an exact nearest-neighbor-distance tie) break
/// on the tied atoms' own coordinates (content), not on their array index
/// -- an index-based tie-break would still vary with input order, since
/// swapping two tied atoms' array positions is exactly the operation that
/// changes which one occupies the "winning" index (see the analogous fix
/// and its doc comment in `usr.rs::extreme_idx`, which found this the hard
/// way empirically).
///
/// This does NOT make the large-molecule fallback fully order-independent,
/// only its tie-breaking: the algorithm is inherently greedy
/// (first-processed mol1 atom claims its nearest mol2 atom, removing it
/// from the pool for every atom processed after), so even with canonical
/// tie-breaks, two mol1 atoms that are NOT tied but are both plausible
/// matches for the same mol2 atom can still get different pairings
/// depending on which one the (now-deterministic, but still a specific
/// choice) processing order visits first. Small and medium molecules use the
/// global assignment below to avoid this order dependence.
fn greedy_correspondence(
    coords1: &[[f64; 3]],
    types1: &[MMFF94Type],
    coords2: &[[f64; 3]],
    types2: &[MMFF94Type],
) -> Vec<(usize, usize)> {
    let c1 = centroid(coords1);
    let mut order: Vec<usize> = (0..coords1.len()).collect();
    order.sort_by(|&a, &b| {
        let da: f64 = (0..3).map(|k| (coords1[a][k] - c1[k]).powi(2)).sum();
        let db: f64 = (0..3).map(|k| (coords1[b][k] - c1[k]).powi(2)).sum();
        db.partial_cmp(&da)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                coords1[a]
                    .partial_cmp(&coords1[b])
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    let mut used = vec![false; coords2.len()];
    let mut pairs = Vec::new();
    for i in order {
        let mut best_j: Option<usize> = None;
        let mut best_d2 = f64::INFINITY;
        for j in 0..coords2.len() {
            if used[j] || types1[i] != types2[j] {
                continue;
            }
            let d2: f64 = (0..3)
                .map(|k| (coords1[i][k] - coords2[j][k]).powi(2))
                .sum();
            let better = d2 < best_d2
                || (d2 == best_d2 && best_j.is_some_and(|bj| coords2[j] < coords2[bj]));
            if better {
                best_d2 = d2;
                best_j = Some(j);
            }
        }
        if let Some(j) = best_j
            && best_d2.sqrt() <= DIST_CUTOFF
        {
            used[j] = true;
            pairs.push((i, j));
        }
    }
    pairs
}

/// Rectangular minimum-cost assignment. There are at least as many columns as
/// rows; one zero-cost dummy column per row permits unmatched atoms. Rows and
/// real columns arrive in coordinate order so exact assignment ties do not
/// depend on the molecule's atom insertion order.
fn assign_rows(costs: &[Vec<f64>]) -> Vec<usize> {
    let rows = costs.len();
    let columns = costs[0].len();
    let mut row_potential = vec![0.0; rows + 1];
    let mut col_potential = vec![0.0; columns + 1];
    let mut matched_row = vec![0usize; columns + 1];
    let mut predecessor = vec![0usize; columns + 1];

    for row in 1..=rows {
        matched_row[0] = row;
        let mut current_col = 0;
        let mut minimum = vec![f64::INFINITY; columns + 1];
        let mut visited = vec![false; columns + 1];
        loop {
            visited[current_col] = true;
            let current_row = matched_row[current_col];
            let mut delta = f64::INFINITY;
            let mut next_col = 0;
            for col in 1..=columns {
                if visited[col] {
                    continue;
                }
                let reduced = costs[current_row - 1][col - 1]
                    - row_potential[current_row]
                    - col_potential[col];
                if reduced < minimum[col] {
                    minimum[col] = reduced;
                    predecessor[col] = current_col;
                }
                if minimum[col] < delta {
                    delta = minimum[col];
                    next_col = col;
                }
            }
            for col in 0..=columns {
                if visited[col] {
                    row_potential[matched_row[col]] += delta;
                    col_potential[col] -= delta;
                } else {
                    minimum[col] -= delta;
                }
            }
            current_col = next_col;
            if matched_row[current_col] == 0 {
                break;
            }
        }
        loop {
            let previous = predecessor[current_col];
            matched_row[current_col] = matched_row[previous];
            current_col = previous;
            if current_col == 0 {
                break;
            }
        }
    }

    let mut assignment = vec![usize::MAX; rows];
    for col in 1..=columns {
        if matched_row[col] != 0 {
            assignment[matched_row[col] - 1] = col - 1;
        }
    }
    assignment
}

fn global_correspondence(
    coords1: &[[f64; 3]],
    types1: &[MMFF94Type],
    coords2: &[[f64; 3]],
    types2: &[MMFF94Type],
) -> Vec<(usize, usize)> {
    if coords1.is_empty() || coords2.is_empty() {
        return Vec::new();
    }
    // Bound cubic assignment work for large molecules. Their former greedy
    // behavior remains available until a sparse matching path is provided.
    if coords1.len().max(coords2.len()) > GLOBAL_ASSIGNMENT_MAX_ATOMS {
        return greedy_correspondence(coords1, types1, coords2, types2);
    }
    let left = canonical_coordinate_order(coords1);
    let right = canonical_coordinate_order(coords2);
    let real_columns = right.len();
    let match_bonus = (left.len().min(real_columns) + 1) as f64;
    let mut costs = vec![vec![0.0; real_columns + left.len()]; left.len()];
    for (row, &i) in left.iter().enumerate() {
        for (col, &j) in right.iter().enumerate() {
            let distance_squared: f64 = (0..3)
                .map(|axis| (coords1[i][axis] - coords2[j][axis]).powi(2))
                .sum();
            costs[row][col] =
                if types1[i] == types2[j] && distance_squared <= DIST_CUTOFF * DIST_CUTOFF {
                    // One additional match outweighs every possible change
                    // in the remaining Gaussian overlap terms.
                    -match_bonus - (-distance_squared / (2.0 * SCORE_SIGMA * SCORE_SIGMA)).dexp()
                } else {
                    1_000_000.0
                };
        }
    }
    assign_rows(&costs)
        .into_iter()
        .enumerate()
        .filter_map(|(row, col)| {
            if col < real_columns && costs[row][col] < 0.0 {
                Some((left[row], right[col]))
            } else {
                None
            }
        })
        .collect()
}

/// Alternate re-fitting the rotation for the current pairing (reusing
/// [`align_coords`]) and re-pairing against the newly rotated coordinates,
/// tracking the best-scoring pairing seen (not just the last — the
/// alternating optimization is not guaranteed monotonic once the pairing
/// set itself changes between iterations).
fn refine(
    coords1: &[[f64; 3]],
    types1: &[MMFF94Type],
    seeded_coords2: &[[f64; 3]],
    types2: &[MMFF94Type],
) -> (Vec<(usize, usize)>, f64) {
    let mut current2 = seeded_coords2.to_vec();
    let mut best_pairs: Vec<(usize, usize)> = Vec::new();
    let mut best_score = f64::NEG_INFINITY;
    let mut prev_score = f64::NEG_INFINITY;

    for _ in 0..MAX_REFINE_ITERS {
        let pairs = global_correspondence(coords1, types1, &current2, types2);
        if pairs.is_empty() {
            break;
        }
        let s = overlap_score(&pairs, coords1, &current2);
        if s > best_score {
            best_score = s;
            best_pairs = pairs.clone();
        }
        if (s - prev_score).abs() < CONVERGENCE_EPS {
            break;
        }
        prev_score = s;

        let sub1: Vec<[f64; 3]> = pairs.iter().map(|&(i, _)| coords1[i]).collect();
        let sub2: Vec<[f64; 3]> = pairs.iter().map(|&(_, j)| current2[j]).collect();
        let result = align_coords(&sub1, &sub2);
        current2 = apply_alignment(&current2, &result);
    }

    (best_pairs, best_score)
}

/// Proper-rotation sign combinations for the 3 principal-axis columns —
/// exactly the 4 (of 8 possible) that preserve det = +1.
const SEED_SIGNS: [[f64; 3]; 4] = [
    [1.0, 1.0, 1.0],
    [1.0, -1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
];

/// Find an atom correspondence between `mol1` and `mol2` given their 3D
/// coordinates, without requiring atom pairs to be known in advance.
///
/// Returns `(mol1_atom_idx, mol2_atom_idx)` pairs. Returns an empty `Vec`
/// when no compatible correspondence is found (e.g. no shared MMFF94 types
/// within the distance cutoff at any seed orientation).
pub(crate) fn correspondence_search(
    mol1: &Molecule,
    coords1: &[[f64; 3]],
    mol2: &Molecule,
    coords2: &[[f64; 3]],
) -> Result<Vec<(usize, usize)>, O3AError> {
    if coords1.len() != mol1.atom_count() || coords2.len() != mol2.atom_count() {
        return Err(O3AError::CoordinateMismatch);
    }
    let types1 = chematic_ff::assign_mmff94_types(mol1).map_err(O3AError::TypeAssignment)?;
    let types2 = chematic_ff::assign_mmff94_types(mol2).map_err(O3AError::TypeAssignment)?;

    let c1 = centroid(coords1);
    let c2 = centroid(coords2);
    let (_, axes1) = jacobi3(covariance_tensor(coords1, c1));
    let (_, axes2) = jacobi3(covariance_tensor(coords2, c2));
    let axes1 = normalize_proper(axes1);
    let axes2 = normalize_proper(axes2);

    let mut best_pairs: Vec<(usize, usize)> = Vec::new();
    let mut best_score = f64::NEG_INFINITY;

    for signs in &SEED_SIGNS {
        let mut axes2_signed = axes2;
        for row in axes2_signed.iter_mut() {
            for (col, sign) in signs.iter().enumerate() {
                row[col] *= sign;
            }
        }
        // Signed rotation matrix mapping mol2's (signed) principal frame
        // onto mol1's: R = axes1 * axes2_signed^T.
        let rot = matmul3(axes1, axes2_signed);
        let seeded = apply_rotation_translation(coords2, c2, rot, c1);

        let (pairs, score) = refine(coords1, &types1, &seeded, &types2);
        if score > best_score {
            best_score = score;
            best_pairs = pairs;
        }
    }

    Ok(best_pairs)
}

// ---------------------------------------------------------------------------
// Full alignment
// ---------------------------------------------------------------------------

/// Result of a full O3A alignment: the atom correspondence found by
/// [`correspondence_search`], its post-alignment Gaussian overlap score, and
/// the Kabsch superposition (via [`crate::align::align_coords`]) of `mol2`
/// onto `mol1`, fit on that correspondence only.
#[derive(Debug, Clone)]
pub struct O3AResult {
    /// `(mol1_atom_idx, mol2_atom_idx)` pairs used for the fit.
    pub pairs: Vec<(usize, usize)>,
    /// Gaussian overlap score of the paired atoms after alignment — higher
    /// means a tighter, more extensive overlap.
    pub score: f64,
    /// Kabsch superposition of `mol2` onto `mol1`, fit on `pairs` only.
    pub alignment: crate::align::AlignResult,
}

/// Find an atom correspondence between `mol1` and `mol2` (via
/// [`correspondence_search`]) and superpose `mol2` onto `mol1` using it.
///
/// A thin wrapper: all algorithmic complexity lives in `correspondence_search`
/// and the already-tested Kabsch alignment in [`crate::align`]. To move all
/// of `mol2`'s coordinates (not just the paired subset) into `mol1`'s frame,
/// apply [`apply_alignment`] to the full `coords2` using `result.alignment`.
pub fn o3a_align(
    mol1: &Molecule,
    coords1: &[[f64; 3]],
    mol2: &Molecule,
    coords2: &[[f64; 3]],
) -> Result<O3AResult, O3AError> {
    let pairs = correspondence_search(mol1, coords1, mol2, coords2)?;
    let sub1: Vec<[f64; 3]> = pairs.iter().map(|&(i, _)| coords1[i]).collect();
    let sub2: Vec<[f64; 3]> = pairs.iter().map(|&(_, j)| coords2[j]).collect();
    let alignment = align_coords(&sub1, &sub2);
    let aligned_sub2 = apply_alignment(&sub2, &alignment);
    let score = overlap_score_paired(&sub1, &aligned_sub2);
    Ok(O3AResult {
        pairs,
        score,
        alignment,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;
    use std::collections::HashSet;

    #[test]
    fn global_assignment_beats_greedy_competition() {
        let atom_type = MMFF94Type::C_sp3;
        let left = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]];
        let right = [[0.05, 0.0, 0.0], [-0.5, 0.0, 0.0]];
        let types = [atom_type; 2];
        let greedy = greedy_correspondence(&left, &types, &right, &types);
        let global = global_correspondence(&left, &types, &right, &types);
        assert_eq!(global, vec![(0, 1), (1, 0)]);
        assert!(overlap_score(&global, &left, &right) > overlap_score(&greedy, &left, &right));
    }

    #[test]
    fn global_assignment_uses_dummy_for_incompatible_atoms() {
        let left = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let right = [[0.1, 0.0, 0.0], [100.0, 0.0, 0.0]];
        let types = [MMFF94Type::C_sp3; 2];
        assert_eq!(
            global_correspondence(&left, &types, &right, &types),
            vec![(0, 0)]
        );
    }

    #[test]
    fn global_assignment_seed_overlap_does_not_retreat_on_symmetry_corpus() {
        for &smiles in ORDER_DEPENDENCE_CORPUS {
            let probe = self_align_baseline_vs_reversed(smiles);
            let types1 = chematic_ff::assign_mmff94_types(&probe.mol1).unwrap();
            let types2 = chematic_ff::assign_mmff94_types(&probe.mol2).unwrap();
            let c1 = centroid(&probe.coords1);
            let c2 = centroid(&probe.coords2);
            let (_, axes1) = jacobi3(covariance_tensor(&probe.coords1, c1));
            let (_, axes2) = jacobi3(covariance_tensor(&probe.coords2, c2));
            let axes1 = normalize_proper(axes1);
            let axes2 = normalize_proper(axes2);
            for signs in SEED_SIGNS {
                let mut signed = axes2;
                for row in &mut signed {
                    for (column, sign) in signs.iter().enumerate() {
                        row[column] *= sign;
                    }
                }
                let seeded =
                    apply_rotation_translation(&probe.coords2, c2, matmul3(axes1, signed), c1);
                let greedy = greedy_correspondence(&probe.coords1, &types1, &seeded, &types2);
                let global = global_correspondence(&probe.coords1, &types1, &seeded, &types2);
                assert!(
                    overlap_score(&global, &probe.coords1, &seeded) + 1e-9
                        >= overlap_score(&greedy, &probe.coords1, &seeded),
                    "{smiles}: global assignment reduced seed overlap"
                );
            }
        }
    }

    #[test]
    fn greedy_correspondence_tie_break_is_content_based() {
        // mol1 has two same-typed atoms (indices 1, 2) confirmed exactly
        // tied for farthest from centroid (500/9 each, by construction).
        // Honest caveat, unlike usr.rs's equivalent test: this specific
        // 3-point configuration passes BOTH before and after the fix (the
        // greedy algorithm still converges to the same final pairing here
        // even when the tied atom is processed in a different order, since
        // there's no actual competition between mol1 atoms for the same
        // mol2 atom in this small a case) -- so this test does not serve as
        // a red-before-green repro the way usr.rs's does. Kept as a
        // regression guard and to document the intent of the fix, not as
        // proof the pre-fix code was reachably wrong.
        let t = MMFF94Type::C_sp3;
        let coords1 = vec![[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 0.0]];
        let types1 = vec![t, t, t];
        let coords2 = vec![[0.1, 0.0, 0.0], [10.0, 0.1, 0.0], [0.1, 10.0, 0.0]];
        let types2 = vec![t, t, t];

        let pairs_a = greedy_correspondence(&coords1, &types1, &coords2, &types2);

        let mut coords1_swapped = coords1.clone();
        coords1_swapped.swap(1, 2);
        let mut types1_swapped = types1.clone();
        types1_swapped.swap(1, 2);
        let mut coords2_swapped = coords2.clone();
        coords2_swapped.swap(1, 2);
        let mut types2_swapped = types2.clone();
        types2_swapped.swap(1, 2);
        let pairs_b = greedy_correspondence(
            &coords1_swapped,
            &types1_swapped,
            &coords2_swapped,
            &types2_swapped,
        );

        // Map indices back through the swap so both results describe pairings
        // in terms of the ORIGINAL (pre-swap) coordinate identities.
        let unswap = |idx: usize| -> usize {
            match idx {
                1 => 2,
                2 => 1,
                other => other,
            }
        };
        let mut normalized_b: Vec<(usize, usize)> = pairs_b
            .iter()
            .map(|&(i, j)| (unswap(i), unswap(j)))
            .collect();
        let mut normalized_a = pairs_a.clone();
        normalized_a.sort_unstable();
        normalized_b.sort_unstable();
        assert_eq!(
            normalized_a, normalized_b,
            "greedy_correspondence's pairing (by coordinate identity) must not depend on array order"
        );
    }

    /// Rotate `coords` by a fixed, non-trivial rotation and translate.
    fn rotate_translate(coords: &[[f64; 3]]) -> Vec<[f64; 3]> {
        // 90° rotation about z, plus a translation.
        coords
            .iter()
            .map(|p| [-p[1] + 5.0, p[0] - 3.0, p[2] + 2.0])
            .collect()
    }

    #[test]
    fn self_alignment_recovers_identity_correspondence() {
        // Asymmetric molecule so the identity mapping is the only sane answer.
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let n = mol.atom_count();
        // Deterministic pseudo-3D layout: spread atoms out non-degenerately.
        let coords1: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let t = i as f64;
                [t * 1.3, (t * 0.7).dsin() * 2.0, (t * 0.3).dcos() * 1.5]
            })
            .collect();
        let coords2 = rotate_translate(&coords1);

        let pairs = correspondence_search(&mol, &coords1, &mol, &coords2).unwrap();
        assert_eq!(pairs.len(), n, "every atom should find its rotated self");
        for (i, j) in &pairs {
            assert_eq!(i, j, "self-alignment must recover the identity mapping");
        }
    }

    #[test]
    fn correspondence_is_injective() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let n = mol.atom_count();
        let coords1: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let t = i as f64;
                [t * 1.3, (t * 0.7).dsin() * 2.0, (t * 0.3).dcos() * 1.5]
            })
            .collect();
        let coords2 = rotate_translate(&coords1);

        let pairs = correspondence_search(&mol, &coords1, &mol, &coords2).unwrap();
        let js: HashSet<usize> = pairs.iter().map(|&(_, j)| j).collect();
        assert_eq!(js.len(), pairs.len(), "no mol2 atom should be used twice");
    }

    #[test]
    fn small_hand_built_system_finds_expected_pairing() {
        // Methanol's two heavy atoms (C, O) have distinct MMFF94 types,
        // giving an unambiguous, hand-verifiable 2-atom pairing.
        let mol = parse("CO").unwrap();
        let coords1: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]];
        let coords2 = rotate_translate(&coords1);

        let pairs = correspondence_search(&mol, &coords1, &mol, &coords2).unwrap();
        assert_eq!(pairs.len(), 2);
        let mut sorted = pairs.clone();
        sorted.sort();
        assert_eq!(sorted, vec![(0, 0), (1, 1)]);
    }

    #[test]
    fn coordinate_mismatch_is_rejected() {
        let mol = parse("CO").unwrap();
        let coords_short = vec![[0.0, 0.0, 0.0]];
        let coords_ok = vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]];
        let err = correspondence_search(&mol, &coords_short, &mol, &coords_ok).unwrap_err();
        assert!(
            err.to_string()
                .contains("coords length does not match atom count")
        );
        assert!(matches!(err, O3AError::CoordinateMismatch));
    }

    #[test]
    fn real_conformer_self_alignment() {
        // Sanity check on a realistically-shaped conformer (not a synthetic
        // sine/cosine layout): generate_coords's actual 3D geometry.
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let coords = crate::dg::generate_coords(&mol);
        let n = mol.atom_count();
        let coords1: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let p = coords.get(chematic_core::AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect();
        let coords2 = rotate_translate(&coords1);

        let pairs = correspondence_search(&mol, &coords1, &mol, &coords2).unwrap();
        assert_eq!(pairs.len(), n, "every atom should find its rotated self");
        for (i, j) in &pairs {
            assert_eq!(
                i, j,
                "self-alignment on a real conformer must recover identity"
            );
        }

        // The correspondence should also be a Kabsch-perfect fit: aligning
        // coords1 onto coords2 via the found pairing must reproduce
        // rotate_translate's rigid transform (near-zero RMSD).
        let sub1: Vec<[f64; 3]> = pairs.iter().map(|&(i, _)| coords1[i]).collect();
        let sub2: Vec<[f64; 3]> = pairs.iter().map(|&(_, j)| coords2[j]).collect();
        let result = align_coords(&sub1, &sub2);
        assert!(result.rmsd < 1e-6, "rmsd={} should be ~0", result.rmsd);
    }

    #[test]
    fn benzene_toluene_shared_ring_atoms_match() {
        // Toluene = benzene + methyl. The 6 aromatic ring atoms should all
        // find compatible-type matches in toluene's own ring.
        let benzene = parse("c1ccccc1").unwrap();
        let toluene = parse("Cc1ccccc1").unwrap();
        let bc = crate::dg::generate_coords(&benzene);
        let tc = crate::dg::generate_coords(&toluene);
        let bn = benzene.atom_count();
        let tn = toluene.atom_count();
        let coords1: Vec<[f64; 3]> = (0..bn)
            .map(|i| {
                let p = bc.get(chematic_core::AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect();
        let coords2: Vec<[f64; 3]> = (0..tn)
            .map(|i| {
                let p = tc.get(chematic_core::AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect();

        let pairs = correspondence_search(&benzene, &coords1, &toluene, &coords2).unwrap();
        assert!(
            pairs.len() >= 4,
            "expected most of benzene's 6 ring atoms to find a match in toluene's ring, got {}",
            pairs.len()
        );
    }

    #[test]
    fn unrelated_molecules_still_return_a_correspondence() {
        // Different molecules can still share some compatible atom types;
        // this just checks the function doesn't panic and returns pairs
        // that are individually within the distance cutoff.
        let mol1 = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let mol2 = parse("c1ccccc1").unwrap();
        let n1 = mol1.atom_count();
        let n2 = mol2.atom_count();
        let coords1: Vec<[f64; 3]> = (0..n1).map(|i| [i as f64 * 1.4, 0.0, 0.0]).collect();
        let coords2: Vec<[f64; 3]> = (0..n2).map(|i| [i as f64 * 1.4, 0.5, 0.0]).collect();

        let pairs = correspondence_search(&mol1, &coords1, &mol2, &coords2).unwrap();
        for &(i, j) in &pairs {
            assert!(i < n1);
            assert!(j < n2);
        }
    }

    fn coords3d_to_vec(coords: &crate::coords::Coords3D, n: usize) -> Vec<[f64; 3]> {
        (0..n)
            .map(|i| {
                let p = coords.get(chematic_core::AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect()
    }

    #[test]
    fn o3a_align_self_rotation_recovers_zero_rmsd() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let n = mol.atom_count();
        let coords1 = coords3d_to_vec(&crate::dg::generate_coords(&mol), n);
        let coords2 = rotate_translate(&coords1);

        let result = o3a_align(&mol, &coords1, &mol, &coords2).unwrap();
        assert_eq!(
            result.pairs.len(),
            n,
            "every atom should find its rotated self"
        );
        for (i, j) in &result.pairs {
            assert_eq!(i, j, "self-alignment must recover the identity mapping");
        }
        assert!(
            result.alignment.rmsd < 1e-6,
            "rmsd={} should be ~0",
            result.alignment.rmsd
        );
    }

    #[test]
    fn o3a_align_rmsd_matches_direct_align_coords_for_known_correspondence() {
        // Methanol's C/O have distinct MMFF94 types, so correspondence_search
        // is guaranteed to recover the identity pairing — giving a known
        // correspondence to compare against a direct align_coords call.
        let mol = parse("CO").unwrap();
        let coords1: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]];
        let coords2 = rotate_translate(&coords1);

        let result = o3a_align(&mol, &coords1, &mol, &coords2).unwrap();
        let direct = align_coords(&coords1, &coords2);
        assert!(
            (result.alignment.rmsd - direct.rmsd).abs() < 1e-9,
            "o3a_align rmsd {} should match direct align_coords rmsd {} for the same correspondence",
            result.alignment.rmsd,
            direct.rmsd
        );
    }

    #[test]
    fn o3a_align_scaffold_shared_pair_scores_higher_than_unrelated() {
        // Toluene shares benzene's exact aromatic ring (same atom types, same
        // flat hexagon geometry). Cyclohexane is a genuine negative control:
        // saturated sp3 carbons (a different MMFF94 type entirely) in a
        // puckered, non-planar ring — not just a substituent swap on the same
        // scaffold. (Aspirin would be a *bad* negative control here since it
        // contains benzene's own aromatic ring as a substructure.)
        let benzene = parse("c1ccccc1").unwrap();
        let toluene = parse("Cc1ccccc1").unwrap();
        let cyclohexane = parse("C1CCCCC1").unwrap();

        let bc = coords3d_to_vec(&crate::dg::generate_coords(&benzene), benzene.atom_count());
        let tc = coords3d_to_vec(&crate::dg::generate_coords(&toluene), toluene.atom_count());
        let cc = coords3d_to_vec(
            &crate::dg::generate_coords(&cyclohexane),
            cyclohexane.atom_count(),
        );

        let related = o3a_align(&benzene, &bc, &toluene, &tc).unwrap();
        let unrelated = o3a_align(&benzene, &bc, &cyclohexane, &cc).unwrap();

        assert!(
            related.score > unrelated.score,
            "benzene/toluene (shared aromatic ring) score {} should exceed benzene/cyclohexane score {}",
            related.score,
            unrelated.score
        );
    }

    #[test]
    fn o3a_align_pairs_are_injective() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let n = mol.atom_count();
        let coords1: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let t = i as f64;
                [t * 1.3, (t * 0.7).dsin() * 2.0, (t * 0.3).dcos() * 1.5]
            })
            .collect();
        let coords2 = rotate_translate(&coords1);

        let result = o3a_align(&mol, &coords1, &mol, &coords2).unwrap();
        let js: HashSet<usize> = result.pairs.iter().map(|&(_, j)| j).collect();
        assert_eq!(
            js.len(),
            result.pairs.len(),
            "no mol2 atom should be used twice"
        );
    }

    /// Rebuild `mol` with atoms relabeled by `perm` (perm[new_idx] = old_idx),
    /// permuting `coords` identically so each new-index atom keeps its own
    /// position -- same technique used repeatedly elsewhere this round for
    /// permutation-invariance probes (see
    /// [[feedback_permutation_invariance_test_template]]).
    fn permute_mol_and_coords(
        mol: &Molecule,
        coords: &[[f64; 3]],
        perm: &[usize],
    ) -> (Molecule, Vec<[f64; 3]>) {
        use chematic_core::{AtomIdx, MoleculeBuilder};
        let mut old_to_new = vec![0u32; perm.len()];
        for (new_idx, &old_idx) in perm.iter().enumerate() {
            old_to_new[old_idx] = new_idx as u32;
        }
        let mut builder = MoleculeBuilder::new();
        for &old_idx in perm {
            builder.add_atom(mol.atom(AtomIdx(old_idx as u32)).clone());
        }
        for (_, bond) in mol.bonds() {
            let a = AtomIdx(old_to_new[bond.atom1.0 as usize]);
            let b = AtomIdx(old_to_new[bond.atom2.0 as usize]);
            let _ = builder.add_bond(a, b, bond.order);
        }
        let permuted_coords: Vec<[f64; 3]> = perm.iter().map(|&old_idx| coords[old_idx]).collect();
        (builder.build(), permuted_coords)
    }

    /// Molecules with internal symmetry likely to expose greedy
    /// order-dependence (`generate_coords` is rule-based, not iteratively
    /// randomized, so it plausibly places symmetric substituents at exactly
    /// equidistant positions), plus one asymmetric control.
    const ORDER_DEPENDENCE_CORPUS: &[&str] = &[
        "CC(=O)Oc1ccccc1C(=O)O", // asymmetric control
        "c1ccccc1",
        "Cc1ccc(C)cc1",
        "CC(=O)CC(=O)C",
        "CC(C)(C)C",
        "CC(C)(C)OC(=O)N",
        "c1ccc2ccccc2c1",
        "CC(C)(C)c1ccc(C(C)(C)C)cc1",
        "O=C1CCC(=O)N1",
        "CC(C)(C)C(=O)C(C)(C)C",
    ];

    /// mol1/coords1 (unpermuted) plus mol2/coords2 in both its original and
    /// atom-order-reversed forms, for `self_align_baseline_vs_reversed`.
    struct OrderDependenceProbe {
        mol1: Molecule,
        coords1: Vec<[f64; 3]>,
        mol2: Molecule,
        coords2: Vec<[f64; 3]>,
        mol2_rev: Molecule,
        coords2_rev: Vec<[f64; 3]>,
        /// `correspondence_search(mol1, coords1, mol2, coords2)`, sorted.
        baseline_pairs: Vec<(usize, usize)>,
        /// `correspondence_search(mol1, coords1, mol2_rev, coords2_rev)`,
        /// mapped back to original mol2 atom identities and sorted.
        reversed_pairs: Vec<(usize, usize)>,
    }

    /// Self-align mol1 against a rotated+translated mol2, once with mol2's
    /// atoms in original order and once reversed (coords permuted to match).
    fn self_align_baseline_vs_reversed(smi: &str) -> OrderDependenceProbe {
        let mol1 = parse(smi).unwrap();
        let n = mol1.atom_count();
        let coords1_map = crate::dg::generate_coords(&mol1);
        let coords1: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let p = coords1_map.get(chematic_core::AtomIdx(i as u32));
                [p.x, p.y, p.z]
            })
            .collect();

        let mol2 = parse(smi).unwrap();
        let coords2: Vec<[f64; 3]> = rotate_translate(&coords1);
        let baseline = correspondence_search(&mol1, &coords1, &mol2, &coords2).unwrap();

        let perm: Vec<usize> = (0..n).rev().collect();
        let (mol2_rev, coords2_rev) = permute_mol_and_coords(&mol2, &coords2, &perm);
        let reversed = correspondence_search(&mol1, &coords1, &mol2_rev, &coords2_rev).unwrap();
        let mut reversed_pairs: Vec<(usize, usize)> =
            reversed.iter().map(|&(i, j)| (i, perm[j])).collect();
        let mut baseline_pairs = baseline.clone();
        baseline_pairs.sort_unstable();
        reversed_pairs.sort_unstable();

        OrderDependenceProbe {
            mol1,
            coords1,
            mol2,
            coords2,
            mol2_rev,
            coords2_rev,
            baseline_pairs,
            reversed_pairs,
        }
    }

    #[test]
    fn correspondence_search_alignment_score_is_order_independent() {
        // Reordering mol2's atoms must preserve the fitted overlap score.
        for &smi in ORDER_DEPENDENCE_CORPUS {
            let p = self_align_baseline_vs_reversed(smi);
            let base_score = o3a_align(&p.mol1, &p.coords1, &p.mol2, &p.coords2)
                .unwrap()
                .score;
            let rev_score = o3a_align(&p.mol1, &p.coords1, &p.mol2_rev, &p.coords2_rev)
                .unwrap()
                .score;
            assert!(
                (base_score - rev_score).abs() < 1e-9,
                "{smi}: alignment score depends on mol2 atom order: {base_score} vs {rev_score}"
            );
        }
    }

    #[test]
    fn correspondence_search_pairs_identity_is_order_independent() {
        for &smi in ORDER_DEPENDENCE_CORPUS {
            let p = self_align_baseline_vs_reversed(smi);
            assert_eq!(
                p.baseline_pairs, p.reversed_pairs,
                "{smi}: correspondence_search's pairing (by atom identity) \
                 changed when mol2's atom order was reversed"
            );
            let n = p.mol2.atom_count();
            let permutations = [
                (1..n).chain(0..1).collect::<Vec<_>>(),
                (0..n).step_by(2).chain((1..n).step_by(2)).collect(),
            ];
            for permutation in permutations {
                let (mol2, coords2) = permute_mol_and_coords(&p.mol2, &p.coords2, &permutation);
                let pairs = correspondence_search(&p.mol1, &p.coords1, &mol2, &coords2).unwrap();
                let mut original_pairs = pairs
                    .into_iter()
                    .map(|(i, j)| (i, permutation[j]))
                    .collect::<Vec<_>>();
                original_pairs.sort_unstable();
                assert_eq!(p.baseline_pairs, original_pairs, "{smi}: {permutation:?}");
            }
        }
    }

    #[test]
    fn rdkit_embedded_pxylene_correspondence_survives_renumbering() {
        // Heavy-atom coordinates from RDKit 2026.09.1 ETKDGv3, seed 42.
        // RDKit GetO3A itself returns different atom maps for this molecule
        // after probe atom reversal; this checks the same fixed input geometry.
        let molecule = parse("Cc1ccc(C)cc1").unwrap();
        let reference = vec![
            [
                2.87688384965826,
                -0.002092989212400923,
                -0.08653065571627763,
            ],
            [
                1.375506197969415,
                0.031376908900133456,
                -0.04181616373635037,
            ],
            [0.6837946804515053, 1.2162483338831502, 0.11090569228195972],
            [-0.6968589029854234, 1.1964285865863726, 0.14633051634528366],
            [
                -1.3672286381292256,
                -0.007579808880645602,
                0.028604200033073776,
            ],
            [
                -2.8659967768237866,
                0.004430820785118485,
                0.07084524663778523,
            ],
            [
                -0.6983511002657733,
                -1.2008939331499826,
                -0.12445273185002001,
            ],
            [
                0.6864109596823279,
                -1.1682896941005392,
                -0.15854404629064228,
            ],
        ];
        let probe = reference
            .iter()
            .map(|&[x, y, z]| [y + 1.0, -x + 2.0, z + 3.0])
            .collect::<Vec<_>>();
        let mut baseline = correspondence_search(&molecule, &reference, &molecule, &probe).unwrap();
        baseline.sort_unstable();
        let permutation = (0..molecule.atom_count()).rev().collect::<Vec<_>>();
        let (reversed_molecule, reversed_probe) =
            permute_mol_and_coords(&molecule, &probe, &permutation);
        let mut reversed =
            correspondence_search(&molecule, &reference, &reversed_molecule, &reversed_probe)
                .unwrap()
                .into_iter()
                .map(|(i, j)| (i, permutation[j]))
                .collect::<Vec<_>>();
        reversed.sort_unstable();
        assert_eq!(baseline, reversed);
    }
}
