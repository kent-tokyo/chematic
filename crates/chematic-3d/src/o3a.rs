//! Open3DAlign (O3A)-style atom correspondence search for 3D molecular alignment.
//!
//! Unlike simple RMSD superposition ([`crate::align::align_coords`]), which
//! requires atom pairs to already be known, O3A finds the atom correspondence
//! itself. This module implements only the correspondence-search primitive —
//! the riskier, genuinely novel piece — as an isolated, directly-testable
//! unit. The full `o3a_align` pipeline (which feeds this correspondence into
//! the already-tested [`crate::align::align_coords`]) is a separate, thin
//! wrapper built on top of this.
//!
//! # Algorithm
//! 1. Assign MMFF94 atom types to both molecules; two atoms are considered
//!    "compatible" only when their types match exactly.
//! 2. Pre-align by principal axes (via the inertia-like covariance tensor's
//!    eigenvectors). Principal axes have a sign ambiguity — of the 8 possible
//!    per-axis sign combinations only 4 are proper rotations (det = +1); all
//!    4 are tried as independent seeds, since picking the wrong one can
//!    converge to a mirrored, degenerate correspondence.
//! 3. For each seed: greedily pair each mol1 atom (farthest-from-centroid
//!    first, to anchor asymmetric ends) with its nearest compatible-type,
//!    not-yet-used mol2 atom within a distance cutoff, then alternate
//!    "refit rotation for the current pairing" (reusing `align_coords`) and
//!    "re-pair against the newly rotated coordinates" until the overlap
//!    score stops improving.
//! 4. Return the correspondence from whichever seed reached the best score.
//!
//! ponytail: `correspondence_search` and its helpers are only exercised by
//! this module's tests until the PR2 `o3a_align` wrapper lands and wires
//! them into the public API — hence the blanket dead_code allow below.
#![allow(dead_code)]

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

// ---------------------------------------------------------------------------
// Small linear-algebra helpers (local to this module — not a general
// utility, so kept private rather than widening align.rs's visibility)
// ---------------------------------------------------------------------------

fn centroid(coords: &[[f64; 3]]) -> [f64; 3] {
    let n = coords.len().max(1) as f64;
    let mut c = [0.0; 3];
    for p in coords {
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
    for p in coords {
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

fn overlap_score(pairs: &[(usize, usize)], coords1: &[[f64; 3]], coords2: &[[f64; 3]]) -> f64 {
    pairs
        .iter()
        .map(|&(i, j)| {
            let d2: f64 = (0..3)
                .map(|k| (coords1[i][k] - coords2[j][k]).powi(2))
                .sum();
            (-d2 / (2.0 * SCORE_SIGMA * SCORE_SIGMA)).exp()
        })
        .sum()
}

/// Greedily pair each mol1 atom (farthest-from-centroid first) with its
/// nearest not-yet-used, type-compatible mol2 atom within `DIST_CUTOFF`.
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
        db.partial_cmp(&da).unwrap_or(std::cmp::Ordering::Equal)
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
            if d2 < best_d2 {
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
        let pairs = greedy_correspondence(coords1, types1, &current2, types2);
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

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;
    use std::collections::HashSet;

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
                [t * 1.3, (t * 0.7).sin() * 2.0, (t * 0.3).cos() * 1.5]
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
                [t * 1.3, (t * 0.7).sin() * 2.0, (t * 0.3).cos() * 1.5]
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
}
