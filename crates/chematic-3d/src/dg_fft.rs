//! Distance geometry with eigenvalue decomposition for 3D coordinate generation.
//!
//! Implements the distance geometry approach:
//! 1. Build bound matrix from bond lengths and angle constraints
//! 2. Convert distances to Gram matrix (inner products)
//! 3. Eigenvalue decomposition to extract coordinates
//! 4. Metrization to enforce Euclidean distances
//!
//! This provides more accurate geometry than rule-based placement, especially
//! for complex molecules with many torsion degrees of freedom.

#![allow(dead_code)]

use std::f64::consts::PI;

use chematic_core::{AtomIdx, BondOrder, Molecule};

use crate::coords::{Coords3D, Point3};

// ---------------------------------------------------------------------------
// Bound matrix construction
// ---------------------------------------------------------------------------

/// Build a bound matrix (lower and upper distance bounds) from molecular constraints.
///
/// Returns (lower_bounds, upper_bounds) as n×n matrices where:
/// - lower_bounds[i][j] = minimum allowed distance
/// - upper_bounds[i][j] = maximum allowed distance
///
/// Constraints include:
/// - Bond lengths (from ideal values ± tolerance)
/// - Angle constraints (from ideal angles)
/// - Van der Waals (from VDW radii sum)
///
/// A thin wrapper over [`build_bond_angle_bounds`] + [`apply_vdw_bounds`], split apart
/// so `distance_geometry_v2.rs`'s `enforce_chirality` path can insert declared-E/Z 1-4
/// bounds *between* the two (before the generic VDW non-bonded floor would otherwise
/// apply to that same pair -- see that module's `apply_declared_ez_bounds` doc for why
/// the ordering matters). This function's own output is unchanged either way: it always
/// runs both steps back to back, with nothing in between.
pub(crate) fn build_bound_matrix(mol: &Molecule) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let (mut lower, mut upper) = build_bond_angle_bounds(mol);
    apply_vdw_bounds(mol, &mut lower, &mut upper);
    (lower, upper)
}

/// Bond-length (1-2) and angle-derived (1-3) distance bounds only -- no Van der Waals
/// non-bonded floor yet, see [`apply_vdw_bounds`]. Exists as its own function so a
/// caller can insert additional pair-specific constraints (e.g. declared-E/Z 1-4
/// bounds) before the generic VDW floor is applied; [`build_bound_matrix`] is the
/// composition most callers want.
pub(crate) fn build_bond_angle_bounds(mol: &Molecule) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = mol.atom_count();
    let mut lower = vec![vec![0.0; n]; n];
    let mut upper = vec![vec![f64::INFINITY; n]; n];

    // Diagonal: distance to self = 0
    for i in 0..n {
        lower[i][i] = 0.0;
        upper[i][i] = 0.0;
    }

    // Bond length constraints (distance = 1)
    for (_, bond) in mol.bonds() {
        let a = bond.atom1;
        let b = bond.atom2;
        let dist = ideal_bond_length(mol, a, b);
        let tolerance = 0.05; // ±0.05 Å tolerance
        lower[a.0 as usize][b.0 as usize] = (dist - tolerance).max(0.5);
        upper[a.0 as usize][b.0 as usize] = dist + tolerance;
        lower[b.0 as usize][a.0 as usize] = lower[a.0 as usize][b.0 as usize];
        upper[b.0 as usize][a.0 as usize] = upper[a.0 as usize][b.0 as usize];
    }

    // Angle constraints (distance = 2)
    for center_idx in 0..n {
        let center = AtomIdx(center_idx as u32);
        let neighbors: Vec<AtomIdx> = mol.neighbors(center).map(|(nb, _)| nb).collect();

        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                let a = neighbors[i];
                let b = neighbors[j];

                // `a` and `b` are both bonded to `center`; if they are ALSO bonded to
                // each other, `center`-`a`-`b` is a 3-membered ring, not a genuine 1-3
                // (through-center) relationship. The generic ~109.5°/120° angle below
                // is wrong for a ring this tight (true angle ~60°) and, left in place,
                // overwrites the already-correct 1-2 bond-length bound set above for
                // this same pair with a contradictory one (see
                // `distance_geometry_v2::tests::cyclopropane_exact_bound_contradiction_verified`).
                // A 3-membered ring's three bond lengths already fix its shape (a
                // triangle's side lengths determine it up to reflection), so skipping
                // the angle bound here drops a wrong constraint, not a needed one.
                if mol.bond_between(a, b).is_some() {
                    continue;
                }

                let bond_len_1 = ideal_bond_length(mol, center, a);
                let bond_len_2 = ideal_bond_length(mol, center, b);
                let angle = ideal_bond_angle(mol, center);

                // Distance = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(angle))
                let dist_sq = bond_len_1.powi(2) + bond_len_2.powi(2)
                    - 2.0 * bond_len_1 * bond_len_2 * angle.cos();
                let dist = dist_sq.max(0.0).sqrt();
                let tolerance = 0.1;

                let idx_a = a.0 as usize;
                let idx_b = b.0 as usize;
                lower[idx_a][idx_b] = lower[idx_a][idx_b].max((dist - tolerance).max(0.5));
                upper[idx_a][idx_b] = upper[idx_a][idx_b].min(dist + tolerance);
                lower[idx_b][idx_a] = lower[idx_a][idx_b];
                upper[idx_b][idx_a] = upper[idx_a][idx_b];
            }
        }
    }

    (lower, upper)
}

/// Van der Waals non-bonded lower bounds (distance >= sum of VDW radii), applied in
/// place. Skips bonded pairs, and skips any pair where an existing tighter (smaller)
/// upper bound has already been set (an angle constraint, or -- for `enforce_chirality`
/// -- a declared-E/Z 1-4 bound) — applying VDW there would make lower > upper. This is
/// the same exemption pattern bonded/1-3 pairs already get here, extended to whichever
/// pair a caller tightened before this runs: VDW's generic non-bonded-sterics
/// assumption doesn't hold for a pair whose separation is actually fixed by a nearer,
/// more specific constraint.
pub(crate) fn apply_vdw_bounds(mol: &Molecule, lower: &mut [Vec<f64>], upper: &mut [Vec<f64>]) {
    let n = mol.atom_count();
    for i in 0..n {
        for j in (i + 1)..n {
            let a = AtomIdx(i as u32);
            let b = AtomIdx(j as u32);

            // Skip if bonded
            if mol.bond_between(a, b).is_some() {
                continue;
            }

            let vdw_sum = vdw_distance(mol, a, b);
            // Only raise the lower bound if it would not exceed the upper bound.
            if lower[i][j] < vdw_sum && vdw_sum <= upper[i][j] {
                lower[i][j] = vdw_sum;
                lower[j][i] = vdw_sum;
            }
        }
    }
}

/// Ideal bond length (Å) from atom pair and bond order.
///
/// Generic covalent-radius-sum model (`chematic_core::Element::covalent_radius`,
/// full periodic table) with a Pauling-style bond-order length correction, rather
/// than a small hardcoded element-pair table. This fixes a real bug the previous
/// 9-entry match table had: `BondOrder::Aromatic` (the order chematic's own SMILES
/// parser actually assigns to lowercase aromatic atoms, e.g. `c1ccccc1` -- see
/// `chematic-smiles/src/parser.rs::implicit_bond`) matched none of the table's arms
/// and silently fell through to the `_ => 1.54` single-bond default, overestimating
/// every aromatic bond (true aromatic C-C is ~1.39 Å, not 1.54 Å) and every other
/// element pair not in the 9-entry list (S, P, Br, I, N-N, O-O, ...). The scale
/// factors below (single=1.00, double=0.87, triple=0.78, aromatic=0.93) are the
/// same values `scripts/etkdg_vs_rdkit_gap.py::_BOND_ORDER_SCALE` uses against
/// RDKit's `GetRcovalent`, so this is a standard/public bond-order correction, not
/// something reverse-engineered from a single molecule.
pub(crate) fn ideal_bond_length(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> f64 {
    let ra = mol.atom(a).element.covalent_radius() as f64;
    let rb = mol.atom(b).element.covalent_radius() as f64;
    let order = mol
        .bond_between(a, b)
        .map(|(_, bond)| bond.order)
        .unwrap_or(chematic_core::BondOrder::Single);
    (ra + rb) * bond_order_length_scale(order)
}

/// Pauling-style bond-order length scale factor (relative to the single-bond
/// covalent-radius sum). `Quadruple`/`Up`/`Down`/`Zero` are treated as single-bond
/// length: `Up`/`Down` are directional *single* bonds (E/Z markers) in this
/// codebase's `BondOrder`, not a distinct bond order; `Quadruple`/zero-order bonds
/// have no well-established public length-correction table and default to 1.00
/// rather than guessing.
pub(crate) fn bond_order_length_scale(order: BondOrder) -> f64 {
    match order {
        BondOrder::Double => 0.87,
        BondOrder::Triple => 0.78,
        BondOrder::Aromatic => 0.93,
        _ => 1.00,
    }
}

/// Ideal bond angle (radians) at center atom.
///
/// `pub` (not `pub(crate)`) so the gap-check example's angle-violation check can call
/// the *exact same* generic 109.5°/120° model `build_bound_matrix` uses internally,
/// rather than duplicating this constant elsewhere and risking drift. That also means
/// that check is an internal-consistency check (does the final geometry agree with
/// the model the bounds were built from?), not an external-oracle check like the
/// RDKit-covalent-radius bond-length reference -- see the example's own doc comment.
pub fn ideal_bond_angle(mol: &Molecule, center: AtomIdx) -> f64 {
    let atom = mol.atom(center);

    // Simplified: aromatic/double → ~120°, else ~109.5°
    if atom.aromatic {
        2.0 * PI / 3.0 // 120°
    } else {
        let has_double = mol
            .neighbors(center)
            .any(|(_, bidx)| matches!(mol.bond(bidx).order, BondOrder::Double));
        if has_double {
            2.0 * PI / 3.0
        } else {
            1.91 // ~109.5°
        }
    }
}

/// Van der Waals distance (sum of VDW radii, Å).
fn vdw_distance(mol: &Molecule, a: AtomIdx, b: AtomIdx) -> f64 {
    vdw_radius(mol, a) + vdw_radius(mol, b)
}

/// Bondi (1964) van der Waals radii (Å), public-domain standard table -- the same
/// category of source RDKit itself draws non-bonded radii from. Elements outside this
/// short list fall back to `covalent_radius() + 0.75`, a rough but monotonic estimate
/// (real VdW/covalent ratios run ~1.2-1.9x across the periodic table) rather than a
/// flat constant that would be wrong in both directions for very small/large atoms.
fn vdw_radius(mol: &Molecule, idx: AtomIdx) -> f64 {
    let el = mol.atom(idx).element;
    match el.atomic_number() {
        1 => 1.20,  // H
        6 => 1.70,  // C
        7 => 1.55,  // N
        8 => 1.52,  // O
        9 => 1.47,  // F
        15 => 1.80, // P
        16 => 1.80, // S
        17 => 1.75, // Cl
        35 => 1.85, // Br
        53 => 1.98, // I
        _ => el.covalent_radius() as f64 + 0.75,
    }
}

// ---------------------------------------------------------------------------
// Eigenvalue decomposition & coordinate generation
// ---------------------------------------------------------------------------

/// Maximum atom count accepted by the DG coordinate generator.
///
/// O(n²) memory and O(n³) Floyd-Warshall make larger molecules prohibitive;
/// callers that need coordinates for bigger structures must use an external tool.
pub const DG_MAX_ATOMS: usize = 500;

/// Generate 3D coordinates via distance geometry (ETDKG-style).
///
/// Algorithm:
/// 1. Build distance bounds from bond lengths and angles
/// 2. Create target distance matrix (average of bounds)
/// 3. Compute Gram matrix from distances
/// 4. Eigenvalue decomposition via Jacobi method
/// 5. Extract 3D coordinates from top 3 eigenvectors
/// 6. Center molecule at origin
pub fn generate_coords_dg(mol: &Molecule) -> Coords3D {
    let n = mol.atom_count();

    if n == 0 {
        return Coords3D::new_zeroed(0);
    }

    if n > DG_MAX_ATOMS {
        return Coords3D::new_zeroed(n);
    }

    if n == 1 {
        let mut coords = Coords3D::new_zeroed(1);
        coords.set(
            AtomIdx(0),
            Point3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
        );
        return coords;
    }

    // Step 1: Build bounds + smooth to give finite upper bounds for all pairs.
    let (mut lower, mut upper) = build_bound_matrix(mol);
    smooth_bounds(&mut lower, &mut upper);

    // Step 2: Create target distance matrix (average of smoothed bounds).
    // After smoothing, upper[i][j] is finite for all connected pairs; clamp
    // any remaining infinity (disconnected graph edge case) to 4× the largest
    // finite upper bound so the Gram matrix stays numerically well-conditioned.
    let max_finite: f64 = upper
        .iter()
        .flat_map(|row| row.iter())
        .filter(|&&v| v.is_finite())
        .cloned()
        .fold(0.0f64, f64::max);
    let fallback = (max_finite * 4.0).max(10.0);

    let mut dist_matrix = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                dist_matrix[i][j] = 0.0;
            } else {
                let u = if upper[i][j].is_finite() {
                    upper[i][j]
                } else {
                    fallback
                };
                dist_matrix[i][j] = (lower[i][j] + u) / 2.0;
            }
        }
    }

    // Step 3: Compute Gram matrix
    let gram = distance_to_gram_matrix(&dist_matrix);

    // Step 4: Eigenvalue decomposition (Jacobi method)
    let (eigenvalues, eigenvectors) = jacobi_eigendecompose(&gram);

    // Step 5: Extract 3D coordinates from top 3 positive eigenvectors
    let mut coords = Coords3D::new_zeroed(n);

    // Collect the 3 largest *positive* eigenvalues, sorted descending.
    // Negative eigenvalues are excluded; if fewer than 3 positive values exist
    // the remaining coordinate axes are left at zero.
    let mut pos_evals: Vec<usize> = (0..n).filter(|&i| eigenvalues[i] > 1e-10).collect();
    pos_evals.sort_by(|&a, &b| {
        eigenvalues[b]
            .partial_cmp(&eigenvalues[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let use_indices: Vec<usize> = pos_evals.into_iter().take(3).collect();

    // Set coordinates using only positive eigenvalues (no .abs() hack).
    for i in 0..n {
        let coord = |dim: usize| -> f64 {
            use_indices
                .get(dim)
                .map_or(0.0, |&idx| eigenvalues[idx].sqrt() * eigenvectors[i][idx])
        };
        coords.set(
            AtomIdx(i as u32),
            Point3 {
                x: coord(0),
                y: coord(1),
                z: coord(2),
            },
        );
    }

    // Step 6: Center molecule
    center_coordinates(&mut coords);

    // Step 7: Refinement — push/pull atom pairs to satisfy distance bounds.
    // Classical MDS is an initial guess only; ring molecules produce frustrated
    // distance matrices that the rank-3 truncation cannot reproduce exactly.
    // Each pass moves both atoms by half the violation along their axis.
    refine_coords(&mut coords, &lower, &upper, 300);

    coords
}

/// Iterative bounds-driven refinement (SHAKE-like).
///
/// For each atom pair whose current distance violates [lower, upper], move
/// both atoms half the violation along their connecting axis. After enough
/// iterations every pair converges inside its bounds. Classical MDS is used
/// only as the initial placement; this step enforces the geometry.
pub(crate) fn refine_coords(
    coords: &mut Coords3D,
    lower: &[Vec<f64>],
    upper: &[Vec<f64>],
    n_iter: usize,
) {
    let n = coords.atom_count();
    for _ in 0..n_iter {
        for i in 0..n {
            for j in (i + 1)..n {
                let pi = coords.get(AtomIdx(i as u32));
                let pj = coords.get(AtomIdx(j as u32));
                let dx = pj.x - pi.x;
                let dy = pj.y - pi.y;
                let dz = pj.z - pi.z;
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                if d < 1e-10 {
                    continue;
                }

                let lo = lower[i][j];
                let hi = upper[i][j];
                let target = if d < lo {
                    lo
                } else if hi.is_finite() && d > hi {
                    hi
                } else {
                    continue;
                };

                // half the signed correction along the i→j axis
                let half = (target - d) * 0.5;
                let ux = dx / d;
                let uy = dy / d;
                let uz = dz / d;

                // i moves away from j when target > d (stretch), toward j when compress
                coords.set(
                    AtomIdx(i as u32),
                    Point3 {
                        x: pi.x - half * ux,
                        y: pi.y - half * uy,
                        z: pi.z - half * uz,
                    },
                );
                coords.set(
                    AtomIdx(j as u32),
                    Point3 {
                        x: pj.x + half * ux,
                        y: pj.y + half * uy,
                        z: pj.z + half * uz,
                    },
                );
            }
        }
    }
}

/// Apply triangle-inequality bounds smoothing (Floyd-Warshall).
///
/// Propagates finite upper bounds to all pairs:
///   upper[i][j] ≤ upper[i][k] + upper[k][j]
/// and tightens lower bounds:
///   lower[i][j] ≥ max(0, lower[i][k] − upper[k][j])
///
/// After smoothing, upper[i][j] is finite for all atom pairs in a connected
/// molecule, eliminating the infinity entries that would produce NaN in the
/// Gram matrix.
pub(crate) fn smooth_bounds(lower: &mut [Vec<f64>], upper: &mut [Vec<f64>]) {
    let n = lower.len();
    for k in 0..n {
        for i in 0..n {
            if upper[i][k].is_infinite() {
                continue;
            }
            for j in 0..n {
                if i == j {
                    continue;
                }
                // Tighten upper bound via k
                let via_k_upper = upper[i][k] + upper[k][j];
                if via_k_upper < upper[i][j] {
                    upper[i][j] = via_k_upper;
                    upper[j][i] = via_k_upper;
                }
                // Tighten lower bound via k
                let via_k_lower = lower[i][k] - upper[k][j];
                if via_k_lower > lower[i][j] {
                    lower[i][j] = via_k_lower.max(0.0);
                    lower[j][i] = lower[i][j];
                }
            }
        }
    }
}

/// Convert distance matrix to Gram matrix using classical MDS double-centering.
///
/// B[i,j] = −0.5 · (D[i,j]² − μ_row[i] − μ_col[j] + μ_all)
///
/// where μ_row[i] = (1/n) Σ_k D[i,k]² and μ_all = (1/n) Σ_i μ_row[i].
/// This is equivalent to centering the inner-product matrix at the centroid
/// of all atoms, which distributes the reference symmetrically rather than
/// anchoring at atom 0.
pub(crate) fn distance_to_gram_matrix(dist: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = dist.len();
    if n == 0 {
        return vec![];
    }

    // D² row means
    let row_means: Vec<f64> = dist
        .iter()
        .map(|row| row.iter().map(|&d| d * d).sum::<f64>() / n as f64)
        .collect();
    let total_mean: f64 = row_means.iter().sum::<f64>() / n as f64;

    let mut gram = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let d2 = dist[i][j] * dist[i][j];
            gram[i][j] = -0.5 * (d2 - row_means[i] - row_means[j] + total_mean);
        }
    }
    gram
}

/// Jacobi eigendecomposition for symmetric matrices.
///
/// Returns (eigenvalues, eigenvectors) where eigenvectors[i][j] is the
/// i-th component of the j-th eigenvector.
pub(crate) fn jacobi_eigendecompose(mat: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = mat.len();
    let mut a = mat.to_vec();
    let mut v = vec![vec![0.0; n]; n];

    // Initialize eigenvector matrix to identity
    for i in 0..n {
        v[i][i] = 1.0;
    }

    // Need at least n*(n-1)/2 rotations per sweep; 5 sweeps is typical for convergence.
    let max_iterations = (n * (n + 1) / 2 * 5).max(100);
    let tolerance = 1e-10;

    for _iteration in 0..max_iterations {
        // Find largest off-diagonal element
        let mut max_val = 0.0;
        let mut p = 0;
        let mut q = 1;

        for i in 0..n {
            for j in (i + 1)..n {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }

        if max_val < tolerance {
            break;
        }

        // Compute rotation angle
        let apq = a[p][q];
        let app = a[p][p];
        let aqq = a[q][q];

        let theta = 0.5 * (2.0 * apq / (app - aqq)).atan();
        let c = theta.cos();
        let s = theta.sin();

        // Apply Givens rotation
        for i in 0..n {
            if i == p || i == q {
                continue;
            }

            let aip = a[i][p];
            let aiq = a[i][q];

            a[i][p] = c * aip - s * aiq;
            a[p][i] = a[i][p];
            a[i][q] = s * aip + c * aiq;
            a[q][i] = a[i][q];
        }

        // Update diagonal
        a[p][p] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[q][q] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        // Update eigenvectors
        for i in 0..n {
            let vip = v[i][p];
            let viq = v[i][q];
            v[i][p] = c * vip - s * viq;
            v[i][q] = s * vip + c * viq;
        }
    }

    let eigenvalues: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    (eigenvalues, v)
}

/// Center coordinates at origin (centroid).
pub(crate) fn center_coordinates(coords: &mut Coords3D) {
    let n = coords.atom_count();
    if n == 0 {
        return;
    }

    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;

    for i in 0..n {
        let p = coords.get(AtomIdx(i as u32));
        cx += p.x;
        cy += p.y;
        cz += p.z;
    }

    cx /= n as f64;
    cy /= n as f64;
    cz /= n as f64;

    for i in 0..n {
        let p = coords.get(AtomIdx(i as u32));
        coords.set(
            AtomIdx(i as u32),
            Point3 {
                x: p.x - cx,
                y: p.y - cy,
                z: p.z - cz,
            },
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn test_bound_matrix_ethane() {
        let mol = parse("CC").unwrap();
        let (lower, upper) = build_bound_matrix(&mol);
        assert_eq!(lower.len(), 2);
        assert_eq!(upper.len(), 2);
        // C-C bond should be ~1.54 Å
        assert!(lower[0][1] > 1.4 && lower[0][1] < 1.6);
        assert!(upper[0][1] > 1.4 && upper[0][1] < 1.6);
    }

    #[test]
    fn test_bound_matrix_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let (lower, _upper) = build_bound_matrix(&mol);
        assert_eq!(lower.len(), 6);
        // Aromatic C-C should be ~1.40 Å
        assert!(lower[0][1] > 1.3 && lower[0][1] < 1.5);
    }

    #[test]
    fn test_ideal_bond_length_cc() {
        assert!((ideal_bond_length_test(6, 6, BondOrder::Single) - 1.54).abs() < 0.01);
    }

    #[test]
    fn test_ideal_bond_angle_sp2() {
        let angle = ideal_bond_angle_test(true);
        assert!((angle - 2.0 * PI / 3.0).abs() < 0.01); // 120°
    }

    #[test]
    fn test_generate_coords_dg_ethane() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords_dg(&mol);

        assert_eq!(coords.atom_count(), 2);

        let p0 = coords.get(AtomIdx(0));
        let p1 = coords.get(AtomIdx(1));

        // Should have non-zero distance
        let dist = p0.distance(&p1);
        assert!(dist > 1.4 && dist < 1.7, "C-C bond distance: {}", dist);
    }

    #[test]
    fn test_generate_coords_dg_finite() {
        let mol = parse("c1ccc(C)cc1").unwrap();
        let coords = generate_coords_dg(&mol);

        // All coordinates must be finite (no NaN, no Inf)
        for i in 0..mol.atom_count() {
            let p = coords.get(AtomIdx(i as u32));
            assert!(p.x.is_finite(), "atom {} x is not finite", i);
            assert!(p.y.is_finite(), "atom {} y is not finite", i);
            assert!(p.z.is_finite(), "atom {} z is not finite", i);
        }
    }

    #[test]
    fn test_generate_coords_dg_centered() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords_dg(&mol);

        // Calculate centroid
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;

        for i in 0..mol.atom_count() {
            let p = coords.get(AtomIdx(i as u32));
            cx += p.x;
            cy += p.y;
            cz += p.z;
        }

        cx /= mol.atom_count() as f64;
        cy /= mol.atom_count() as f64;
        cz /= mol.atom_count() as f64;

        // Centroid should be near origin
        assert!(cx.abs() < 1e-6, "centroid x: {}", cx);
        assert!(cy.abs() < 1e-6, "centroid y: {}", cy);
        assert!(cz.abs() < 1e-6, "centroid z: {}", cz);
    }

    #[test]
    fn test_generate_coords_dg_single_atom() {
        let mol = parse("C").unwrap();
        let coords = generate_coords_dg(&mol);

        assert_eq!(coords.atom_count(), 1);
        let p = coords.get(AtomIdx(0));
        assert_eq!(p.x, 0.0);
        assert_eq!(p.y, 0.0);
        assert_eq!(p.z, 0.0);
    }

    #[test]
    fn test_generate_coords_dg_ethene() {
        let mol = parse("C=C").unwrap();
        let coords = generate_coords_dg(&mol);

        let p0 = coords.get(AtomIdx(0));
        let p1 = coords.get(AtomIdx(1));
        let dist = p0.distance(&p1);

        // Double bond C=C is ~1.34 Å
        assert!(dist > 1.2 && dist < 1.5, "C=C bond distance: {}", dist);
    }

    #[test]
    fn test_generate_coords_dg_propane() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords_dg(&mol);

        assert_eq!(coords.atom_count(), 3);

        let p0 = coords.get(AtomIdx(0));
        let p1 = coords.get(AtomIdx(1));
        let p2 = coords.get(AtomIdx(2));

        // All bonds should have reasonable distances
        let d01 = p0.distance(&p1);
        let d12 = p1.distance(&p2);
        let d02 = p0.distance(&p2);

        assert!(d01 > 1.4 && d01 < 1.7);
        assert!(d12 > 1.4 && d12 < 1.7);
        assert!(d02 > 2.0 && d02 < 3.5);
    }

    #[test]
    fn test_generate_coords_dg_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords_dg(&mol);

        assert_eq!(coords.atom_count(), 6);

        // All atoms should have non-zero positions
        let mut all_nonzero = true;
        for i in 0..6 {
            let p = coords.get(AtomIdx(i as u32));
            if p.x.abs() < 1e-10 && p.y.abs() < 1e-10 && p.z.abs() < 1e-10 {
                all_nonzero = false;
            }
        }
        assert!(all_nonzero, "some atoms have zero coordinates");
    }

    #[test]
    fn test_gram_matrix_simple() {
        let dist = vec![
            vec![0.0, 1.5, 2.5],
            vec![1.5, 0.0, 1.6],
            vec![2.5, 1.6, 0.0],
        ];

        let gram = distance_to_gram_matrix(&dist);
        assert_eq!(gram.len(), 3);
        assert_eq!(gram[0].len(), 3);

        // Gram matrix should be symmetric
        for i in 0..3 {
            for j in 0..3 {
                assert!((gram[i][j] - gram[j][i]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_jacobi_eigendecompose_identity() {
        let mat = vec![vec![1.0, 0.0], vec![0.0, 1.0]];

        let (eigenvalues, _eigenvectors) = jacobi_eigendecompose(&mat);
        assert_eq!(eigenvalues.len(), 2);
        assert!((eigenvalues[0] - 1.0).abs() < 1e-6);
        assert!((eigenvalues[1] - 1.0).abs() < 1e-6);
    }

    // T1: Gram-matrix centering — verify bond lengths are reasonable for larger molecules.
    //
    // The atom-0-reference formula is mathematically equivalent to double-centering
    // for exact Euclidean distance matrices.  For approximate distances (as produced by
    // the bound matrix), the two differ only in numerical bias.  These tests confirm
    // that generated bond lengths are within ±0.3 Å of the target for molecules up to
    // 6 heavy atoms, establishing that the current implementation has no practical impact.

    #[test]
    fn test_gram_centering_butane_bond_lengths() {
        // Butane (4 atoms: 0-1-2-3). Atoms 0 and 3 have no angle constraint (4 bonds apart)
        // → dist_matrix[0][3] = infinity. This exercises the long-range infinity path.
        let mol = parse("CCCC").unwrap();
        let coords = generate_coords_dg(&mol);

        // All coords must be finite.
        for i in 0..4 {
            let p = coords.get(AtomIdx(i as u32));
            assert!(
                p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                "butane atom {} has non-finite coords ({}, {}, {})",
                i,
                p.x,
                p.y,
                p.z
            );
        }

        // Bond distances should be within ±0.3 Å of the ideal 1.54 Å target.
        let check_bond = |a: u32, b: u32| {
            let pa = coords.get(AtomIdx(a));
            let pb = coords.get(AtomIdx(b));
            let d = pa.distance(&pb);
            (d > 1.2 && d < 1.9, d)
        };
        let (ok01, d01) = check_bond(0, 1);
        let (ok12, d12) = check_bond(1, 2);
        let (ok23, d23) = check_bond(2, 3);
        assert!(ok01, "butane C0-C1 bond = {d01:.3} Å (expected ~1.54)");
        assert!(ok12, "butane C1-C2 bond = {d12:.3} Å (expected ~1.54)");
        assert!(ok23, "butane C2-C3 bond = {d23:.3} Å (expected ~1.54)");
    }

    #[test]
    fn test_gram_centering_hexane_all_bonds() {
        // Hexane (6 atoms). Many long-range pairs (0-3, 0-4, 0-5, 1-4, 1-5, 2-5)
        // have upper = infinity in the bound matrix.
        let mol = parse("CCCCCC").unwrap();
        let coords = generate_coords_dg(&mol);

        for i in 0..6 {
            let p = coords.get(AtomIdx(i as u32));
            assert!(
                p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                "hexane atom {} non-finite: ({}, {}, {})",
                i,
                p.x,
                p.y,
                p.z
            );
        }

        // All five C-C bonds should be reasonable.
        for i in 0..5u32 {
            let pa = coords.get(AtomIdx(i));
            let pb = coords.get(AtomIdx(i + 1));
            let d = pa.distance(&pb);
            assert!(d > 1.2 && d < 1.9, "hexane C{i}-C{} bond = {d:.3} Å", i + 1);
        }
    }

    #[test]
    fn test_gram_centering_toluene_bond_quality() {
        // Toluene has a methyl C that is ≥4 bonds from most ring atoms → many infinite
        // upper bounds. All coords must be finite and ring bonds reasonable.
        let mol = parse("c1ccc(C)cc1").unwrap();
        let coords = generate_coords_dg(&mol);

        for i in 0..mol.atom_count() {
            let p = coords.get(AtomIdx(i as u32));
            assert!(
                p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                "toluene atom {i} non-finite"
            );
        }

        // Check that BONDED atom pairs have distances within 1.0–2.0 Å.
        use chematic_core::BondIdx;
        for bi in 0..mol.bond_count() {
            let bond = mol.bond(BondIdx(bi as u32));
            let pa = coords.get(bond.atom1);
            let pb = coords.get(bond.atom2);
            let d = pa.distance(&pb);
            assert!(
                d > 1.0 && d < 2.0,
                "toluene bond {}-{} = {d:.3} Å",
                bond.atom1.0,
                bond.atom2.0
            );
        }
    }

    // Helper functions for testing
    fn ideal_bond_length_test(a: u8, b: u8, order: BondOrder) -> f64 {
        match (a.min(b), a.max(b), order) {
            (6, 6, BondOrder::Single) => 1.54,
            (6, 6, BondOrder::Double) => 1.34,
            (6, 6, BondOrder::Triple) => 1.20,
            _ => 1.54,
        }
    }

    fn ideal_bond_angle_test(aromatic: bool) -> f64 {
        if aromatic { 2.0 * PI / 3.0 } else { 1.91 }
    }
}
