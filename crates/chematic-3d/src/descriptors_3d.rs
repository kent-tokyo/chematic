//! 3D molecular descriptors for ML pipelines.
//!
//! Implements WHIM, GETAWAY, RDF, and AutoCorr 3D descriptors:
//! - WHIM (22 dims): Weighted Holistic Invariant Molecular descriptors, two weight schemes
//! - GETAWAY (19 dims): leverage-based geometric/topologic descriptors
//! - RDF (20 dims): Radial Distribution Function descriptors
//! - AutoCorr3D (8 dims): Moreau-Broto self-correlation at Euclidean distance lags

use crate::coords::Coords3D;
use crate::shape_descriptors::jacobi3;
use chematic_core::Molecule;

// ---------------------------------------------------------------------------
// WHIM internals
// ---------------------------------------------------------------------------

/// Compute 11 WHIM descriptors for one weight scheme.
///
/// Returns `[λ₁, λ₂, λ₃, ν₁, ν₂, ν₃, T, A, V, K, D]` where:
/// - λₖ: eigenvalues of the weighted covariance matrix (descending)
/// - νₖ: skewness of projection onto k-th principal axis
/// - T = Σλ, A = Σᵢ<ⱼ λᵢλⱼ, V = λ₁λ₂λ₃, K = Σᵢ<ⱼ(λᵢ−λⱼ)²/T², D = T/3
fn whim_11(xs: &[[f64; 3]], weights: &[f64]) -> [f64; 11] {
    let n = xs.len();
    let total_w: f64 = weights.iter().sum();
    if total_w < 1e-10 {
        return [0.0; 11];
    }

    // Weighted centroid
    let mut com = [0.0f64; 3];
    for i in 0..n {
        for d in 0..3 {
            com[d] += weights[i] * xs[i][d];
        }
    }
    for d in 0..3 {
        com[d] /= total_w;
    }

    // Weighted covariance matrix (normalised by total weight)
    let mut cov = [[0.0f64; 3]; 3];
    for i in 0..n {
        let dx = [xs[i][0] - com[0], xs[i][1] - com[1], xs[i][2] - com[2]];
        for a in 0..3 {
            for b in 0..3 {
                cov[a][b] += weights[i] * dx[a] * dx[b];
            }
        }
    }
    for a in 0..3 {
        for b in 0..3 {
            cov[a][b] /= total_w;
        }
    }

    // Eigendecompose — jacobi3 returns ascending eigenvalues
    // evecs[row][col] = row-th component of the col-th eigenvector
    let (evals_asc, evecs) = jacobi3(cov);

    // Reorder to descending: λ₁ ≥ λ₂ ≥ λ₃
    // ascending index 2 → largest; 1 → middle; 0 → smallest
    let lam = [
        evals_asc[2].max(0.0),
        evals_asc[1].max(0.0),
        evals_asc[0].max(0.0),
    ];
    // Skewness-based symmetry coefficients νₖ
    // Both numerator (weighted 3rd moment) and denominator (λ^(3/2)) must use
    // the same weighting scheme; dividing by total_w matches the covariance.
    let mut nu = [0.0f64; 3];
    for k in 0..3 {
        let lambda = lam[k];
        if lambda < 1e-10 {
            continue;
        }
        let col = 2 - k; // ascending index 0→smallest, 2→largest; descend by inverting
        let mut sum_cube = 0.0f64;
        for i in 0..n {
            let proj = (xs[i][0] - com[0]) * evecs[0][col]
                + (xs[i][1] - com[1]) * evecs[1][col]
                + (xs[i][2] - com[2]) * evecs[2][col];
            sum_cube += weights[i] * proj.powi(3);
        }
        nu[k] = (sum_cube / total_w) / lambda.powf(1.5);
    }

    // Global 3D indices
    let (l1, l2, l3) = (lam[0], lam[1], lam[2]);
    let t = l1 + l2 + l3;
    let a = l1 * l2 + l1 * l3 + l2 * l3;
    let v = l1 * l2 * l3;
    let k = if t > 1e-10 {
        ((l1 - l2).powi(2) + (l1 - l3).powi(2) + (l2 - l3).powi(2)) / t.powi(2)
    } else {
        0.0
    };
    let d = t / 3.0;

    [l1, l2, l3, nu[0], nu[1], nu[2], t, a, v, k, d]
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Compute WHIM (Weighted Holistic Invariant Molecular) 3D descriptors.
///
/// Returns a 22-element vector: 11 descriptors for unit weights followed by
/// 11 descriptors for atomic-mass weights.
///
/// Each 11-element block is `[λ₁, λ₂, λ₃, ν₁, ν₂, ν₃, T, A, V, K, D]`:
/// - `λ₁ ≥ λ₂ ≥ λ₃`: eigenvalues of the weighted covariance matrix (Å²)
/// - `ν₁, ν₂, ν₃`: skewness-based symmetry coefficients for each principal axis
/// - `T = λ₁+λ₂+λ₃`, `A = Σᵢ<ⱼ λᵢλⱼ`, `V = λ₁λ₂λ₃`: global extent indices
/// - `K = Σᵢ<ⱼ(λᵢ−λⱼ)²/T²`: anisotropy index (0 = sphere, 1 = rod/disk)
/// - `D = T/3`: average spread
pub fn whim_descriptors(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    let n = mol.atom_count();
    if n < 2 {
        return vec![0.0; 22];
    }

    let xs: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let p = coords.get(chematic_core::AtomIdx(i as u32));
            [p.x, p.y, p.z]
        })
        .collect();

    // W1: unit weights
    let w1 = vec![1.0f64; n];
    // W2: atomic-mass weights
    let w2: Vec<f64> = (0..n)
        .map(|i| {
            mol.atom(chematic_core::AtomIdx(i as u32))
                .element
                .atomic_mass()
        })
        .collect();

    let mut result = Vec::with_capacity(22);
    result.extend_from_slice(&whim_11(&xs, &w1));
    result.extend_from_slice(&whim_11(&xs, &w2));
    result
}

/// GETAWAY (GEometry, Topology and Atom-Weights AssemblY) descriptors.
///
/// Returns a 19-element vector:
/// - `H[1..8]` — leverage autocorrelation at topological lags 1–8:
///   `H[k] = Σ_{d(i,j)=k} √(h_i · h_j)` (sum over heavy-atom pairs)
/// - `R[1..8]` — normalised: `H[k] / W_k` (W_k = number of pairs at lag k)
/// - `Hmax`, `Hmean`, `Htot` — leverage statistics
///
/// The per-atom *leverage* h_i is the diagonal of the hat matrix
/// `H = X(X^T X)^{-1} X^T`, where X is the centred 3D coordinate matrix
/// of heavy atoms.  Leverage measures each atom's geometric influence.
pub fn getaway_descriptors(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    // Collect heavy-atom indices (exclude H).
    let heavy: Vec<usize> = (0..mol.atom_count())
        .filter(|&i| {
            mol.atom(chematic_core::AtomIdx(i as u32))
                .element
                .atomic_number()
                != 1
        })
        .collect();
    let hn = heavy.len();
    if hn < 2 {
        return vec![0.0; 19];
    }

    // ── Step 1: centred coordinate matrix X (hn×3) ──────────────────────────
    let mut cx = 0.0f64;
    let mut cy = 0.0f64;
    let mut cz = 0.0f64;
    for &h in &heavy {
        let p = coords.get(chematic_core::AtomIdx(h as u32));
        cx += p.x;
        cy += p.y;
        cz += p.z;
    }
    let hn_f = hn as f64;
    cx /= hn_f;
    cy /= hn_f;
    cz /= hn_f;

    let xs: Vec<[f64; 3]> = heavy
        .iter()
        .map(|&h| {
            let p = coords.get(chematic_core::AtomIdx(h as u32));
            [p.x - cx, p.y - cy, p.z - cz]
        })
        .collect();

    // ── Step 2: X^T X (3×3) ──────────────────────────────────────────────────
    let mut xtx = [[0.0f64; 3]; 3];
    for row in &xs {
        for a in 0..3 {
            for b in 0..3 {
                xtx[a][b] += row[a] * row[b];
            }
        }
    }

    // ── Step 3: (X^T X)^{-1} — analytical 3×3 inverse ───────────────────────
    let inv = mat3_inv(&xtx);

    // ── Step 4: leverage h_i = X_i^T (X^T X)^{-1} X_i ───────────────────────
    let leverage: Vec<f64> = xs
        .iter()
        .map(|xi| {
            let mut h = 0.0f64;
            for a in 0..3 {
                for b in 0..3 {
                    h += xi[a] * inv[a][b] * xi[b];
                }
            }
            h.max(0.0) // numerical safety against tiny negatives
        })
        .collect();

    // ── Step 5: topological distance matrix (BFS) ────────────────────────────
    let topo = heavy_topo_dist_local(mol, &heavy);

    // ── Step 6: GETAWAY H and R descriptors ─────────────────────────────────
    const MAX_LAG: usize = 8;
    let mut h_lags = vec![0.0f64; MAX_LAG];
    let mut w_lags = [0usize; MAX_LAG];

    for i in 0..hn {
        for j in (i + 1)..hn {
            let d = topo[i][j];
            if d == 0 || d as usize > MAX_LAG {
                continue;
            }
            let k = (d - 1) as usize;
            h_lags[k] += (leverage[i] * leverage[j]).sqrt();
            w_lags[k] += 1;
        }
    }

    let r_lags: Vec<f64> = h_lags
        .iter()
        .zip(w_lags.iter())
        .map(|(&h, &w)| if w == 0 { 0.0 } else { h / w as f64 })
        .collect();

    // ── Step 7: leverage statistics ─────────────────────────────────────────
    let hmax = leverage.iter().cloned().fold(0.0f64, f64::max);
    let hmean = leverage.iter().sum::<f64>() / hn_f;
    let htot = leverage.iter().sum::<f64>();

    let mut out = h_lags;
    out.extend(r_lags);
    out.extend([hmax, hmean, htot]);
    out // 8 + 8 + 3 = 19 elements
}

/// Inverse of a 3×3 matrix.  Returns the identity when the determinant is
/// near-zero (degenerate / planar coordinate set).
fn mat3_inv(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() < 1e-10 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    let d = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * d,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * d,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * d,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * d,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * d,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * d,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * d,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * d,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * d,
        ],
    ]
}

/// BFS topological distance matrix for a subset of atoms.
/// Avoids a cross-crate dependency on chematic-chem.
fn heavy_topo_dist_local(mol: &Molecule, heavy: &[usize]) -> Vec<Vec<u32>> {
    use std::collections::{HashSet, VecDeque};
    let heavy_set: HashSet<usize> = heavy.iter().copied().collect();
    let hn = heavy.len();
    let mut matrix = vec![vec![u32::MAX; hn]; hn];
    for (p, &start) in heavy.iter().enumerate() {
        matrix[p][p] = 0;
        let n_atoms = mol.atom_count();
        let mut dist = vec![usize::MAX; n_atoms];
        dist[start] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(start);
        while let Some(cur) = queue.pop_front() {
            let d = dist[cur];
            for (nb, _) in mol.neighbors(chematic_core::AtomIdx(cur as u32)) {
                let ni = nb.0 as usize;
                if heavy_set.contains(&ni) && dist[ni] == usize::MAX {
                    dist[ni] = d + 1;
                    queue.push_back(ni);
                }
            }
        }
        for (q, &h) in heavy.iter().enumerate() {
            let d = dist[h];
            if d != usize::MAX {
                matrix[p][q] = d as u32;
            }
        }
    }
    matrix
}

/// Combined WHIM + GETAWAY descriptor vector for ML.
///
/// Returns a 41-element feature vector: WHIM (22) followed by GETAWAY (19).
pub fn whim_getaway_combined(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    let mut result = whim_descriptors(mol, coords);
    result.extend(getaway_descriptors(mol, coords));
    result
}

/// RDF (Radial Distribution Function) 3D descriptors.
///
/// Computes 20 mass-weighted RDF values at shells r_k = 0.5, 1.0, ..., 10.0 Å:
///
/// `g(r_k) = Σᵢ Σⱼ>ᵢ mᵢmⱼ · exp(−β·(r_k − rᵢⱼ)²)`
///
/// with β = 100 and mᵢ = atomic mass.  Returns a 20-element vector.
pub fn rdf_descriptors(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    const N_SHELLS: usize = 20;
    const BETA: f64 = 100.0;

    let n = mol.atom_count();
    if n < 2 {
        return vec![0.0; N_SHELLS];
    }

    let mut g = vec![0.0f64; N_SHELLS];

    for i in 0..n {
        let idx_i = chematic_core::AtomIdx(i as u32);
        let mi = mol.atom(idx_i).element.atomic_mass();
        let pi = coords.get(idx_i);

        for j in (i + 1)..n {
            let idx_j = chematic_core::AtomIdx(j as u32);
            let mj = mol.atom(idx_j).element.atomic_mass();
            let pj = coords.get(idx_j);

            let rij = pi.distance(&pj);
            let weight = mi * mj;

            for k in 0..N_SHELLS {
                let r_k = (k as f64 + 1.0) * 0.5; // 0.5, 1.0, ..., 10.0
                let diff = r_k - rij;
                g[k] += weight * (-BETA * diff * diff).exp();
            }
        }
    }

    g
}

/// AutoCorr3D: Moreau-Broto Self-Correlation (Euclidean Distance).
///
/// Compute self-correlation for 3D coordinates using binned Euclidean distances.
/// Each lag corresponds to a distance bin (k * 1Å):
/// - lag 1: 0-1 Å
/// - lag 2: 1-2 Å
/// - lag 3: 2-3 Å
/// - ... lag 8: 7-8 Å
///
/// For each lag, sum over all atom pairs (i,j) with distance in bin of v(i) * v(j),
/// where v(i) is the atomic mass (simplified feature for 3D self-correlation).
pub fn autocorr_3d(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    if mol.atom_count() < 2 {
        return vec![0.0; 8];
    }

    let n = mol.atom_count();
    let mut result = vec![0.0; 8];

    for lag in 1..=8 {
        let lower = (lag - 1) as f64;
        let upper = lag as f64;
        let mut sum = 0.0;

        for i in 0..n {
            let idx_i = chematic_core::AtomIdx(i as u32);
            let atom_i = mol.atom(idx_i);
            let mass_i = atom_i.element.atomic_mass();
            let p_i = coords.get(idx_i);

            for j in (i + 1)..n {
                let idx_j = chematic_core::AtomIdx(j as u32);
                let atom_j = mol.atom(idx_j);
                let mass_j = atom_j.element.atomic_mass();
                let p_j = coords.get(idx_j);

                let dist = p_i.distance(&p_j);
                if dist >= lower && dist < upper {
                    sum += mass_i * mass_j;
                }
            }
        }
        result[lag - 1] = sum;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dg::generate_coords;
    use chematic_smiles::parse;

    #[test]
    fn test_whim_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let desc = whim_descriptors(&mol, &coords);
        assert_eq!(desc.len(), 22, "WHIM should be 22 elements (11×W1 + 11×W2)");
        assert!(
            desc.iter().all(|&d| d.is_finite()),
            "all WHIM descriptors should be finite"
        );
        // λ₁ (W1 block[0]) should be non-negative
        assert!(desc[0] >= 0.0, "λ₁(W1) should be ≥ 0: {}", desc[0]);
    }

    #[test]
    fn test_whim_heteroatom_w1_w2_differ() {
        // Ethanol has C, C, O → mass-weight scheme gives more weight to O
        // → centroid shifts → W1 ≠ W2
        let mol = parse("CCO").unwrap();
        let coords = generate_coords(&mol);
        let desc = whim_descriptors(&mol, &coords);
        assert_eq!(desc.len(), 22);
        let w1 = &desc[..11];
        let w2 = &desc[11..];
        assert!(
            w1.iter().zip(w2.iter()).any(|(a, b)| (a - b).abs() > 1e-6),
            "W1 and W2 blocks should differ for a heteroatom molecule"
        );
    }

    #[test]
    fn test_getaway_propane() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords(&mol);
        let desc = getaway_descriptors(&mol, &coords);
        // HATs-based GETAWAY: 8(H) + 8(R) + 3(Hmax/Hmean/Htot) = 19
        assert_eq!(desc.len(), 19);
        assert!(
            desc.iter().all(|&d| d.is_finite()),
            "all GETAWAY descriptors should be finite"
        );
    }

    #[test]
    fn test_getaway_single_atom() {
        // Molecule with fewer than 2 heavy atoms → zero vector of length 19
        let mol = parse("[H]").unwrap();
        let coords = generate_coords(&mol);
        let desc = getaway_descriptors(&mol, &coords);
        assert_eq!(desc.len(), 19);
        assert!(desc.iter().all(|&d| d == 0.0));
    }

    #[test]
    fn test_getaway_leverage_nonnegative() {
        // Htot (index 18) should be non-negative (sum of leverages)
        let mol = parse("CCO").unwrap();
        let coords = generate_coords(&mol);
        let desc = getaway_descriptors(&mol, &coords);
        assert!(desc[18] >= 0.0, "Htot must be non-negative: {}", desc[18]);
        assert!(desc[17] >= 0.0, "Hmean must be non-negative");
        assert!(desc[16] >= 0.0, "Hmax must be non-negative");
    }

    #[test]
    fn test_combined_aspirin() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let coords = generate_coords(&mol);
        let desc = whim_getaway_combined(&mol, &coords);
        // WHIM (22) + GETAWAY (19) = 41
        assert_eq!(desc.len(), 41);
        assert!(
            desc.iter().all(|&d| d.is_finite()),
            "all combined descriptors should be finite"
        );
    }

    #[test]
    fn test_autocorr_3d_single_atom() {
        let mol = parse("C").unwrap();
        let coords = generate_coords(&mol);
        let ac = autocorr_3d(&mol, &coords);
        assert_eq!(ac.len(), 8);
        // Single atom: no pairs → all zeros
        for val in ac {
            assert!((val - 0.0).abs() < 1e-9);
        }
    }

    #[test]
    fn test_autocorr_3d_ethane() {
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);
        let ac = autocorr_3d(&mol, &coords);
        assert_eq!(ac.len(), 8);
        // Ethane: C-C distance ≈ 1.54 Å → lag 2 (1-2 Å)
        // Mass of C ≈ 12.0, so product ≈ 144
        assert!(ac[0] < 1.0, "lag 1 (0-1Å) should be minimal: {}", ac[0]);
        assert!(ac[1] > 100.0, "lag 2 (1-2Å) should be ~144: {}", ac[1]);
    }

    #[test]
    fn test_autocorr_3d_propane() {
        let mol = parse("CCC").unwrap();
        let coords = generate_coords(&mol);
        let ac = autocorr_3d(&mol, &coords);
        assert_eq!(ac.len(), 8);
        // Should have non-zero values in appropriate distance bins
        assert!(
            ac.iter().any(|&x| x > 0.0),
            "should have non-zero autocorr values"
        );
        assert!(
            ac.iter().all(|&x| x.is_finite()),
            "all values should be finite"
        );
    }

    #[test]
    fn test_autocorr_3d_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let ac = autocorr_3d(&mol, &coords);
        assert_eq!(ac.len(), 8);
        // Benzene: ring structure with various distances
        assert!(
            ac.iter().any(|&x| x > 0.0),
            "benzene should have non-zero autocorr"
        );
        assert!(
            ac.iter().all(|&x| x.is_finite()),
            "all values should be finite"
        );
    }

    // -----------------------------------------------------------------------
    // RDF descriptor tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_rdf_length_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let coords = generate_coords(&mol);
        let rdf = rdf_descriptors(&mol, &coords);
        assert_eq!(rdf.len(), 20, "RDF should have 20 shells");
        assert!(rdf.iter().all(|&v| v.is_finite()), "all RDF values finite");
        assert!(
            rdf.iter().any(|&v| v > 0.0),
            "RDF should be non-zero for benzene"
        );
    }

    #[test]
    fn test_rdf_single_atom_zeros() {
        let mol = parse("C").unwrap();
        let coords = generate_coords(&mol);
        let rdf = rdf_descriptors(&mol, &coords);
        assert_eq!(rdf.len(), 20);
        assert!(rdf.iter().all(|&v| v == 0.0), "single atom → all zero RDF");
    }

    #[test]
    fn test_rdf_ethane_nonzero_near_bond() {
        // C-C distance ≈ 1.54 Å → peak near shell 3 (r=1.5) or 4 (r=2.0)
        let mol = parse("CC").unwrap();
        let coords = generate_coords(&mol);
        let rdf = rdf_descriptors(&mol, &coords);
        assert_eq!(rdf.len(), 20);
        // At least one of shells 2-4 (r=1.0..2.0 Å) should be non-negligible
        assert!(
            rdf[2] + rdf[3] > 1e-5,
            "RDF near C-C bond length should be non-zero"
        );
    }

    // -----------------------------------------------------------------------
    // WHIM property tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_whim_single_atom_zeros() {
        let mol = parse("C").unwrap();
        let coords = generate_coords(&mol);
        let desc = whim_descriptors(&mol, &coords);
        assert_eq!(desc.len(), 22);
        // Single atom: covariance is zero → all descriptors zero
        assert!(
            desc.iter().all(|&v| v == 0.0),
            "single atom WHIM should be all zeros"
        );
    }

    #[test]
    fn test_whim_eigenvalues_nonneg() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let coords = generate_coords(&mol);
        let desc = whim_descriptors(&mol, &coords);
        // λ₁ ≥ λ₂ ≥ λ₃ ≥ 0 in both W1 and W2 blocks
        for block_start in [0, 11] {
            let (l1, l2, l3) = (
                desc[block_start],
                desc[block_start + 1],
                desc[block_start + 2],
            );
            assert!(l1 >= l2 - 1e-9, "W block {block_start}: λ₁ ≥ λ₂");
            assert!(l2 >= l3 - 1e-9, "W block {block_start}: λ₂ ≥ λ₃");
            assert!(l3 >= -1e-9, "W block {block_start}: λ₃ ≥ 0");
        }
    }
}
