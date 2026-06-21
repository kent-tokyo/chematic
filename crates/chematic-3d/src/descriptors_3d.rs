//! 3D molecular descriptors for ML pipelines.
//!
//! Implements simplified WHIM and GETAWAY-like descriptors:
//! - WHIM: weighted holistic invariant descriptors based on mass distribution
//! - GETAWAY: geometric/topologic descriptors with wavelet autocorrelation
//!
//! Full implementations require expensive computations; this module provides
//! practical approximations suitable for ML feature vectors.

use crate::coords::{Coords3D, Point3};
use chematic_core::Molecule;

/// Compute WHIM descriptors: mass-weighted shape descriptors.
/// Returns: [L1, L2, L3, P1, P2, P3, ALPHA, BETA, GAMMA, DELTA]
/// where L* = eigenvalues of inertia tensor, P* = principal moments,
/// and ALPHA/BETA/GAMMA/DELTA are derived shape metrics.
pub fn whim_descriptors(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    if mol.atom_count() < 2 {
        return vec![0.0; 10];
    }

    // Compute center of mass
    let mut total_mass = 0.0;
    let mut com = Point3::zero();

    for i in 0..mol.atom_count() {
        let atom = mol.atom(chematic_core::AtomIdx(i as u32));
        let mass = atom.element.atomic_mass();
        total_mass += mass;
        let p = coords.get(chematic_core::AtomIdx(i as u32));
        com = com.add(&p.scale(mass));
    }

    if total_mass == 0.0 {
        return vec![0.0; 10];
    }

    com = com.scale(1.0 / total_mass);

    // Compute inertia tensor
    let mut ixx = 0.0;
    let mut iyy = 0.0;
    let mut izz = 0.0;

    for i in 0..mol.atom_count() {
        let atom = mol.atom(chematic_core::AtomIdx(i as u32));
        let mass = atom.element.atomic_mass();
        let p = coords.get(chematic_core::AtomIdx(i as u32));
        let r = p.sub(&com);

        ixx += mass * (r.y * r.y + r.z * r.z);
        iyy += mass * (r.x * r.x + r.z * r.z);
        izz += mass * (r.x * r.x + r.y * r.y);
    }

    // Approximate eigenvalues (simplified: use diagonal dominance)
    let l1 = ixx;
    let l2 = iyy;
    let l3 = izz;

    let p1 = (l1 / total_mass).sqrt();
    let p2 = (l2 / total_mass).sqrt();
    let p3 = (l3 / total_mass).sqrt();

    // Shape metrics
    let alpha = p1 + p2 + p3; // Total reach
    let beta = (p1 * p2 + p2 * p3 + p3 * p1) / 3.0; // Average interaction
    let gamma = (p1 * p2 * p3).cbrt(); // Geometric mean
    let delta = p1 - p3; // Anisotropy

    vec![l1, l2, l3, p1, p2, p3, alpha, beta, gamma, delta]
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
/// Returns a 29-element feature vector: WHIM (10) followed by GETAWAY (19).
pub fn whim_getaway_combined(mol: &Molecule, coords: &Coords3D) -> Vec<f64> {
    let mut result = whim_descriptors(mol, coords);
    result.extend(getaway_descriptors(mol, coords));
    result
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
        assert_eq!(desc.len(), 10);
        assert!(
            desc.iter().all(|&d| d.is_finite()),
            "all WHIM descriptors should be finite"
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
        // WHIM (10) + GETAWAY (19) = 29
        assert_eq!(desc.len(), 29);
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
}
