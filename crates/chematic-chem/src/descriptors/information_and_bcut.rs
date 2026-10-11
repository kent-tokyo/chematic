// ---------------------------------------------------------------------------
// Information Content descriptors (Mordred IC family)
// ---------------------------------------------------------------------------

/// Shannon-entropy-based information content descriptors.
///
/// Based on the vertex equivalence class partition of the molecular graph:
/// atoms are grouped by `(element, heavy_degree)` as their Morgan-0 invariant.
///
/// - **IC**  = Shannon entropy of the partition (bits)
/// - **TIC** = N × IC  (total information content)
/// - **SIC** = IC / log₂(N)  (structural information content; 0 for fully symmetric)
/// - **BIC** = IC × (N-1) / N  (bond information content)
/// - **CIC** = log₂(N) − IC  (complementary information content)
///
/// All values are 0 for a single-atom molecule.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct InformationContent {
    pub ic: f64,
    pub tic: f64,
    pub sic: f64,
    pub bic: f64,
    pub cic: f64,
}

/// Compute the Information Content descriptor family (IC, TIC, SIC, BIC, CIC).
pub fn information_content(mol: &Molecule) -> InformationContent {
    let heavy: Vec<_> = mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() != 1)
        .collect();
    let n = heavy.len();
    if n < 2 {
        return InformationContent::default();
    }
    // Vertex invariant: (atomic_number, heavy_degree)
    let mut counts: FxHashMap<(u8, u32), usize> = FxHashMap::default();
    for (idx, atom) in &heavy {
        let deg = mol
            .neighbors(*idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() != 1)
            .count() as u32;
        *counts
            .entry((atom.element.atomic_number(), deg))
            .or_insert(0) += 1;
    }
    let n_f = n as f64;
    let ic: f64 = counts
        .values()
        .map(|&c| {
            let p = c as f64 / n_f;
            -p * p.log2()
        })
        .sum();
    let log_n = n_f.log2();
    InformationContent {
        ic,
        tic: n_f * ic,
        sic: if log_n > 0.0 { ic / log_n } else { 0.0 },
        bic: ic * (n_f - 1.0) / n_f,
        cic: (log_n - ic).max(0.0),
    }
}

// ---------------------------------------------------------------------------
// Molecular Distance Edge (MDE) descriptors
// ---------------------------------------------------------------------------

/// Molecular Distance Edge Carbon descriptor (MDEC-*xy*).
///
/// Σ 1/√d(i,j) over all ordered pairs (i,j) where *i* is a carbon with
/// *x* carbon neighbours and *j* is a carbon with *y* carbon neighbours.
/// Here *x,y* ∈ {1,2,3,4} (primary, secondary, tertiary, quaternary C).
///
/// Returns 10 values in the order:
/// MDEC11, MDEC12, MDEC13, MDEC14, MDEC22, MDEC23, MDEC24, MDEC33, MDEC34, MDEC44.
pub fn mde_carbon(mol: &Molecule) -> [f64; 10] {
    let heavy: Vec<_> = mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() != 1)
        .map(|(idx, _)| idx)
        .collect();
    let heavy_count = heavy.len();
    if heavy_count == 0 {
        return [0.0; 10];
    }
    // Build carbon-neighbor-count for each heavy atom index.
    let mut c_neighbors: Vec<Option<u8>> = vec![None; mol.atom_count()];
    for &idx in &heavy {
        if mol.atom(idx).element.atomic_number() == 6 {
            let cn = mol
                .neighbors(idx)
                .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 6)
                .count()
                .min(4) as u8;
            c_neighbors[idx.0 as usize] = Some(cn);
        }
    }
    let dist = crate::topo_descriptors::topological_distance_matrix(mol);
    // Map original AtomIdx → row in distance matrix (heavy-atom order).
    let heavy_pos: FxHashMap<u32, usize> = heavy
        .iter()
        .enumerate()
        .map(|(p, &idx)| (idx.0, p))
        .collect();
    let mut out = [0.0f64; 10];
    let n = heavy.len();
    for pi in 0..n {
        let ai = heavy[pi];
        let xi = match c_neighbors[ai.0 as usize] {
            Some(v) if v >= 1 => v,
            _ => continue,
        };
        for pj in (pi + 1)..n {
            let aj = heavy[pj];
            let xj = match c_neighbors[aj.0 as usize] {
                Some(v) if v >= 1 => v,
                _ => continue,
            };
            let d = dist[pi][pj];
            if d == 0 || d == u32::MAX {
                continue;
            }
            let contrib = 1.0 / (d as f64).sqrt();
            let (lo, hi) = if xi <= xj { (xi, xj) } else { (xj, xi) };
            // Map (lo, hi) from (1,1)..(4,4) to index 0..9
            let slot = match (lo, hi) {
                (1, 1) => 0,
                (1, 2) => 1,
                (1, 3) => 2,
                (1, 4) => 3,
                (2, 2) => 4,
                (2, 3) => 5,
                (2, 4) => 6,
                (3, 3) => 7,
                (3, 4) => 8,
                (4, 4) => 9,
                _ => continue,
            };
            out[slot] += contrib;
        }
    }
    let _ = heavy_pos; // used indirectly via dist which was built from heavy-atom ordering
    out
}

// ---------------------------------------------------------------------------
// BCUT2D — Burden Connectivity Matrix eigenvalue descriptors
// ---------------------------------------------------------------------------

/// BCUT2D descriptors: largest and smallest eigenvalues of the Burden
/// connectivity matrix weighted by four atomic properties.
///
/// The Burden matrix B (n×n, n = heavy atoms) is defined as:
/// - `B[i,i]` = per-atom property value
/// - `B[i,j]` = √(Z_i × Z_j) × bond_weight / 8  for bonded pairs
/// - `B[i,j]` = 0.001  for non-bonded pairs
///
/// Four atomic properties are used, each giving a HI/LO eigenvalue pair:
/// - **CHG** formal charge per atom
/// - **LOGP** Crippen LogP atomic contribution
/// - **MR** molar refractivity atomic contribution
/// - **MW** atomic mass
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Bcut2D {
    pub chghi: f64,
    pub chglo: f64,
    pub logphi: f64,
    pub logplo: f64,
    pub mrhi: f64,
    pub mrlo: f64,
    pub mwhi: f64,
    pub mwlo: f64,
}

/// Compute BCUT2D descriptors (8 values: 4 properties × HI/LO eigenvalue).
pub fn bcut2d(mol: &Molecule) -> Bcut2D {
    let heavy: Vec<AtomIdx> = mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() != 1)
        .map(|(idx, _)| idx)
        .collect();
    let n = heavy.len();
    if n == 0 {
        return Bcut2D::default();
    }
    // Per-atom properties.
    let logp_v = logp_crippen_per_atom(mol);
    let mr_v = mr_per_atom(mol);
    let charge_v: Vec<f64> = mol.atoms().map(|(_, a)| a.charge as f64).collect();
    let mw_v: Vec<f64> = mol.atoms().map(|(_, a)| a.element.atomic_mass()).collect();
    let props: [Vec<f64>; 4] = [
        heavy.iter().map(|&i| charge_v[i.0 as usize]).collect(),
        heavy.iter().map(|&i| logp_v[i.0 as usize]).collect(),
        heavy.iter().map(|&i| mr_v[i.0 as usize]).collect(),
        heavy.iter().map(|&i| mw_v[i.0 as usize]).collect(),
    ];
    // Connectivity (bonded) lookup.
    let mut bonded = vec![vec![0.0f64; n]; n];
    for (pi, &ai) in heavy.iter().enumerate() {
        for (pj, &aj) in heavy.iter().enumerate().skip(pi + 1) {
            if let Some((bidx, _)) = mol.bond_between(ai, aj) {
                let bo = match mol.bond(bidx).order {
                    chematic_core::BondOrder::Single => 1.0,
                    chematic_core::BondOrder::Double => 2.0,
                    chematic_core::BondOrder::Triple => 3.0,
                    chematic_core::BondOrder::Aromatic => 1.5,
                    chematic_core::BondOrder::Up | chematic_core::BondOrder::Down => 1.0,
                    _ => 1.0,
                };
                let zi = mol.atom(ai).element.atomic_number() as f64;
                let zj = mol.atom(aj).element.atomic_number() as f64;
                let w = (zi * zj).sqrt() * bo / 8.0;
                bonded[pi][pj] = w;
                bonded[pj][pi] = w;
            }
        }
    }
    let compute = |prop: &[f64]| -> (f64, f64) {
        let mut mat = vec![vec![0.001f64; n]; n];
        for i in 0..n {
            mat[i][i] = prop[i];
            for j in (i + 1)..n {
                if bonded[i][j] > 0.0 {
                    mat[i][j] = bonded[i][j];
                    mat[j][i] = bonded[i][j];
                }
            }
        }
        let eigs = jacobi_eigenvalues(&mat);
        let hi = eigs.iter().copied().fold(f64::MIN, f64::max);
        let lo = eigs.iter().copied().fold(f64::MAX, f64::min);
        (hi, lo)
    };
    let (chghi, chglo) = compute(&props[0]);
    let (logphi, logplo) = compute(&props[1]);
    let (mrhi, mrlo) = compute(&props[2]);
    let (mwhi, mwlo) = compute(&props[3]);
    Bcut2D {
        chghi,
        chglo,
        logphi,
        logplo,
        mrhi,
        mrlo,
        mwhi,
        mwlo,
    }
}

/// Jacobi eigenvalue algorithm: all eigenvalues of a real symmetric matrix.
///
/// Same classical rotation technique as `chematic_3d::shape_descriptors::jacobi3`
/// (3×3-specific, and not importable here since `chematic-chem` is a dependency
/// *of* `chematic-3d`, not the reverse), generalized to arbitrary N. Burden
/// matrices here are small (heavy-atom count, typically well under 100), so a
/// full eigendecomposition is simple and avoids power iteration's convergence
/// pitfall: unshifted power iteration converges to the eigenvalue of largest
/// *magnitude*, not the algebraic max — which made an earlier min/max split
/// via `min_eig(A) = -max_eig(-A)` always return the same value for both.
#[allow(clippy::needless_range_loop)] // direct (i, j) indexing matches the textbook algorithm
fn jacobi_eigenvalues(mat: &[Vec<f64>]) -> Vec<f64> {
    let n = mat.len();
    if n == 0 {
        return Vec::new();
    }
    let mut a: Vec<Vec<f64>> = mat.to_vec();
    for _ in 0..100 {
        // Find the largest-magnitude off-diagonal element.
        let (mut p, mut q, mut max_val) = (0, 1, 0.0f64);
        for i in 0..n {
            for j in (i + 1)..n {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if max_val < 1e-12 {
            break;
        }
        let (app, aqq, apq) = (a[p][p], a[q][q], a[p][q]);
        let tau = (aqq - app) / (2.0 * apq);
        let t = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
        } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        a[p][p] = app - t * apq;
        a[q][q] = aqq + t * apq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
        for i in 0..n {
            if i != p && i != q {
                let (aip, aiq) = (a[i][p], a[i][q]);
                a[i][p] = c * aip - s * aiq;
                a[p][i] = a[i][p];
                a[i][q] = s * aip + c * aiq;
                a[q][i] = a[i][q];
            }
        }
    }
    (0..n).map(|i| a[i][i]).collect()
}

