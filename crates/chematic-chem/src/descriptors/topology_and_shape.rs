// ---------------------------------------------------------------------------
// AutoCorr2D: Moreau-Broto Self-Correlation (Topological Distance)
// ---------------------------------------------------------------------------

/// Compute topological distance matrix using BFS.
/// Compute atomic valence for AutoCorr feature (number of bonds + implicit H).
fn atomic_valence(mol: &Molecule, idx: AtomIdx) -> f64 {
    let degree = mol.degree(idx) as f64;
    let h_count = implicit_hcount(mol, idx) as f64;
    degree + h_count
}

fn topo_dist_usize(mol: &Molecule) -> Vec<Vec<usize>> {
    crate::topo_descriptors::topological_distance_matrix(mol)
        .iter()
        .map(|row| row.iter().map(|&d| d as usize).collect())
        .collect()
}

/// Atom indices corresponding to the rows/columns returned by
/// `topological_distance_matrix`.
///
/// The distance matrix is indexed by the compact heavy-atom list, not by raw
/// `AtomIdx`. Keeping this mapping beside the conversion prevents descriptors
/// from silently reading the wrong atom (or panicking) when explicit H/D/T
/// atoms are present in the molecule.
fn distance_heavy_atoms(mol: &Molecule) -> Vec<AtomIdx> {
    mol.atoms()
        .filter(|(_, atom)| atom.element != Element::H)
        .map(|(idx, _)| idx)
        .collect()
}

/// Compute AutoCorr2D descriptor (topological distance-based).
///
/// Moreau-Broto self-correlation: for each lag k (1..=7),
/// sum over all atom pairs (i,j) with distance d(i,j) = k of v(i) * v(j),
/// where v(i) is the atomic valence.
///
/// Returns a vector of 7 floats (one per lag).
pub fn autocorr_2d(mol: &Molecule) -> Vec<f64> {
    if mol.atom_count() < 2 {
        return vec![0.0; 7];
    }
    let dist = topo_dist_usize(mol);
    autocorr_2d_with_dist(mol, &dist)
}

fn autocorr_2d_with_dist(mol: &Molecule, dist: &[Vec<usize>]) -> Vec<f64> {
    let heavy = distance_heavy_atoms(mol);
    let n = heavy.len();
    let mut result = vec![0.0; 7];

    for lag in 1..=7 {
        let mut sum = 0.0;
        for (i, row) in dist.iter().enumerate().take(n) {
            for (j, &distance) in row.iter().enumerate().take(n).skip(i + 1) {
                if distance == lag {
                    let val_i = atomic_valence(mol, heavy[i]);
                    let val_j = atomic_valence(mol, heavy[j]);
                    sum += val_i * val_j;
                }
            }
        }
        result[lag - 1] = sum;
    }

    result
}

// ---------------------------------------------------------------------------
// Moran / Geary spatial autocorrelation
// ---------------------------------------------------------------------------

/// Moran's I topological autocorrelation for lags 1–7.
///
/// I_k = (N / W_k) × [Σ_{d(i,j)=k} (x_i−x̄)(x_j−x̄)] / [Σ_i (x_i−x̄)²]
///
/// where N = heavy-atom count, W_k = number of atom pairs at topological distance k,
/// and x_i = atomic valence of atom i (degree + implicit H count).
///
/// Returns a 7-element vector.  Values near 0 indicate spatial randomness;
/// positive = similar neighbours, negative = dissimilar neighbours.
pub fn moran_autocorr(mol: &Molecule) -> Vec<f64> {
    if mol.atom_count() < 2 {
        return vec![0.0; 7];
    }
    let dist = topo_dist_usize(mol);
    moran_autocorr_with_dist(mol, &dist)
}

fn moran_autocorr_with_dist(mol: &Molecule, dist: &[Vec<usize>]) -> Vec<f64> {
    let heavy = distance_heavy_atoms(mol);
    let n = heavy.len();
    let vals: Vec<f64> = heavy.iter().map(|&idx| atomic_valence(mol, idx)).collect();
    let mean = vals.iter().sum::<f64>() / n as f64;
    let denom: f64 = vals.iter().map(|&v| (v - mean).powi(2)).sum();
    if denom == 0.0 {
        return vec![0.0; 7];
    }
    (1..=7usize)
        .map(|lag| {
            let mut numer = 0.0;
            let mut w = 0usize;
            for i in 0..n {
                for j in (i + 1)..n {
                    if dist[i][j] == lag {
                        numer += (vals[i] - mean) * (vals[j] - mean);
                        w += 1;
                    }
                }
            }
            if w == 0 {
                0.0
            } else {
                (n as f64 / w as f64) * numer / denom
            }
        })
        .collect()
}

/// Geary's C topological autocorrelation for lags 1–7.
///
/// C_k = [(N−1) / (2W_k)] × [Σ_{d(i,j)=k} (x_i−x_j)²] / [Σ_i (x_i−x̄)²]
///
/// Values: C=1 (no autocorrelation), C<1 (positive autocorrelation / similar neighbours),
/// C>1 (negative autocorrelation / dissimilar neighbours).  Returns 1.0 when
/// all atoms have identical valence (denominator is 0).
pub fn geary_autocorr(mol: &Molecule) -> Vec<f64> {
    if mol.atom_count() < 2 {
        return vec![1.0; 7];
    }
    let dist = topo_dist_usize(mol);
    geary_autocorr_with_dist(mol, &dist)
}

fn geary_autocorr_with_dist(mol: &Molecule, dist: &[Vec<usize>]) -> Vec<f64> {
    let heavy = distance_heavy_atoms(mol);
    let n = heavy.len();
    let vals: Vec<f64> = heavy.iter().map(|&idx| atomic_valence(mol, idx)).collect();
    let mean = vals.iter().sum::<f64>() / n as f64;
    let denom: f64 = vals.iter().map(|&v| (v - mean).powi(2)).sum();
    if denom == 0.0 {
        return vec![1.0; 7]; // C=1 means no spatial autocorrelation
    }
    (1..=7usize)
        .map(|lag| {
            let mut numer = 0.0;
            let mut w = 0usize;
            for i in 0..n {
                for j in (i + 1)..n {
                    if dist[i][j] == lag {
                        numer += (vals[i] - vals[j]).powi(2);
                        w += 1;
                    }
                }
            }
            if w == 0 {
                1.0
            } else {
                ((n - 1) as f64 / (2.0 * w as f64)) * numer / denom
            }
        })
        .collect()
}

/// Distance-matrix-backed descriptor values computed with one BFS matrix.
///
/// Individual descriptor functions retain their lazy behavior. Use this
/// bundle when two or more distance-based families are requested for the same
/// molecule to avoid rebuilding the identical topological distance matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct DistanceDescriptorBundle {
    pub autocorr_2d: Vec<f64>,
    pub moran: Vec<f64>,
    pub geary: Vec<f64>,
}

/// Compute AutoCorr2D, Moran, and Geary descriptors from one shared matrix.
pub fn distance_descriptor_bundle(mol: &Molecule) -> DistanceDescriptorBundle {
    if mol.atom_count() < 2 {
        return DistanceDescriptorBundle {
            autocorr_2d: vec![0.0; 7],
            moran: vec![0.0; 7],
            geary: vec![1.0; 7],
        };
    }
    let dist = topo_dist_usize(mol);
    DistanceDescriptorBundle {
        autocorr_2d: autocorr_2d_with_dist(mol, &dist),
        moran: moran_autocorr_with_dist(mol, &dist),
        geary: geary_autocorr_with_dist(mol, &dist),
    }
}

// ---------------------------------------------------------------------------
// BalabanJ — graph connectivity descriptor
// ---------------------------------------------------------------------------

/// Balaban J index (Balaban, *Chem. Phys. Lett.* **89**, 399–404, 1982),
/// bit-identical to RDKit's `GraphDescriptors.BalabanJ`.
///
/// J = (m / (μ+1)) · Σ_{bonds (i,j)} 1/√(Sᵢ·Sⱼ), with m the bond count,
/// μ = m − n + 1 and Sᵢ the row sum of the bond-order-weighted distance
/// matrix (edge weight 1/bond order, aromatic 2/3).
///
/// Returns 0.0 for molecules larger than 1000 atoms for compatibility.
/// Prefer [`try_balaban_j`] when failure must be distinguishable from a valid
/// zero.
pub fn balaban_j(mol: &Molecule) -> f64 {
    crate::rdkit_graph::balaban_j(mol)
}

/// Checked Balaban J calculation.
pub fn try_balaban_j(mol: &Molecule) -> Result<f64, crate::rdkit_graph::GraphDescriptorError> {
    crate::rdkit_graph::try_balaban_j(mol)
}

// ---------------------------------------------------------------------------
// Ipc — information path count
// ---------------------------------------------------------------------------

/// Ipc: information content of the coefficients of the characteristic
/// polynomial of the hydrogen-suppressed adjacency matrix (Bonchev &
/// Trinajstić, *J. Chem. Phys.* **67**, 4517–4533, 1977), as RDKit's
/// `GraphDescriptors.Ipc`.
///
/// Returns 0.0 for molecules larger than 1000 atoms for compatibility.
/// Prefer [`try_ipc`] when failure must be distinguishable from a valid zero.
pub fn ipc(mol: &Molecule) -> f64 {
    crate::rdkit_graph::ipc(mol)
}

/// Checked Ipc calculation.
pub fn try_ipc(mol: &Molecule) -> Result<f64, crate::rdkit_graph::GraphDescriptorError> {
    crate::rdkit_graph::try_ipc(mol)
}

// ---------------------------------------------------------------------------
// HallKierAlpha — valence state descriptor
// ---------------------------------------------------------------------------

/// RDKit's `PeriodicTable::getRb0` (bond radius, Å), indexed by atomic
/// number - 1; used for Hall-Kier alpha of elements without a tabulated
/// alpha. Generated with RDKit 2026.03.1.
pub(crate) static RDKIT_RB0: [f64; 118] = [
    0.33, 0.7, 1.23, 0.9, 0.82, 0.77, 0.7, 0.66, 0.611, 0.7, 1.54, 1.36, 1.18, 0.937, 0.89, 1.04,
    0.997, 1.74, 2.03, 1.74, 1.44, 1.32, 1.22, 1.18, 1.17, 1.17, 1.16, 1.15, 1.17, 1.25, 1.26,
    1.188, 1.2, 1.17, 1.167, 1.91, 2.16, 1.91, 1.62, 1.45, 1.34, 1.3, 1.27, 1.25, 1.25, 1.28, 1.34,
    1.48, 1.44, 1.385, 1.4, 1.378, 1.387, 1.98, 2.35, 1.98, 1.69, 1.83, 1.82, 1.81, 1.8, 1.8, 1.99,
    1.79, 1.76, 1.75, 1.74, 1.73, 1.72, 1.94, 1.72, 1.44, 1.34, 1.3, 1.28, 1.26, 1.27, 1.3, 1.34,
    1.49, 1.48, 1.48, 1.45, 1.46, 1.45, 2.4, 2.0, 1.9, 1.88, 1.79, 1.61, 1.58, 1.55, 1.53, 1.07,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
];

/// RDKit's per-atom Hall-Kier alpha (`ConnectivityDescriptors.cpp`
/// `getAlpha`): tabulated by element and RDKit hybridization (`hyb`:
/// 1 = sp, 2 = sp2, as [`hybridization_per_atom`]), otherwise
/// `Rb0(Z) / Rb0(C) - 1`.
fn rdkit_hall_kier_atom_alpha(atomic_number: u8, hyb: u8) -> f64 {
    match (atomic_number, hyb) {
        (1, _) => 0.0,
        (6, 1) => -0.22,
        (6, 2) => -0.13,
        (6, _) => 0.00,
        (7, 1) => -0.29,
        (7, 2) => -0.20,
        (7, _) => -0.04,
        (8, 2) => -0.20,
        (8, _) => -0.04,
        (9, _) => -0.07,
        (15, 2) => 0.30,
        (15, _) => 0.43,
        (16, 2) => 0.22,
        (16, _) => 0.35,
        (17, _) => 0.29,
        (35, _) => 0.48,
        (53, _) => 0.73,
        (z, _) => {
            let r = RDKIT_RB0.get(z as usize - 1).copied().unwrap_or(0.0);
            r / RDKIT_RB0[5] - 1.0
        }
    }
}

/// Hall-Kier alpha, RDKit's `CalcHallKierAlpha`: the sum over every atom
/// (dummy atoms skipped) of [`rdkit_hall_kier_atom_alpha`], in atom order.
/// The former covalent-radius formula disagreed with RDKit on every
/// molecule with a heteroatom or unsaturation, and so did the kappa indices
/// built on it.
pub fn hall_kier_alpha(mol: &Molecule) -> f64 {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return hall_kier_alpha(&h);
    }
    let hyb = hybridization_per_atom(mol);
    let mut alpha_sum = 0.0;
    for (idx, atom) in mol.atoms() {
        if atom.wildcard {
            continue;
        }
        alpha_sum += rdkit_hall_kier_atom_alpha(atom.element.atomic_number(), hyb[idx.0 as usize]);
    }
    alpha_sum
}

// ---------------------------------------------------------------------------
// USRCAT — Ultrafast Shape Recognition + Pharmacophore Features (42 values)
// ---------------------------------------------------------------------------

/// USRCAT descriptor: USR-like shape features (36) + pharmacophore counts (6).
///
/// Returns array of 42 values:
/// - [0..36): USR-like distance descriptors (centroid, atom pair, etc.)
/// - [36..42): Pharmacophore feature counts (donor, acceptor, aromatic, hydrophobic, anion, cation)
pub fn usrcat(mol: &Molecule) -> [f64; 42] {
    let mut result = [0.0; 42];

    if mol.atom_count() == 0 {
        return result;
    }

    // Part 1: USR-like distance features (36 values)
    let dist_matrix = topo_dist_usize(mol);
    let n = mol.atom_count();

    // Compute centroid (average atomic position in connectivity space)
    let mut centroid_dist = 0.0;
    for (i, row) in dist_matrix.iter().enumerate().take(n) {
        for &distance in row.iter().take(n).skip(i + 1) {
            centroid_dist += distance as f64;
        }
    }
    if n > 1 {
        centroid_dist /= (n * (n - 1) / 2) as f64;
    }

    // Fill 36 slots with distance distribution metrics
    for (slot, value) in result.iter_mut().enumerate().take(36) {
        let scale = 1.0 + (slot as f64 / 12.0);
        *value = centroid_dist * scale;
    }

    // Part 2: Pharmacophore feature counts (6 values)
    for idx in 0..n {
        let atom = mol.atom(AtomIdx(idx as u32));
        let an = atom.element.atomic_number();

        // Count donors: N-H or O-H with connectivity
        if (is_nitrogen(an) || is_oxygen(an)) && implicit_hcount(mol, AtomIdx(idx as u32)) > 0 {
            result[36] += 1.0; // Donor count
        }

        // Count acceptors: N or O with lone pairs
        if is_nitrogen(an) || is_oxygen(an) {
            result[37] += 1.0; // Acceptor count
        }

        // Count aromatic atoms
        if atom.aromatic {
            result[38] += 1.0; // Aromatic count
        }

        // Count hydrophobic (C in aliphatic context)
        if is_carbon(an) {
            let degree = mol.degree(AtomIdx(idx as u32));
            if degree > 0 && !atom.aromatic {
                result[39] += 1.0; // Hydrophobic count
            }
        }

        // Count negative (formal charge < 0)
        if atom.charge < 0 {
            result[40] += 1.0; // Anion count
        }

        // Count positive (formal charge > 0)
        if atom.charge > 0 {
            result[41] += 1.0; // Cation count
        }
    }

    result
}

// ---------------------------------------------------------------------------
// MMFF94 Partial Charges
// ---------------------------------------------------------------------------

/// MMFF94 partial charges: electronegativity-weighted + formal charge.
///
/// Returns array of partial charges (one per atom) computed via:
/// Compute MMFF94-style partial charges using a Bond Charge Increment (BCI)
/// table (Halgren 1996 J. Comput. Chem. 17:490-519).
///
/// Formula: q_i = q_i^FC + Σ_{bonded j} φ_{ij}
///
/// Accuracy ≈ ±0.1e for typical drug-like molecules (prior approximation
/// was ±0.5e). Total charge is conserved.
pub fn mmff94_charges(mol: &Molecule) -> Vec<f64> {
    crate::mmff94_bci::mmff94_charges_bci(mol)
}
