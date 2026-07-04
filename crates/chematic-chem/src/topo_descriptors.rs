//! Topological molecular descriptors.
//!
//! Implements the Wiener topological index, Hall-Kier Kappa shape indices
//! (κ1/κ2/κ3), Kier-Hall Chi connectivity indices (χ0–χ4 and χ0v–χ4v),
//! Bertz complexity, and Labute approximate surface area (LabuteASA).
//!
//! All descriptors except LabuteASA operate on the heavy-atom subgraph
//! (hydrogen atoms excluded from path/distance calculations).

use rustc_hash::{FxHashMap, FxHashSet};
use std::collections::VecDeque;
use std::f64::consts::PI;

use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};

// ─── Helpers ────────────────────────────────────────────────────────────────

/// Indices of all heavy (non-hydrogen) atoms.
fn heavy_indices(mol: &Molecule) -> Vec<usize> {
    mol.atoms()
        .filter(|(_, a)| a.element.atomic_number() != 1)
        .map(|(idx, _)| idx.0 as usize)
        .collect()
}

/// Heavy-atom degree for each atom index (0 for H atoms).
///
/// Avoids repeating the neighbor-filter loop inside bond-level loops.
fn heavy_degrees(mol: &Molecule) -> Vec<u32> {
    let n = mol.atom_count();
    let is_heavy: Vec<bool> = (0..n)
        .map(|i| mol.atom(AtomIdx(i as u32)).element.atomic_number() != 1)
        .collect();
    (0..n)
        .map(|i| {
            if !is_heavy[i] {
                return 0;
            }
            mol.neighbors(AtomIdx(i as u32))
                .filter(|(nb, _)| is_heavy[nb.0 as usize])
                .count() as u32
        })
        .collect()
}

/// BFS shortest-path distances from `start` in the heavy-atom subgraph.
/// Returns `usize::MAX` for disconnected pairs or hydrogen atoms.
fn bfs_from(mol: &Molecule, start: usize, heavy_set: &FxHashSet<usize>) -> Vec<usize> {
    let n = mol.atom_count();
    let mut dist = vec![usize::MAX; n];
    dist[start] = 0;
    let mut queue = VecDeque::new();
    queue.push_back(start);
    while let Some(cur) = queue.pop_front() {
        let d = dist[cur];
        for (nb, _) in mol.neighbors(AtomIdx(cur as u32)) {
            let ni = nb.0 as usize;
            if heavy_set.contains(&ni) && dist[ni] == usize::MAX {
                dist[ni] = d + 1;
                queue.push_back(ni);
            }
        }
    }
    dist
}

/// Connectivity delta (degree in heavy-atom graph) for atom `idx`.
fn delta(mol: &Molecule, idx: AtomIdx, heavy_set: &FxHashSet<usize>) -> f64 {
    mol.neighbors(idx)
        .filter(|(nb, _)| heavy_set.contains(&(nb.0 as usize)))
        .count() as f64
}

/// Valence-corrected delta: δᵥ = (Zᵥ − H) / (Z − Zᵥ − 1).
fn delta_v(mol: &Molecule, idx: AtomIdx) -> f64 {
    let atom = mol.atom(idx);
    let z = atom.element.atomic_number() as i32;
    let zv = valence_electrons(atom.element.atomic_number()) as i32;
    let h = implicit_hcount(mol, idx) as i32 + atom.hydrogen_count.unwrap_or(0) as i32;
    let denom = z - zv - 1;
    if denom <= 0 {
        return (zv - h).max(1) as f64;
    }
    ((zv - h).max(0)) as f64 / denom as f64
}

fn valence_electrons(z: u8) -> u8 {
    match z {
        1 => 1,
        2 => 2,
        3 => 1,
        4 => 2,
        5 => 3,
        6 => 4,
        7 => 5,
        8 => 6,
        9 => 7,
        10 => 8,
        11 => 1,
        12 => 2,
        13 => 3,
        14 => 4,
        15 => 5,
        16 => 6,
        17 => 7,
        18 => 8,
        35 => 7,
        53 => 7,
        _ => z.min(8),
    }
}

/// DFS: accumulate chi contributions for paths of exactly `target_len` bonds
/// starting from `cur`.  `running_product` is the product of delta values of
/// atoms visited so far (including `cur`).
#[allow(clippy::too_many_arguments)]
fn chi_dfs(
    mol: &Molecule,
    cur: usize,
    target_len: usize,
    cur_len: usize,
    running_product: f64,
    visited: &mut Vec<bool>,
    heavy_set: &FxHashSet<usize>,
    use_valence: bool,
) -> f64 {
    if cur_len == target_len {
        return running_product.powf(-0.5);
    }
    let mut sum = 0.0f64;
    for (nb, _) in mol.neighbors(AtomIdx(cur as u32)) {
        let ni = nb.0 as usize;
        if heavy_set.contains(&ni) && !visited[ni] {
            let d_nb = if use_valence {
                delta_v(mol, AtomIdx(ni as u32))
            } else {
                delta(mol, AtomIdx(ni as u32), heavy_set)
            };
            if d_nb > 0.0 {
                visited[ni] = true;
                sum += chi_dfs(
                    mol,
                    ni,
                    target_len,
                    cur_len + 1,
                    running_product * d_nb,
                    visited,
                    heavy_set,
                    use_valence,
                );
                visited[ni] = false;
            }
        }
    }
    sum
}

/// Count simple paths of exactly `length` bonds in the heavy-atom subgraph.
/// Returns undirected path count (each path counted once).
fn count_paths(mol: &Molecule, heavy: &[usize], length: usize) -> usize {
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();
    let mut total = 0usize;
    for &start in heavy {
        let mut visited = vec![false; mol.atom_count()];
        visited[start] = true;
        total += count_paths_dfs(mol, start, length, 0, &mut visited, &heavy_set);
    }
    total / 2
}

fn count_paths_dfs(
    mol: &Molecule,
    cur: usize,
    target_len: usize,
    cur_len: usize,
    visited: &mut Vec<bool>,
    heavy_set: &FxHashSet<usize>,
) -> usize {
    if cur_len == target_len {
        return 1;
    }
    let mut count = 0;
    for (nb, _) in mol.neighbors(AtomIdx(cur as u32)) {
        let ni = nb.0 as usize;
        if heavy_set.contains(&ni) && !visited[ni] {
            visited[ni] = true;
            count += count_paths_dfs(mol, ni, target_len, cur_len + 1, visited, heavy_set);
            visited[ni] = false;
        }
    }
    count
}

/// Compute chi sum for paths of exactly `n` bonds (n ≥ 1).
/// Each undirected path is counted once (sum divided by 2).
fn chi_n_with(
    mol: &Molecule,
    heavy: &[usize],
    heavy_set: &FxHashSet<usize>,
    n: usize,
    use_valence: bool,
) -> f64 {
    let mut total = 0.0f64;
    for &start in heavy {
        let d_start = if use_valence {
            delta_v(mol, AtomIdx(start as u32))
        } else {
            delta(mol, AtomIdx(start as u32), heavy_set)
        };
        if d_start <= 0.0 {
            continue;
        }
        let mut visited = vec![false; mol.atom_count()];
        visited[start] = true;
        total += chi_dfs(
            mol,
            start,
            n,
            0,
            d_start,
            &mut visited,
            heavy_set,
            use_valence,
        );
    }
    total / 2.0
}

fn chi_n(mol: &Molecule, n: usize, use_valence: bool) -> f64 {
    let heavy = heavy_indices(mol);
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();
    chi_n_with(mol, &heavy, &heavy_set, n, use_valence)
}

// ─── Wiener Index ────────────────────────────────────────────────────────────

/// Wiener topological index.
///
/// Sum of all pairwise shortest-path distances between heavy atoms.
/// Computed on the hydrogen-depleted graph.
pub fn wiener_index(mol: &Molecule) -> f64 {
    let heavy = heavy_indices(mol);
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();
    let mut sum = 0u64;
    for i in 0..heavy.len() {
        let row = bfs_from(mol, heavy[i], &heavy_set);
        for j in (i + 1)..heavy.len() {
            let d = row[heavy[j]];
            if d != usize::MAX {
                sum += d as u64;
            }
        }
    }
    sum as f64
}

// ─── Padmakar-Ivan (PI) Index ────────────────────────────────────────────────

/// Padmakar-Ivan (PI) topological index (Khadikar et al. 2001).
///
/// For each bond e = (u, v) in the heavy-atom graph, let:
/// - n_u(e) = number of heavy atoms strictly closer to u than to v
/// - n_v(e) = number of heavy atoms strictly closer to v than to u
///
/// PI = Σ_e [n_u(e) + n_v(e)]
///
/// Reference values: ethane = 2, propane = 6, butane = 12, benzene = 36.
pub fn padmakar_ivan_index(mol: &Molecule) -> u64 {
    let heavy = heavy_indices(mol);
    let n = heavy.len();
    if n < 2 {
        return 0;
    }
    // Guard: O(n²) distance matrix would OOM for large molecules (same pattern as hosoya_index).
    if n > 1000 {
        return u64::MAX;
    }
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();

    // Map original atom index → compressed heavy-atom position 0..n
    let mut pos: FxHashMap<usize, usize> =
        FxHashMap::with_capacity_and_hasher(n, Default::default());
    for (p, &h) in heavy.iter().enumerate() {
        pos.insert(h, p);
    }

    // Full BFS distance matrix for heavy atoms: dist[p][q] = d(heavy[p], heavy[q])
    let mut dist = vec![vec![usize::MAX; n]; n];
    for p in 0..n {
        dist[p][p] = 0;
        let row = bfs_from(mol, heavy[p], &heavy_set);
        for q in 0..n {
            let d = row[heavy[q]];
            if d != usize::MAX {
                dist[p][q] = d;
            }
        }
    }

    // Sum n_u + n_v over each heavy-atom bond
    let mut pi_val = 0u64;
    for (_, bond) in mol.bonds() {
        let u = bond.atom1.0 as usize;
        let v = bond.atom2.0 as usize;
        if !heavy_set.contains(&u) || !heavy_set.contains(&v) {
            continue;
        }
        let pu = pos[&u];
        let pv = pos[&v];

        let mut n_u = 0u64;
        let mut n_v = 0u64;
        for (du, dv) in dist[pu].iter().zip(dist[pv].iter()) {
            if du < dv {
                n_u += 1;
            } else if dv < du {
                n_v += 1;
            }
            // equidistant vertices contribute 0
        }
        pi_val += n_u + n_v;
    }
    pi_val
}

// ─── Kappa Shape Indices ─────────────────────────────────────────────────────

/// Hall-Kier κ1 shape index (alpha-corrected, matches RDKit `CalcKappa1`).
///
/// κ1 = (A+α)·(A+α−1)² / (P1+α)²  where A = heavy atom count, P1 = bond
/// count, α = [`hall_kier_alpha`](crate::descriptors::hall_kier_alpha).
/// A larger value indicates a more linear graph.
pub fn kappa1(mol: &Molecule) -> f64 {
    let heavy = heavy_indices(mol);
    let n = heavy.len();
    if n < 2 {
        return 0.0;
    }
    let p1 = count_paths(mol, &heavy, 1);
    if p1 == 0 {
        return 0.0;
    }
    let alpha = crate::descriptors::hall_kier_alpha(mol);
    let a = n as f64 + alpha;
    let p1 = p1 as f64 + alpha;
    a * (a - 1.0).powi(2) / p1.powi(2)
}

/// Hall-Kier κ2 shape index (alpha-corrected, matches RDKit `CalcKappa2`).
///
/// κ2 = (A+α−1)·(A+α−2)² / (P2+α)²  where P2 = count of 2-bond paths.
pub fn kappa2(mol: &Molecule) -> f64 {
    let heavy = heavy_indices(mol);
    let n = heavy.len();
    if n < 3 {
        return 0.0;
    }
    let p2 = count_paths(mol, &heavy, 2);
    if p2 == 0 {
        return 0.0;
    }
    let alpha = crate::descriptors::hall_kier_alpha(mol);
    let a = n as f64 + alpha;
    let p2 = p2 as f64 + alpha;
    (a - 1.0) * (a - 2.0).powi(2) / p2.powi(2)
}

/// Hall-Kier κ3 shape index (alpha-corrected, matches RDKit `CalcKappa3`).
///
/// Formula depends on parity of heavy-atom count:
/// - odd n:  κ3 = (A+α−1)·(A+α−3)² / (P3+α)²
/// - even n: κ3 = (A+α−2)·(A+α−3)² / (P3+α)²
///
/// Returns 0.0 when fewer than 4 heavy atoms or no 3-bond paths exist.
pub fn kappa3(mol: &Molecule) -> f64 {
    let heavy = heavy_indices(mol);
    let n = heavy.len();
    if n < 4 {
        return 0.0;
    }
    let p3 = count_paths(mol, &heavy, 3);
    if p3 == 0 {
        return 0.0;
    }
    let alpha = crate::descriptors::hall_kier_alpha(mol);
    let a = n as f64 + alpha;
    let p3 = p3 as f64 + alpha;
    let factor = if n % 2 == 1 { a - 1.0 } else { a - 2.0 };
    factor * (a - 3.0).powi(2) / p3.powi(2)
}

/// Compute κ1, κ2, κ3 in a single `heavy_indices` pass.
///
/// Returns `(κ1, κ2, κ3)`. Use when all three are needed to avoid
/// three redundant `heavy_indices` computations.
pub fn kappa_all(mol: &Molecule) -> (f64, f64, f64) {
    let heavy = heavy_indices(mol);
    let n = heavy.len();
    let alpha = crate::descriptors::hall_kier_alpha(mol);
    let a = n as f64 + alpha;

    let k1 = if n >= 2 {
        let p1 = count_paths(mol, &heavy, 1);
        if p1 == 0 {
            0.0
        } else {
            let p1 = p1 as f64 + alpha;
            a * (a - 1.0).powi(2) / p1.powi(2)
        }
    } else {
        0.0
    };

    let k2 = if n >= 3 {
        let p2 = count_paths(mol, &heavy, 2);
        if p2 == 0 {
            0.0
        } else {
            let p2 = p2 as f64 + alpha;
            (a - 1.0) * (a - 2.0).powi(2) / p2.powi(2)
        }
    } else {
        0.0
    };

    let k3 = if n >= 4 {
        let p3 = count_paths(mol, &heavy, 3);
        if p3 == 0 {
            0.0
        } else {
            let p3 = p3 as f64 + alpha;
            let factor = if n % 2 == 1 { a - 1.0 } else { a - 2.0 };
            factor * (a - 3.0).powi(2) / p3.powi(2)
        }
    } else {
        0.0
    };

    (k1, k2, k3)
}

// ─── Chi Connectivity Indices ────────────────────────────────────────────────

/// Kier-Hall χ0 connectivity index.
///
/// χ0 = Σᵢ δᵢ^(−0.5) over all heavy atoms, where δᵢ = heavy-atom degree.
/// Atoms with δ = 0 contribute 0.
pub fn chi0(mol: &Molecule) -> f64 {
    let heavy = heavy_indices(mol);
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();
    heavy
        .iter()
        .map(|&i| {
            let d = delta(mol, AtomIdx(i as u32), &heavy_set);
            if d > 0.0 { d.powf(-0.5) } else { 0.0 }
        })
        .sum()
}

/// Kier-Hall χ1 connectivity index (bond-path sum).
pub fn chi1(mol: &Molecule) -> f64 {
    chi_n(mol, 1, false)
}

/// Kier-Hall χ2 connectivity index (2-bond path sum).
pub fn chi2(mol: &Molecule) -> f64 {
    chi_n(mol, 2, false)
}

/// Kier-Hall χ3 connectivity index (3-bond path sum).
pub fn chi3(mol: &Molecule) -> f64 {
    chi_n(mol, 3, false)
}

/// Kier-Hall χ4 connectivity index (4-bond path sum).
pub fn chi4(mol: &Molecule) -> f64 {
    chi_n(mol, 4, false)
}

/// Valence-corrected χ0v connectivity index.
///
/// Uses δᵥ = (Zᵥ − H) / (Z − Zᵥ − 1) instead of the simple degree.
pub fn chi0v(mol: &Molecule) -> f64 {
    let heavy = heavy_indices(mol);
    heavy
        .iter()
        .map(|&i| {
            let d = delta_v(mol, AtomIdx(i as u32));
            if d > 0.0 { d.powf(-0.5) } else { 0.0 }
        })
        .sum()
}

/// Valence-corrected χ1v connectivity index.
pub fn chi1v(mol: &Molecule) -> f64 {
    chi_n(mol, 1, true)
}

/// Valence-corrected χ2v connectivity index.
pub fn chi2v(mol: &Molecule) -> f64 {
    chi_n(mol, 2, true)
}

/// Valence-corrected χ3v connectivity index.
pub fn chi3v(mol: &Molecule) -> f64 {
    chi_n(mol, 3, true)
}

/// Valence-corrected χ4v connectivity index.
pub fn chi4v(mol: &Molecule) -> f64 {
    chi_n(mol, 4, true)
}

/// Compute all 10 Hall-Kier connectivity indices in a single pass.
///
/// Returns `(χ0, χ1, χ2, χ3, χ4, χ0v, χ1v, χ2v, χ3v, χ4v)`.
/// Use when all indices are needed to avoid 10 redundant `heavy_indices` computations.
pub fn chi_all(mol: &Molecule) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let heavy = heavy_indices(mol);
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();

    let c0 = heavy
        .iter()
        .map(|&i| {
            let d = delta(mol, AtomIdx(i as u32), &heavy_set);
            if d > 0.0 { d.powf(-0.5) } else { 0.0 }
        })
        .sum();
    let c0v = heavy
        .iter()
        .map(|&i| {
            let d = delta_v(mol, AtomIdx(i as u32));
            if d > 0.0 { d.powf(-0.5) } else { 0.0 }
        })
        .sum();

    (
        c0,
        chi_n_with(mol, &heavy, &heavy_set, 1, false),
        chi_n_with(mol, &heavy, &heavy_set, 2, false),
        chi_n_with(mol, &heavy, &heavy_set, 3, false),
        chi_n_with(mol, &heavy, &heavy_set, 4, false),
        c0v,
        chi_n_with(mol, &heavy, &heavy_set, 1, true),
        chi_n_with(mol, &heavy, &heavy_set, 2, true),
        chi_n_with(mol, &heavy, &heavy_set, 3, true),
        chi_n_with(mol, &heavy, &heavy_set, 4, true),
    )
}

// ─── Bertz Complexity ────────────────────────────────────────────────────────

/// Simplified Bertz CT molecular complexity index.
///
/// CT = m_total + Σᵢ C(deg_total_i, 2)
///
/// where m_total = total bond count including implicit C-H bonds,
/// deg_total_i = heavy-atom degree + implicit H count for atom i, and
/// C(n, 2) = n·(n−1)/2.  This is the additive topology formula from
/// Bertz (1981) JACS 103, 3599 without logarithmic weighting.
pub fn bertz_ct(mol: &Molecule) -> f64 {
    let mut total_h_bonds = 0u64;
    let mut complexity = 0.0f64;
    for (idx, _) in mol.atoms() {
        let heavy_deg = mol.degree(idx);
        let h = implicit_hcount(mol, idx) as usize;
        total_h_bonds += h as u64;
        let total_deg = heavy_deg + h;
        complexity += (total_deg * total_deg.saturating_sub(1) / 2) as f64;
    }
    let heavy_bonds = mol.bond_count() as u64;
    let m_total = heavy_bonds + total_h_bonds;
    complexity + m_total as f64
}

// ─── Labute ASA ─────────────────────────────────────────────────────────────

/// Covalent (Rb0) radius of an element (Å), from RDKit's ptable.GetRb0.
/// Returns 0.0 for unrecognized elements (they contribute no surface area).
fn rb0(atomic_number: u8) -> f64 {
    match atomic_number {
        1 => 0.33,   // H
        6 => 0.77,   // C
        7 => 0.70,   // N
        8 => 0.66,   // O
        9 => 0.611,  // F
        14 => 1.04,  // Si
        15 => 0.89,  // P
        16 => 1.04,  // S
        17 => 0.997, // Cl
        33 => 1.21,  // As
        34 => 1.20,  // Se
        35 => 1.167, // Br
        53 => 1.387, // I
        _ => 0.0,
    }
}

/// Bond-type scale factor used in the Labute formula (Å subtracted from Ri+Rj).
///
/// Shorter bonds (double, triple, aromatic) bring atoms closer, increasing
/// surface overlap.  Single bonds have scale 0 (spheres just touching, no overlap).
fn bond_scale(order: BondOrder) -> f64 {
    match order {
        BondOrder::Aromatic => 0.1,
        BondOrder::Single
        | BondOrder::Up
        | BondOrder::Down
        | BondOrder::Zero
        | BondOrder::Dative
        | BondOrder::QueryAny
        | BondOrder::QuerySingleOrDouble
        | BondOrder::QuerySingleOrAromatic
        | BondOrder::QueryDoubleOrAromatic => 0.0,
        BondOrder::Double => 0.2,
        BondOrder::Triple | BondOrder::Quadruple => 0.3,
    }
}

/// Per-atom Labute approximate surface area contributions (Å²) plus the
/// pooled implicit-hydrogen area term, computed together in one pass.
///
/// Implements: P. Labute, 2000, *J. Mol. Graph. Mod.* **18**, 464–477, matching
/// RDKit's `_CalcLabuteASAContribs`: implicit hydrogens are **not** counted
/// per atom. Instead each heavy atom contributes exactly once (regardless of
/// its actual implicit H count) to a single molecule-wide pooled hydrogen
/// term, which is excluded from the per-atom values (RDKit's `ats`) and only
/// added into the whole-molecule total ([`labute_asa`]). This is a faithful,
/// numerically-verified port of RDKit's behavior, not a simplification.
fn labute_asa_parts(mol: &Molecule) -> (Vec<f64>, f64) {
    let n = mol.atom_count();
    if n == 0 {
        return (Vec::new(), 0.0);
    }

    const R_H: f64 = 0.33;
    let mut v: Vec<f64> = vec![0.0; n];
    let mut h_pool = 0.0f64;
    let radii: Vec<f64> = (0..n)
        .map(|i| rb0(mol.atom(AtomIdx(i as u32)).element.atomic_number()))
        .collect();

    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        let ri = radii[i];
        let rj = radii[j];
        if ri < 1e-10 || rj < 1e-10 {
            continue;
        }
        let scale = bond_scale(bond.order);
        let bij = ri + rj - scale;
        let dij = (ri - rj).abs().max(bij).min(ri + rj);
        v[i] += rj * rj - (ri - dij) * (ri - dij) / dij;
        v[j] += ri * ri - (rj - dij) * (rj - dij) / dij;
    }

    for i in 0..n {
        let ri = radii[i];
        if ri < 1e-10 {
            continue;
        }
        // Runs once per heavy atom regardless of its actual implicit H
        // count — see doc comment above.
        let dij = ri + R_H;
        v[i] += R_H * R_H - (ri - dij) * (ri - dij) / dij;
        h_pool += ri * ri - (R_H - dij) * (R_H - dij) / dij;
    }

    let per_atom = (0..n)
        .map(|i| {
            let ri = radii[i];
            if ri < 1e-10 {
                return 0.0;
            }
            (4.0 * PI * ri * ri - PI * ri * v[i]).max(0.0)
        })
        .collect();
    let h_pool_area = (4.0 * PI * R_H * R_H - PI * R_H * h_pool).max(0.0);

    (per_atom, h_pool_area)
}

/// Per-atom Labute approximate surface area contributions (Å²), excluding
/// the pooled implicit-hydrogen term (see [`labute_asa_parts`]). This is the
/// per-atom weight used by the VSA descriptor families
/// ([`crate::vsa`]), matching RDKit's `ats` output.
pub fn labute_asa_per_atom(mol: &Molecule) -> Vec<f64> {
    labute_asa_parts(mol).0
}

/// Pooled implicit-hydrogen area term (Å²) excluded from
/// [`labute_asa_per_atom`] but included in [`labute_asa`]'s total.
/// Only used by `vsa.rs` tests that check the VSA-sum-vs-total invariant.
#[cfg(test)]
pub(crate) fn labute_h_pool_area(mol: &Molecule) -> f64 {
    labute_asa_parts(mol).1
}

/// Labute approximate surface area (Å²).
///
/// Implements: P. Labute, 2000, *J. Mol. Graph. Mod.* **18**, 464–477.
///
/// Formula per atom i:
/// ```text
/// V_i  = Σ_j max(0, (Rj² − (Ri − dij)²) / dij)
/// A_i  = max(0, 4π Ri² − π Ri V_i)
/// ASA  = Σ A_i
/// ```
///
/// Bond distance: `dij = clamp(|Ri−Rj|, Ri+Rj−scale, Ri+Rj)`.
/// Implicit H atoms (radius 0.33 Å, single-bond scale 0) are included.
/// Randić connectivity index (χ).
///
/// χ = Σ_{bonds} 1 / √(deg(u) × deg(v))
///
/// Measures branching: lower values = more branched.
pub fn randic_index(mol: &Molecule) -> f64 {
    let deg = heavy_degrees(mol);
    let mut sum = 0.0f64;
    for i in 0..mol.bond_count() {
        let bond = mol.bond(chematic_core::BondIdx(i as u32));
        let da = deg[bond.atom1.0 as usize] as f64;
        let db = deg[bond.atom2.0 as usize] as f64;
        if da > 0.0 && db > 0.0 {
            sum += 1.0 / (da * db).sqrt();
        }
    }
    sum
}

/// Zagreb topological index M1.
///
/// M1 = Σ_{atoms} deg(v)²  (heavy-atom graph only).
pub fn zagreb_index_m1(mol: &Molecule) -> u32 {
    heavy_degrees(mol).iter().map(|&d| d * d).sum()
}

/// Zagreb index M2 (second Zagreb index).
///
/// M2 = Σ (deg(a) × deg(b)) for each heavy-atom bond (a, b).
///
/// Complements [`zagreb_index_m1`] (Σ deg(v)²); both quantify molecular branching.
/// Higher M2 indicates more branching or denser connectivity.
pub fn zagreb_index_m2(mol: &Molecule) -> u32 {
    let deg = heavy_degrees(mol);
    let mut sum = 0u32;
    for bidx in 0..mol.bond_count() {
        let bond = mol.bond(chematic_core::BondIdx(bidx as u32));
        let da = deg[bond.atom1.0 as usize];
        let db = deg[bond.atom2.0 as usize];
        // Skip bonds involving H atoms (degree == 0).
        if da > 0 && db > 0 {
            sum += da * db;
        }
    }
    sum
}

// ─── Eccentricity-based descriptors ─────────────────────────────────────────

/// Eccentricity of each heavy atom: max shortest-path distance to any other heavy atom.
///
/// Returns a vector of length `n_heavy`.  For a disconnected molecule the eccentricity
/// of isolated atoms is 0.
pub fn graph_eccentricities(mol: &Molecule) -> Vec<u32> {
    let heavy = heavy_indices(mol);
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();
    heavy
        .iter()
        .map(|&h| {
            let row = bfs_from(mol, h, &heavy_set);
            heavy
                .iter()
                .filter_map(|&h2| {
                    let d = row[h2];
                    if d == usize::MAX {
                        None
                    } else {
                        Some(d as u32)
                    }
                })
                .max()
                .unwrap_or(0)
        })
        .collect()
}

/// Graph diameter: maximum eccentricity over all heavy atoms.
///
/// Equals the longest shortest path in the heavy-atom graph.
pub fn graph_diameter(mol: &Molecule) -> u32 {
    graph_eccentricities(mol).into_iter().max().unwrap_or(0)
}

/// Graph radius: minimum eccentricity over all heavy atoms.
pub fn graph_radius(mol: &Molecule) -> u32 {
    graph_eccentricities(mol).into_iter().min().unwrap_or(0)
}

/// Petitjean topological index: `(diameter - radius) / diameter`.
///
/// Ranges from 0 (perfectly symmetric / linear chain) to 1 (highly asymmetric).
/// Returns 0 for single-atom molecules or when the diameter is 0.
pub fn petitjean_index(mol: &Molecule) -> f64 {
    let ecc = graph_eccentricities(mol);
    if ecc.is_empty() {
        return 0.0;
    }
    let d = ecc.iter().copied().max().unwrap_or(0);
    let r = ecc.iter().copied().min().unwrap_or(0);
    if d == 0 {
        0.0
    } else {
        (d - r) as f64 / d as f64
    }
}

/// Eccentric Connectivity Index: Σ_v [deg(v) × ecc(v)] over heavy atoms.
///
/// Introduced by Sharma et al. (1997).  Higher values indicate more
/// branched or elongated structures.
pub fn eccentric_connectivity_index(mol: &Molecule) -> u64 {
    let heavy = heavy_indices(mol);
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();
    // Reuse graph_eccentricities to avoid a separate O(n²) BFS pass.
    let ecc = graph_eccentricities(mol);
    heavy
        .iter()
        .zip(ecc.iter())
        .map(|(&h, &e)| {
            let deg = mol
                .neighbors(AtomIdx(h as u32))
                .filter(|(nb, _)| heavy_set.contains(&(nb.0 as usize)))
                .count() as u64;
            deg * e as u64
        })
        .sum()
}

// ─── Hosoya Index ────────────────────────────────────────────────────────────

/// Hosoya topological index Z: total number of matchings (including the empty matching).
///
/// Z(G) = Σ_k p(G, k) where p(G, k) is the number of k-matchings.
/// Computed via the vertex-removal recurrence: Z(G) = Z(G−v) + Σ_{u∈N(v)} Z(G−v−u).
///
/// Practical for drug-like molecules (< 60 heavy atoms).  For larger graphs
/// the exponential worst-case may be slow.
pub fn hosoya_index(mol: &Molecule) -> u64 {
    let heavy = heavy_indices(mol);
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();
    let n = heavy.len();
    if n == 0 {
        return 1; // empty graph has one matching (the empty one)
    }
    // The Hosoya index grows as the Fibonacci sequence for path graphs (Fib(n+1)).
    // Fib(50) ≈ 1.3×10¹⁰ calls — this already takes seconds on modern hardware.
    // Cap at 40 heavy atoms to avoid unbounded CPU spin on large/pathological inputs.
    if n > 40 {
        return 0; // sentinel: too large to compute efficiently
    }
    let pos_of: FxHashMap<usize, usize> = heavy.iter().enumerate().map(|(i, &h)| (h, i)).collect();
    let mut adj = vec![vec![false; n]; n];
    for (_, bond) in mol.bonds() {
        let a = bond.atom1.0 as usize;
        let b = bond.atom2.0 as usize;
        if let (Some(&pa), Some(&pb)) = (pos_of.get(&a), pos_of.get(&b))
            && heavy_set.contains(&a)
            && heavy_set.contains(&b)
        {
            adj[pa][pb] = true;
            adj[pb][pa] = true;
        }
    }
    let mut available = vec![true; n];
    count_matchings_hosoya(&adj, &mut available, n)
}

fn count_matchings_hosoya(adj: &[Vec<bool>], available: &mut Vec<bool>, n: usize) -> u64 {
    // Find first available vertex.
    let v = match available.iter().position(|&a| a) {
        Some(v) => v,
        None => return 1, // empty graph → exactly one matching (the empty matching)
    };
    // Case 1: v is left unmatched.
    available[v] = false;
    let mut z = count_matchings_hosoya(adj, available, n);
    // Case 2: v is matched to each available neighbor u.
    for u in 0..n {
        if adj[v][u] && available[u] {
            available[u] = false;
            z += count_matchings_hosoya(adj, available, n);
            available[u] = true;
        }
    }
    available[v] = true;
    z
}

// ─── Topological distance matrix ─────────────────────────────────────────────

/// Topological distance matrix for heavy atoms.
///
/// Entry `[i][j]` is the length of the shortest path (in bonds) between
/// heavy atom `i` and heavy atom `j`.  Diagonal entries are 0.
/// Disconnected atoms get `u32::MAX`.
///
/// The row/column index matches the atom's position in the heavy-atom list
/// (atoms sorted by their original `AtomIdx`).
pub fn topological_distance_matrix(mol: &Molecule) -> Vec<Vec<u32>> {
    let heavy = heavy_indices(mol);
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();
    let n = heavy.len();
    // Map original index → heavy-atom position.
    let mut pos_of: FxHashMap<usize, usize> = FxHashMap::default();
    for (p, &h) in heavy.iter().enumerate() {
        pos_of.insert(h, p);
    }
    let mut matrix = vec![vec![u32::MAX; n]; n];
    for p in 0..n {
        matrix[p][p] = 0;
        let row = bfs_from(mol, heavy[p], &heavy_set);
        for q in 0..n {
            let d = row[heavy[q]];
            if d != usize::MAX {
                matrix[p][q] = d as u32;
            }
        }
    }
    matrix
}

pub fn labute_asa(mol: &Molecule) -> f64 {
    let (per_atom, h_pool_area) = labute_asa_parts(mol);
    per_atom.iter().sum::<f64>() + h_pool_area
}

// ─── VABC — van der Waals atomic bonded-contribution volume ──────────────────

/// Bondi van der Waals radius (Å).
fn r_vdw_bondi(z: u8) -> f64 {
    match z {
        1 => 1.20,
        5 => 1.80,
        6 => 1.70,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        14 => 2.10,
        15 => 1.80,
        16 => 1.80,
        17 => 1.75,
        33 => 1.85,
        34 => 1.90,
        35 => 1.85,
        53 => 1.98,
        _ => 2.00,
    }
}

/// Intersection volume of two spheres with radii r1, r2 and center distance d (Å³).
fn sphere_intersection(r1: f64, r2: f64, d: f64) -> f64 {
    if d <= 0.0 || d >= r1 + r2 {
        return 0.0;
    }
    if d <= (r1 - r2).abs() {
        let r_min = r1.min(r2);
        return 4.0 / 3.0 * PI * r_min * r_min * r_min;
    }
    let h1 = r1 - (d * d + r1 * r1 - r2 * r2) / (2.0 * d);
    let h2 = r2 - (d * d + r2 * r2 - r1 * r1) / (2.0 * d);
    let cap = |r: f64, h: f64| {
        if h <= 0.0 {
            0.0
        } else {
            PI / 3.0 * h * h * (3.0 * r - h)
        }
    };
    cap(r1, h1) + cap(r2, h2)
}

/// VABC — van der Waals atomic bonded-contribution volume approximation (Å³).
///
/// Estimates the molecular van der Waals volume using Bondi radii for each atom
/// and subtracting spherical-cap overlap volumes for each bond (heavy–heavy and
/// heavy–implicit-H).  Bond lengths are estimated from covalent (Rb0) radii.
/// Does not require 3D coordinates.
///
/// Ref: Zhao, Y. H. et al., *J. Org. Chem.* **2003**, *68*, 7368–7373.
pub fn vabc(mol: &Molecule) -> f64 {
    let n = mol.atom_count();
    if n == 0 {
        return 0.0;
    }
    const R_H_VDW: f64 = 1.20;
    const R_H_COV: f64 = 0.33;
    let sphere = |r: f64| 4.0 / 3.0 * PI * r * r * r;

    let mut v = 0.0f64;

    for (idx, atom) in mol.atoms() {
        let z = atom.element.atomic_number();
        v += sphere(r_vdw_bondi(z));
        let nh = implicit_hcount(mol, idx) as usize;
        v += nh as f64 * sphere(R_H_VDW);
    }

    for (_, bond) in mol.bonds() {
        let z1 = mol.atom(bond.atom1).element.atomic_number();
        let z2 = mol.atom(bond.atom2).element.atomic_number();
        let rb1 = rb0(z1);
        let rb2 = rb0(z2);
        if rb1 > 1e-10 && rb2 > 1e-10 {
            v -= sphere_intersection(r_vdw_bondi(z1), r_vdw_bondi(z2), rb1 + rb2);
        }
    }

    for (idx, atom) in mol.atoms() {
        let z = atom.element.atomic_number();
        let rb_heavy = rb0(z);
        if rb_heavy < 1e-10 {
            continue;
        }
        let ov = sphere_intersection(r_vdw_bondi(z), R_H_VDW, rb_heavy + R_H_COV);
        let nh = implicit_hcount(mol, idx) as usize;
        v -= nh as f64 * ov;
    }

    v.max(0.0)
}

// ─── Gravitational index ─────────────────────────────────────────────────────

fn avg_mass_for_grav(z: u8) -> f64 {
    match z {
        1 => 1.008,
        6 => 12.011,
        7 => 14.007,
        8 => 15.999,
        9 => 18.998,
        14 => 28.086,
        15 => 30.974,
        16 => 32.065,
        17 => 35.453,
        35 => 79.904,
        53 => 126.904,
        n => n as f64,
    }
}

/// Gravitational topological index.
///
/// G = Σ_{i<j} mᵢ × mⱼ / dᵢⱼ²
///
/// where mᵢ is the average atomic mass and dᵢⱼ is the topological (bond-graph)
/// distance between heavy atoms i and j.
///
/// Ref: Karelson, M. *Molecular Descriptors in QSAR/QSPR*; Wiley, 2000.
pub fn gravitational_index(mol: &Molecule) -> f64 {
    let heavy = heavy_indices(mol);
    let n = heavy.len();
    if n == 0 {
        return 0.0;
    }
    let masses: Vec<f64> = heavy
        .iter()
        .map(|&a| avg_mass_for_grav(mol.atom(AtomIdx(a as u32)).element.atomic_number()))
        .collect();
    let dist = topological_distance_matrix(mol);
    let mut g = 0.0f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let d = dist[i][j];
            if d != u32::MAX && d > 0 {
                g += masses[i] * masses[j] / (d as f64 * d as f64);
            }
        }
    }
    g
}

// ─── Schultz & Gutman MTI ────────────────────────────────────────────────────

/// Schultz molecular topological index (MTI).
///
/// MTI = Σ_{i<j} (δᵢ + δⱼ) × dᵢⱼ
///
/// where δᵢ is the heavy-atom degree of atom i and dᵢⱼ is the topological distance.
///
/// Ref: Schultz, H. P. *J. Chem. Inf. Comput. Sci.* **1989**, *29*, 227–228.
pub fn schultz_mti(mol: &Molecule) -> u64 {
    let heavy = heavy_indices(mol);
    let n = heavy.len();
    if n == 0 {
        return 0;
    }
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();
    let deg: Vec<u64> = heavy
        .iter()
        .map(|&a| {
            mol.neighbors(AtomIdx(a as u32))
                .filter(|(nb, _)| heavy_set.contains(&(nb.0 as usize)))
                .count() as u64
        })
        .collect();
    let dist = topological_distance_matrix(mol);
    let mut s = 0u64;
    for i in 0..n {
        for j in (i + 1)..n {
            let d = dist[i][j];
            if d < u32::MAX {
                s += (deg[i] + deg[j]) * d as u64;
            }
        }
    }
    s
}

/// Gutman molecular topological index (MTI*).
///
/// MTI* = Σ_{i<j} δᵢ × δⱼ × dᵢⱼ
///
/// where δᵢ is the heavy-atom degree and dᵢⱼ is the topological distance.
///
/// Ref: Gutman, I. *J. Serb. Chem. Soc.* **1994**, *59*, 619–626.
pub fn gutman_mti(mol: &Molecule) -> u64 {
    let heavy = heavy_indices(mol);
    let n = heavy.len();
    if n == 0 {
        return 0;
    }
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();
    let deg: Vec<u64> = heavy
        .iter()
        .map(|&a| {
            mol.neighbors(AtomIdx(a as u32))
                .filter(|(nb, _)| heavy_set.contains(&(nb.0 as usize)))
                .count() as u64
        })
        .collect();
    let dist = topological_distance_matrix(mol);
    let mut s = 0u64;
    for i in 0..n {
        for j in (i + 1)..n {
            let d = dist[i][j];
            if d < u32::MAX {
                s += deg[i] * deg[j] * d as u64;
            }
        }
    }
    s
}

/// Total number of valence electrons (heavy atoms + implicit H).
///
/// Equivalent to RDKit `NumValenceElectrons`.
pub fn num_valence_electrons(mol: &Molecule) -> u32 {
    mol.atoms()
        .filter(|(_, a)| a.element.atomic_number() > 1)
        .map(|(idx, a)| {
            let heavy = u32::from(valence_electrons(a.element.atomic_number()));
            let h = u32::from(implicit_hcount(mol, idx));
            heavy + h // H contributes 1 valence electron each
        })
        .sum()
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap_or_else(|e| panic!("parse '{s}': {e}"))
    }

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    // ── Wiener Index ──────────────────────────────────────────────────────────

    #[test]
    fn wiener_ethane() {
        // CC: only 1 pair at distance 1 → W = 1
        assert_eq!(wiener_index(&mol("CC")) as u32, 1);
    }

    #[test]
    fn wiener_propane() {
        // CCC: pairs (C1,C2)=1, (C1,C3)=2, (C2,C3)=1 → W = 4
        assert_eq!(wiener_index(&mol("CCC")) as u32, 4);
    }

    #[test]
    fn wiener_benzene() {
        // c1ccccc1: 6 atoms in ring, W = 1+2+3+2+1 + 1+2+3+2 + 1+2+3 + 1+2 + 1 = 27
        assert_eq!(wiener_index(&mol("c1ccccc1")) as u32, 27);
    }

    #[test]
    fn wiener_increases_with_chain_length() {
        let w2 = wiener_index(&mol("CC")); // ethane
        let w3 = wiener_index(&mol("CCC")); // propane
        let w4 = wiener_index(&mol("CCCC")); // butane
        assert!(w2 < w3 && w3 < w4);
    }

    #[test]
    fn wiener_single_atom_zero() {
        assert_eq!(wiener_index(&mol("C")) as u32, 0);
    }

    // ── Kappa Indices ─────────────────────────────────────────────────────────

    #[test]
    fn kappa1_propane() {
        // n=3, p1=2 (C-C, C-C): κ1 = 3·4/4 = 3.0
        assert!(close(kappa1(&mol("CCC")), 3.0, 0.01));
    }

    #[test]
    fn kappa1_benzene() {
        // Aromatic C: alpha = 6·(0.67/0.77 - 1) = -0.78 ≠ 0 (RDKit CalcKappa1
        // verified: 3.4116), unlike sp3-only alkanes where alpha = 0.
        assert!(close(kappa1(&mol("c1ccccc1")), 3.4116, 0.001));
    }

    #[test]
    fn kappa2_propane() {
        // n=3, p2=1: κ2 = 2·1/1 = 2.0
        assert!(close(kappa2(&mol("CCC")), 2.0, 0.01));
    }

    #[test]
    fn kappa2_benzene() {
        // RDKit CalcKappa2(benzene) verified: 1.6058 (alpha-corrected).
        assert!(close(kappa2(&mol("c1ccccc1")), 1.6058, 0.001));
    }

    #[test]
    fn kappa3_propane_zero() {
        // n=3, no 3-bond paths: κ3 = 0
        assert_eq!(kappa3(&mol("CCC")), 0.0);
    }

    #[test]
    fn kappa3_benzene() {
        // RDKit CalcKappa3(benzene) verified: 0.5824 (alpha-corrected).
        assert!(close(kappa3(&mol("c1ccccc1")), 0.5824, 0.001));
    }

    #[test]
    fn kappa1_single_atom_zero() {
        assert_eq!(kappa1(&mol("C")), 0.0);
    }

    #[test]
    fn kappa_alpha_corrected_matches_rdkit_aspirin() {
        // RDKit CalcKappa1/2/3("CC(=O)Oc1ccccc1C(=O)O") verified: 9.2496 /
        // 3.7093 / 2.2974. Aspirin mixes sp3/sp2 C, O, and aromatic C, so a
        // nonzero Hall-Kier alpha correction is exercised on all three.
        // Tolerance is 0.1 (not 0.001 like the benzene/aspirin-free cases)
        // because `hall_kier_alpha`'s covalent-radius table has its own
        // separately-tracked precision gap for O (chematic's alpha for
        // aspirin is -1.766 vs RDKit's -1.840) — see
        // tasks/descriptor_validation_coverage.md. This test only guards the
        // kappa formula's alpha-wiring, not that residual table precision.
        let m = mol("CC(=O)Oc1ccccc1C(=O)O");
        assert!(close(kappa1(&m), 9.2496, 0.1), "kappa1 = {}", kappa1(&m));
        assert!(close(kappa2(&m), 3.7093, 0.1), "kappa2 = {}", kappa2(&m));
        assert!(close(kappa3(&m), 2.2974, 0.1), "kappa3 = {}", kappa3(&m));
    }

    #[test]
    fn kappa_all_matches_individual() {
        for smi in ["CC", "CCC", "CCCC", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"] {
            let m = mol(smi);
            let (k1, k2, k3) = kappa_all(&m);
            assert!((k1 - kappa1(&m)).abs() < 1e-10, "{smi}: kappa1 mismatch");
            assert!((k2 - kappa2(&m)).abs() < 1e-10, "{smi}: kappa2 mismatch");
            assert!((k3 - kappa3(&m)).abs() < 1e-10, "{smi}: kappa3 mismatch");
        }
    }

    // ── Chi Connectivity ─────────────────────────────────────────────────────

    #[test]
    fn chi0_benzene() {
        // 6 aromatic C each with δ=2: χ0 = 6 · 2^(-0.5) ≈ 4.243
        let c = chi0(&mol("c1ccccc1"));
        assert!(close(c, 4.243, 0.01), "chi0(benzene) = {c}");
    }

    #[test]
    fn chi1_benzene() {
        // 6 bonds, each (2·2)^(-0.5) = 0.5: χ1 = 3.0
        let c = chi1(&mol("c1ccccc1"));
        assert!(close(c, 3.0, 0.01), "chi1(benzene) = {c}");
    }

    #[test]
    fn chi0_propane() {
        // δ(C1)=1, δ(C2)=2, δ(C3)=1: χ0 = 1 + 1/√2 + 1 ≈ 2.707
        let c = chi0(&mol("CCC"));
        assert!(close(c, 2.707, 0.01), "chi0(propane) = {c}");
    }

    #[test]
    fn chi1_propane() {
        // bonds (1·2) and (2·1): χ1 = 2·(2)^(-0.5) ≈ 1.414
        let c = chi1(&mol("CCC"));
        assert!(close(c, 1.414, 0.01), "chi1(propane) = {c}");
    }

    #[test]
    fn chi_increases_with_chain() {
        // Longer chain → larger χ0
        assert!(chi0(&mol("CCC")) < chi0(&mol("CCCC")));
    }

    #[test]
    fn chi0v_benzene() {
        // Each aromatic C: Zv=4, H=1, denom=1 → δv=3: χ0v = 6/√3 ≈ 3.464
        let c = chi0v(&mol("c1ccccc1"));
        assert!(close(c, 3.464, 0.01), "chi0v(benzene) = {c}");
    }

    #[test]
    fn chi1v_benzene() {
        // bonds (3·3): χ1v = 6·(9)^(-0.5) = 6/3 = 2.0
        let c = chi1v(&mol("c1ccccc1"));
        assert!(close(c, 2.0, 0.01), "chi1v(benzene) = {c}");
    }

    // ── Bertz CT ──────────────────────────────────────────────────────────────

    #[test]
    fn bertz_ct_increases_with_complexity() {
        let bz = bertz_ct(&mol("c1ccccc1"));
        let asp = bertz_ct(&mol("CC(=O)Oc1ccccc1C(=O)O"));
        assert!(bz < asp, "benzene BertzCT {bz} should be < aspirin {asp}");
    }

    #[test]
    fn bertz_ct_ethane_less_than_propane() {
        assert!(bertz_ct(&mol("CC")) < bertz_ct(&mol("CCC")));
    }

    #[test]
    fn bertz_ct_methane() {
        // C: deg=0, h=4, total=4, C(4,2)=6; m=4 H bonds → CT = 6+4 = 10
        assert!(close(bertz_ct(&mol("C")), 10.0, 0.01));
    }

    #[test]
    fn bertz_ct_benzene() {
        // 6 C with total_deg=3: 6·C(3,2)=18; bonds=6+6=12 → CT = 18+12 = 30
        assert!(close(bertz_ct(&mol("c1ccccc1")), 30.0, 0.01));
    }

    // ── LabuteASA ─────────────────────────────────────────────────────────────

    #[test]
    fn labute_asa_positive() {
        assert!(labute_asa(&mol("c1ccccc1")) > 0.0, "benzene ASA > 0");
    }

    #[test]
    fn labute_asa_single_oxygen() {
        // RDKit CalcLabuteASA("O") verified: 6.8492. The implicit-H
        // contribution is a single pooled term per heavy atom (not scaled by
        // the atom's actual H count) — see `labute_asa_parts` doc comment.
        let asa = labute_asa(&mol("O"));
        assert!(
            (asa - 6.8492).abs() < 0.001,
            "water ASA {asa:.4} != RDKit 6.8492"
        );
    }

    #[test]
    fn labute_asa_matches_rdkit_aspirin() {
        // RDKit CalcLabuteASA("CC(=O)Oc1ccccc1C(=O)O") verified: 74.7571.
        let asa = labute_asa(&mol("CC(=O)Oc1ccccc1C(=O)O"));
        assert!((asa - 74.7571).abs() < 0.001, "aspirin ASA {asa:.4}");
    }

    #[test]
    fn labute_asa_per_atom_excludes_pooled_h_term() {
        // Quaternary carbon (CC(C)(C)C) has an atom with zero implicit H —
        // RDKit still runs the pooled-H pass once for it (see
        // `labute_asa_parts`). RDKit `_CalcLabuteASAContribs` ats verified:
        // [6.9237, 5.4150, 6.9237, 6.9237, 6.9237], hs=1.0891.
        let m = mol("CC(C)(C)C");
        let per_atom = labute_asa_per_atom(&m);
        let expected = [6.9237, 5.4150, 6.9237, 6.9237, 6.9237];
        for (got, want) in per_atom.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 0.001, "{per_atom:?} vs {expected:?}");
        }
        let total = labute_asa(&m);
        assert!((total - (expected.iter().sum::<f64>() + 1.0891)).abs() < 0.001);
    }

    #[test]
    fn labute_asa_monotone_with_size() {
        // Larger molecule → larger ASA.
        let bz = labute_asa(&mol("c1ccccc1"));
        let asp = labute_asa(&mol("CC(=O)Oc1ccccc1C(=O)O"));
        assert!(
            bz < asp,
            "benzene ASA {bz:.2} should be < aspirin ASA {asp:.2}"
        );
    }

    #[test]
    fn labute_asa_aromatic_reduces_vs_saturated() {
        // Aromatic C-C bonds (scale=0.1) create more surface overlap than
        // single C-C bonds (scale=0), so benzene ASA < cyclohexane ASA
        // (same atom count, but less overlap in cyclohexane).
        let bz = labute_asa(&mol("c1ccccc1")); // aromatic
        let ch = labute_asa(&mol("C1CCCCC1")); // saturated
        assert!(bz < ch, "benzene ASA {bz:.2} < cyclohexane ASA {ch:.2}");
    }

    #[test]
    fn randic_index_ethane() {
        // Ethane: 1 edge between two degree-1 nodes → Randic = 1/sqrt(1*1) = 1.0
        let m = mol("CC");
        assert!((randic_index(&m) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn zagreb_m1_ethane() {
        // Ethane: 2 atoms each with degree 1 → Σ d² = 1+1 = 2
        let m = mol("CC");
        assert_eq!(zagreb_index_m1(&m), 2);
    }

    // ── Zagreb M2 ───────────────────────────────────────────────────────────

    #[test]
    fn zagreb_m2_ethane() {
        // Ethane C-C: one bond, deg(C)=1, deg(C)=1 → M2 = 1*1 = 1
        assert_eq!(zagreb_index_m2(&mol("CC")), 1);
    }

    #[test]
    fn zagreb_m2_propane() {
        // Propane C-C-C: bonds (C1-C2) and (C2-C3)
        // deg: C1=1, C2=2, C3=1
        // M2 = 1*2 + 2*1 = 4
        assert_eq!(zagreb_index_m2(&mol("CCC")), 4);
    }

    #[test]
    fn zagreb_m2_benzene() {
        // Benzene: 6 bonds, each between degree-2 atoms → M2 = 6 * (2*2) = 24
        assert_eq!(zagreb_index_m2(&mol("c1ccccc1")), 24);
    }

    #[test]
    fn zagreb_m2_ge_m1_for_branched() {
        // For any graph, M2 captures edge-based branching; M1 is vertex-based.
        // Both should be positive for any non-trivial molecule.
        for smi in ["CC(C)C", "CC(=O)O", "c1ccccc1"] {
            let m = mol(smi);
            assert!(zagreb_index_m2(&m) > 0, "M2 should be > 0 for {smi}");
        }
    }

    // ── Eccentricity / Petitjean / ECI ───────────────────────────────────────

    #[test]
    fn eccentricities_propane() {
        // CCC: terminal Cs have ecc=2, middle C has ecc=1
        let m = mol("CCC");
        let ecc = graph_eccentricities(&m);
        assert_eq!(ecc.len(), 3);
        assert_eq!(*ecc.iter().max().unwrap(), 2);
        assert_eq!(*ecc.iter().min().unwrap(), 1);
    }

    #[test]
    fn graph_diameter_propane() {
        assert_eq!(graph_diameter(&mol("CCC")), 2);
    }

    #[test]
    fn graph_radius_propane() {
        assert_eq!(graph_radius(&mol("CCC")), 1);
    }

    #[test]
    fn petitjean_propane() {
        // (2 - 1) / 2 = 0.5
        let v = petitjean_index(&mol("CCC"));
        assert!((v - 0.5).abs() < 1e-9, "expected 0.5 got {v}");
    }

    #[test]
    fn petitjean_benzene() {
        // c1ccccc1: 6-cycle — all atoms have eccentricity 3, so diameter=radius=3 → Petitjean=0
        let v = petitjean_index(&mol("c1ccccc1"));
        assert!(
            (v - 0.0).abs() < 1e-9,
            "expected 0.0 (symmetric ring) got {v}"
        );
    }

    #[test]
    fn petitjean_single_atom() {
        assert_eq!(petitjean_index(&mol("C")), 0.0);
    }

    #[test]
    fn eci_propane() {
        // CCC: middle C deg=2 ecc=1 → 2; each terminal C deg=1 ecc=2 → 2 each → total 6
        assert_eq!(eccentric_connectivity_index(&mol("CCC")), 6);
    }

    // ── Hosoya Index ─────────────────────────────────────────────────────────

    #[test]
    fn hosoya_methane() {
        // Single atom: 0 edges → Z = 1 (empty matching only)
        assert_eq!(hosoya_index(&mol("C")), 1);
    }

    #[test]
    fn hosoya_ethane() {
        // CC: 1 edge → Z = 2 (empty + {e1})
        assert_eq!(hosoya_index(&mol("CC")), 2);
    }

    #[test]
    fn hosoya_propane() {
        // CCC: 2 edges, no two share a vertex → Z = 1 + 2 = 3
        assert_eq!(hosoya_index(&mol("CCC")), 3);
    }

    #[test]
    fn hosoya_butane() {
        // CCCC: 3 edges; matchings: empty(1) + each edge(3) + {e1,e3}(1) = 5
        assert_eq!(hosoya_index(&mol("CCCC")), 5);
    }

    #[test]
    fn hosoya_benzene() {
        // c1ccccc1: C6 cycle — Z(C6) = M0+M1+M2+M3 = 1+6+9+2 = 18
        // M2 = C(6,2) - 6 adjacent pairs = 15 - 6 = 9; M3 = 2 perfect matchings
        assert_eq!(hosoya_index(&mol("c1ccccc1")), 18);
    }

    #[test]
    fn distance_matrix_ethane() {
        let m = mol("CC");
        let dm = topological_distance_matrix(&m);
        assert_eq!(dm.len(), 2);
        assert_eq!(dm[0][0], 0);
        assert_eq!(dm[0][1], 1);
        assert_eq!(dm[1][0], 1);
        assert_eq!(dm[1][1], 0);
    }

    #[test]
    fn distance_matrix_propane() {
        // C-C-C: d(0,1)=1, d(0,2)=2, d(1,2)=1
        let m = mol("CCC");
        let dm = topological_distance_matrix(&m);
        assert_eq!(dm[0][2], 2);
        assert_eq!(dm[1][2], 1);
    }

    // ── Padmakar-Ivan (PI) Index ─────────────────────────────────────────────

    #[test]
    fn pi_single_atom() {
        assert_eq!(padmakar_ivan_index(&mol("C")), 0);
    }

    #[test]
    fn pi_ethane() {
        // CC: 1 bond, n_u=1 (C1), n_v=1 (C2) → PI = 2
        assert_eq!(padmakar_ivan_index(&mol("CC")), 2);
    }

    #[test]
    fn pi_propane() {
        // CCC: edge(1-2) → n_u=1, n_v=2 → 3; edge(2-3) → n_u=2, n_v=1 → 3; PI = 6
        assert_eq!(padmakar_ivan_index(&mol("CCC")), 6);
    }

    #[test]
    fn pi_butane() {
        // CCCC: each of 3 edges contributes 4 → PI = 12
        assert_eq!(padmakar_ivan_index(&mol("CCCC")), 12);
    }

    #[test]
    fn pi_benzene() {
        // c1ccccc1: 6-ring, each of 6 edges contributes 6 → PI = 36
        assert_eq!(padmakar_ivan_index(&mol("c1ccccc1")), 36);
    }

    #[test]
    fn pi_linear_chain_formula() {
        // For a linear chain of n heavy atoms: PI = n(n-1)
        for n in 2..=6 {
            let smiles = "C".repeat(n);
            let expected = (n as u64) * (n as u64 - 1);
            assert_eq!(padmakar_ivan_index(&mol(&smiles)), expected, "chain n={n}");
        }
    }

    // ── Schultz / Gutman MTI ─────────────────────────────────────────────────

    #[test]
    fn schultz_mti_ethane() {
        // CC: 2 atoms, deg=[1,1], d=1 → MTI = (1+1)*1 = 2
        assert_eq!(schultz_mti(&mol("CC")), 2);
    }

    #[test]
    fn gutman_mti_ethane() {
        // CC: 2 atoms, deg=[1,1], d=1 → MTI* = 1*1*1 = 1
        assert_eq!(gutman_mti(&mol("CC")), 1);
    }

    #[test]
    fn schultz_mti_propane() {
        // CCC: atoms 0,1,2 with deg=[1,2,1], distances: d(0,1)=1, d(0,2)=2, d(1,2)=1
        // (1+2)*1 + (1+1)*2 + (2+1)*1 = 3 + 4 + 3 = 10
        assert_eq!(schultz_mti(&mol("CCC")), 10);
    }

    #[test]
    fn gutman_mti_propane() {
        // CCC: atoms 0,1,2 with deg=[1,2,1], distances: d(0,1)=1, d(0,2)=2, d(1,2)=1
        // 1*2*1 + 1*1*2 + 2*1*1 = 2 + 2 + 2 = 6
        assert_eq!(gutman_mti(&mol("CCC")), 6);
    }

    #[test]
    fn schultz_mti_empty() {
        // edge case: single atom
        let m = chematic_smiles::parse("[C]").unwrap();
        assert_eq!(schultz_mti(&m), 0);
    }

    // ── VABC ─────────────────────────────────────────────────────────────────

    #[test]
    fn vabc_methane() {
        // Single carbon + 4 implicit H — result should be positive and roughly 30–50 Å³
        let v = vabc(&mol("C"));
        assert!(v > 10.0, "VABC(methane) should be > 10 Å³, got {v}");
        assert!(v < 100.0, "VABC(methane) should be < 100 Å³, got {v}");
    }

    #[test]
    fn vabc_water() {
        let v = vabc(&mol("O"));
        assert!(v > 5.0, "VABC(water) should be > 5 Å³, got {v}");
        assert!(v < 50.0, "VABC(water) should be < 50 Å³, got {v}");
    }

    #[test]
    fn vabc_increases_with_size() {
        // Larger molecule should have larger VABC
        let v_ethane = vabc(&mol("CC"));
        let v_propane = vabc(&mol("CCC"));
        let v_butane = vabc(&mol("CCCC"));
        assert!(v_propane > v_ethane, "propane > ethane");
        assert!(v_butane > v_propane, "butane > propane");
    }

    // ── Gravitational index ──────────────────────────────────────────────────

    #[test]
    fn gravitational_index_single_atom() {
        let m = chematic_smiles::parse("[C]").unwrap();
        assert_eq!(gravitational_index(&m), 0.0);
    }

    #[test]
    fn gravitational_index_ethane() {
        // CC: 2 C, mass=12.011 each, d=1 → G = 12.011*12.011/1 ≈ 144.26
        let g = gravitational_index(&mol("CC"));
        assert!((g - 12.011_f64 * 12.011).abs() < 0.01, "got {g}");
    }

    #[test]
    fn gravitational_index_positive() {
        let g = gravitational_index(&mol("c1ccccc1"));
        assert!(g > 0.0, "benzene gravitational index should be positive");
    }

    #[test]
    fn chi_all_matches_individual() {
        // Verify chi_all returns identical values to individual chi0-chi4/v functions.
        for smi in [
            "CC",
            "CCC",
            "c1ccccc1",
            "CC(=O)Oc1ccccc1C(=O)O",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
        ] {
            let m = mol(smi);
            let (c0, c1, c2, c3, c4, c0v, c1v, c2v, c3v, c4v) = chi_all(&m);
            assert!((c0 - chi0(&m)).abs() < 1e-10, "{smi}: chi0 mismatch");
            assert!((c1 - chi1(&m)).abs() < 1e-10, "{smi}: chi1 mismatch");
            assert!((c2 - chi2(&m)).abs() < 1e-10, "{smi}: chi2 mismatch");
            assert!((c3 - chi3(&m)).abs() < 1e-10, "{smi}: chi3 mismatch");
            assert!((c4 - chi4(&m)).abs() < 1e-10, "{smi}: chi4 mismatch");
            assert!((c0v - chi0v(&m)).abs() < 1e-10, "{smi}: chi0v mismatch");
            assert!((c1v - chi1v(&m)).abs() < 1e-10, "{smi}: chi1v mismatch");
            assert!((c2v - chi2v(&m)).abs() < 1e-10, "{smi}: chi2v mismatch");
            assert!((c3v - chi3v(&m)).abs() < 1e-10, "{smi}: chi3v mismatch");
            assert!((c4v - chi4v(&m)).abs() < 1e-10, "{smi}: chi4v mismatch");
        }
    }
}
