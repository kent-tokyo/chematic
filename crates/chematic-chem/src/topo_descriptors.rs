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
/// atoms visited so far (including `cur`). `delta[i]` is atom `i`'s delta,
/// 0 for atoms outside the heavy set, so atoms with no positive delta
/// (including every non-heavy atom) are never entered.
fn chi_dfs(
    mol: &Molecule,
    cur: usize,
    target_len: usize,
    cur_len: usize,
    running_product: f64,
    visited: &mut [bool],
    delta: &[f64],
) -> f64 {
    if cur_len == target_len {
        return running_product.powf(-0.5);
    }
    let mut sum = 0.0f64;
    for &(nb, _) in mol.neighbor_slice(AtomIdx(cur as u32)) {
        let ni = nb.0 as usize;
        if !visited[ni] {
            let d_nb = delta[ni];
            if d_nb > 0.0 {
                visited[ni] = true;
                sum += chi_dfs(
                    mol,
                    ni,
                    target_len,
                    cur_len + 1,
                    running_product * d_nb,
                    visited,
                    delta,
                );
                visited[ni] = false;
            }
        }
    }
    sum
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
    // Per-atom deltas and heavy flags, computed once; the DFS visits atoms
    // and sums contributions in exactly the historical order.
    let atom_count = mol.atom_count();
    let mut delta_of = vec![0.0f64; atom_count];
    for &h in heavy {
        delta_of[h] = if use_valence {
            delta_v(mol, AtomIdx(h as u32))
        } else {
            delta(mol, AtomIdx(h as u32), heavy_set)
        };
    }
    let mut visited = vec![false; atom_count];
    let mut total = 0.0f64;
    for &start in heavy {
        let d_start = delta_of[start];
        if d_start <= 0.0 {
            continue;
        }
        visited[start] = true;
        total += chi_dfs(mol, start, n, 0, d_start, &mut visited, &delta_of);
        visited[start] = false;
    }
    total / 2.0
}

fn chi_n(mol: &Molecule, n: usize, use_valence: bool) -> f64 {
    let heavy = heavy_indices(mol);
    // The heavy-atom set only feeds the plain (non-valence) delta.
    let heavy_set: FxHashSet<usize> = if use_valence {
        FxHashSet::default()
    } else {
        heavy.iter().copied().collect()
    };
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
    wiener_index_with_topology(mol, &heavy, &heavy_set)
}

fn wiener_index_with_topology(
    mol: &Molecule,
    heavy: &[usize],
    heavy_set: &FxHashSet<usize>,
) -> f64 {
    let mut sum = 0u64;
    for i in 0..heavy.len() {
        let row = bfs_from(mol, heavy[i], heavy_set);
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

/// RDKit's kappa helpers (`ConnectivityDescriptors.cpp`), operation for
/// operation: `p` is the path count, `a` the heavy-atom count, and a zero
/// denominator gives 0.
fn rdkit_kappa1(p1: f64, a: f64, alpha: f64) -> f64 {
    let denom = p1 + alpha;
    if denom != 0.0 {
        (a + alpha) * (a + alpha - 1.0) * (a + alpha - 1.0) / (denom * denom)
    } else {
        0.0
    }
}

fn rdkit_kappa2(p2: f64, a: f64, alpha: f64) -> f64 {
    let denom = (p2 + alpha) * (p2 + alpha);
    if denom != 0.0 {
        (a + alpha - 1.0) * (a + alpha - 2.0) * (a + alpha - 2.0) / denom
    } else {
        0.0
    }
}

fn rdkit_kappa3(p3: f64, a: usize, alpha: f64) -> f64 {
    let denom = (p3 + alpha) * (p3 + alpha);
    if denom == 0.0 {
        return 0.0;
    }
    let a_f = a as f64;
    if a % 2 == 1 {
        (a_f + alpha - 1.0) * (a_f + alpha - 3.0) * (a_f + alpha - 3.0) / denom
    } else {
        (a_f + alpha - 2.0) * (a_f + alpha - 3.0) * (a_f + alpha - 3.0) / denom
    }
}

/// Atoms RDKit's `getNumHeavyAtoms` counts (atomic number > 1).
fn rdkit_heavy_atom_count(mol: &Molecule) -> usize {
    mol.atoms()
        .filter(|(_, a)| !a.wildcard && a.element.atomic_number() > 1)
        .count()
}

/// Hall-Kier κ1 shape index, RDKit's `CalcKappa1`:
/// `(A+α)(A+α−1)² / (P1+α)²` with P1 the bond count.
pub fn kappa1(mol: &Molecule) -> f64 {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return kappa1(&h);
    }
    rdkit_kappa1(
        mol.bond_count() as f64,
        rdkit_heavy_atom_count(mol) as f64,
        crate::descriptors::hall_kier_alpha(mol),
    )
}

/// Hall-Kier κ2 shape index, RDKit's `CalcKappa2`:
/// `(A+α−1)(A+α−2)² / (P2+α)²` with P2 the count of 2-bond paths.
pub fn kappa2(mol: &Molecule) -> f64 {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return kappa2(&h);
    }
    rdkit_kappa2(
        rdkit_paths_of_length(mol, 2, true).len() as f64,
        rdkit_heavy_atom_count(mol) as f64,
        crate::descriptors::hall_kier_alpha(mol),
    )
}

/// Hall-Kier κ3 shape index, RDKit's `CalcKappa3` (the numerator depends on
/// the parity of the heavy-atom count).
pub fn kappa3(mol: &Molecule) -> f64 {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return kappa3(&h);
    }
    rdkit_kappa3(
        rdkit_paths_of_length(mol, 3, true).len() as f64,
        rdkit_heavy_atom_count(mol),
        crate::descriptors::hall_kier_alpha(mol),
    )
}

/// κ1, κ2 and κ3 with one Hall-Kier alpha: `(κ1, κ2, κ3)`.
pub fn kappa_all(mol: &Molecule) -> (f64, f64, f64) {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return kappa_all(&h);
    }
    let alpha = crate::descriptors::hall_kier_alpha(mol);
    let a = rdkit_heavy_atom_count(mol);
    (
        rdkit_kappa1(mol.bond_count() as f64, a as f64, alpha),
        rdkit_kappa2(
            rdkit_paths_of_length(mol, 2, true).len() as f64,
            a as f64,
            alpha,
        ),
        rdkit_kappa3(rdkit_paths_of_length(mol, 3, true).len() as f64, a, alpha),
    )
}

// ─── Chi Connectivity Indices ────────────────────────────────────────────────

/// Kier-Hall χ0 connectivity index.
///
/// χ0 = Σᵢ δᵢ^(−0.5) over all heavy atoms, where δᵢ = heavy-atom degree.
/// Atoms with δ = 0 contribute 0.
pub fn chi0(mol: &Molecule) -> f64 {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return chi0(&h);
    }
    rdkit_chi0(mol)
}

/// Kier-Hall χ1 connectivity index (bond-path sum).
pub fn chi1(mol: &Molecule) -> f64 {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return chi1(&h);
    }
    rdkit_chi1(mol)
}

/// Kier-Hall χ2 connectivity index (2-bond path sum).
pub fn chi2(mol: &Molecule) -> f64 {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return chi2(&h);
    }
    chi_n(mol, 2, false)
}

/// Kier-Hall χ3 connectivity index (3-bond path sum).
pub fn chi3(mol: &Molecule) -> f64 {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return chi3(&h);
    }
    chi_n(mol, 3, false)
}

/// Kier-Hall χ4 connectivity index (4-bond path sum).
pub fn chi4(mol: &Molecule) -> f64 {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return chi4(&h);
    }
    chi_n(mol, 4, false)
}

/// Valence-corrected χ0v connectivity index (RDKit `CalcChi0v`).
pub fn chi0v(mol: &Molecule) -> f64 {
    chematic_smiles::with_rdkit_model_molecule(mol, chi0v_impl)
}

fn chi0v_impl(mol: &Molecule) -> f64 {
    rdkit_chi_v_order(mol, 0)
}

/// Valence-corrected χ1v connectivity index (RDKit `CalcChi1v`).
pub fn chi1v(mol: &Molecule) -> f64 {
    chematic_smiles::with_rdkit_model_molecule(mol, chi1v_impl)
}

fn chi1v_impl(mol: &Molecule) -> f64 {
    rdkit_chi_v_order(mol, 1)
}

/// Valence-corrected χ2v connectivity index (RDKit `CalcChi2v`).
pub fn chi2v(mol: &Molecule) -> f64 {
    chematic_smiles::with_rdkit_model_molecule(mol, chi2v_impl)
}

fn chi2v_impl(mol: &Molecule) -> f64 {
    rdkit_chi_v_order(mol, 2)
}

/// Valence-corrected χ3v connectivity index (RDKit `CalcChi3v`).
pub fn chi3v(mol: &Molecule) -> f64 {
    chematic_smiles::with_rdkit_model_molecule(mol, chi3v_impl)
}

fn chi3v_impl(mol: &Molecule) -> f64 {
    rdkit_chi_v_order(mol, 3)
}

/// Valence-corrected χ4v connectivity index (RDKit `CalcChi4v`).
pub fn chi4v(mol: &Molecule) -> f64 {
    chematic_smiles::with_rdkit_model_molecule(mol, chi4v_impl)
}

fn chi4v_impl(mol: &Molecule) -> f64 {
    rdkit_chi_v_order(mol, 4)
}

/// Compute all 10 Hall-Kier connectivity indices in a single pass.
///
/// Returns `(χ0, χ1, χ2, χ3, χ4, χ0v, χ1v, χ2v, χ3v, χ4v)`.
/// Use when all indices are needed to avoid 10 redundant `heavy_indices` computations.
pub fn chi_all(mol: &Molecule) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    chematic_smiles::with_rdkit_model_molecule(mol, chi_all_impl)
}

fn chi_all_impl(mol: &Molecule) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let heavy = heavy_indices(mol);
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();
    chi_all_with_topology(mol, &heavy, &heavy_set)
}

fn chi_all_with_topology(
    mol: &Molecule,
    heavy: &[usize],
    heavy_set: &FxHashSet<usize>,
) -> (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64) {
    let c0 = rdkit_chi0(mol);
    let [c0v, c1v, c2v, c3v, c4v] = rdkit_chi_v(mol);

    (
        c0,
        rdkit_chi1(mol),
        chi_n_with(mol, heavy, heavy_set, 2, false),
        chi_n_with(mol, heavy, heavy_set, 3, false),
        chi_n_with(mol, heavy, heavy_set, 4, false),
        c0v,
        c1v,
        c2v,
        c3v,
        c4v,
    )
}

/// Shared topology-derived descriptor values for one molecule.
///
/// Bulk descriptor callers should compute this once when requesting multiple
/// topology descriptors. The contained values preserve the existing scalar
/// function results while sharing heavy-atom extraction, topology membership,
/// and path-count preparation.
#[derive(Debug, Clone, PartialEq)]
pub struct TopologyBundle {
    pub wiener: f64,
    pub kappa: (f64, f64, f64),
    pub chi: (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64),
}

/// Compute the shared topology descriptor bundle for one molecule.
pub fn topology_bundle(mol: &Molecule) -> TopologyBundle {
    let heavy = heavy_indices(mol);
    let heavy_set: FxHashSet<usize> = heavy.iter().copied().collect();
    TopologyBundle {
        wiener: wiener_index_with_topology(mol, &heavy, &heavy_set),
        kappa: kappa_all(mol),
        chi: chi_all_with_topology(mol, &heavy, &heavy_set),
    }
}

// ─── Bertz Complexity ────────────────────────────────────────────────────────

/// Bertz CT molecular complexity index (Bertz, *J. Am. Chem. Soc.* **103**,
/// 3599–3601, 1981), bit-identical to RDKit's `GraphDescriptors.BertzCT`
/// (default `cutoff=100`; aromatic bonds as order 1.5).
///
/// Returns 0.0 for fewer than two atoms or more than 1000 atoms.
pub fn bertz_ct(mol: &Molecule) -> f64 {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return bertz_ct(&h);
    }
    crate::rdkit_graph::bertz_ct(mol)
}

// ─── Labute ASA ─────────────────────────────────────────────────────────────

/// Covalent (Rb0) radius of an element (Å), from RDKit's ptable.GetRb0.
/// Returns 0.0 for unrecognized elements (they contribute no surface area).
/// Bond radii of the van der Waals volume estimate (unchanged; this
/// descriptor has no RDKit counterpart).
fn vdw_volume_rb0(atomic_number: u8) -> f64 {
    match atomic_number {
        1 => 0.33,
        6 => 0.77,
        7 => 0.70,
        8 => 0.66,
        9 => 0.611,
        14 => 1.04,
        15 => 0.89,
        16 => 1.04,
        17 => 0.997,
        33 => 1.21,
        34 => 1.20,
        35 => 1.167,
        53 => 1.387,
        _ => 0.0,
    }
}

fn rb0(atom: &chematic_core::Atom) -> f64 {
    // RDKit's getRb0; a dummy atom (atomic number 0) has none.
    if atom.wildcard {
        return 0.0;
    }
    crate::descriptors::RDKIT_RB0
        .get(atom.element.atomic_number() as usize - 1)
        .copied()
        .unwrap_or(0.0)
}

/// Bond-type scale factor used in the Labute formula (Å subtracted from Ri+Rj).
///
/// Shorter bonds (double, triple, aromatic) bring atoms closer, increasing
/// surface overlap.  Single bonds have scale 0 (spheres just touching, no overlap).
fn bond_scale(view: &Molecule, bond: chematic_core::BondIdx) -> f64 {
    // RDKit: aromatic bonds 0.1; otherwise by bond type for types < 4
    // (single 0, double 0.2, triple 0.3), nothing for any other type.
    let b = view.bond(bond);
    // Aromaticity belongs to the bond type here. In particular, a dative
    // bond between two aromatic atoms remains dative in RDKit and must not
    // receive the 0.1 aromatic overlap adjustment.
    if b.order == BondOrder::Aromatic {
        return 0.1;
    }
    match b.order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => 0.0,
        BondOrder::Double => 0.2,
        BondOrder::Triple => 0.3,
        _ => 0.0,
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
fn labute_asa_core(mol: &Molecule, collect_per_atom: bool) -> (Vec<f64>, f64, f64) {
    // RDKit's getLabuteAtomContribs (MolSurf.cpp) with includeHs, operation
    // for operation; bond aromaticity as RDKit perceives it.
    let n = mol.atom_count();
    if n == 0 {
        return (Vec::new(), 0.0, 0.0);
    }
    let view = crate::descriptors::descriptor_aromaticity(mol);
    let view: &Molecule = &view;
    // Most drug-like molecules fit inline.  Keeping radius and overlap in one
    // buffer avoids two heap allocations and improves locality without
    // changing RDKit's accumulation order.
    let mut atom_data: smallvec::SmallVec<[(f64, f64); 64]> =
        mol.atoms().map(|(_, a)| (rb0(a), 0.0)).collect();
    let mut accumulate_bond = |bidx| {
        let bond = view.bond(bidx);
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        let (ri, rj) = (atom_data[i].0, atom_data[j].0);
        let mut bij = ri + rj;
        bij -= bond_scale(view, bidx);
        let dij = (ri - rj).abs().max(bij).min(ri + rj);
        atom_data[i].1 += rj * rj - (ri - dij) * (ri - dij) / dij;
        atom_data[j].1 += ri * ri - (rj - dij) * (rj - dij) / dij;
    };
    // Bonds in RDKit's numbering (ring closures last), without materializing
    // the complete order vector for every molecule.
    mol.for_each_bond_in_rdkit_order(&mut accumulate_bond);
    // One hydrogen per atom, pooled into a single term (as RDKit does).
    let rh = crate::descriptors::RDKIT_RB0[0];
    let mut h_contrib = 0.0f64;
    for i in 0..n {
        let ri = atom_data[i].0;
        let bij = ri + rh;
        let dij = (ri - rh).abs().max(bij).min(ri + rh);
        atom_data[i].1 += rh * rh - (ri - dij) * (ri - dij) / dij;
        h_contrib += ri * ri - (rh - dij) * (rh - dij) / dij;
    }
    let mut per_atom = if collect_per_atom {
        Vec::with_capacity(n)
    } else {
        Vec::new()
    };
    let mut total = 0.0;
    for i in 0..n {
        let (ri, overlap) = atom_data[i];
        let contribution = PI * ri * (4.0 * ri - overlap);
        total += contribution;
        if collect_per_atom {
            per_atom.push(contribution);
        }
    }
    let h_pool_area = if h_contrib.abs() > 1e-4 {
        PI * rh * (4.0 * rh - h_contrib)
    } else {
        0.0
    };
    (per_atom, h_pool_area, total + h_pool_area)
}

fn labute_asa_parts(mol: &Molecule) -> (Vec<f64>, f64) {
    let (per_atom, h_pool_area, _) = labute_asa_core(mol, true);
    (per_atom, h_pool_area)
}

/// Per-atom Labute approximate surface area contributions (Å²), excluding
/// the pooled implicit-hydrogen term (see [`labute_asa_parts`]). This is the
/// per-atom weight used by the VSA descriptor families
/// ([`crate::vsa`]), matching RDKit's `ats` output.
pub fn labute_asa_per_atom(mol: &Molecule) -> Vec<f64> {
    labute_asa_parts(mol).0
}

// ─── RDKit-exact path enumeration and valence connectivity ───────────────────

/// RDKit's `PeriodicTable::getNouterElecs`, indexed by atomic number - 1.
/// Generated with RDKit 2026.03.1.
const RDKIT_N_OUTER_ELECS: [u8; 118] = [
    1, 2, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3,
    4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 3, 4, 5, 6,
    7, 8, 9, 10, 11, 12, 13, 14, 15, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 3,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
];

/// RDKit's `findAllPathsOfLengthN(mol, target_len, use_bonds, useHs=false)`
/// (`Subgraphs.cpp`), in RDKit's order: atom paths grown one atom at a time
/// from every start atom (neighbours in ascending index), a ring may close
/// on the last step, and paths with the same bond set keep their first
/// occurrence. With `use_bonds` a path of `target_len` bonds is returned as
/// its bond indices, otherwise a path of `target_len` atoms as atoms.
/// Sums over these paths add in RDKit's order, so they are bit-identical.
pub(crate) fn rdkit_paths_of_length(
    mol: &Molecule,
    target_len: usize,
    use_bonds: bool,
) -> Vec<Vec<u32>> {
    let dim = mol.atom_count();
    let upper = if use_bonds {
        target_len + 1
    } else {
        target_len
    };
    if dim == 0 || upper == 0 {
        return Vec::new();
    }
    let is_h = |i: usize| {
        let a = mol.atom(AtomIdx(i as u32));
        !a.wildcard && a.element.atomic_number() == 1
    };
    // Sorted neighbour lists stand in for RDKit's adjacency-matrix scan.
    let mut adj: Vec<Vec<u32>> = vec![Vec::new(); dim];
    for (_, bond) in mol.bonds() {
        let (a, b) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        if a == b || is_h(a) || is_h(b) {
            continue;
        }
        adj[a].push(b as u32);
        adj[b].push(a as u32);
    }
    for list in &mut adj {
        list.sort_unstable();
        list.dedup();
    }
    let mut paths: Vec<Vec<u32>> = (0..dim as u32).map(|i| vec![i]).collect();
    for _ in 1..upper {
        let mut next = Vec::new();
        for path in &paths {
            let end = *path.last().expect("paths are never empty") as usize;
            for &other in &adj[end] {
                if !path.contains(&other)
                    || (upper > 2 && path.len() == upper - 1 && path[path.len() - 2] != other)
                {
                    let mut p = path.clone();
                    p.push(other);
                    next.push(p);
                }
            }
        }
        paths = next;
    }
    if !use_bonds && upper == 1 {
        return paths;
    }
    let mut seen: FxHashSet<Vec<u32>> = FxHashSet::default();
    let mut out = Vec::new();
    for path in paths {
        let mut bonds: Vec<u32> = path
            .windows(2)
            .map(|w| {
                mol.bond_between(AtomIdx(w[0]), AtomIdx(w[1]))
                    .expect("consecutive path atoms are bonded")
                    .0
                    .0
            })
            .collect();
        let ordered = bonds.clone();
        bonds.sort_unstable();
        if seen.insert(bonds) {
            out.push(if use_bonds { ordered } else { path });
        }
    }
    out
}

/// RDKit's Hall-Kier valence deltas (`hkDeltas`), already as `1/sqrt(δv)`
/// (0 where δv is 0): `nOuter - H` up to neon, otherwise
/// `(nOuter - H) / (Z - nOuter - 1)`; hydrogen and dummy atoms give 0.
fn rdkit_hk_deltas(mol: &Molecule) -> Vec<f64> {
    mol.atoms()
        .map(|(idx, atom)| {
            let n = if atom.wildcard {
                0
            } else {
                atom.element.atomic_number() as i32
            };
            let d = if n <= 1 {
                0.0
            } else {
                let outer = RDKIT_N_OUTER_ELECS[n as usize - 1] as i32;
                let num = outer - i32::from(implicit_hcount(mol, idx));
                if n <= 10 {
                    num as f64
                } else {
                    num as f64 / (n - outer - 1) as f64
                }
            };
            if d != 0.0 { 1.0 / d.sqrt() } else { d }
        })
        .collect()
}

/// RDKit's `calcChiNv` for n >= 2 over [`rdkit_paths_of_length`].
fn rdkit_chi_nv(mol: &Molecule, hk: &[f64], n: usize) -> f64 {
    let mut res = 0.0;
    for p in rdkit_paths_of_length(mol, n + 1, false) {
        let mut accum = 1.0;
        for &a in &p[..n] {
            accum *= hk[a as usize];
        }
        // A closed ring path counts its first atom once (RDKit github #463).
        if p[n] != p[0] {
            accum *= hk[p[n] as usize];
        }
        res += accum;
    }
    res
}

/// RDKit `CalcChi0v` .. `CalcChi4v`, bit-identical.
fn rdkit_chi_v(mol: &Molecule) -> [f64; 5] {
    let hk = rdkit_hk_deltas(mol);
    [
        rdkit_chi0v_from(&hk),
        rdkit_chi1v_from(mol, &hk),
        rdkit_chi_nv(mol, &hk, 2),
        rdkit_chi_nv(mol, &hk, 3),
        rdkit_chi_nv(mol, &hk, 4),
    ]
}

/// One of [`rdkit_chi_v`]'s values, without computing the others.
fn rdkit_chi_v_order(mol: &Molecule, n: usize) -> f64 {
    let hk = rdkit_hk_deltas(mol);
    match n {
        0 => rdkit_chi0v_from(&hk),
        1 => rdkit_chi1v_from(mol, &hk),
        _ => rdkit_chi_nv(mol, &hk, n),
    }
}

fn rdkit_chi0v_from(hk: &[f64]) -> f64 {
    hk.iter().fold(0.0, |acc, x| acc + x)
}

/// Bonds in RDKit's numbering (ring closures last), as RDKit sums them.
fn rdkit_chi1v_from(mol: &Molecule, hk: &[f64]) -> f64 {
    let mut chi1v = 0.0;
    for b in mol.rdkit_bond_order() {
        let bond = mol.bond(b);
        chi1v += hk[bond.atom1.0 as usize] * hk[bond.atom2.0 as usize];
    }
    chi1v
}

/// RDKit's Python `GraphDescriptors.Chi0`: `sqrt(1 / degree)` over every
/// atom with a neighbour (graph hydrogens included), summed in atom order.
fn rdkit_chi0(mol: &Molecule) -> f64 {
    let mut res = 0.0;
    for (idx, _) in mol.atoms() {
        let d = mol.degree(idx) as f64;
        if d != 0.0 {
            res += (1.0 / d).sqrt();
        }
    }
    res
}

/// RDKit's Python `GraphDescriptors.Chi1`: `sqrt(1 / (deg_i * deg_j))` over
/// the bonds with nonzero product, summed left to right.
fn rdkit_chi1(mol: &Molecule) -> f64 {
    let mut res = 0.0;
    for b in mol.rdkit_bond_order() {
        let bond = mol.bond(b);
        let c = (mol.degree(bond.atom1) * mol.degree(bond.atom2)) as f64;
        if c != 0.0 {
            res += (1.0 / c).sqrt();
        }
    }
    res
}

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
    chematic_smiles::with_rdkit_model_scalar(
        mol,
        chematic_core::derived_cache::DerivedScalar::LabuteAsa,
        labute_asa_model,
    )
}

fn labute_asa_model(mol: &Molecule) -> f64 {
    if let Some(h) = chematic_smiles::rdkit_hydrogen_suppressed(mol) {
        return labute_asa_core(&h, false).2;
    }
    labute_asa_core(mol, false).2
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
        let rb1 = vdw_volume_rb0(z1);
        let rb2 = vdw_volume_rb0(z2);
        if rb1 > 1e-10 && rb2 > 1e-10 {
            v -= sphere_intersection(r_vdw_bondi(z1), r_vdw_bondi(z2), rb1 + rb2);
        }
    }

    for (idx, atom) in mol.atoms() {
        let z = atom.element.atomic_number();
        let rb_heavy = vdw_volume_rb0(z);
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
include!("topo_descriptors_tests.rs");
