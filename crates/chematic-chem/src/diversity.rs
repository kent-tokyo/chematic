//! Diversity picking and clustering algorithms.
//!
//! Functions accept a `sim_fn` closure so callers can supply any Tanimoto
//! variant (e.g. `chematic_fp::tanimoto_topo_path`).

use chematic_core::Molecule;

/// MaxMin diversity picking.
///
/// Returns the indices of `n` maximally diverse molecules from `mols`.
/// The first picked molecule is the one with the highest internal diversity
/// (maximum minimum distance to all others); subsequent picks maximise the
/// minimum distance to already-selected molecules.
///
/// `sim_fn(a, b) -> f64` must return a value in [0, 1] (Tanimoto-like).
/// Distance is `1 - sim`.
///
/// If `n >= mols.len()`, all indices are returned in selection order.
pub fn maxmin_picks<F>(mols: &[Molecule], n: usize, sim_fn: F) -> Vec<usize>
where
    F: Fn(&Molecule, &Molecule) -> f64,
{
    if mols.is_empty() || n == 0 {
        return Vec::new();
    }
    let n = n.min(mols.len());

    // min_dist[i] = minimum distance from mol i to any already-selected mol.
    let mut min_dist: Vec<f64> = vec![f64::MAX; mols.len()];

    // Seed: pick the molecule with highest average distance to all others
    // (most "isolated"). A cheap O(n²) seed scan.
    let seed = (0..mols.len())
        .max_by(|&a, &b| {
            let da: f64 = mols.iter().map(|m| 1.0 - sim_fn(&mols[a], m)).sum();
            let db: f64 = mols.iter().map(|m| 1.0 - sim_fn(&mols[b], m)).sum();
            da.partial_cmp(&db).unwrap()
        })
        .unwrap_or(0);

    let mut picked = vec![seed];

    // Update min_dist for seed.
    for i in 0..mols.len() {
        let d = 1.0 - sim_fn(&mols[seed], &mols[i]);
        if d < min_dist[i] {
            min_dist[i] = d;
        }
    }

    while picked.len() < n {
        // Next pick: mol with maximum min_dist (most distant from all selected).
        let next = (0..mols.len())
            .filter(|i| !picked.contains(i))
            .max_by(|&a, &b| min_dist[a].partial_cmp(&min_dist[b]).unwrap())
            .unwrap();

        picked.push(next);

        // Update min_dist.
        for i in 0..mols.len() {
            let d = 1.0 - sim_fn(&mols[next], &mols[i]);
            if d < min_dist[i] {
                min_dist[i] = d;
            }
        }
    }

    picked
}

/// Butina clustering (Taylor-Butina).
///
/// Returns clusters as `Vec<Vec<usize>>` where each inner Vec contains the
/// indices of molecules in that cluster (centroid first). Clusters are ordered
/// by size (largest first).
///
/// `cutoff` is a similarity threshold in [0, 1]; molecules within `cutoff`
/// Tanimoto similarity of the centroid belong to the same cluster.
/// Common default: 0.65 (Tanimoto similarity ≥ 0.65 → same cluster).
pub fn butina_cluster<F>(mols: &[Molecule], cutoff: f64, sim_fn: F) -> Vec<Vec<usize>>
where
    F: Fn(&Molecule, &Molecule) -> f64,
{
    if mols.is_empty() {
        return Vec::new();
    }

    // Build neighbour lists: for each mol, which others are ≥ cutoff similar.
    let n = mols.len();
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            if sim_fn(&mols[i], &mols[j]) >= cutoff {
                neighbours[i].push(j);
                neighbours[j].push(i);
            }
        }
    }

    // Butina: repeatedly pick the mol with most neighbours as centroid,
    // remove it and its neighbours from the pool.
    let mut assigned = vec![false; n];
    let mut clusters: Vec<Vec<usize>> = Vec::new();

    loop {
        // Find unassigned mol with most unassigned neighbours.
        let centroid = (0..n)
            .filter(|&i| !assigned[i])
            .max_by_key(|&i| neighbours[i].iter().filter(|&&j| !assigned[j]).count());
        let centroid = match centroid {
            Some(c) => c,
            None => break,
        };

        let mut cluster = vec![centroid];
        assigned[centroid] = true;

        for &nb in &neighbours[centroid].clone() {
            if !assigned[nb] {
                cluster.push(nb);
                assigned[nb] = true;
            }
        }

        clusters.push(cluster);
    }

    clusters.sort_by_key(|b| std::cmp::Reverse(b.len()));
    clusters
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn sim(a: &Molecule, b: &Molecule) -> f64 {
        chematic_fp::tanimoto_topo_path(a, b)
    }

    #[test]
    fn maxmin_picks_correct_count() {
        let smiles = ["c1ccccc1", "CC(=O)O", "CCN", "c1ccncc1", "CCCO"];
        let mols: Vec<Molecule> = smiles.iter().map(|s| parse(s).unwrap()).collect();
        let picked = maxmin_picks(&mols, 3, sim);
        assert_eq!(picked.len(), 3);
        // No duplicates.
        let mut sorted = picked.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(sorted.len(), 3);
    }

    #[test]
    fn maxmin_clamps_to_mol_count() {
        let mols: Vec<Molecule> = ["C", "CC"].iter().map(|s| parse(s).unwrap()).collect();
        assert_eq!(maxmin_picks(&mols, 100, sim).len(), 2);
    }

    #[test]
    fn butina_identical_molecules_one_cluster() {
        let mols: Vec<Molecule> = (0..3).map(|_| parse("c1ccccc1").unwrap()).collect();
        let clusters = butina_cluster(&mols, 0.5, sim);
        // All identical → one cluster of size 3.
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].len(), 3);
    }

    #[test]
    fn butina_dissimilar_molecules_separate_clusters() {
        let smiles = ["c1ccccc1", "CCCCCCC", "c1ccncc1", "NCCN"];
        let mols: Vec<Molecule> = smiles.iter().map(|s| parse(s).unwrap()).collect();
        let clusters = butina_cluster(&mols, 0.8, sim);
        // At high cutoff, each mol should be its own cluster.
        assert_eq!(clusters.len(), smiles.len());
    }

    #[test]
    fn butina_all_molecules_assigned() {
        let smiles = ["c1ccccc1", "CC(=O)O", "CCN", "c1ccncc1", "CCCO"];
        let mols: Vec<Molecule> = smiles.iter().map(|s| parse(s).unwrap()).collect();
        let clusters = butina_cluster(&mols, 0.5, sim);
        let total: usize = clusters.iter().map(|c| c.len()).sum();
        assert_eq!(total, mols.len());
    }
}
