//! EState indices (Hall & Kier, 1991).
//!
//! Computes per-atom electrotopological state values that encode both
//! electronic and topological environment.  Only heavy atoms receive
//! meaningful values; hydrogen positions receive 0.0.

use std::collections::{HashSet, VecDeque};

use chematic_core::{AtomIdx, Molecule, implicit_hcount};

// ─── Atom-level helpers ──────────────────────────────────────────────────────

fn principal_quantum_number(z: u8) -> f64 {
    if z <= 2 {
        1.0
    } else if z <= 10 {
        2.0
    } else if z <= 18 {
        3.0
    } else if z <= 36 {
        4.0
    } else if z <= 54 {
        5.0
    } else {
        6.0
    }
}

fn valence_electrons(z: u8) -> u8 {
    match z {
        1 => 1,
        2 => 2,
        3..=10 => z - 2,
        11..=18 => z - 10,
        19 | 20 => 1 + (z - 18), // K, Ca
        // period 4 d-block: use group valence for common elements
        35 => 7, // Br
        36 => 8, // Kr
        53 => 7, // I
        54 => 8, // Xe
        // fall-back: number of electrons beyond last noble gas shell (period 4+)
        _ => {
            let rem = z % 18;
            if rem == 0 { 8 } else { rem.min(8) }
        }
    }
}

/// BFS shortest-path distances from `start` in the heavy-atom subgraph.
/// Returns `usize::MAX` for unreachable atoms.
fn bfs_heavy(mol: &Molecule, start: usize, heavy_set: &HashSet<usize>) -> Vec<usize> {
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

// ─── Public API ──────────────────────────────────────────────────────────────

/// Per-atom EState values indexed by atom position in `mol.atoms()`.
///
/// Hydrogen atoms and atoms with no heavy-atom neighbours (δ = 0) that
/// are isolated receive `0.0`.
pub fn estate_indices(mol: &Molecule) -> Vec<f64> {
    let n = mol.atom_count();
    if n == 0 {
        return vec![];
    }

    let heavy: Vec<usize> = mol
        .atoms()
        .filter(|(_, a)| a.element.atomic_number() != 1)
        .map(|(idx, _)| idx.0 as usize)
        .collect();
    let heavy_set: HashSet<usize> = heavy.iter().copied().collect();

    // ── Intrinsic state I_i ──────────────────────────────────────────────────
    // I_i = ((2/N)² · δᵥ + 1) / δ   when δ > 0
    //      = (2/N)² · δᵥ + 1         when δ = 0  (isolated heavy atom)
    let mut intrinsic = vec![0.0f64; n];
    for &i in &heavy {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        let z = atom.element.atomic_number();
        let pqn = principal_quantum_number(z);
        let zv = valence_electrons(z) as f64;
        let h_total = implicit_hcount(mol, idx) as f64 + atom.hydrogen_count.unwrap_or(0) as f64;
        let delta_v = (zv - h_total).max(0.0);
        let delta = mol
            .neighbors(idx)
            .filter(|(nb, _)| heavy_set.contains(&(nb.0 as usize)))
            .count() as f64;

        let base = (2.0 / pqn).powi(2) * delta_v + 1.0;
        intrinsic[i] = if delta > 0.0 { base / delta } else { base };
    }

    // ── EState perturbation ──────────────────────────────────────────────────
    // S_i = I_i + Σ_{j≠i} (I_i − I_j) / r_ij²
    // r_ij = topological distance + 1 (Kier & Hall, 1991) — the "+1" is not
    // an off-by-one, it's part of the definition (adjacent atoms get r_ij=2).
    let mut estate = intrinsic.clone();
    for &i in &heavy {
        let dists = bfs_heavy(mol, i, &heavy_set);
        for &j in &heavy {
            if i == j {
                continue;
            }
            let r = dists[j];
            if r == usize::MAX {
                continue;
            }
            let rij = (r + 1) as f64;
            estate[i] += (intrinsic[i] - intrinsic[j]) / (rij * rij);
        }
    }

    estate
}

/// Maximum EState index over all heavy atoms.
pub fn max_estate(mol: &Molecule) -> f64 {
    estate_indices(mol)
        .into_iter()
        .filter(|v| *v != 0.0)
        .fold(f64::NEG_INFINITY, f64::max)
}

/// Minimum EState index over all heavy atoms.
pub fn min_estate(mol: &Molecule) -> f64 {
    estate_indices(mol)
        .into_iter()
        .filter(|v| *v != 0.0)
        .fold(f64::INFINITY, f64::min)
}

/// Sum of EState indices over all heavy atoms.
pub fn sum_estate(mol: &Molecule) -> f64 {
    estate_indices(mol).into_iter().sum()
}

/// Compute sum, max, and min EState in a single `estate_indices` pass.
///
/// Use when all three values are needed to avoid three redundant computations.
/// Returns `(sum, max, min)`.  `max` and `min` ignore zero-valued atoms.
pub fn estate_all(mol: &Molecule) -> (f64, f64, f64) {
    let indices = estate_indices(mol);
    let sum = indices.iter().copied().sum();
    let max = indices
        .iter()
        .copied()
        .filter(|v| *v != 0.0)
        .fold(f64::NEG_INFINITY, f64::max);
    let min = indices
        .iter()
        .copied()
        .filter(|v| *v != 0.0)
        .fold(f64::INFINITY, f64::min);
    (sum, max, min)
}

// ─── Tests ───────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn ethane_symmetric_estate() {
        let mol = parse("CC").unwrap();
        let es = estate_indices(&mol);
        // Both carbons in ethane are symmetric — EState values must be equal.
        assert!(
            (es[0] - es[1]).abs() < 1e-9,
            "ethane carbons should have equal EState: {} vs {}",
            es[0],
            es[1]
        );
    }

    #[test]
    fn acetic_acid_oxygen_higher_than_carbon() {
        // Oxygen is more electronegative → higher intrinsic state.
        let mol = parse("CC(=O)O").unwrap(); // acetic acid: C, C(=O), O(=C), O-H
        let es = estate_indices(&mol);
        // At least one O should have higher EState than both Cs.
        let max_o = mol
            .atoms()
            .filter(|(_, a)| a.element == chematic_core::Element::O)
            .map(|(idx, _)| es[idx.0 as usize])
            .fold(f64::NEG_INFINITY, f64::max);
        let max_c = mol
            .atoms()
            .filter(|(_, a)| a.element == chematic_core::Element::C)
            .map(|(idx, _)| es[idx.0 as usize])
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            max_o > max_c,
            "O EState ({max_o:.3}) should exceed C EState ({max_c:.3})"
        );
    }

    #[test]
    fn pyridine_nitrogen_highest_estate() {
        let mol = parse("c1ccncc1").unwrap();
        let es = estate_indices(&mol);
        let n_idx = mol
            .atoms()
            .find(|(_, a)| a.element == chematic_core::Element::N)
            .map(|(idx, _)| idx.0 as usize)
            .unwrap();
        let max_c = mol
            .atoms()
            .filter(|(_, a)| a.element == chematic_core::Element::C)
            .map(|(idx, _)| es[idx.0 as usize])
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            es[n_idx] > max_c,
            "pyridine N ({:.3}) should have higher EState than max C ({:.3})",
            es[n_idx],
            max_c
        );
    }

    #[test]
    fn sum_estate_positive() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap(); // aspirin
        assert!(
            sum_estate(&mol) > 0.0,
            "aspirin sum_estate should be positive"
        );
    }

    #[test]
    fn max_min_estate_ordering() {
        let mol = parse("CC(=O)O").unwrap();
        assert!(
            max_estate(&mol) >= min_estate(&mol),
            "max_estate should be >= min_estate"
        );
    }

    #[test]
    fn estate_indices_length_matches_atom_count() {
        let mol = parse("c1ccccc1").unwrap();
        assert_eq!(estate_indices(&mol).len(), mol.atom_count());
    }

    #[test]
    fn estate_indices_match_rdkit_ethanol() {
        // RDKit EState.EStateIndices("CCO") verified: [1.6806, 0.25, 7.5694].
        // r_ij must be (topological distance + 1), not the raw distance —
        // that's the specific bug this test guards against.
        let mol = parse("CCO").unwrap();
        let es = estate_indices(&mol);
        let expected = [1.6806, 0.25, 7.5694];
        for (got, want) in es.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 0.001, "{es:?} vs {expected:?}");
        }
    }

    #[test]
    fn estate_all_matches_individual() {
        for smi in ["CC", "CC(=O)O", "c1ccccc1", "CN1C=NC2=C1C(=O)N(C)C(=O)N2C"] {
            let mol = parse(smi).unwrap();
            let (s, mx, mn) = estate_all(&mol);
            assert!((s - sum_estate(&mol)).abs() < 1e-10, "{smi}: sum mismatch");
            assert!((mx - max_estate(&mol)).abs() < 1e-10, "{smi}: max mismatch");
            assert!((mn - min_estate(&mol)).abs() < 1e-10, "{smi}: min mismatch");
        }
    }
}
