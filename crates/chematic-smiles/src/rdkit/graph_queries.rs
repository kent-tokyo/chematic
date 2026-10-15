//! Checked graph queries following RDKit Release_2026_09_1 Matrices.cpp.
//! Upstream copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC.
//! Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors.
//! BSD-3-Clause; see THIRD_PARTY_NOTICES.md. The existing sanitized model
//! supplies atom indices, bond ordering and endpoint valence contributions.
use super::{RdkitSmilesError, matrices, mol::Mol};
use chematic_core::Molecule;

fn unsupported(detail: impl Into<String>) -> RdkitSmilesError {
    RdkitSmilesError::Unsupported(detail.into())
}

pub(super) fn prepare(mol: &Molecule) -> Result<Mol, RdkitSmilesError> {
    if mol.atom_count() > 1024 || mol.bond_count() > 16_384 {
        return Err(unsupported(format!(
            "graph query resource limit: {} atoms / {} bonds; maxima are 1024 / 16384",
            mol.atom_count(),
            mol.bond_count()
        )));
    }
    if mol.bonds().any(|(_, b)| b.atom1 == b.atom2) {
        return Err(unsupported("graph queries do not support self bonds"));
    }
    let model = super::rdkit_mol_for_writing(mol)?;
    matrices::checked_len(model.atoms.len())?;
    Ok(model)
}

/// RDKit 2026.09.1 `Chem.GetAdjacencyMatrix(m, useBO, emptyVal=0)`.
/// Rows use the same sanitized atom order as `rdkit_distance_matrix`,
/// including graph-H removal (or retained explicitly added H atoms).
/// Unweighted entries are 1 for each bond. Weighted entries are endpoint
/// valence contributions: aromatic 1.5, dative donor→acceptor 0 / reverse 1.
///
/// Only `empty_value=0` is supported. RDKit's nonzero `emptyVal` uses a byte
/// fill and can produce NaNs or platform-dependent integer conversions.
/// Inputs are capped before preparation at 1024 atoms and 16384 bonds;
/// dense storage uses checked dimensions and fallible allocation.
pub fn rdkit_adjacency_matrix(
    mol: &Molecule,
    use_bo: bool,
    empty_value: i32,
) -> Result<Vec<Vec<f64>>, RdkitSmilesError> {
    if empty_value != 0 {
        return Err(unsupported(
            "nonzero adjacency empty_value is unsupported; use 0",
        ));
    }
    let model = prepare(mol)?;
    let n = model.atoms.len();
    let mut matrix = matrices::filled(matrices::checked_len(n)?, 0.0)?;
    for bond in &model.bonds {
        matrix[bond.begin * n + bond.end] = if use_bo {
            bond.valence_contrib(bond.begin)
        } else {
            1.0
        };
        matrix[bond.end * n + bond.begin] = if use_bo {
            bond.valence_contrib(bond.end)
        } else {
            1.0
        };
    }
    matrices::into_rows(matrix, n)
}

fn indices(capacity: usize) -> Result<Vec<usize>, RdkitSmilesError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(capacity)
        .map_err(|e| unsupported(format!("graph query allocation failed: {e}")))?;
    Ok(result)
}

/// RDKit 2026.09.1 `Chem.GetShortestPath(m, start, end)`, including both
/// endpoints. Unweighted, undirected BFS in RDKit's bond-index neighbor order.
/// Tied paths therefore need not be reversals of the opposite-direction call.
/// Different fragments return an empty vector; equal or out-of-range endpoints
/// return `Unsupported` instead of RDKit's invariant/range exceptions.
/// Endpoints refer to the sanitized RDKit model, after graph-H removal.
/// The same input caps as `rdkit_adjacency_matrix` apply; BFS uses O(n) storage.
pub fn rdkit_shortest_path(
    mol: &Molecule,
    start: usize,
    end: usize,
) -> Result<Vec<usize>, RdkitSmilesError> {
    let model = prepare(mol)?;
    let n = model.atoms.len();
    if start >= n || end >= n || start == end {
        return Err(unsupported(format!(
            "shortest path requires distinct atom indices below {n}; got {start}, {end}"
        )));
    }
    let mut parent = indices(n)?;
    parent.resize(n, usize::MAX);
    parent[start] = start;
    let mut queue = indices(n)?;
    queue.push(start);
    let mut cursor = 0;
    while cursor < queue.len() && parent[end] == usize::MAX {
        let current = queue[cursor];
        cursor += 1;
        for &bond in &model.atom_bonds[current] {
            let neighbor = model.bonds[bond].other(current);
            if parent[neighbor] == usize::MAX {
                parent[neighbor] = current;
                queue.push(neighbor);
                if neighbor == end {
                    break;
                }
            }
        }
    }
    if parent[end] == usize::MAX {
        return Ok(Vec::new());
    }
    let mut path = indices(n)?;
    let mut current = end;
    path.push(current);
    while current != start {
        current = parent[current];
        path.push(current);
    }
    path.reverse();
    Ok(path)
}
