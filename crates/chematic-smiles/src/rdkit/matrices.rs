//! `MolOps::getDistanceMat` and `MolOps::get3DDistanceMat` (RDKit
//! 2026.03.1 `Matrices.cpp`).

use super::{RdkitSmilesError, mol::Mol};
use chematic_core::Molecule;

/// `LOCAL_INF`: the distance of atoms in different fragments.
const LOCAL_INF: f64 = 1e8;

// Bound both input preparation and quadratic storage. The topological path
// holds three flat f64 buffers during Floyd-Warshall (24 MiB at this cap).
const MAX_MATRIX_ENTRIES: usize = 1024 * 1024;

pub(super) fn checked_len(n: usize) -> Result<usize, RdkitSmilesError> {
    let entries = n.checked_mul(n).ok_or_else(|| {
        RdkitSmilesError::Unsupported(format!("distance matrix dimensions overflow for {n} atoms"))
    })?;
    if entries > MAX_MATRIX_ENTRIES {
        return Err(RdkitSmilesError::Unsupported(format!(
            "distance matrix resource limit: {n} atoms require {entries} entries; maximum is {MAX_MATRIX_ENTRIES} entries (1024 atoms)"
        )));
    }
    Ok(entries)
}

fn filled(len: usize, value: f64) -> Result<Vec<f64>, RdkitSmilesError> {
    let mut values = Vec::new();
    values.try_reserve_exact(len).map_err(|err| {
        RdkitSmilesError::Unsupported(format!("distance matrix allocation failed: {err}"))
    })?;
    values.resize(len, value);
    Ok(values)
}

pub(super) fn into_rows(flat: Vec<f64>, n: usize) -> Result<Vec<Vec<f64>>, RdkitSmilesError> {
    let mut rows = Vec::new();
    rows.try_reserve_exact(n).map_err(|err| {
        RdkitSmilesError::Unsupported(format!("distance matrix row allocation failed: {err}"))
    })?;
    for chunk in flat.chunks(n.max(1)).take(n) {
        let mut row = filled(chunk.len(), 0.0)?;
        row.copy_from_slice(chunk);
        rows.push(row);
    }
    Ok(rows)
}

/// `getDistanceMat(mol, useBO, useAtomWts)`: Floyd-Warshall over bond
/// weights 1 (`useBO`: `1 / order`, aromatic bonds `2/3`), then with
/// `useAtomWts` the diagonal set to `6 / atomic number`.
pub(crate) fn distance_mat(
    mol: &Mol,
    use_bo: bool,
    use_atom_wts: bool,
) -> Result<Vec<f64>, RdkitSmilesError> {
    let n = mol.atoms.len();
    let mut d = filled(checked_len(n)?, LOCAL_INF)?;
    for i in 0..n {
        d[i * n + i] = 0.0;
    }
    for bond in &mol.bonds {
        let contrib = if !use_bo {
            1.0
        } else if bond.aromatic {
            2.0 / 3.0
        } else {
            1.0 / bond.bt.as_double()
        };
        d[bond.begin * n + bond.end] = contrib;
        d[bond.end * n + bond.begin] = contrib;
    }
    floyd_warshall(n, &mut d)?;
    if use_atom_wts {
        for i in 0..n {
            d[i * n + i] = 6.0 / f64::from(mol.atoms[i].anum);
        }
    }
    Ok(d)
}

/// RDKit's `FloydWarshall` (`v1 <= v2` keeps the old distance).
fn floyd_warshall(n: usize, d: &mut [f64]) -> Result<(), RdkitSmilesError> {
    let mut last = filled(d.len(), 0.0)?;
    last.copy_from_slice(d);
    let mut curr = filled(d.len(), 0.0)?;
    for k in 0..n {
        for i in 0..n {
            let dik = last[i * n + k];
            for j in 0..n {
                let v1 = last[i * n + j];
                let v2 = dik + last[k * n + j];
                curr[i * n + j] = if v1 <= v2 { v1 } else { v2 };
            }
        }
        std::mem::swap(&mut curr, &mut last);
    }
    d.copy_from_slice(&last);
    Ok(())
}

/// `get3DDistanceMat(mol, confId, useAtomWts)` for the conformer `coords`.
pub(crate) fn distance_mat_3d(
    mol: &Mol,
    coords: &[[f64; 3]],
    use_atom_wts: bool,
) -> Result<Vec<f64>, RdkitSmilesError> {
    let n = mol.atoms.len();
    let mut d = filled(checked_len(n)?, 0.0)?;
    for i in 0..n {
        d[i * n + i] = if use_atom_wts {
            6.0 / f64::from(mol.atoms[i].anum)
        } else {
            0.0
        };
        for j in i + 1..n {
            let (p, q) = (coords[i], coords[j]);
            let (dx, dy, dz) = (p[0] - q[0], p[1] - q[1], p[2] - q[2]);
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            d[i * n + j] = dist;
            d[j * n + i] = dist;
        }
    }
    Ok(d)
}

/// `Chem.GetDistanceMatrix(m, useBO, useAtomWts)` for
/// `m = Chem.MolFromSmiles(s)` (RDKit 2026.03.1; the `Chem.AddHs` molecule
/// when `mol` carries added hydrogens): topological distances by
/// Floyd-Warshall, `1e8` between fragments; `use_bo` weights bonds by
/// `1 / bond order` (aromatic `2/3`), `use_atom_wts` sets the diagonal to
/// `6 / atomic number`. One row per atom in RDKit's atom order.
///
/// Inputs and prepared models are limited to 1024 atoms (1,048,576 entries).
/// The input check precedes RDKit-model preparation, including H removal.
/// Size/overflow/allocation errors return [`RdkitSmilesError::Unsupported`]
/// rather than attempting an unbounded dense allocation (RDKit #9647).
pub fn rdkit_distance_matrix(
    mol: &Molecule,
    use_bo: bool,
    use_atom_wts: bool,
) -> Result<Vec<Vec<f64>>, RdkitSmilesError> {
    checked_len(mol.atom_count())?;
    let m = super::rdkit_mol_for_writing(mol)?;
    let n = m.atoms.len();
    let d = distance_mat(&m, use_bo, use_atom_wts)?;
    into_rows(d, n)
}

/// `Chem.Get3DDistanceMatrix(m, useAtomWts=...)` for RDKit's molecule of
/// `mol` (as in [`rdkit_distance_matrix`]) with the conformer `coords`
/// (one `[x, y, z]` per atom in RDKit's atom order).
/// The same 1024-input/prepared-atom cap and fallible-allocation policy as
/// [`rdkit_distance_matrix`] apply, before model preparation.
pub fn rdkit_distance_matrix_3d(
    mol: &Molecule,
    coords: &[[f64; 3]],
    use_atom_wts: bool,
) -> Result<Vec<Vec<f64>>, RdkitSmilesError> {
    checked_len(mol.atom_count())?;
    let m = super::rdkit_mol_for_writing(mol)?;
    let n = m.atoms.len();
    if coords.len() != n {
        return Err(RdkitSmilesError::Unsupported(format!(
            "{} coordinates for {n} atoms",
            coords.len()
        )));
    }
    let d = distance_mat_3d(&m, coords, use_atom_wts)?;
    into_rows(d, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checked_dimensions_cover_the_rdkit_9647_boundary_without_allocating() {
        assert_eq!(checked_len(0).unwrap(), 0);
        assert_eq!(checked_len(1024).unwrap(), MAX_MATRIX_ENTRIES);
        for n in [1025, 46_341, 65_535, 65_536, 65_537] {
            assert!(
                checked_len(n)
                    .unwrap_err()
                    .to_string()
                    .contains("resource limit")
            );
        }
        assert!(
            checked_len(usize::MAX)
                .unwrap_err()
                .to_string()
                .contains("overflow")
        );
    }

    #[test]
    fn nested_rows_handle_empty_and_small_matrices() {
        assert!(into_rows(Vec::new(), 0).unwrap().is_empty());
        assert_eq!(
            into_rows(vec![0., 1., 1., 0.], 2).unwrap(),
            [[0., 1.], [1., 0.]]
        );
    }
}
