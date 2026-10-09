//! `MolOps::getDistanceMat` and `MolOps::get3DDistanceMat` (RDKit
//! 2026.03.1 `Matrices.cpp`).

use super::mol::Mol;

/// `LOCAL_INF`: the distance of atoms in different fragments.
const LOCAL_INF: f64 = 1e8;

/// `getDistanceMat(mol, useBO, useAtomWts)`: Floyd-Warshall over bond
/// weights 1 (`useBO`: `1 / order`, aromatic bonds `2/3`), then with
/// `useAtomWts` the diagonal set to `6 / atomic number`.
pub(crate) fn distance_mat(mol: &Mol, use_bo: bool, use_atom_wts: bool) -> Vec<f64> {
    let n = mol.atoms.len();
    let mut d = vec![LOCAL_INF; n * n];
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
    floyd_warshall(n, &mut d);
    if use_atom_wts {
        for i in 0..n {
            d[i * n + i] = 6.0 / f64::from(mol.atoms[i].anum);
        }
    }
    d
}

/// RDKit's `FloydWarshall` (`v1 <= v2` keeps the old distance).
fn floyd_warshall(n: usize, d: &mut [f64]) {
    let mut last = d.to_vec();
    let mut curr = vec![0.0; n * n];
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
}

/// `get3DDistanceMat(mol, confId, useAtomWts)` for the conformer `coords`.
pub(crate) fn distance_mat_3d(mol: &Mol, coords: &[[f64; 3]], use_atom_wts: bool) -> Vec<f64> {
    let n = mol.atoms.len();
    let mut d = vec![0.0; n * n];
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
    d
}
