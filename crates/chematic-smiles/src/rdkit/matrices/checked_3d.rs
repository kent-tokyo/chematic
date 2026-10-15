//! Opt-in numeric checks extending the RDKit #9647 allocation boundary.
use super::{Mol, RdkitSmilesError, checked_len, filled, into_rows};
use chematic_core::Molecule;

fn invalid(message: impl Into<String>) -> RdkitSmilesError {
    RdkitSmilesError::Unsupported(message.into())
}

/// Bounded 3D distances with a finite-result contract.
///
/// Atom preparation/order, H removal and `6 / Z` diagonals follow
/// [`super::rdkit_distance_matrix_3d`]. Inputs and prepared models are capped at
/// 1024 atoms; coordinates must match the prepared model exactly. Nonfinite
/// coordinates, undefined weighted dummy-atom diagonals and distances outside
/// finite f64 range return an error, without a partial result.
///
/// Ordinary finite calculations retain the legacy sum-of-squares expression.
/// If squaring overflows or underflows, a scaled `hypot` norm recovers a finite
/// representable result. The legacy compatibility API is unchanged.
pub fn rdkit_distance_matrix_3d_checked(
    mol: &Molecule,
    coords: &[[f64; 3]],
    use_atom_wts: bool,
) -> Result<Vec<Vec<f64>>, RdkitSmilesError> {
    checked_len(mol.atom_count())?;
    if mol.atom_count() == 0 {
        return if coords.is_empty() {
            Ok(Vec::new())
        } else {
            Err(invalid(format!("{} coordinates for 0 atoms", coords.len())))
        };
    }
    let prepared = super::super::rdkit_mol_for_writing(mol)?;
    let n = prepared.atoms.len();
    let entries = checked_len(n)?;
    if coords.len() != n {
        return Err(invalid(format!(
            "{} coordinates for {n} atoms",
            coords.len()
        )));
    }
    for (i, xyz) in coords.iter().enumerate() {
        if xyz.iter().any(|value| !value.is_finite()) {
            return Err(invalid(format!("non-finite coordinate at atom {i}")));
        }
    }
    if use_atom_wts {
        for (i, atom) in prepared.atoms.iter().enumerate() {
            if atom.anum == 0 {
                return Err(invalid(format!(
                    "atomic-weight diagonal is undefined for atomic number zero at atom {i}"
                )));
            }
        }
    }
    let flat = fill_checked(&prepared, coords, use_atom_wts, entries)?;
    into_rows(flat, n)
}

fn fill_checked(
    mol: &Mol,
    coords: &[[f64; 3]],
    use_atom_wts: bool,
    entries: usize,
) -> Result<Vec<f64>, RdkitSmilesError> {
    let n = mol.atoms.len();
    let mut result = filled(entries, 0.0)?;
    for i in 0..n {
        if use_atom_wts {
            result[i * n + i] = 6.0 / f64::from(mol.atoms[i].anum);
        }
        for j in i + 1..n {
            let p = coords[i];
            let q = coords[j];
            let (dx, dy, dz) = (p[0] - q[0], p[1] - q[1], p[2] - q[2]);
            let squared = dx * dx + dy * dy + dz * dz;
            let norm = if squared.is_infinite()
                || squared.is_subnormal()
                || (squared == 0.0 && (dx != 0.0 || dy != 0.0 || dz != 0.0))
            {
                dx.hypot(dy).hypot(dz)
            } else {
                squared.sqrt()
            };
            if !norm.is_finite() {
                return Err(invalid(format!(
                    "non-finite 3D distance between atoms {i} and {j}"
                )));
            }
            result[i * n + j] = norm;
            result[j * n + i] = norm;
        }
    }
    Ok(result)
}
