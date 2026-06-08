//! Gasteiger-Marsili PEOE partial charges.
//!
//! Implements the Partial Equalization of Orbital Electronegativities (PEOE)
//! method: Gasteiger & Marsili, Tetrahedron 1980, 36, 3219–3228.
//!
//! Parameters from RDKit GasteigerParams.cpp. Runs on an H-explicit molecule
//! (via `add_hydrogens`) and returns charges for heavy atoms only.

use crate::hydrogen::add_hydrogens;
use chematic_core::{AtomIdx, BondOrder, Molecule};

/// (a, b, c) electronegativity polynomial: χ(q) = a + b·q + c·q²
type AbC = (f64, f64, f64);

fn hybridisation_sp(mol: &Molecule, idx: AtomIdx) -> u8 {
    if mol.atom(idx).aromatic {
        return 2;
    }
    let mut has_double = false;
    let mut has_triple = false;
    for (_, bidx) in mol.neighbors(idx) {
        match mol.bond(bidx).order {
            BondOrder::Double => has_double = true,
            BondOrder::Triple => has_triple = true,
            _ => {}
        }
    }
    if has_triple {
        1
    } else if has_double {
        2
    } else {
        3
    }
}

/// Return (a, b, c) for the given atom; None for unsupported elements.
fn params(mol: &Molecule, idx: AtomIdx) -> Option<AbC> {
    let atom = mol.atom(idx);
    let an = atom.element.atomic_number();
    let hyb = hybridisation_sp(mol, idx);
    let ar = atom.aromatic;

    Some(match an {
        1 => (7.17, 6.24, -0.56),
        6 => match hyb {
            3 => (7.98, 9.18, 1.88),
            2 => (8.79, 9.32, 1.51),
            _ => (10.39, 9.45, 0.73),
        },
        7 => match hyb {
            3 => (11.54, 10.82, 1.36),
            2 => (12.87, 11.15, 0.85),
            _ => (15.68, 11.70, -0.27),
        },
        8 => {
            if ar || hyb == 2 {
                (17.07, 13.79, 0.47)
            } else {
                (14.18, 12.92, 1.39)
            }
        }
        16 => {
            if ar {
                (10.88, 9.49, 1.33)
            } else {
                (10.14, 9.13, 1.38)
            }
        }
        9 => (14.66, 13.85, 2.31),
        17 => (11.00, 9.69, 1.35),
        35 => (10.08, 8.47, 1.16),
        53 => (9.90, 7.96, 0.96),
        15 => {
            if hyb == 2 {
                (9.665, 8.530, 0.735)
            } else {
                (8.90, 8.24, 0.96)
            }
        }
        _ => return None,
    })
}

/// Compute Gasteiger-Marsili PEOE partial charges.
///
/// Returns a `Vec<f64>` indexed by atom order (same as `mol.atoms()`).
/// Implicit hydrogens are added internally for accurate charge propagation,
/// then stripped. Atoms with unsupported elements receive charge 0.0.
pub fn gasteiger_charges(mol: &Molecule) -> Vec<f64> {
    const N_ITER: usize = 12;
    const DAMP0: f64 = 0.5;

    let hmol = add_hydrogens(mol);
    let n = hmol.atom_count();
    let mut q: Vec<f64> = vec![0.0; n];

    let pars: Vec<Option<AbC>> = (0..n).map(|i| params(&hmol, AtomIdx(i as u32))).collect();

    let chi = |par: AbC, qi: f64| par.0 + par.1 * qi + par.2 * qi * qi;

    let mut damp = DAMP0;
    for _ in 0..N_ITER {
        let ens: Vec<f64> = (0..n)
            .map(|i| pars[i].map_or(0.0, |p| chi(p, q[i])))
            .collect();

        let mut dq: Vec<f64> = vec![0.0; n];
        for (_, bond) in hmol.bonds() {
            let i = bond.atom1.0 as usize;
            let j = bond.atom2.0 as usize;
            let (pi, pj) = match (pars[i], pars[j]) {
                (Some(a), Some(b)) => (a, b),
                _ => continue,
            };
            let ei = ens[i];
            let ej = ens[j];
            if (ei - ej).abs() < 1e-10 {
                continue;
            }

            // Electron density flows from lower EN (donor) to higher EN (acceptor).
            // Donor becomes more positive → q[donor] increases.
            // Acceptor becomes more negative → q[acceptor] decreases.
            let (donor, acceptor, e_donor, e_acceptor, p_donor) = if ei < ej {
                (i, j, ei, ej, pi)
            } else {
                (j, i, ej, ei, pj)
            };

            // Cation EN of donor = a+b+c (q=+1 state).
            let cation_en = p_donor.0 + p_donor.1 + p_donor.2;
            if cation_en.abs() < 1e-10 {
                continue;
            }

            let delta = (e_acceptor - e_donor) / cation_en * damp;
            // Donor loses electron density → more positive
            dq[donor] += delta;
            // Acceptor gains electron density → more negative
            dq[acceptor] -= delta;
        }

        for i in 0..n {
            q[i] += dq[i];
        }
        damp *= DAMP0;
    }

    // Return only heavy-atom charges (first `mol.atom_count()` atoms in hmol).
    q.truncate(mol.atom_count());
    q
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    #[test]
    fn methanol_oxygen_more_negative_than_carbon() {
        let mol = parse("CO").unwrap();
        let q = gasteiger_charges(&mol);
        let o_idx = mol
            .atoms()
            .find(|(_, a)| a.element.atomic_number() == 8)
            .map(|(i, _)| i.0 as usize)
            .unwrap();
        let c_idx = mol
            .atoms()
            .find(|(_, a)| a.element.atomic_number() == 6)
            .map(|(i, _)| i.0 as usize)
            .unwrap();
        assert!(
            q[o_idx] < q[c_idx],
            "O charge {:.4} should be < C charge {:.4}",
            q[o_idx],
            q[c_idx]
        );
    }

    #[test]
    fn water_oxygen_negative() {
        let mol = parse("O").unwrap();
        let q = gasteiger_charges(&mol);
        assert!(
            q[0] < 0.0,
            "water O should have negative charge, got {:.4}",
            q[0]
        );
    }

    #[test]
    fn charge_sum_near_zero_neutral_molecule() {
        // Only heavy atoms returned, but the heavy-atom sum should be small.
        let mol = parse("CC(=O)O").unwrap();
        let q = gasteiger_charges(&mol);
        let sum: f64 = q.iter().sum();
        // Heavy atoms alone won't sum to exactly 0 (H charges excluded).
        // Just verify it's not wildly off.
        assert!(sum.abs() < 2.0, "heavy-atom charge sum = {sum:.4}");
    }

    #[test]
    fn aspirin_charges_vector_length() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let q = gasteiger_charges(&mol);
        assert_eq!(q.len(), mol.atom_count());
    }

    #[test]
    fn electronegative_atoms_negative() {
        // In acetic acid, both oxygens should be negative.
        let mol = parse("CC(=O)O").unwrap();
        let q = gasteiger_charges(&mol);
        for (idx, atom) in mol.atoms() {
            if atom.element.atomic_number() == 8 {
                assert!(
                    q[idx.0 as usize] < 0.0,
                    "O charge should be negative, got {:.4}",
                    q[idx.0 as usize]
                );
            }
        }
    }
}
