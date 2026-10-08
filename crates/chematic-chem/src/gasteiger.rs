//! Gasteiger-Marsili PEOE partial charges.
//!
//! Implements the Partial Equalization of Orbital Electronegativities (PEOE)
//! method: Gasteiger & Marsili, Tetrahedron 1980, 36, 3219–3228.
//!
//! A port of RDKit's `computeGasteigerCharges` (`GasteigerCharges.cpp`,
//! default 12 iterations), operation for operation, so the charges are
//! bit-identical to RDKit's `_GasteigerCharge`: hydrogens stay implicit and
//! each heavy atom carries one pooled hydrogen charge; parameters are looked
//! up by element and RDKit hybridization (`GasteigerParams.cpp`); formal
//! charges are first split over same-element atoms of a conjugated system.

use chematic_core::{AtomIdx, BondIdx, Molecule, implicit_hcount};

/// (a, b, c) electronegativity polynomial: χ(q) = a + b·q + c·q²
type AbC = (f64, f64, f64);

/// RDKit `IONXH`.
const IONXH: f64 = 20.02;
/// RDKit `DAMP` and `DAMP_SCALE`.
const DAMP: f64 = 0.5;
const DAMP_SCALE: f64 = 0.5;
const N_ITER: usize = 12;

/// RDKit's parameter mode for an atom (`computeGasteigerCharges`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    Sp3,
    Sp2,
    Sp,
    Star,
    So,
    So2,
    Unknown,
}

/// `GasteigerParams::getParams(elem, mode)`; unknown pairs fall back to
/// `X *` = (0, 0, 0).
fn params(z: u8, mode: Mode) -> AbC {
    use Mode::*;
    match (z, mode) {
        (1, Star) => (7.17, 6.24, -0.56),
        (6, Sp3) => (7.98, 9.18, 1.88),
        (6, Sp2) => (8.79, 9.32, 1.51),
        (6, Sp) => (10.39, 9.45, 0.73),
        (7, Sp3) => (11.54, 10.82, 1.36),
        (7, Sp2) => (12.87, 11.15, 0.85),
        (7, Sp) => (15.68, 11.7, -0.27),
        (8, Sp3) => (14.18, 12.92, 1.39),
        (8, Sp2) => (17.07, 13.79, 0.47),
        (9, Sp3) => (14.66, 13.85, 2.31),
        (17, Sp3) => (11.00, 9.69, 1.35),
        (35, Sp3) => (10.08, 8.47, 1.16),
        (53, Sp3) => (9.9, 7.96, 0.96),
        (16, Sp3) | (16, So) => (10.14, 9.13, 1.38),
        (16, So2) => (12.00, 10.81, 1.20),
        (16, Sp2) => (10.88, 9.49, 1.33),
        (15, Sp3) => (8.90, 8.24, 0.96),
        (15, Sp2) => (9.665, 8.530, 0.735),
        (14, Sp3) => (7.300, 6.567, 0.657),
        (14, Sp2) => (7.905, 6.748, 0.443),
        (14, Sp) => (9.065, 7.027, -0.002),
        (5, Sp3) => (5.980, 6.820, 1.605),
        (5, Sp2) => (6.420, 6.807, 1.322),
        (4, Sp3) => (3.845, 6.755, 3.165),
        (4, Sp2) => (4.005, 6.725, 3.035),
        (12, Sp2) => (3.565, 5.572, 2.197),
        (12, Sp3) => (3.300, 5.587, 2.447),
        (12, Sp) => (4.040, 5.472, 1.823),
        (13, Sp3) => (5.375, 4.953, 0.867),
        (13, Sp2) => (5.795, 5.020, 0.695),
        _ => (0.0, 0.0, 0.0),
    }
}

fn atomic_number(mol: &Molecule, idx: AtomIdx) -> u8 {
    let a = mol.atom(idx);
    if a.wildcard {
        0
    } else {
        a.element.atomic_number()
    }
}

/// Compute Gasteiger-Marsili PEOE partial charges, bit-identical to RDKit's
/// `ComputeGasteigerCharges` (`_GasteigerCharge`; 12 iterations).
///
/// Returns a `Vec<f64>` indexed by atom order (same as `mol.atoms()`).
/// Atoms without parameters use RDKit's all-zero fallback, which can make
/// charges non-finite exactly as in RDKit.
pub fn gasteiger_charges(mol: &Molecule) -> Vec<f64> {
    let n = mol.atom_count();
    if n == 0 {
        return Vec::new();
    }
    let view = crate::descriptors::descriptor_aromaticity(mol);
    let view: &Molecule = &view;

    // Each atom's bonds in RDKit's adjacency order (bond creation order).
    let mut atom_bonds: Vec<Vec<(usize, BondIdx)>> = vec![Vec::new(); n];
    for b in mol.rdkit_bond_order() {
        let bond = view.bond(b);
        let (i, j) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        atom_bonds[i].push((j, b));
        atom_bonds[j].push((i, b));
    }
    let z: Vec<u8> = (0..n)
        .map(|i| atomic_number(view, AtomIdx(i as u32)))
        .collect();

    let mut charges = vec![0.0f64; n];
    split_charge_conjugated(view, &atom_bonds, &z, &mut charges);

    let mut atm_ps: Vec<AbC> = Vec::with_capacity(n);
    let mut ion_x = vec![0.0f64; n];
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let mode = match chematic_smarts::rdkit_hybridization(view, idx) {
            Some(3) => Mode::Sp3,
            Some(2) => Mode::Sp2,
            Some(1) => Mode::Sp,
            _ if z[i] == 1 => Mode::Star,
            _ if z[i] == 16 => {
                let no = atom_bonds[i].iter().filter(|(nb, _)| z[*nb] == 8).count();
                match no {
                    2 => Mode::So2,
                    1 => Mode::So,
                    _ => Mode::Sp3,
                }
            }
            _ => Mode::Unknown,
        };
        let p = params(z[i], mode);
        atm_ps.push(p);
        ion_x[i] = if z[i] == 1 { IONXH } else { p.0 + p.1 + p.2 };
    }
    let h_params = params(1, Mode::Star);
    let n_hs: Vec<u32> = (0..n)
        .map(|i| u32::from(implicit_hcount(view, AtomIdx(i as u32))))
        .collect();

    let mut h_chrg = vec![0.0f64; n];
    let mut energ = vec![0.0f64; n];
    let mut damp = DAMP;
    for _ in 0..N_ITER {
        for aix in 0..n {
            let p = atm_ps[aix];
            energ[aix] = p.0 + charges[aix] * (p.1 + p.2 * charges[aix]);
        }
        for aix in 0..n {
            let mut dq = 0.0f64;
            for &(nb, _) in &atom_bonds[aix] {
                let dx = energ[nb] - energ[aix];
                let sgn = if dx < 0.0 { 0.0 } else { 1.0 };
                dq += dx / ((sgn * (ion_x[aix] - ion_x[nb])) + ion_x[nb]);
            }
            let ni_hs = n_hs[aix];
            if ni_hs > 0 {
                let ni = f64::from(ni_hs);
                let q_hs = h_chrg[aix] / ni;
                let enr = h_params.0 + q_hs * (h_params.1 + h_params.2 * q_hs);
                let dx = enr - energ[aix];
                let sgn = if dx < 0.0 { 0.0 } else { 1.0 };
                let dq_h = dx / ((sgn * (ion_x[aix] - IONXH)) + IONXH);
                dq += ni * dq_h;
                h_chrg[aix] -= ni * dq_h * damp;
            }
            charges[aix] += damp * dq;
        }
        damp *= DAMP_SCALE;
    }
    charges
}

/// RDKit `Gasteiger::splitChargeConjugated`: a formal charge on an atom is
/// shared equally with the same-element atoms two conjugated bonds away
/// (the two nitrogens of an amidinium start at +0.5 each).
fn split_charge_conjugated(
    view: &Molecule,
    atom_bonds: &[Vec<(usize, BondIdx)>],
    z: &[u8],
    charges: &mut [f64],
) {
    let n = charges.len();
    for aix in 0..n {
        let mut formal = f64::from(view.atom(AtomIdx(aix as u32)).charge);
        if !(formal.abs() > f64::EPSILON && charges[aix].abs() < f64::EPSILON) {
            continue;
        }
        let mut marker = vec![aix];
        for &(aax, b1) in &atom_bonds[aix] {
            if !chematic_smarts::rdkit_bond_is_conjugated(view, b1) {
                continue;
            }
            for &(yax, b2) in &atom_bonds[aax] {
                if b1 != b2
                    && chematic_smarts::rdkit_bond_is_conjugated(view, b2)
                    && z[aix] == z[yax]
                {
                    formal += f64::from(view.atom(AtomIdx(yax as u32)).charge);
                    marker.push(yax);
                }
            }
        }
        let share = formal / marker.len() as f64;
        for &m in &marker {
            charges[m] = share;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    /// Reference: RDKit 2026.03.1 `ComputeGasteigerCharges` (`_GasteigerCharge`),
    /// compared bit for bit. Covers charge splitting over a conjugated
    /// amidinium and carboxylate, sulfonyl `so2` mode and a nitro group.
    #[test]
    fn charges_match_rdkit_exactly() {
        let cases: [(&str, &[f64]); 4] = [
            (
                "C[NH+]=C(N)c1ccccc1",
                &[
                    0.07215551602567336,
                    -0.2771338393004865,
                    0.2720715966029391,
                    -0.2869310585025476,
                    0.06286713146056176,
                    -0.046773743793832666,
                    -0.061366422175682056,
                    -0.06221591245877897,
                    -0.061366422175682056,
                    -0.046773743793832666,
                ],
            ),
            (
                "CC(=O)[O-]",
                &[
                    -0.02496536893238246,
                    0.03827859892113979,
                    -0.5504770280073037,
                    -0.5504770280073037,
                ],
            ),
            (
                "c1ccsc1S(=O)(=O)N",
                &[
                    -0.031117254920057964,
                    -0.050102663706403726,
                    -0.007900599221896128,
                    -0.13100828577209192,
                    0.13414824857655913,
                    0.2470289930063984,
                    -0.20647805897111213,
                    -0.20647805897111213,
                    -0.22392952188652163,
                ],
            ),
            (
                "O=[N+]([O-])c1ccccc1",
                &[
                    -0.2583096384435056,
                    0.06218339391776295,
                    -0.2583096384435056,
                    0.2689216270439583,
                    0.015874619206005438,
                    -0.05543595891821206,
                    -0.061769815628152806,
                    -0.05543595891821206,
                    0.015874619206005438,
                ],
            ),
        ];
        for (smiles, want) in cases {
            let got = gasteiger_charges(&parse(smiles).unwrap());
            assert_eq!(got, want, "{smiles}");
        }
    }

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
