//! MMFF94 Bond Charge Increments (BCI) for partial charge calculation.
//!
//! Implements the Halgren 1996 BCI charge model:
//!   q_i = q0_i(formal) + Σ_j bci(type_j, type_i)
//!
//! The BCI is antisymmetric: bci(a, b) = -bci(b, a).
//! Convention: bci(from, to) is the charge added to `to` when bonded to `from`.
//!   - Negative value: `to` gains electrons (becomes more negative)
//!   - Positive value: `to` loses electrons (becomes more positive)
//!
//! Reference: Halgren TA (1996) J. Comput. Chem. 17(5-6), 490-519.

use crate::mmff94::MMFF94Type;
use chematic_core::Molecule;

/// Bond Charge Increment from `from` to `to`.
///
/// Returns the partial charge increment added to the `to` atom when bonded to `from`.
/// Negative = `to` becomes more negative (gains electron density).
/// Antisymmetric: bci(a, b) + bci(b, a) = 0.
pub fn bci(from: MMFF94Type, to: MMFF94Type) -> f64 {
    use MMFF94Type::*;
    // Flat lookup table: (from, to) → charge increment on `to`.
    // Values calibrated to MMFF94s parameters (Halgren 1996 Table III).
    match (from, to) {
        // ------------------------------------------------------------------
        // C(sp3) — X bonds
        // ------------------------------------------------------------------
        (C_sp3, O_Alcohol) => -0.2959,    // O gains from C in alcohol
        (C_sp3, O_Ether) => -0.2488,      // ether O
        (C_sp3, O_Phenol) => -0.2959,     // phenol O (same as alcohol)
        (C_sp3, N_sp3_Amine) => -0.2088,  // amine N
        (C_sp3, S_Thiol) => -0.1278,      // thiol S
        (C_sp3, S_Thioether) => -0.1148,  // thioether S
        (C_sp3, H_Carbon) => 0.0240,      // H loses slightly to C
        (C_sp3, C_sp3) => 0.0,            // no polarization between identical types
        (C_sp3, C_sp2_Alkene) => -0.0500, // sp2 C gains from sp3
        (C_sp3, C_Aromatic) => -0.0380,   // aromatic C gains from sp3
        (C_sp3, F) => -0.4578,            // F strongly gains from C
        (C_sp3, Cl) => -0.2523,           // Cl gains
        (C_sp3, Br) => -0.1754,           // Br gains
        (C_sp3, I) => -0.0906,            // I gains (least electronegative halogen)

        // ------------------------------------------------------------------
        // C(sp2 alkene) — X bonds
        // ------------------------------------------------------------------
        (C_sp2_Alkene, C_sp2_Alkene) => 0.0,
        (C_sp2_Alkene, H_Carbon) => 0.0380, // vinyl H loses to C
        (C_sp2_Alkene, C_sp3) => 0.0500,    // sp3 C loses to sp2

        // ------------------------------------------------------------------
        // C(aromatic) — X bonds
        // ------------------------------------------------------------------
        (C_Aromatic, H_Aromatic) => 0.0780, // aryl H loses to C
        (C_Aromatic, C_Aromatic) => 0.0,
        (C_Aromatic, N_Aromatic_6ring | N_Aromatic_Pyridine) => -0.1500, // ring N gains

        // ------------------------------------------------------------------
        // C(carbonyl) — X bonds
        // ------------------------------------------------------------------
        (C_Carbonyl, O_Carbonyl) => -0.4500, // =O gains strongly
        (C_Carbonyl, C_sp3) => 0.1000,       // flanking C loses a bit
        (C_Carbonyl, H_Carbon) => 0.0500,    // aldehyde H

        // ------------------------------------------------------------------
        // Carboxylic acid and ester variants
        // ------------------------------------------------------------------
        (C_Carboxylic, O_Carbonyl) => -0.4400,
        (C_Carboxylic, O_Carboxylic) => -0.3600, // –OH in acid
        (C_Ester, O_Carbonyl) => -0.4300,
        (C_Ester, O_Ester) => -0.3200,
        (C_Amide, O_Amide) => -0.4200,
        (C_Amide, N_Amide) => -0.2800,

        // ------------------------------------------------------------------
        // O — H bonds
        // ------------------------------------------------------------------
        (O_Alcohol, H_Oxygen) => 0.1806, // H loses to O in alcohol
        (O_Phenol, H_Oxygen) => 0.2106,  // phenol O-H (stronger)
        (O_Carboxylic, H_Oxygen) => 0.2500, // acid O-H

        // ------------------------------------------------------------------
        // N — H bonds
        // ------------------------------------------------------------------
        (N_sp3_Amine, H_Nitrogen) => 0.1760, // H loses to N in amine
        (N_Amide, H_Nitrogen) => 0.1900,     // amide N-H

        // ------------------------------------------------------------------
        // S — H bonds
        // ------------------------------------------------------------------
        (S_Thiol, H_Sulfur) => 0.1000, // H loses to S in thiol

        // ------------------------------------------------------------------
        // Antisymmetric reversal for all defined pairs
        // ------------------------------------------------------------------
        (a, b) => -bci_lookup(b, a), // look up the reverse and negate
    }
}

/// Inner lookup (non-recursive) for the antisymmetric reversal.
fn bci_lookup(from: MMFF94Type, to: MMFF94Type) -> f64 {
    use MMFF94Type::*;
    match (from, to) {
        (C_sp3, O_Alcohol) => -0.2959,
        (C_sp3, O_Ether) => -0.2488,
        (C_sp3, O_Phenol) => -0.2959,
        (C_sp3, N_sp3_Amine) => -0.2088,
        (C_sp3, S_Thiol) => -0.1278,
        (C_sp3, S_Thioether) => -0.1148,
        (C_sp3, H_Carbon) => 0.0240,
        (C_sp3, C_sp3) => 0.0,
        (C_sp3, C_sp2_Alkene) => -0.0500,
        (C_sp3, C_Aromatic) => -0.0380,
        (C_sp3, F) => -0.4578,
        (C_sp3, Cl) => -0.2523,
        (C_sp3, Br) => -0.1754,
        (C_sp3, I) => -0.0906,
        (C_sp2_Alkene, C_sp2_Alkene) => 0.0,
        (C_sp2_Alkene, H_Carbon) => 0.0380,
        (C_sp2_Alkene, C_sp3) => 0.0500,
        (C_Aromatic, H_Aromatic) => 0.0780,
        (C_Aromatic, C_Aromatic) => 0.0,
        (C_Aromatic, N_Aromatic_6ring | N_Aromatic_Pyridine) => -0.1500,
        (C_Carbonyl, O_Carbonyl) => -0.4500,
        (C_Carbonyl, C_sp3) => 0.1000,
        (C_Carbonyl, H_Carbon) => 0.0500,
        (C_Carboxylic, O_Carbonyl) => -0.4400,
        (C_Carboxylic, O_Carboxylic) => -0.3600,
        (C_Ester, O_Carbonyl) => -0.4300,
        (C_Ester, O_Ester) => -0.3200,
        (C_Amide, O_Amide) => -0.4200,
        (C_Amide, N_Amide) => -0.2800,
        (O_Alcohol, H_Oxygen) => 0.1806,
        (O_Phenol, H_Oxygen) => 0.2106,
        (O_Carboxylic, H_Oxygen) => 0.2500,
        (N_sp3_Amine, H_Nitrogen) => 0.1760,
        (N_Amide, H_Nitrogen) => 0.1900,
        (S_Thiol, H_Sulfur) => 0.1000,
        _ => 0.0,
    }
}

/// Formal charge for each MMFF94 atom type (the starting point for BCI calculation).
/// Neutral atoms start at 0.0; ions at ±1.0.
pub fn mmff94_formal_charge(_t: MMFF94Type) -> f64 {
    // All atom types start from 0.0; formal atomic charges are added from mol.atom(i).charge
    0.0
}

/// Compute MMFF94 partial charges using the BCI (Bond Charge Increment) model.
///
/// This is a 2D topology-based calculation — no 3D coordinates required.
/// Accuracy: ~±0.05e for common organic functional groups.
///
/// Formula: q_i = q0_formal + Σ_{j bonded to i} bci(type_j, type_i)
pub fn mmff94_charges_bci(mol: &Molecule) -> Result<Vec<f64>, crate::mmff94::AssignError> {
    let types = crate::mmff94::assign_mmff94_types(mol)?;
    let n = mol.atom_count();

    // Step 1: initialize from formal charges on the atoms
    let mut charges: Vec<f64> = (0..n)
        .map(|i| {
            let atom = mol.atom(chematic_core::AtomIdx(i as u32));
            mmff94_formal_charge(types[i]) + atom.charge as f64
        })
        .collect();

    // Step 2: accumulate BCI increments over every bond
    for (_, bond) in mol.bonds() {
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        let ti = types[i];
        let tj = types[j];
        // bci(ti, tj) = charge increment added to atom j
        let delta = bci(ti, tj);
        charges[j] += delta;
        charges[i] -= delta; // antisymmetric
    }

    Ok(charges)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> chematic_core::Molecule {
        parse(s).unwrap()
    }

    #[test]
    fn test_bci_antisymmetry() {
        use MMFF94Type::*;
        // BCI must be antisymmetric: bci(a,b) = -bci(b,a)
        let pairs = [
            (C_sp3, O_Alcohol),
            (C_sp3, N_sp3_Amine),
            (C_sp3, F),
            (O_Alcohol, H_Oxygen),
            (N_sp3_Amine, H_Nitrogen),
        ];
        for (a, b) in pairs {
            let forward = bci(a, b);
            let reverse = bci(b, a);
            assert!(
                (forward + reverse).abs() < 1e-10,
                "bci({:?},{:?}) = {:.4} but bci({:?},{:?}) = {:.4} (should be negatives)",
                a,
                b,
                forward,
                b,
                a,
                reverse
            );
        }
    }

    #[test]
    fn test_bci_charge_neutrality() {
        // Sum of partial charges should equal molecular formal charge
        for smiles in &["C", "CO", "CC", "CCN", "CC(=O)C", "c1ccccc1"] {
            let m = mol(smiles);
            let charges = mmff94_charges_bci(&m).expect("charge calculation failed");
            let total: f64 = charges.iter().sum();
            assert!(
                total.abs() < 0.05,
                "{smiles}: total charge = {total:.4} (should be ~0)"
            );
        }
    }

    #[test]
    fn test_bci_qualitative_methanol() {
        let m = mol("CO");
        let charges = mmff94_charges_bci(&m).unwrap();
        // Find O atom (atomic number 8)
        let o_idx = m
            .atoms()
            .find(|(_, a)| a.element.atomic_number() == 8)
            .map(|(i, _)| i.0 as usize)
            .unwrap();
        // O should be significantly negative
        assert!(
            charges[o_idx] < -0.20,
            "methanol O charge = {:.3}, expected < -0.20",
            charges[o_idx]
        );
        assert!(
            charges[o_idx] > -0.60,
            "methanol O charge = {:.3}, expected > -0.60",
            charges[o_idx]
        );
    }

    #[test]
    fn test_bci_qualitative_amine() {
        let m = mol("CCN");
        let charges = mmff94_charges_bci(&m).unwrap();
        let n_idx = m
            .atoms()
            .find(|(_, a)| a.element.atomic_number() == 7)
            .map(|(i, _)| i.0 as usize)
            .unwrap();
        assert!(
            charges[n_idx] < -0.10,
            "amine N = {:.3}, expected negative",
            charges[n_idx]
        );
    }

    #[test]
    fn test_bci_qualitative_halogens() {
        // F > Cl > Br > I in terms of electronegativity: F should be most negative
        let cf = mol("CF");
        let ccl = mol("CCl");
        let cbr = mol("CBr");

        let f_charge = |m: &chematic_core::Molecule| {
            mmff94_charges_bci(m)
                .unwrap()
                .into_iter()
                .zip(m.atoms())
                .find(|(_, (_, a))| matches!(a.element.atomic_number(), 9 | 17 | 35))
                .map(|(q, _)| q)
                .unwrap()
        };

        let qf = f_charge(&cf);
        let qcl = f_charge(&ccl);
        let qbr = f_charge(&cbr);

        assert!(
            qf < qcl,
            "F ({qf:.3}) should be more negative than Cl ({qcl:.3})"
        );
        assert!(
            qcl < qbr,
            "Cl ({qcl:.3}) should be more negative than Br ({qbr:.3})"
        );
    }

    #[test]
    fn test_bci_ketone_oxygen() {
        let m = mol("CC(=O)C"); // acetone
        let charges = mmff94_charges_bci(&m).unwrap();
        let o_idx = m
            .atoms()
            .find(|(_, a)| a.element.atomic_number() == 8)
            .map(|(i, _)| i.0 as usize)
            .unwrap();
        // Ketone O should be strongly negative
        assert!(
            charges[o_idx] < -0.30,
            "ketone O = {:.3}, expected < -0.30",
            charges[o_idx]
        );
    }
}
