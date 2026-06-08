//! XLogP3 lipophilicity prediction.
//!
//! Implements the atom-contribution model from:
//!   Cheng T, Zhao Y, Li X, Lin F, Xu Y, Zhang X, Li Y, Wang R, Lai L.
//!   "Computation of Octanol-Water Partition Coefficients by Guiding an Additive
//!    Model with Knowledge." J. Chem. Inf. Model. 2007, 47, 2140–2148.
//!
//! Atom type contributions are summed over all heavy atoms.  Correction factors
//! (hydrophobic pairs, internal H-bonds) are not applied in this implementation;
//! for most drug-like molecules the atom-type sum alone gives a good approximation.
//!
//! For a comparable but different method see [`crate::logp_crippen`].

use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};
use chematic_perception::find_sssr;

/// Predict the octanol/water logP using the XLogP3 atom-contribution model.
///
/// Returns the sum of per-atom contributions; see [`xlogp3_per_atom`] for
/// individual values.
pub fn xlogp3(mol: &Molecule) -> f64 {
    xlogp3_per_atom(mol).iter().sum()
}

/// Return the XLogP3 contribution of each atom in `mol` (one value per atom,
/// indexed by atom position).
///
/// Explicit hydrogen atoms in the molecule contribute `0.1230` each.
/// All other element types not covered by the model return `0.0`.
pub fn xlogp3_per_atom(mol: &Molecule) -> Vec<f64> {
    let rings = find_sssr(mol);

    mol.atoms()
        .map(|(idx, atom)| {
            let an = atom.element.atomic_number();
            let h = implicit_hcount(mol, idx);
            let arom = atom.aromatic;
            let in_ring = rings.contains_atom(idx);
            let charge = atom.charge;

            let has_double = mol
                .neighbors(idx)
                .any(|(_, bi)| mol.bond(bi).order == BondOrder::Double);
            let has_triple = mol
                .neighbors(idx)
                .any(|(_, bi)| mol.bond(bi).order == BondOrder::Triple);

            match an {
                6 => carbon(mol, idx, h, arom, in_ring, has_double, has_triple),
                7 => nitrogen(mol, idx, h, arom, in_ring, has_double, charge),
                8 => oxygen(mol, idx, h, arom, in_ring, has_double),
                16 => sulfur(mol, idx, arom, in_ring, has_double),
                15 => phosphorus(mol, idx),
                9 => 0.4202,  // F
                17 => 0.6437, // Cl
                35 => 0.6923, // Br
                53 => 0.7937, // I
                1 => 0.1230,  // explicit H
                _ => 0.0,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Per-element helpers
// ---------------------------------------------------------------------------

fn carbon(
    mol: &Molecule,
    idx: AtomIdx,
    h: u8,
    arom: bool,
    in_ring: bool,
    has_double: bool,
    has_triple: bool,
) -> f64 {
    if has_triple {
        return if h > 0 { 0.2100 } else { 0.1537 };
    }
    if arom {
        if h > 0 {
            return 0.3395; // aromatic CH
        }
        // Aromatic C without H: bonded to heteroatom → higher contribution
        let bonded_heteroatom = mol.neighbors(idx).any(|(nb, _)| {
            let a = mol.atom(nb).element.atomic_number();
            a != 6 && a != 1
        });
        return if bonded_heteroatom { 0.2445 } else { 0.1840 };
    }
    if has_double {
        // sp2 carbon
        return if in_ring {
            if h > 0 { 0.3015 } else { 0.1255 }
        } else {
            match h {
                2.. => 0.3035,
                1 => 0.2756,
                _ => 0.0918,
            }
        };
    }
    // sp3 carbon
    if in_ring {
        match h {
            0 => -0.0460,
            1 => 0.0460,
            _ => 0.1059,
        }
    } else {
        match h {
            0 => 0.0516,
            1 => 0.1551,
            2 => 0.3741,
            _ => 0.5441,
        }
    }
}

fn nitrogen(
    mol: &Molecule,
    idx: AtomIdx,
    h: u8,
    arom: bool,
    in_ring: bool,
    has_double: bool,
    charge: i8,
) -> f64 {
    // Charged nitrogen (ammonium, etc.)
    if charge > 0 {
        return -2.1789;
    }
    // Nitro / nitroso: N bonded to two O (one double, one single or both single)
    let o_neighbor_count = mol
        .neighbors(idx)
        .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() == 8)
        .count();
    if o_neighbor_count >= 2 {
        return -1.0379;
    }
    if arom {
        return if h > 0 { -0.9139 } else { -0.3239 };
    }
    if has_double {
        return if in_ring {
            -0.3239 // pyridine-like
        } else if h > 0 {
            -0.7322 // imino NH
        } else {
            -0.3456 // imino N
        };
    }
    // sp3 nitrogen
    if in_ring {
        -0.5247
    } else {
        match h {
            2.. => -0.9923,
            1 => -0.5893,
            _ => -0.8406,
        }
    }
}

fn oxygen(
    _mol: &Molecule,
    _idx: AtomIdx,
    h: u8,
    arom: bool,
    _in_ring: bool,
    has_double: bool,
) -> f64 {
    if arom {
        return 0.1105; // furan-type aromatic O
    }
    if has_double {
        return -0.0908; // carbonyl =O
    }
    if h > 0 {
        return -0.2677; // alcohol/phenol OH
    }
    0.0623 // ether O (ring or non-ring)
}

fn sulfur(mol: &Molecule, idx: AtomIdx, arom: bool, _in_ring: bool, has_double: bool) -> f64 {
    if arom {
        return 0.3765;
    }
    // Sulfoxide: S with exactly one =O
    let double_o = mol
        .neighbors(idx)
        .filter(|(nb, bi)| {
            mol.atom(*nb).element.atomic_number() == 8 && mol.bond(*bi).order == BondOrder::Double
        })
        .count();
    if double_o == 1 {
        return -0.0093; // sulfoxide
    }
    // Sulfone: S with two =O
    if double_o >= 2 {
        return -0.1441;
    }
    let h = implicit_hcount(mol, idx);
    if has_double {
        return 0.3649; // thione =S
    }
    if h > 0 {
        return 0.6482; // thiol -SH
    }
    0.3649 // thioether (ring or chain)
}

fn phosphorus(mol: &Molecule, idx: AtomIdx) -> f64 {
    // P=O context
    let has_oxo = mol.neighbors(idx).any(|(nb, bi)| {
        mol.atom(nb).element.atomic_number() == 8 && mol.bond(bi).order == BondOrder::Double
    });
    let _has_oxo = has_oxo;
    0.1543
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::parse;

    fn mol(s: &str) -> Molecule {
        parse(s).unwrap()
    }

    #[test]
    fn test_xlogp3_methane_positive() {
        // Methane should have positive logP (lipophilic carbon)
        let v = xlogp3(&mol("C"));
        assert!(v > 0.0, "methane xlogp3 should be positive, got {v}");
    }

    #[test]
    fn test_xlogp3_water_negative() {
        // Water (single O) should be negative (hydrophilic)
        let v = xlogp3(&mol("O"));
        assert!(v < 0.0, "water xlogp3 should be negative, got {v}");
    }

    #[test]
    fn test_xlogp3_benzene_positive() {
        let v = xlogp3(&mol("c1ccccc1"));
        assert!(v > 0.0, "benzene xlogp3 should be positive, got {v}");
    }

    #[test]
    fn test_xlogp3_ethanol_less_than_hexane() {
        // Ethanol is more hydrophilic than hexane
        let ethanol = xlogp3(&mol("CCO"));
        let hexane = xlogp3(&mol("CCCCCC"));
        assert!(
            ethanol < hexane,
            "ethanol ({ethanol:.3}) should be less lipophilic than hexane ({hexane:.3})"
        );
    }

    #[test]
    fn test_xlogp3_per_atom_len() {
        let m = mol("CC(=O)O"); // acetic acid, 4 heavy atoms
        assert_eq!(xlogp3_per_atom(&m).len(), 4);
    }
}
