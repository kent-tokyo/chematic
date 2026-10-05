//! RDKit's hybridization model for the SMARTS `^n` primitive.
//!
//! Follows RDKit 2026.03.6 `ConjugHybrid.cpp`: the number of orbitals is the
//! total degree plus lone pairs, and an atom with four orbitals, fewer than
//! four neighbours and a conjugated bond is SP2 (amide N, ester O,
//! carboxylate O⁻). Conjugation follows `setConjugation`: aromatic bonds,
//! plus `markConjAtomBonds` around each 2- or 3-coordinate candidate atom.
//! Plain terminal H atoms report no hybridization, as for H atoms added by
//! RDKit's `AddHs`.

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule, implicit_hcount};

/// Outer-shell electrons of a main-group element (`None` otherwise).
fn outer_electrons(z: u8) -> Option<i32> {
    let z = i32::from(z);
    Some(match z {
        1 => 1,
        2 => 2,
        3..=10 => z - 2,
        11..=18 => z - 10,
        19 | 37 | 55 | 87 => 1,
        20 | 38 | 56 | 88 => 2,
        31..=36 => z - 28,
        49..=54 => z - 46,
        81..=86 => z - 78,
        // Transition metals, as RDKit's periodic table counts them
        // (group 3–11 → 3–11, group 12 → 2).
        21..=29 => z - 18,
        39..=47 => z - 36,
        72..=79 => z - 68,
        57 => 3,
        30 | 48 | 80 => 2,
        _ => return None,
    })
}

/// RDKit's default (first) valence for the elements conjugation looks at.
fn default_valence(mol: &Molecule, idx: AtomIdx) -> i32 {
    let element = mol.atom(idx).element;
    element
        .normal_valences()
        .first()
        .map_or(-1, |&v| i32::from(v))
}

fn total_hs(mol: &Molecule, idx: AtomIdx) -> i32 {
    i32::from(implicit_hcount(mol, idx))
}

fn heavy_degree(mol: &Molecule, idx: AtomIdx) -> i32 {
    mol.neighbors(idx).count() as i32
}

/// RDKit `countAtomElec`: electrons available to a pi system (`-1`: none).
fn count_atom_elec(mol: &Molecule, idx: AtomIdx) -> i32 {
    let atom = mol.atom(idx);
    let dv = default_valence(mol, idx);
    if dv <= 1 {
        return -1;
    }
    let degree = heavy_degree(mol, idx) + total_hs(mol, idx);
    if degree > 3 {
        return -1;
    }
    let Some(nouter) = outer_electrons(atom.element.atomic_number()) else {
        return -1;
    };
    let nlp = (nouter - dv - i32::from(atom.charge)).max(0);
    let mut res = (dv - degree) + nlp;
    if res > 1 {
        let explicit_valence =
            crate::match_vf2::total_valence(mol, idx) as i32 - total_hs(mol, idx);
        if explicit_valence - heavy_degree(mol, idx) > 1 {
            res = 1;
        }
    }
    res
}

/// RDKit `isAtomConjugCand`.
fn conjugation_candidate(mol: &Molecule, idx: AtomIdx) -> bool {
    let z = mol.atom(idx).element.atomic_number();
    let nouter = outer_electrons(z).unwrap_or(0);
    // A group-16 atom beyond the first row only with one substituent
    // counting H (`C=C[S-]` conjugates, the thiol in `C=CS` does not).
    (z <= 10 || (nouter != 5 && nouter != 6) || (nouter == 6 && substituents(mol, idx) < 2))
        && count_atom_elec(mol, idx) > 0
}

fn valence_contrib(order: BondOrder) -> f32 {
    match order {
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Aromatic => 1.5,
        BondOrder::Quadruple => 4.0,
        _ => 1.0,
    }
}

fn substituents(mol: &Molecule, idx: AtomIdx) -> i32 {
    heavy_degree(mol, idx) + total_hs(mol, idx)
}

/// Whether `markConjAtomBonds` run at `center` marks `bond`.
fn marked_at(mol: &Molecule, center: AtomIdx, bond: BondIdx) -> bool {
    if !conjugation_candidate(mol, center) {
        return false;
    }
    let sbo = substituents(mol, center);
    if !(2..=3).contains(&sbo) {
        return false;
    }
    let other_ok = |b: BondIdx| {
        let e = mol.bond(b);
        let other = if e.atom1 == center { e.atom2 } else { e.atom1 };
        substituents(mol, other) <= 3 && conjugation_candidate(mol, other)
    };
    let multiple = |b: BondIdx| valence_contrib(mol.bond(b).order) >= 1.5;
    let bonds: Vec<BondIdx> = mol.neighbors(center).map(|(_, b)| b).collect();
    // `bond` as the multiple bond, paired with any other candidate bond …
    (multiple(bond) && bonds.iter().any(|&b2| b2 != bond && other_ok(b2)))
        // … or as the partner of another multiple bond.
        || (other_ok(bond) && bonds.iter().any(|&b1| b1 != bond && multiple(b1)))
}

fn bond_is_conjugated(mol: &Molecule, bond: BondIdx) -> bool {
    let e = mol.bond(bond);
    e.order == BondOrder::Aromatic || marked_at(mol, e.atom1, bond) || marked_at(mol, e.atom2, bond)
}

/// RDKit hybridization as the SMARTS `^n` code (0 = S, 1 = SP, 2 = SP2,
/// 3 = SP3, 4 = SP3D, 5 = SP3D2); `None` = unspecified.
///
/// Follows RDKit's `ConjugHybrid` model: orbital count = total degree plus
/// lone pairs, and a four-orbital atom with fewer than four neighbours and
/// a conjugated bond (amide N, aryl ether O, enamine N) is SP2.
pub fn rdkit_hybridization(mol: &Molecule, idx: AtomIdx) -> Option<u8> {
    let atom = mol.atom(idx);
    let z = atom.element.atomic_number();
    if atom.wildcard || z == 0 {
        return None;
    }
    // A plain terminal H is what RDKit's `AddHs` adds, with no hybridization;
    // an H RDKit keeps from SMILES (`[2H]`, `[H+]`, `[H][H]`) is S.
    if z == 1 {
        let plain = atom.isotope.is_none() && atom.charge == 0 && mol.degree(idx) == 1 && {
            let (nb, _) = mol.neighbors(idx).next().unwrap();
            mol.atom(nb).element.atomic_number() != 1
        };
        return (!plain).then_some(0);
    }
    let total_degree = heavy_degree(mol, idx) + total_hs(mol, idx);
    // An aromatic atom with two or three connections has three orbitals or
    // a conjugated lone pair: SP2 in RDKit's model. Answering here skips the
    // Kekulé valence lookup the general rule needs for aromatic atoms.
    if atom.aromatic && (2..=3).contains(&total_degree) {
        return Some(2);
    }
    let norbs = match outer_electrons(z) {
        Some(nouter) if z < 89 => {
            let total_valence = crate::match_vf2::total_valence(mol, idx) as i32;
            // Not clamped: RDKit lets a negative lone-pair count lower the
            // orbital count (`[Zn++]` with four bonds is SP).
            let free = nouter - (total_valence + i32::from(atom.charge));
            total_degree + free / 2
        }
        _ => total_degree,
    };
    match norbs {
        0 | 1 => Some(0),
        2 => Some(1),
        3 => Some(2),
        4 => {
            let conjugated = mol.neighbors(idx).any(|(_, b)| bond_is_conjugated(mol, b));
            Some(if total_degree < 4 && conjugated { 2 } else { 3 })
        }
        5 => Some(4),
        6 => Some(5),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hybridization_matches_rdkit() {
        // (SMILES, per-atom RDKit 2026.03.6 `^n` code; `-` = unspecified).
        for (smiles, want) in [
            ("CC(=O)[O-]", "3222"),
            ("CC(=O)NC", "32223" /* amide N is SP2 */),
            ("COc1ccccc1", "32222222"),
            ("CC(S)=N", "3232"),
            ("C=C[S-]", "222"),
            ("C=CCl", "223"),
            ("CC#N", "311"),
            ("C[N+](C)(C)C", "33333"),
            ("[2H]C", "03"),
        ] {
            let mol = chematic_smiles::parse(smiles).unwrap();
            let got: String = (0..mol.atom_count())
                .map(|i| match rdkit_hybridization(&mol, AtomIdx(i as u32)) {
                    Some(h) => char::from(b'0' + h),
                    None => '-',
                })
                .collect();
            assert_eq!(got, want, "{smiles}");
        }
    }
}
