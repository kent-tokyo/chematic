//! `rdkit.Chem.EnumerateStereoisomers.EnumerateStereoisomers` with its
//! default options (`onlyUnassigned=True`, `unique=True`,
//! `maxIsomers=1024`, `tryEmbedding=False`) on the RDKit-model molecule.
//!
//! RDKit flips every unspecified centre and double bond that
//! `FindPotentialStereo` reports, re-runs the (legacy) stereo perception on
//! each isomer and keeps the distinct canonical SMILES. Here the flipped
//! centres are RDKit's own candidates (`isAtomPotentialTetrahedralCenter`,
//! stereo-capable double bonds) that the legacy perception keeps when they
//! are all specified; flipping a candidate it would drop only produces
//! duplicates, which uniqueness removes.

use super::RdkitSmilesError;
use super::mol::{BondStereo, BondType, ChiralTag, Mol};
use super::stereo::{
    is_atom_potential_tetrahedral_center, legacy_stereo_perception,
    should_detect_double_bond_stereo,
};
use super::write::mol_to_smiles;

/// Unspecified candidate centres and double bonds of `m` (after perception).
fn candidates(m: &Mol) -> (Vec<usize>, Vec<(usize, usize, usize)>) {
    let atoms: Vec<usize> = (0..m.atoms.len())
        .filter(|&a| {
            m.atoms[a].chiral == ChiralTag::Unspecified
                && is_atom_potential_tetrahedral_center(m, a)
        })
        .collect();
    let mut bonds = Vec::new();
    for b in 0..m.bonds.len() {
        let bond = &m.bonds[b];
        if bond.bt != BondType::Double
            || bond.aromatic
            || !matches!(bond.stereo, BondStereo::None | BondStereo::Any)
            || !should_detect_double_bond_stereo(m, b)
        {
            continue;
        }
        let (beg, end) = (bond.begin, bond.end);
        if !matches!(m.degree(beg), 2 | 3) || !matches!(m.degree(end), 2 | 3) {
            continue;
        }
        let first = |x: usize| {
            m.atom_bonds[x]
                .iter()
                .find(|&&nb| nb != b)
                .map(|&nb| m.bonds[nb].other(x))
        };
        if let (Some(sa), Some(sb)) = (first(beg), first(end)) {
            bonds.push((b, sa, sb));
        }
    }
    (atoms, bonds)
}

/// Apply one flip assignment: bit `i` of `bits` set = `CW` / cis.
fn assign(
    base: &Mol,
    atoms: &[usize],
    bonds: &[(usize, usize, usize)],
    bits: impl Fn(usize) -> bool,
) -> Mol {
    let mut m = base.clone();
    for (i, &a) in atoms.iter().enumerate() {
        m.atoms[a].chiral = if bits(i) {
            ChiralTag::Cw
        } else {
            ChiralTag::Ccw
        };
        m.atoms[a].cip_code = None;
    }
    for (j, &(b, sa, sb)) in bonds.iter().enumerate() {
        m.bonds[b].requested = Some((sa, sb, !bits(atoms.len() + j)));
    }
    for a in &mut m.atoms {
        a.cip_code = None;
    }
    legacy_stereo_perception(&mut m, true, true);
    m
}

/// The distinct canonical SMILES `EnumerateStereoisomers` yields for the
/// sanitized, stereo-perceived `base`, or `Err` where RDKit would sample at
/// random (more than `max_isomers` flip combinations).
pub(crate) fn enumerate(base: &Mol, max_isomers: usize) -> Result<Vec<String>, RdkitSmilesError> {
    let (atoms, bonds) = candidates(base);
    let n = atoms.len() + bonds.len();
    if n == 0 {
        return Ok(vec![mol_to_smiles(base)?]);
    }
    // Keep the candidates the perception retains under a few assignments
    // (a pseudo-asymmetric centre survives only for some of them).
    let mut keep_atom = vec![false; atoms.len()];
    let mut keep_bond = vec![false; bonds.len()];
    let probes: [fn(usize) -> bool; 4] = [|_| true, |_| false, |i| i % 2 == 0, |i| i % 3 == 0];
    for probe in probes {
        let m = assign(base, &atoms, &bonds, probe);
        for (i, &a) in atoms.iter().enumerate() {
            keep_atom[i] |= m.atoms[a].chiral != ChiralTag::Unspecified;
        }
        for (j, &(b, ..)) in bonds.iter().enumerate() {
            keep_bond[j] |= matches!(m.bonds[b].stereo, BondStereo::E | BondStereo::Z);
        }
    }
    let atoms: Vec<usize> = atoms
        .iter()
        .zip(&keep_atom)
        .filter(|&(_, &k)| k)
        .map(|(&a, _)| a)
        .collect();
    let bonds: Vec<(usize, usize, usize)> = bonds
        .iter()
        .zip(&keep_bond)
        .filter(|&(_, &k)| k)
        .map(|(&b, _)| b)
        .collect();
    let n = atoms.len() + bonds.len();
    if n == 0 {
        let m = assign(base, &[], &[], |_| false);
        return Ok(vec![mol_to_smiles(&m)?]);
    }
    // RDKit enumerates every flip combination when there are at most
    // `max_isomers`; otherwise it draws combinations at random until it has
    // `max_isomers` distinct isomers or has seen every combination. Either
    // way it yields every distinct isomer when there are at most
    // `max_isomers` of them, which full enumeration reproduces; beyond that
    // its random sample is not reproduced.
    const MAX_FULL_ENUMERATION_FLIPS: usize = 16;
    if n > MAX_FULL_ENUMERATION_FLIPS {
        return Err(RdkitSmilesError::Unsupported(format!(
            "{n} stereo flips: RDKit samples {max_isomers} of 2^{n} at random"
        )));
    }
    let mut seen = std::collections::BTreeSet::new();
    for bitflag in 0..(1usize << n) {
        let m = assign(base, &atoms, &bonds, |i| bitflag & (1 << i) != 0);
        seen.insert(mol_to_smiles(&m)?);
    }
    if seen.len() > max_isomers {
        return Err(RdkitSmilesError::Unsupported(format!(
            "{} distinct stereoisomers: RDKit yields a random sample of {max_isomers}",
            seen.len()
        )));
    }
    Ok(seen.into_iter().collect())
}

#[cfg(test)]
mod tests {
    use crate::rdkit_stereoisomer_smiles;

    fn isomers(smiles: &str) -> Vec<String> {
        rdkit_stereoisomer_smiles(&crate::parse(smiles).unwrap(), 1024).unwrap()
    }

    // The examples of RDKit 2026.03.1's `EnumerateStereoisomers` docstring.
    #[test]
    fn rdkit_docstring_examples() {
        let all = isomers("BrC=CC1OC(C2)(F)C2(Cl)C1");
        assert_eq!(all.len(), 16);
        assert_eq!(all[0], "F[C@@]12C[C@@]1(Cl)C[C@@H](/C=C/Br)O2");
        assert_eq!(all[15], "F[C@]12C[C@]1(Cl)C[C@H](/C=C\\Br)O2");
        assert_eq!(
            isomers("FC(Cl)C=CC=CC(F)Cl"),
            [
                "F[C@@H](Cl)/C=C/C=C/[C@@H](F)Cl",
                "F[C@@H](Cl)/C=C\\C=C/[C@@H](F)Cl",
                "F[C@@H](Cl)/C=C\\C=C\\[C@@H](F)Cl",
                "F[C@H](Cl)/C=C/C=C/[C@@H](F)Cl",
                "F[C@H](Cl)/C=C/C=C/[C@H](F)Cl",
                "F[C@H](Cl)/C=C/C=C\\[C@@H](F)Cl",
                "F[C@H](Cl)/C=C\\C=C/[C@@H](F)Cl",
                "F[C@H](Cl)/C=C\\C=C/[C@H](F)Cl",
                "F[C@H](Cl)/C=C\\C=C\\[C@@H](F)Cl",
                "F[C@H](Cl)/C=C\\C=C\\[C@H](F)Cl",
            ]
        );
        // Only unassigned centres are expanded.
        let some = isomers("BrC=C[C@H]1OC(C2)(F)C2(Cl)C1");
        assert_eq!(some.len(), 8);
        assert!(some.iter().all(|s| s.contains("[C@@H](/C=C")));
    }

    #[test]
    fn salts_and_conjugated_bonds_are_enumerated() {
        assert_eq!(isomers("CC=CC.C"), ["C.C/C=C/C", "C.C/C=C\\C"]);
        assert_eq!(isomers("F/C=C/C=CC"), ["C/C=C/C=C/F", "C/C=C\\C=C\\F"]);
        assert_eq!(isomers("CCO"), ["CCO"]);
    }

    #[test]
    fn random_sampling_is_refused() {
        let mol = crate::parse(&format!("Br{}F", "[CH](Cl)".repeat(20))).unwrap();
        assert!(rdkit_stereoisomer_smiles(&mol, 1024).is_err());
    }
}
