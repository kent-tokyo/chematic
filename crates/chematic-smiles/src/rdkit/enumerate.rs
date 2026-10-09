//! `rdkit.Chem.EnumerateStereoisomers.EnumerateStereoisomers` with its
//! default options (`onlyUnassigned=True`, `unique=True`,
//! `maxIsomers=1024`, `tryEmbedding=False`) on the RDKit-model molecule.
//!
//! RDKit flips every unspecified centre and double bond that
//! `FindPotentialStereo` reports (ported in `findstereo`), re-runs the
//! (legacy) stereo perception on each isomer and keeps the distinct
//! canonical SMILES.

use super::findstereo::{NOATOM, Specified, StereoType, find_potential_stereo};
use super::mol::{BondStereo, ChiralTag, Mol};
use super::pyrandom::{PyRandom, hash_pair_tuple};
use super::stereo::legacy_stereo_perception;
use super::write::mol_to_smiles;
use super::{RdkitSmilesError, RdkitSmilesParams};

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
    // `_BondFlipper.flip`: STEREOCIS / STEREOTRANS relative to the stereo
    // atoms; then `Chem.SetDoubleBondNeighborDirections(isomer)`, whose
    // directions (one per single bond, so conjugated double bonds sharing
    // one can override each other) are what the perception below reads.
    for (j, &(b, sa, sb)) in bonds.iter().enumerate() {
        m.bonds[b].stereo = if bits(atoms.len() + j) {
            BondStereo::Z
        } else {
            BondStereo::E
        };
        m.bonds[b].stereo_atoms = vec![sa, sb];
    }
    if !bonds.is_empty()
        && super::mol2_read::set_double_bond_neighbor_directions(&mut m, None).is_err()
    {
        // Fall back to the requested configurations.
        let mut fallback = base.clone();
        for (j, &(b, sa, sb)) in bonds.iter().enumerate() {
            fallback.bonds[b].requested = Some((sa, sb, !bits(atoms.len() + j)));
        }
        for (i, &a) in atoms.iter().enumerate() {
            fallback.atoms[a].chiral = m.atoms[a].chiral;
            let _ = i;
        }
        m = fallback;
    }
    for a in &mut m.atoms {
        a.cip_code = None;
    }
    legacy_stereo_perception(&mut m, true, true);
    m
}

/// The centres and double bonds `EnumerateStereoisomers` flips
/// (`_getFlippers` with `onlyUnassigned=True`): the unspecified (or
/// unknown) tetrahedral atoms and double bonds `Chem.FindPotentialStereo`
/// reports, in its order. A double bond's flip is relative to its first
/// controlling atom on each side (ordered by RDKit's symmetry ranks).
pub(crate) fn flippers(base: &Mol) -> (Vec<usize>, Vec<(usize, usize, usize)>) {
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    for info in find_potential_stereo(base) {
        if !matches!(info.specified, Specified::Unspecified | Specified::Unknown) {
            continue;
        }
        match info.kind {
            StereoType::AtomTetrahedral => atoms.push(info.centered_on),
            StereoType::BondDouble => {
                let b = info.centered_on;
                let sa = &base.bonds[b].stereo_atoms;
                if sa.len() == 2 {
                    bonds.push((b, sa[0], sa[1]));
                } else if info.controlling[0] != NOATOM && info.controlling[2] != NOATOM {
                    bonds.push((b, info.controlling[0], info.controlling[2]));
                }
            }
            _ => {}
        }
    }
    (atoms, bonds)
}

/// The distinct canonical SMILES `EnumerateStereoisomers` yields for the
/// sanitized, stereo-perceived `base` (sorted), including RDKit's
/// deterministic default random sample when there are more than
/// `max_isomers` flip combinations; `Err` only for `max_isomers == 0` with
/// more than 16 flips.
pub(crate) fn enumerate(base: &Mol, max_isomers: usize) -> Result<Vec<String>, RdkitSmilesError> {
    let (atoms, bonds) = flippers(base);
    let n = atoms.len() + bonds.len();
    if n == 0 {
        return Ok(vec![mol_to_smiles(base, &RdkitSmilesParams::default())?]);
    }
    // RDKit enumerates every flip combination when there are at most
    // `max_isomers` of them (or `max_isomers == 0`); otherwise it draws
    // combinations with `random.Random(hash(tuple(sorted((degree, Z)))))`
    // `.getrandbits(n)`, skipping repeats, until it has yielded
    // `max_isomers` distinct isomers or seen every combination.
    let full = max_isomers == 0 || (n < 64 && (1u64 << n) <= max_isomers as u64);
    if full {
        const MAX_FULL_ENUMERATION_FLIPS: usize = 16;
        if n > MAX_FULL_ENUMERATION_FLIPS {
            return Err(RdkitSmilesError::Unsupported(format!(
                "{n} stereo flips: 2^{n} isomers requested"
            )));
        }
        let mut seen = std::collections::BTreeSet::new();
        for bitflag in 0..(1usize << n) {
            let m = assign(base, &atoms, &bonds, |i| bitflag & (1 << i) != 0);
            seen.insert(mol_to_smiles(&m, &RdkitSmilesParams::default())?);
        }
        return Ok(seen.into_iter().collect());
    }
    let mut key: Vec<(i64, i64)> = (0..base.atoms.len())
        .map(|a| (base.degree(a) as i64, base.atoms[a].anum as i64))
        .collect();
    key.sort_unstable();
    let mut rng = PyRandom::from_int_seed(hash_pair_tuple(&key));
    let combos: Option<u128> = (n < 128).then(|| 1u128 << n);
    let mut tried = std::collections::HashSet::new();
    let mut seen = std::collections::BTreeSet::new();
    // RDKit keeps drawing until every combination has been seen; refuse
    // rather than loop for an astronomically long time.
    const MAX_DRAWS: usize = 1 << 20;
    while combos.is_none_or(|c| (tried.len() as u128) < c) {
        if tried.len() >= MAX_DRAWS {
            return Err(RdkitSmilesError::Unsupported(format!(
                "{n} stereo flips: fewer than {max_isomers} distinct isomers in {MAX_DRAWS} draws"
            )));
        }
        let bits = rng.getrandbits(n);
        if !tried.insert(bits.clone()) {
            continue;
        }
        let m = assign(base, &atoms, &bonds, |i| bits[i / 32] >> (i % 32) & 1 != 0);
        if seen.insert(mol_to_smiles(&m, &RdkitSmilesParams::default())?)
            && seen.len() >= max_isomers
        {
            break;
        }
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
    fn random_sampling_follows_cpython_random() {
        let mol = crate::parse(&format!("Br{}F", "[CH](Cl)".repeat(20))).unwrap();
        let got = rdkit_stereoisomer_smiles(&mol, 1024).unwrap();
        assert_eq!(got.len(), 1024);
        let unlimited = rdkit_stereoisomer_smiles(&mol, 0);
        assert!(unlimited.is_err());
    }
}
