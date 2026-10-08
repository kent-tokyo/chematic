//! `Canon::rankMolAtoms` (RDKit 2026.03.1 `new_canon`) for a [`Mol`], through
//! chematic-perception's port of the ranking.

use chematic_perception::{RdkitRankAtom, RdkitRankBond, rdkit_rank_mol_atoms};

use super::mol::{BondType, ChiralTag, Mol};

/// `rankMolAtoms(mol, ranks, breakTies=true, includeChirality=true,
/// includeIsotopes=true, includeAtomMaps=true, includeChiralPresence=false,
/// includeStereoGroups=true, useNonStereoRanks=false)` on a molecule
/// without stereo groups.
pub(crate) fn rank_mol_atoms(mol: &Mol) -> Vec<u32> {
    rank_mol_atoms_with(mol, true)
}

/// [`rank_mol_atoms`] with `includeChirality` and `includeIsotopes` both
/// set to `isomeric` (`MolToSmiles` passes `doIsomericSmiles` for both).
pub(crate) fn rank_mol_atoms_with(mol: &Mol, isomeric: bool) -> Vec<u32> {
    let atoms: Vec<RdkitRankAtom> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let chiral_tag = match a.chiral {
                _ if !isomeric => 0,
                ChiralTag::Unspecified => 0,
                ChiralTag::Cw => 1,
                ChiralTag::Ccw => 2,
            };
            RdkitRankAtom {
                atomic_num: a.anum,
                isotope: if isomeric { a.isotope } else { 0 },
                formal_charge: a.charge,
                atom_map: a.map.map_or(0, |m| m as i32),
                chiral_tag,
                total_num_hs: mol.total_num_hs(i),
                num_rings: mol.num_atom_rings(i) as u32,
                ring_stereo: isomeric
                    && a.chiral != ChiralTag::Unspecified
                    && a.ring_stereo_atoms.is_some(),
            }
        })
        .collect();
    let bonds: Vec<RdkitRankBond> = mol
        .bonds
        .iter()
        .map(|b| RdkitRankBond {
            begin: b.begin as u32,
            end: b.end as u32,
            bond_type: if b.aromatic {
                BondType::Aromatic as u32
            } else {
                b.bt as u32
            },
            stereo: if isomeric { b.stereo as u32 } else { 0 },
            stereo_atoms: if isomeric && b.stereo_atoms.len() >= 2 {
                Some((b.stereo_atoms[0] as u32, b.stereo_atoms[1] as u32))
            } else {
                None
            },
        })
        .collect();
    rdkit_rank_mol_atoms(&atoms, &bonds)
}
