//! `Canon::rankMolAtoms` (RDKit 2026.03.1 `new_canon`) for a [`Mol`], through
//! chematic-perception's port of the ranking.

use chematic_perception::{
    RdkitRankAtom, RdkitRankBond, rdkit_rank_fragment_atoms, rdkit_rank_mol_atoms,
};

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
    let (atoms, bonds) = rank_input(mol, isomeric, |i| mol.num_atom_rings(i) as u32);
    rdkit_rank_mol_atoms(&atoms, &bonds)
}

/// `rankFragmentAtoms(mol, ranks, atomsInPlay, bondsInPlay, nullptr,
/// nullptr, breakTies=true, includeChirality=isomeric,
/// includeIsotopes=isomeric, includeAtomMaps=true,
/// includeChiralPresence=false, includeRingStereo=true)` as
/// `MolFragmentToSmiles` calls it: the copied molecule's ring information
/// is not "fast or better", so the ranking counts the rings
/// `MolOps::fastFindRings` finds in the whole molecule.
pub(crate) fn rank_fragment_atoms(
    mol: &Mol,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    isomeric: bool,
) -> Vec<u32> {
    let counts = fast_find_ring_counts(mol);
    let (atoms, bonds) = rank_input(mol, isomeric, |i| counts[i]);
    rdkit_rank_fragment_atoms(&atoms, &bonds, atoms_in_play, bonds_in_play, isomeric)
}

/// Per atom, the number of rings `MolOps::fastFindRings` (`_DFS`) stores.
pub(crate) fn fast_find_ring_counts(mol: &Mol) -> Vec<u32> {
    let n = mol.atoms.len();
    let mut colors = vec![0u8; n];
    let mut counts = vec![0u32; n];
    let mut order: Vec<usize> = Vec::new();
    // Explicit stack of (atom, parent, next neighbour position).
    for root in 0..n {
        if colors[root] != 0 {
            continue;
        }
        if mol.degree(root) < 2 {
            colors[root] = 2;
            continue;
        }
        let mut stack: Vec<(usize, Option<usize>, usize)> = vec![(root, None, 0)];
        colors[root] = 1;
        order.push(root);
        while let Some(&mut (atom, from, ref mut pos)) = stack.last_mut() {
            let bonds = &mol.atom_bonds[atom];
            if *pos == bonds.len() {
                colors[atom] = 2;
                order.pop();
                stack.pop();
                continue;
            }
            let nbr = mol.bonds[bonds[*pos]].other(atom);
            *pos += 1;
            if colors[nbr] == 0 {
                if mol.degree(nbr) < 2 {
                    colors[nbr] = 2;
                } else {
                    colors[nbr] = 1;
                    order.push(nbr);
                    stack.push((nbr, Some(atom), 0));
                }
            } else if colors[nbr] == 1
                && let Some(f) = from
                && nbr != f
            {
                // The ring: `atom` back along the traversal to `nbr`.
                let last = order.iter().rposition(|&x| x == atom).expect("on path");
                for k in (0..=last).rev() {
                    if order[k] == nbr {
                        break;
                    }
                    counts[order[k]] += 1;
                }
                counts[nbr] += 1;
            }
        }
    }
    counts
}

/// The ranking's view of `mol`'s atoms and bonds (`num_rings` per atom
/// from `num_rings`).
fn rank_input(
    mol: &Mol,
    isomeric: bool,
    num_rings: impl Fn(usize) -> u32,
) -> (Vec<RdkitRankAtom>, Vec<RdkitRankBond>) {
    let atoms: Vec<RdkitRankAtom> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let chiral_tag = if isomeric { a.chiral.rdkit_value() } else { 0 };
            RdkitRankAtom {
                atomic_num: a.anum,
                isotope: if isomeric { a.isotope } else { 0 },
                formal_charge: a.charge,
                atom_map: a.map.map_or(0, |m| m as i32),
                chiral_tag,
                total_num_hs: mol.total_num_hs(i),
                num_rings: num_rings(i),
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
    (atoms, bonds)
}
