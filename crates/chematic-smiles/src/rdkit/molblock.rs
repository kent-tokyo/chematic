//! The connection table `Chem.MolToMolBlock(Chem.MolFromSmiles(s))` writes
//! (RDKit 2026.03.1 `MolFileWriter.cpp`: `prepareMol`,
//! `AtomGetMolFileSymbol`, `BondGetMolFileSymbol`, `GetMolFileChargeInfo`),
//! without coordinates, stereo flags or the text layout.

use super::kekulize::kekulize_ranked;
use super::mol::{BondType, Mol};
use super::rank::rank_mol_atoms;
use super::{RdkitSmilesError, parse, periodic, sanitize, stereo};
use chematic_core::Molecule;

/// One atom of the MOL block.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RdkitMolBlockAtom {
    /// Atomic number (0 for a dummy atom).
    pub atomic_num: u32,
    /// The atom symbol written (`"R"` for a plain dummy atom).
    pub symbol: &'static str,
    /// Formal charge (`M  CHG` / `CHG=`).
    pub formal_charge: i32,
    /// MDL radical code (`M  RAD` / `RAD=`): 0 none, 2 doublet, 3 triplet.
    pub radical: i32,
}

/// One bond of the MOL block.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RdkitMolBlockBond {
    /// 0-based index of the first atom.
    pub begin: usize,
    /// 0-based index of the second atom.
    pub end: usize,
    /// MDL bond type (1, 2, 3, 4; 9 for a dative bond).
    pub bond_type: i32,
}

/// The connection table of RDKit's MOL block for a molecule.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RdkitMolBlock {
    /// Atoms in RDKit's order.
    pub atoms: Vec<RdkitMolBlockAtom>,
    /// Bonds in RDKit's order.
    pub bonds: Vec<RdkitMolBlockBond>,
}

/// The atoms and bonds (symbols, charges, radicals, Kekulé bond types) of
/// `Chem.MolToMolBlock(Chem.MolFromSmiles(s))` for the SMILES `s` chematic
/// parsed `mol` from.
///
/// `MolFromSmiles` is modelled as in [`super::rdkit_canonical_smiles`]
/// (`removeHs`, `sanitizeMol`, legacy stereo perception); `prepareMol` then
/// kekulizes with `canonical=true`. The canonical atom ranks of that call
/// (`Canon::rankFragmentAtoms` over the whole molecule) are taken from
/// `Canon::rankMolAtoms`.
pub fn rdkit_mol_block(mol: &Molecule) -> Result<RdkitMolBlock, RdkitSmilesError> {
    let mut m = parse::from_chematic(mol)?;
    sanitize::remove_hs_and_sanitize(&mut m)?;
    stereo::legacy_stereo_perception(&mut m, true, true);
    mol_block_of(&mut m)
}

fn mol_block_of(m: &mut Mol) -> Result<RdkitMolBlock, RdkitSmilesError> {
    if !m.bonds.is_empty() {
        let ranks = rank_mol_atoms(m);
        kekulize_ranked(m, Some(&ranks))?;
    }
    let atoms = (0..m.atoms.len())
        .map(|a| {
            let atom = &m.atoms[a];
            let symbol = if atom.anum == 0 {
                "R"
            } else {
                periodic::symbol(atom.anum)
            };
            let radical = if atom.radicals != 0 && m.total_degree(a) != 0 {
                if atom.radicals % 2 == 1 { 2 } else { 3 }
            } else {
                0
            };
            RdkitMolBlockAtom {
                atomic_num: atom.anum,
                symbol,
                formal_charge: atom.charge,
                radical,
            }
        })
        .collect();
    let bonds = m
        .bonds
        .iter()
        .map(|b| RdkitMolBlockBond {
            begin: b.begin,
            end: b.end,
            bond_type: match b.bt {
                BondType::Single | BondType::Double if b.aromatic => 4,
                BondType::Single => 1,
                BondType::Double => 2,
                BondType::Triple => 3,
                BondType::Aromatic => 4,
                BondType::Dative => 9,
                BondType::Quadruple => 0,
            },
        })
        .collect();
    Ok(RdkitMolBlock { atoms, bonds })
}
