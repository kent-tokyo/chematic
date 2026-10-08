//! RDKit's `molToReaccs`: `MolToMolBlock(mol, includeStereo=true)` read back
//! with `MolStr2Mol`, for a chematic molecule.

use chematic_core::Molecule;
use chematic_smiles::{RdkitMolBlock, RdkitSmilesError, rdkit_mol_block};

use super::{AvalonAtom, AvalonBond, AvalonMolecule};

/// Why [`super::rdkit_avalon_fp`] produced no fingerprint.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RdkitAvalonError {
    /// The molecule is outside what the RDKit model handles, or RDKit
    /// would reject it (`Chem.MolFromSmiles` returns `None`).
    Rdkit(RdkitSmilesError),
}

impl core::fmt::Display for RdkitAvalonError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Rdkit(e) => write!(f, "RDKit Avalon fingerprint: {e}"),
        }
    }
}

impl std::error::Error for RdkitAvalonError {}

impl From<RdkitSmilesError> for RdkitAvalonError {
    fn from(e: RdkitSmilesError) -> Self {
        Self::Rdkit(e)
    }
}

/// The Avalon-toolkit molecule `MolStr2Mol` reads from RDKit's MOL block.
pub fn avalon_molecule_from_rdkit_view(block: &RdkitMolBlock) -> AvalonMolecule {
    AvalonMolecule {
        atoms: block
            .atoms
            .iter()
            .map(|a| AvalonAtom {
                symbol: a.symbol.to_string(),
                charge: a.formal_charge,
                radical: a.radical,
                atext: String::new(),
            })
            .collect(),
        bonds: block
            .bonds
            .iter()
            .map(|b| AvalonBond {
                atoms: [b.begin as i32 + 1, b.end as i32 + 1],
                bond_type: b.bond_type,
            })
            .collect(),
    }
}

pub(crate) fn avalon_molecule_from_molecule(
    mol: &Molecule,
) -> Result<AvalonMolecule, RdkitAvalonError> {
    let block = rdkit_mol_block(mol)?;
    Ok(avalon_molecule_from_rdkit_view(&block))
}
