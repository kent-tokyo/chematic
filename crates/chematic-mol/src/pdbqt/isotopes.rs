//! Explicit isotope projection for Open Babel #3005's labelled PDBQT workflow.
use super::{PdbqtWriteError, checked::validate_pdbqt, write_pdbqt};
use chematic_core::Molecule;
use serde::Serialize;

/// An isotope annotation omitted from PDBQT, linked to the output atom serial.
/// The sidecar records source metadata; it is not part of the docking input.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct PdbqtIsotopeLoss {
    pub source_atom_index: usize,
    pub pdbqt_serial: usize,
    pub atomic_number: u8,
    pub mass_number: u16,
}

/// Rigid PDBQT plus every omitted isotope annotation in source atom order.
/// This is not a lossless graph round trip: bonds, formal charges, atom maps,
/// stereo and implicit H counts are not stored by this rigid PDBQT profile.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct PdbqtIsotopeProjection {
    pub pdbqt: String,
    pub isotopes_omitted: Vec<PdbqtIsotopeLoss>,
}

/// Explicitly project isotope-labelled atoms into ordinary AutoDock types.
///
/// Retains every explicit graph atom and its order, coordinates and supplied
/// partial charge (rounded to three decimals). No H merging, charge folding,
/// ligand torsion tree or chemical preparation is performed. Isotope identity
/// survives only in `isotopes_omitted`, never in PDBQT itself. The default
/// `write_pdbqt_checked` refuses labelled atoms instead. All other checked
/// writer limits still apply and the source molecule is never modified.
pub fn project_pdbqt_isotopes(
    mol: &Molecule,
    coords: &[(f64, f64, f64)],
    charges: &[f64],
    ligand_name: &str,
) -> Result<PdbqtIsotopeProjection, PdbqtWriteError> {
    validate_pdbqt(mol, coords, charges, ligand_name, false)?;
    let isotopes_omitted = mol
        .atoms()
        .filter_map(|(idx, atom)| {
            atom.isotope.map(|mass_number| PdbqtIsotopeLoss {
                source_atom_index: idx.0 as usize,
                pdbqt_serial: idx.0 as usize + 1,
                atomic_number: atom.element.atomic_number(),
                mass_number,
            })
        })
        .collect();
    Ok(PdbqtIsotopeProjection {
        pdbqt: write_pdbqt(mol, coords, charges, ligand_name),
        isotopes_omitted,
    })
}
