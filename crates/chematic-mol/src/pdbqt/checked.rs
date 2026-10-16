//! Fixed-column output validation motivated by Open Babel issue #2101.
use super::{autodock_atom_type, write_pdbqt};
use chematic_core::Molecule;

/// Checked writer failures. Atom indices are zero based; fields name the
/// coordinate component or partial charge that could not be represented.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum PdbqtWriteError {
    EmptyMolecule,
    WildcardAtom {
        atom: usize,
    },
    AtomCount {
        count: usize,
        limit: usize,
    },
    CountMismatch {
        field: &'static str,
        atoms: usize,
        values: usize,
    },
    InvalidName,
    NonFinite {
        atom: usize,
        field: &'static str,
    },
    FieldOverflow {
        atom: usize,
        field: &'static str,
        width: usize,
    },
    UnsupportedAtom {
        atom: usize,
        atomic_number: u8,
    },
}

impl core::fmt::Display for PdbqtWriteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyMolecule => write!(f, "PDBQT writer requires at least one atom"),
            Self::WildcardAtom { atom } => write!(
                f,
                "PDBQT atom {atom}: wildcard has no concrete AutoDock type"
            ),
            Self::AtomCount { count, limit } => {
                write!(f, "PDBQT atom count {count} exceeds serial limit {limit}")
            }
            Self::CountMismatch {
                field,
                atoms,
                values,
            } => write!(f, "PDBQT expected {atoms} {field} values, got {values}"),
            Self::InvalidName => write!(
                f,
                "PDBQT name must contain at most 1024 printable ASCII bytes"
            ),
            Self::NonFinite { atom, field } => {
                write!(f, "PDBQT atom {atom}: {field} is not finite")
            }
            Self::FieldOverflow { atom, field, width } => write!(
                f,
                "PDBQT atom {atom}: {field} does not fit {width} columns after rounding to three decimals"
            ),
            Self::UnsupportedAtom {
                atom,
                atomic_number,
            } => write!(
                f,
                "PDBQT atom {atom}: no AutoDock type for atomic number {atomic_number}"
            ),
        }
    }
}
impl std::error::Error for PdbqtWriteError {}

fn check_numeric(
    value: f64,
    width: usize,
    atom: usize,
    field: &'static str,
) -> Result<(), PdbqtWriteError> {
    if !value.is_finite() {
        return Err(PdbqtWriteError::NonFinite { atom, field });
    }
    // Rust format widths are minimum widths. Inspect the actual rounded text,
    // including its sign, before delegating to the legacy formatter.
    if format!("{value:.3}").len() > width {
        return Err(PdbqtWriteError::FieldOverflow { atom, field, width });
    }
    Ok(())
}

/// Write rigid PDBQT after validating every input before producing output.
///
/// Requires exactly one coordinate (Å) and partial charge (e) per graph atom,
/// including explicit hydrogen atoms. No missing-value defaults are applied.
/// Coordinates use 8 columns and charges 6 columns, each with 3 decimals;
/// rounded values that exceed these widths fail. Serial numbers are 1..99999.
/// Names may be empty and otherwise contain printable ASCII including spaces.
///
/// AutoDock types use the existing heuristic; atoms mapped to `X` fail. This
/// validates serialization, not docking chemistry, protonation, stereochemistry
/// or ligand preparation. No atoms are removed/merged and no charges are
/// calculated. Output has synthetic LIG/A/1 residue identity and no torsion tree.
/// Ordinary accepted input produces exactly the legacy writer's bytes.
pub fn write_pdbqt_checked(
    mol: &Molecule,
    coords: &[(f64, f64, f64)],
    charges: &[f64],
    ligand_name: &str,
) -> Result<String, PdbqtWriteError> {
    let atoms = mol.atom_count();
    if atoms == 0 {
        return Err(PdbqtWriteError::EmptyMolecule);
    }
    if atoms > 99_999 {
        return Err(PdbqtWriteError::AtomCount {
            count: atoms,
            limit: 99_999,
        });
    }
    for (field, values) in [("coordinates", coords.len()), ("charges", charges.len())] {
        if values != atoms {
            return Err(PdbqtWriteError::CountMismatch {
                field,
                atoms,
                values,
            });
        }
    }
    if ligand_name.len() > 1024 || !ligand_name.bytes().all(|b| (0x20..=0x7e).contains(&b)) {
        return Err(PdbqtWriteError::InvalidName);
    }
    for (idx, a) in mol.atoms() {
        let i = idx.0 as usize;
        if a.wildcard {
            return Err(PdbqtWriteError::WildcardAtom { atom: i });
        }
        if autodock_atom_type(mol, idx) == "X" {
            return Err(PdbqtWriteError::UnsupportedAtom {
                atom: i,
                atomic_number: a.element.atomic_number(),
            });
        }
        let (x, y, z) = coords[i];
        for (value, width, field) in [
            (x, 8, "x"),
            (y, 8, "y"),
            (z, 8, "z"),
            (charges[i], 6, "charge"),
        ] {
            check_numeric(value, width, i, field)?;
        }
    }
    Ok(write_pdbqt(mol, coords, charges, ligand_name))
}
