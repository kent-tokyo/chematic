//! Explicit molecular bond geometry, without bond perception or periodic images.
use super::{CifError, UnitCell, write_cif_checked};
use chematic_core::{BondOrder, Molecule};

/// Failure to write the bounded molecular CIF bond profile.
/// Separate from the established exhaustive [`CifError`] API.
#[derive(Debug)]
#[non_exhaustive]
pub enum CifBondWriteError {
    /// Invalid atomic coordinates or cell parameters.
    InvalidInput(CifError),
    /// This profile cannot represent the bond order without changing its meaning.
    UnsupportedBond { bond: usize, order: BondOrder },
    /// A bond's Cartesian distance must be finite and strictly positive.
    InvalidDistance { bond: usize, distance: f64 },
}

impl core::fmt::Display for CifBondWriteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput(error) => error.fmt(f),
            Self::UnsupportedBond { bond, order } => {
                write!(f, "CIF bond {bond} has unsupported order {order:?}")
            }
            Self::InvalidDistance { bond, distance } => write!(
                f,
                "CIF bond {bond} requires a finite positive Cartesian distance, got {distance}"
            ),
        }
    }
}

impl std::error::Error for CifBondWriteError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidInput(error) => Some(error),
            _ => None,
        }
    }
}

/// Write checked P1 atomic positions plus the molecule's explicit bond table.
///
/// Coordinates are Cartesian Å in exact molecule atom order. Applies all
/// [`write_cif_checked`] validations. Atom labels match that writer's labels;
/// bond rows retain graph order and endpoint order. `_geom_bond_distance` is
/// the Euclidean distance in the supplied Cartesian frame, even with a cell:
/// no minimum-image wrapping, bond perception or symmetry expansion occurs.
/// Both symmetry columns are `.` (the listed sites, without transformation).
/// `_ccdc_geom_bond_type` uses S/D/T/A for single/double/triple/aromatic bonds.
/// Other bond orders and zero/non-finite distances fail without returning output.
///
/// This is a coordinate and connectivity projection, not a full molecular
/// serialization: charge, isotopes, implicit H, stereochemical annotations and
/// metadata are not exported. [`super::parse_cif`] still reads positions only.
/// An edgeless molecule produces no bond loop.
pub fn write_cif_with_bonds_checked(
    mol: &Molecule,
    coords: &[(f64, f64, f64)],
    cell: Option<&UnitCell>,
) -> Result<String, CifBondWriteError> {
    let mut out = write_cif_checked(mol, coords, cell).map_err(CifBondWriteError::InvalidInput)?;
    if mol.bond_count() == 0 {
        return Ok(out);
    }
    out.push_str(
        "\nloop_\n_geom_bond_atom_site_label_1\n_geom_bond_atom_site_label_2\n\
         _geom_bond_distance\n_geom_bond_site_symmetry_1\n\
         _geom_bond_site_symmetry_2\n_ccdc_geom_bond_type\n",
    );
    for (idx, bond) in mol.bonds() {
        let kind = match bond.order {
            BondOrder::Single => 'S',
            BondOrder::Double => 'D',
            BondOrder::Triple => 'T',
            BondOrder::Aromatic => 'A',
            order => {
                return Err(CifBondWriteError::UnsupportedBond {
                    bond: idx.0 as usize,
                    order,
                });
            }
        };
        let (x1, y1, z1) = coords[bond.atom1.0 as usize];
        let (x2, y2, z2) = coords[bond.atom2.0 as usize];
        let distance = (x1 - x2).hypot(y1 - y2).hypot(z1 - z2);
        if !distance.is_finite() || distance <= 0.0 {
            return Err(CifBondWriteError::InvalidDistance {
                bond: idx.0 as usize,
                distance,
            });
        }
        let label1 = format!(
            "{}{}",
            mol.atom(bond.atom1).element.symbol(),
            bond.atom1.0 + 1
        );
        let label2 = format!(
            "{}{}",
            mol.atom(bond.atom2).element.symbol(),
            bond.atom2.0 + 1
        );
        out.push_str(&format!("{label1} {label2} {distance} . . {kind}\n"));
    }
    Ok(out)
}
