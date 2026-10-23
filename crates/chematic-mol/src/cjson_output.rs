//! Checked projection to the existing ChemicalJSON graph/Cartesian profile.
use chematic_core::{Atom, BondOrder, Chirality, Molecule};

/// Information that the local CJSON graph writer cannot retain.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CjsonWriteError {
    CoordinateCount {
        expected: usize,
        actual: usize,
    },
    NonfiniteCoordinate {
        atom_index: usize,
        axis: usize,
    },
    AtomProperty {
        atom_index: usize,
        property: &'static str,
    },
    BondOrder {
        bond_index: usize,
        order: BondOrder,
    },
    BondDirection {
        bond_index: usize,
    },
    SelfBond {
        bond_index: usize,
        atom_index: usize,
    },
    StereoGroups,
}

impl std::fmt::Display for CjsonWriteError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::CoordinateCount { expected, actual } => write!(
                f,
                "CJSON: expected {expected} coordinate rows or none, got {actual}"
            ),
            Self::NonfiniteCoordinate { atom_index, axis } => write!(
                f,
                "CJSON: nonfinite coordinate at atom {atom_index}, axis {axis}"
            ),
            Self::AtomProperty {
                atom_index,
                property,
            } => write!(f, "CJSON: cannot preserve {property} at atom {atom_index}"),
            Self::BondOrder { bond_index, order } => write!(
                f,
                "CJSON: cannot preserve bond order {order:?} at bond {bond_index}"
            ),
            Self::BondDirection { bond_index } => write!(
                f,
                "CJSON: cannot preserve stereo direction at bond {bond_index}"
            ),
            Self::SelfBond {
                bond_index,
                atom_index,
            } => write!(
                f,
                "CJSON: self bond {bond_index} at atom {atom_index} is outside the reader profile"
            ),
            Self::StereoGroups => write!(f, "CJSON: cannot preserve enhanced stereo groups"),
        }
    }
}

impl std::error::Error for CjsonWriteError {}

fn unsupported_atom_property(atom: &Atom) -> Option<&'static str> {
    if atom.chirality != Chirality::None {
        Some("atom chirality")
    } else if atom.wildcard {
        Some("wildcard identity")
    } else if atom.hydrogen_count.is_some() {
        Some("explicit hydrogen count")
    } else if atom.aromatic {
        Some("atom aromaticity")
    } else if atom.atom_map.is_some() {
        Some("atom map")
    } else if atom.cip_code.is_some() {
        Some("CIP descriptor")
    } else if atom.isotope == Some(0) {
        Some("isotope mass zero")
    } else {
        None
    }
}

fn validate_graph(mol: &Molecule) -> Result<(), CjsonWriteError> {
    if !mol.stereo_groups().is_empty() {
        return Err(CjsonWriteError::StereoGroups);
    }
    for (index, atom) in mol.atoms() {
        let property = unsupported_atom_property(atom)
            .or_else(|| mol.r_group_label(index).map(|_| "R-group label"))
            .or_else(|| mol.atom_tag(index).map(|_| "atom tag"));
        if let Some(property) = property {
            return Err(CjsonWriteError::AtomProperty {
                atom_index: index.0 as usize,
                property,
            });
        }
    }
    for (index, bond) in mol.bonds() {
        if bond.atom1 == bond.atom2 {
            return Err(CjsonWriteError::SelfBond {
                bond_index: index.0 as usize,
                atom_index: bond.atom1.0 as usize,
            });
        }
        if !matches!(
            bond.order,
            BondOrder::Single | BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
        ) {
            return Err(CjsonWriteError::BondOrder {
                bond_index: index.0 as usize,
                order: bond.order,
            });
        }
        if mol.bond_direction(index).is_some() {
            return Err(CjsonWriteError::BondDirection {
                bond_index: index.0 as usize,
            });
        }
    }
    Ok(())
}

/// Serialize only if the existing CJSON reader can retain the in-memory graph.
///
/// Retains atom/bond indices, elements, formal charges, positive isotope masses,
/// single/double/triple/aromatic bond orders, and optional finite Cartesian rows
/// in Angstrom. Rejects all other atom properties, bond orders/directions, tags,
/// R-group labels and enhanced stereo groups. In particular, quadruple bonds
/// are refused because this reader treats order 4 as a legacy aromatic alias.
/// Aromatic *atom flags* are not encoded, even when aromatic bonds are present.
/// No perception, hydrogen expansion, kekulization or stereo inference is run.
///
/// This checks an already parsed molecule, not uninterpreted source metadata or
/// properties discarded by a source reader. Empty coordinates mean absent 3D;
/// otherwise exactly one finite row per atom is required. Errors return no text.
/// Self bonds are refused because the checked reader does not accept them.
/// Use [`crate::write_cjson_record`] for exact retention of a CJSON source, or
/// [`crate::write_cjson`] for the legacy lossy graph projection.
pub fn write_cjson_checked(
    mol: &Molecule,
    coords: &[(f64, f64, f64)],
) -> Result<String, CjsonWriteError> {
    if !coords.is_empty() && coords.len() != mol.atom_count() {
        return Err(CjsonWriteError::CoordinateCount {
            expected: mol.atom_count(),
            actual: coords.len(),
        });
    }
    for (atom_index, &(x, y, z)) in coords.iter().enumerate() {
        for (axis, value) in [x, y, z].into_iter().enumerate() {
            if !value.is_finite() {
                return Err(CjsonWriteError::NonfiniteCoordinate { atom_index, axis });
            }
        }
    }
    validate_graph(mol)?;
    Ok(crate::write_cjson(mol, coords))
}
