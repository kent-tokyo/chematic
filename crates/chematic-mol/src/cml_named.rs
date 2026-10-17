//! Named CML graph export, motivated by Open Babel #2120.

use std::collections::HashSet;
use std::fmt::Write;

use chematic_core::{AtomIdx, BondOrder, Chirality, Element, Molecule, implicit_hcount};

/// Invalid input or graph semantics outside the named writer's profile.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CmlNamedError {
    SizeLimit,
    IdCount,
    InvalidId { index: usize },
    DuplicateId { index: usize },
    CoordinateCount,
    NonFiniteCoordinate { index: usize },
    Unsupported(String),
}

impl std::fmt::Display for CmlNamedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SizeLimit => write!(
                f,
                "named CML requires 1..=4096 atoms and at most 8192 bonds"
            ),
            Self::IdCount => write!(f, "atom_ids must have one ID per atom"),
            Self::InvalidId { index } => write!(
                f,
                "invalid CML atom ID at index {index}: expected [A-Za-z][A-Za-z0-9_-]*, at most 128 bytes"
            ),
            Self::DuplicateId { index } => write!(f, "duplicate CML atom ID at index {index}"),
            Self::CoordinateCount => write!(f, "coordinates must have one 3D row per atom"),
            Self::NonFiniteCoordinate { index } => {
                write!(f, "non-finite coordinate at atom {index}")
            }
            Self::Unsupported(message) => write!(f, "named CML unsupported: {message}"),
        }
    }
}

impl std::error::Error for CmlNamedError {}

fn checked_order(order: BondOrder) -> Result<&'static str, CmlNamedError> {
    match order {
        BondOrder::Single => Ok("1"),
        BondOrder::Double => Ok("2"),
        BondOrder::Triple => Ok("3"),
        BondOrder::Aromatic => Ok("A"),
        _ => Err(CmlNamedError::Unsupported(format!("bond order {order:?}"))),
    }
}

fn validate(
    mol: &Molecule,
    ids: &[String],
    coords: Option<&[[f64; 3]]>,
) -> Result<(), CmlNamedError> {
    if mol.atom_count() == 0 || mol.atom_count() > 4096 || mol.bond_count() > 8192 {
        return Err(CmlNamedError::SizeLimit);
    }
    if ids.len() != mol.atom_count() {
        return Err(CmlNamedError::IdCount);
    }
    let mut seen = HashSet::new();
    for (index, id) in ids.iter().enumerate() {
        if id.len() > 128
            || !id.as_bytes().first().is_some_and(u8::is_ascii_alphabetic)
            || !id
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || b"_-".contains(&c))
        {
            return Err(CmlNamedError::InvalidId { index });
        }
        if !seen.insert(id) {
            return Err(CmlNamedError::DuplicateId { index });
        }
    }
    if let Some(coords) = coords {
        if coords.len() != mol.atom_count() {
            return Err(CmlNamedError::CoordinateCount);
        }
        if let Some(index) = coords
            .iter()
            .position(|row| row.iter().any(|v| !v.is_finite()))
        {
            return Err(CmlNamedError::NonFiniteCoordinate { index });
        }
    }
    if !mol.stereo_groups().is_empty() || mol.has_bond_directions() {
        return Err(CmlNamedError::Unsupported("stereo metadata".into()));
    }
    for (idx, atom) in mol.atoms() {
        if atom.wildcard
            || atom.atom_map.is_some()
            || mol.atom_tag(idx).is_some()
            || atom.chirality != Chirality::None
            || atom.cip_code.is_some()
            || atom.isotope == Some(0)
        {
            return Err(CmlNamedError::Unsupported(format!(
                "atom metadata at index {}",
                idx.0
            )));
        }
        let aromatic_bond = mol
            .neighbors(idx)
            .any(|(_, b)| mol.bond(b).order == BondOrder::Aromatic);
        if aromatic_bond != atom.aromatic {
            return Err(CmlNamedError::Unsupported(
                "inconsistent aromatic atom/bond flags".into(),
            ));
        }
    }
    for (_, bond) in mol.bonds() {
        checked_order(bond.order)?;
    }
    Ok(())
}

fn write_atom(out: &mut String, mol: &Molecule, idx: AtomIdx, id: &str, coord: Option<[f64; 3]>) {
    let atom = mol.atom(idx);
    let h = usize::from(implicit_hcount(mol, idx))
        + mol
            .neighbors(idx)
            .filter(|(n, _)| mol.atom(*n).element == Element::H)
            .count();
    write!(
        out,
        "    <atom id=\"{id}\" elementType=\"{}\" formalCharge=\"{}\" hydrogenCount=\"{h}\"",
        atom.element.symbol(),
        atom.charge
    )
    .unwrap();
    if let Some(isotope) = atom.isotope {
        write!(out, " isotopeNumber=\"{isotope}\"").unwrap();
    }
    if let Some([x, y, z]) = coord {
        write!(out, " x3=\"{x:.17e}\" y3=\"{y:.17e}\" z3=\"{z:.17e}\"").unwrap();
    }
    out.push_str("/>\n");
}

/// Write a named molecular graph using exactly the supplied atom IDs.
///
/// IDs are in atom insertion order, unique, ASCII IDs in the intersection of
/// the CML molecular-convention and Schema 3 atomRefs2 grammars,
/// with at most 128 bytes. Optional Cartesian coordinates are in Angstrom and
/// must contain one finite 3D row per atom. Names are never repaired, reordered
/// or replaced with generated IDs; all bond references use the supplied IDs.
///
/// Preserves elements, formal charges, isotopes, total H counts and single,
/// double, triple or aromatic bonds. Refuses stereo, maps, tags, wildcard atoms,
/// unsupported bonds and inconsistent aromatic flags. Does not encode source
/// file metadata, spin multiplicity, force-field types or partial charges.
/// Existing CML readers return graph/2D data and do not expose these IDs or 3D
/// coordinates; consumers needing them must inspect the exported CML.
pub fn write_cml_named_checked(
    mol: &Molecule,
    atom_ids: &[String],
    coordinates: Option<&[[f64; 3]]>,
) -> Result<String, CmlNamedError> {
    validate(mol, atom_ids, coordinates)?;
    let mut out = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<molecule xmlns=\"http://www.xml-cml.org/schema\" id=\"m1\">\n  <atomArray>\n",
    );
    for (idx, _) in mol.atoms() {
        let i = idx.0 as usize;
        write_atom(&mut out, mol, idx, &atom_ids[i], coordinates.map(|c| c[i]));
    }
    out.push_str("  </atomArray>\n  <bondArray>\n");
    for (_, bond) in mol.bonds() {
        writeln!(
            out,
            "    <bond atomRefs2=\"{} {}\" order=\"{}\"/>",
            atom_ids[bond.atom1.0 as usize],
            atom_ids[bond.atom2.0 as usize],
            checked_order(bond.order)?
        )
        .unwrap();
    }
    out.push_str("  </bondArray>\n</molecule>\n");
    Ok(out)
}
