//! XYZ file format parser and writer.
//!
//! The XYZ format is a simple line-oriented format:
//! - Line 1: number of atoms (integer)
//! - Line 2: comment line (arbitrary text)
//! - Lines 3+: `<symbol>  <x>  <y>  <z>` (space-separated floats)

use chematic_core::{Atom, AtomIdx, Element, Molecule, MoleculeBuilder};

use crate::coords::{Coords3D, Point3};

/// Errors that can occur when parsing XYZ format.
#[derive(Debug, Clone, PartialEq)]
pub enum XyzError {
    /// The first line did not parse as a valid positive integer.
    InvalidAtomCount,
    /// A coordinate line (1-indexed, including header lines) could not be parsed.
    InvalidLine(usize),
    /// An element symbol from the file was not recognised.
    UnknownElement(String),
}

impl core::fmt::Display for XyzError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidAtomCount => write!(f, "invalid atom count in XYZ header"),
            Self::InvalidLine(n) => write!(f, "invalid XYZ coordinate line {n}"),
            Self::UnknownElement(s) => write!(f, "unknown element symbol '{s}' in XYZ file"),
        }
    }
}

/// Parse an XYZ format string into a [`Molecule`] and [`Coords3D`].
///
/// The returned `Molecule` contains only heavy atoms with no bonds; XYZ files
/// do not encode connectivity.
pub fn parse_xyz(input: &str) -> Result<(Molecule, Coords3D), XyzError> {
    let mut lines = input.lines();

    // Line 1: atom count.
    let count_line = lines.next().unwrap_or("").trim();
    let n: usize = count_line.parse().map_err(|_| XyzError::InvalidAtomCount)?;

    // Line 2: comment — consumed and discarded.
    lines.next();

    let mut builder = MoleculeBuilder::new();
    let mut points: Vec<Point3> = Vec::with_capacity(n);

    for i in 0..n {
        // Line index in the file is i + 3 (1-indexed), but we just use i for clarity.
        let line = lines.next().ok_or(XyzError::InvalidLine(i + 3))?.trim();

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(XyzError::InvalidLine(i + 3));
        }

        let symbol = parts[0];
        let x: f64 = parts[1].parse().map_err(|_| XyzError::InvalidLine(i + 3))?;
        let y: f64 = parts[2].parse().map_err(|_| XyzError::InvalidLine(i + 3))?;
        let z: f64 = parts[3].parse().map_err(|_| XyzError::InvalidLine(i + 3))?;

        let element = Element::from_symbol(symbol)
            .ok_or_else(|| XyzError::UnknownElement(symbol.to_string()))?;

        builder.add_atom(Atom::new(element));
        points.push(Point3::new(x, y, z));
    }

    // XYZ carries no bond information — return molecule with atoms only.
    let mol = builder.build();
    let coords = Coords3D { points };
    Ok((mol, coords))
}

/// Write a molecule and its coordinates as an XYZ format string.
///
/// `comment` is placed on the second line; it must not contain a newline.
pub fn write_xyz(mol: &Molecule, coords: &Coords3D, comment: &str) -> String {
    let n = mol.atom_count();
    let mut out = String::new();

    // Line 1: atom count.
    out.push_str(&n.to_string());
    out.push('\n');

    // Line 2: comment.
    out.push_str(comment);
    out.push('\n');

    // Atom lines.
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        let p = coords.get(idx);
        out.push_str(&format!(
            "{:<3} {:12.6} {:12.6} {:12.6}\n",
            atom.element.symbol(),
            p.x,
            p.y,
            p.z
        ));
    }

    out
}
