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
    /// The complete input exceeded the configured byte limit.
    InputTooLarge { limit: usize },
    /// The declared atom count exceeded the configured limit.
    TooManyAtoms { count: usize, limit: usize },
    /// A physical input line exceeded the configured byte limit.
    LineTooLong { line: usize, limit: usize },
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
            Self::InputTooLarge { limit } => write!(f, "XYZ input exceeds {limit}-byte limit"),
            Self::TooManyAtoms { count, limit } => {
                write!(f, "XYZ atom count {count} exceeds {limit}-atom limit")
            }
            Self::LineTooLong { line, limit } => {
                write!(f, "XYZ line {line} exceeds {limit}-byte limit")
            }
            Self::InvalidAtomCount => write!(f, "invalid atom count in XYZ header"),
            Self::InvalidLine(n) => write!(f, "invalid XYZ coordinate line {n}"),
            Self::UnknownElement(s) => write!(f, "unknown element symbol '{s}' in XYZ file"),
        }
    }
}

/// Resource limits applied by the 3D XYZ parser.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct XyzParseLimits {
    pub max_input_bytes: usize,
    pub max_atoms: usize,
    pub max_line_bytes: usize,
}

impl Default for XyzParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_atoms: 2_000_000,
            max_line_bytes: 1024,
        }
    }
}

/// Parse an XYZ format string into a [`Molecule`] and [`Coords3D`].
///
/// The returned `Molecule` contains only heavy atoms with no bonds; XYZ files
/// do not encode connectivity.
pub fn parse_xyz(input: &str) -> Result<(Molecule, Coords3D), XyzError> {
    parse_xyz_with_limits(input, XyzParseLimits::default())
}

/// Parse XYZ with explicit resource limits.
pub fn parse_xyz_with_limits(
    input: &str,
    limits: XyzParseLimits,
) -> Result<(Molecule, Coords3D), XyzError> {
    if input.len() > limits.max_input_bytes {
        return Err(XyzError::InputTooLarge {
            limit: limits.max_input_bytes,
        });
    }

    let mut lines = input.lines();

    // Line 1: atom count.
    let count_line = lines.next().unwrap_or("");
    if count_line.len() > limits.max_line_bytes {
        return Err(XyzError::LineTooLong {
            line: 1,
            limit: limits.max_line_bytes,
        });
    }
    let count_line = count_line.trim();
    let n: usize = count_line.parse().map_err(|_| XyzError::InvalidAtomCount)?;
    if n > limits.max_atoms {
        return Err(XyzError::TooManyAtoms {
            count: n,
            limit: limits.max_atoms,
        });
    }

    // Line 2: comment — consumed and discarded.
    if let Some(comment) = lines.next()
        && comment.len() > limits.max_line_bytes
    {
        return Err(XyzError::LineTooLong {
            line: 2,
            limit: limits.max_line_bytes,
        });
    }

    let mut builder = MoleculeBuilder::new();
    let mut points: Vec<Point3> = Vec::with_capacity(n);

    for i in 0..n {
        // Line index in the file is i + 3 (1-indexed), but we just use i for clarity.
        let raw_line = lines.next().ok_or(XyzError::InvalidLine(i + 3))?;
        if raw_line.len() > limits.max_line_bytes {
            return Err(XyzError::LineTooLong {
                line: i + 3,
                limit: limits.max_line_bytes,
            });
        }
        let line = raw_line.trim();

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

#[cfg(test)]
mod shared_contract_tests {
    use super::{parse_xyz, write_xyz};
    use chematic_smiles::canonical_smiles;
    use serde_json::Value;

    const CONTRACT: &str = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/cross_binding_contract.json"
    ));

    #[test]
    fn shared_xyz_contract_matches() {
        let document: Value = serde_json::from_str(CONTRACT).expect("contract JSON");
        let contract = &document["xyz_contract"];
        assert_eq!(contract["schema_version"], 1);
        let (mol, coords) = parse_xyz(contract["input"].as_str().unwrap()).unwrap();
        let expected = &contract["expected"];
        assert_eq!(
            mol.atom_count(),
            expected["atom_count"].as_u64().unwrap() as usize
        );
        assert_eq!(
            canonical_smiles(&mol),
            expected["canonical_smiles"].as_str().unwrap()
        );
        for (point, expected_point) in coords
            .points
            .iter()
            .zip(expected["coords"].as_array().unwrap())
        {
            let values = expected_point.as_array().unwrap();
            assert!((point.x - values[0].as_f64().unwrap()).abs() < 1e-9);
            assert!((point.y - values[1].as_f64().unwrap()).abs() < 1e-9);
            assert!((point.z - values[2].as_f64().unwrap()).abs() < 1e-9);
        }
    }

    #[test]
    fn shared_xyz_roundtrip_contract_matches() {
        let document: Value = serde_json::from_str(CONTRACT).expect("contract JSON");
        let contract = &document["xyz_contract"];
        let (mol, coords) = parse_xyz(contract["input"].as_str().unwrap()).unwrap();
        let serialized = write_xyz(&mol, &coords, "shared-contract-roundtrip");
        let (roundtripped, roundtripped_coords) = parse_xyz(&serialized).unwrap();
        let expected = &contract["expected"];
        assert_eq!(roundtripped.atom_count(), expected["atom_count"]);
        assert_eq!(
            canonical_smiles(&roundtripped),
            expected["canonical_smiles"].as_str().unwrap()
        );
        assert_eq!(roundtripped_coords.points.len(), coords.points.len());
        for (actual, original) in roundtripped_coords.points.iter().zip(coords.points.iter()) {
            assert!((actual.x - original.x).abs() < 1e-6);
            assert!((actual.y - original.y).abs() < 1e-6);
            assert!((actual.z - original.z).abs() < 1e-6);
        }
    }
}
