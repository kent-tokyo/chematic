//! Tinker XYZ Cartesian records, implemented from the public Tinker manual.
//! Names and force-field types remain opaque: this format declares connectivity,
//! not chemical bond orders, isotopes, charges, or an element mapping.

use serde::{Deserialize, Serialize};
use std::fmt::Write;

mod reader;
pub use reader::{parse_tinker_xyz, parse_tinker_xyz_with_limits};

/// Resource limits for one XYZ frame (not a concatenated Tinker ARC archive).
#[derive(Debug, Clone, Copy)]
pub struct TinkerParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
    pub max_neighbors: usize,
}

impl Default for TinkerParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 65_536,
            max_atoms: 100_000,
            max_neighbors: 64,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TinkerError {
    /// One-based input line, or zero for structured-record validation.
    pub line: usize,
    pub detail: String,
}

impl core::fmt::Display for TinkerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Tinker XYZ line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for TinkerError {}

fn error(line: usize, detail: impl Into<String>) -> TinkerError {
    TinkerError {
        line,
        detail: detail.into(),
    }
}

/// Six cell parameters; Cartesian atom positions are never transformed.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TinkerCell {
    pub lengths_angstrom: [f64; 3],
    pub angles_degrees: [f64; 3],
}

/// One-based IDs and neighbor order are retained exactly.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TinkerAtom {
    pub index: usize,
    /// Opaque ASCII atom name (1–3 characters, first alphabetic), not necessarily an element.
    pub name: String,
    pub position_angstrom: [f64; 3],
    /// Force-field-local identifier, including zero for an unassigned type.
    pub atom_type: u32,
    pub neighbors: Vec<usize>,
}

/// Editable single-frame record. No force-field parameters are inferred.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TinkerXyz {
    pub title: String,
    pub cell: Option<TinkerCell>,
    pub atoms: Vec<TinkerAtom>,
}

impl TinkerXyz {
    pub fn validate(&self) -> Result<(), TinkerError> {
        self.validate_with_limits(&TinkerParseLimits::default())
    }

    pub(crate) fn validate_with_limits(
        &self,
        limits: &TinkerParseLimits,
    ) -> Result<(), TinkerError> {
        if self.title.len() > limits.max_line_len
            || self.title.chars().any(char::is_control)
            || self.title.starts_with(char::is_whitespace)
        {
            return Err(error(
                0,
                "title must be one line without leading whitespace or controls",
            ));
        }
        if self.atoms.is_empty() || self.atoms.len() > limits.max_atoms {
            return Err(error(0, "atom count outside supported limits"));
        }
        if let Some(cell) = &self.cell {
            validate_cell(cell)?;
        }
        for (i, atom) in self.atoms.iter().enumerate() {
            validate_atom(atom, i, self.atoms.len(), limits)?;
            for &other in &atom.neighbors {
                if !self.atoms[other - 1].neighbors.contains(&atom.index) {
                    return Err(error(
                        0,
                        format!("non-reciprocal connection {} -> {other}", atom.index),
                    ));
                }
            }
        }
        Ok(())
    }
}

fn validate_atom(
    atom: &TinkerAtom,
    i: usize,
    count: usize,
    limits: &TinkerParseLimits,
) -> Result<(), TinkerError> {
    if atom.index != i + 1 {
        return Err(error(0, "atom IDs must be sequential, starting at 1"));
    }
    if atom.name.is_empty()
        || atom.name.len() > 3
        || !atom.name.as_bytes()[0].is_ascii_alphabetic()
        || !atom.name.bytes().all(|b| b.is_ascii_graphic())
    {
        return Err(error(
            0,
            "atom names need 1–3 printable ASCII characters, starting with a letter",
        ));
    }
    if atom.position_angstrom.iter().any(|v| !v.is_finite()) {
        return Err(error(0, "non-finite atom coordinate"));
    }
    if atom.neighbors.len() > limits.max_neighbors {
        return Err(error(0, "neighbor count exceeds limit"));
    }
    for (j, &other) in atom.neighbors.iter().enumerate() {
        if other == 0
            || other > count
            || other == atom.index
            || atom.neighbors[..j].contains(&other)
        {
            return Err(error(0, "dangling, self, or duplicate connection"));
        }
    }
    Ok(())
}

fn validate_cell(cell: &TinkerCell) -> Result<(), TinkerError> {
    if cell
        .lengths_angstrom
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || cell
            .angles_degrees
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0 || *v >= 180.0)
    {
        return Err(error(
            0,
            "cell needs finite positive lengths and angles strictly between 0 and 180 degrees",
        ));
    }
    let [a, b, c] = cell.angles_degrees.map(|v| v.to_radians().cos());
    let determinant = 1.0 + 2.0 * a * b * c - a * a - b * b - c * c;
    if determinant <= 1e-12 {
        return Err(error(
            0,
            "singular, near-singular, or impossible cell metric",
        ));
    }
    Ok(())
}

/// Serialize current structured data with round-trip f64 precision.
/// Whitespace is canonicalized; no source-text replay or chemical perception.
pub fn write_tinker_xyz(record: &TinkerXyz) -> Result<String, TinkerError> {
    record.validate()?;
    let limits = TinkerParseLimits::default();
    let mut out = format!("{}  {}\n", record.atoms.len(), record.title);
    if let Some(cell) = &record.cell {
        for (i, v) in cell
            .lengths_angstrom
            .iter()
            .chain(&cell.angles_degrees)
            .enumerate()
        {
            if i > 0 {
                out.push(' ');
            }
            write!(out, "{v:.17e}").expect("writing a String");
        }
        out.push('\n');
    }
    for atom in &record.atoms {
        let start = out.len();
        let [x, y, z] = atom.position_angstrom;
        write!(
            out,
            "{} {} {x:.17e} {y:.17e} {z:.17e} {}",
            atom.index, atom.name, atom.atom_type
        )
        .expect("writing a String");
        for id in &atom.neighbors {
            write!(out, " {id}").expect("writing a String");
        }
        out.push('\n');
        if out.len() - start - 1 > limits.max_line_len || out.len() > limits.max_input_bytes {
            return Err(error(0, "serialized record exceeds line or byte limit"));
        }
    }
    // Header includes a count and separator in addition to the public title.
    if out
        .lines()
        .next()
        .is_some_and(|line| line.len() > limits.max_line_len)
    {
        return Err(error(0, "serialized header exceeds line limit"));
    }
    Ok(out)
}

/// JSON helpers shared by Python, WASM and CLI; unknown fields are refused.
pub fn tinker_xyz_from_json(json: &str) -> Result<TinkerXyz, TinkerError> {
    if json.len() > TinkerParseLimits::default().max_input_bytes {
        return Err(error(0, "JSON exceeds byte limit"));
    }
    let record: TinkerXyz = serde_json::from_str(json).map_err(|e| error(0, e.to_string()))?;
    record.validate()?;
    Ok(record)
}
