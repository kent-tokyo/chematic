//! ORCA-AICCM EXYZ geometry, distinct from ASE extended XYZ.
//! Implemented from the public Open Babel format description, without copying
//! Open Babel source. All coordinates, lattice rows and offsets are in angstroms.
use chematic_core::Element;

mod reader;
pub use reader::{parse_exyz, parse_exyz_with_limits};

#[derive(Debug, Clone, Copy)]
pub struct ExyzParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
}
impl Default for ExyzParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 65_536,
            max_atoms: 1_000_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExyzError {
    /// One-based source line; zero means validation of constructed geometry.
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for ExyzError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "EXYZ line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for ExyzError {}
fn error(line: usize, detail: impl Into<String>) -> ExyzError {
    ExyzError {
        line,
        detail: detail.into(),
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ExyzAtom {
    pub element: Element,
    pub position_angstrom: [f64; 3],
    pub virtual_atom: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ExyzCell {
    /// Cartesian row vectors, preserving basis orientation and handedness.
    pub lattice_angstrom: [[f64; 3]; 3],
    /// Preserved verbatim as geometry; never added to atom positions.
    pub offset_angstrom: [f64; 3],
}

#[derive(Debug, Clone, PartialEq)]
pub struct ExyzGeometry {
    /// Complete second line, including exact `%PBC` / `%VIRTUAL` tokens.
    pub header: String,
    pub atoms: Vec<ExyzAtom>,
    pub cell: Option<ExyzCell>,
    /// Original input; the writer always serializes current structured data.
    pub source_text: String,
}

fn header_flags(header: &str) -> Result<(bool, bool), ExyzError> {
    if header.chars().any(|c| c.is_control() && c != '\t') {
        return Err(error(
            2,
            "header must be one line without control characters",
        ));
    }
    let (mut pbc, mut virtual_enabled) = (false, false);
    for token in header.split_whitespace() {
        let flag = match token {
            "%PBC" => &mut pbc,
            "%VIRTUAL" => &mut virtual_enabled,
            _ if token.starts_with('%') => return Err(error(2, "unsupported header directive")),
            _ if token.starts_with("Lattice=") || token.starts_with("Properties=") => {
                return Err(error(
                    2,
                    "ASE ExtXYZ is a different format; use parse_extxyz",
                ));
            }
            _ => continue,
        };
        if *flag {
            return Err(error(2, "duplicate header directive"));
        }
        *flag = true;
    }
    Ok((pbc, virtual_enabled))
}

impl ExyzGeometry {
    /// Validate public mutable data, including consistency with header flags.
    pub fn validate(&self) -> Result<(), ExyzError> {
        let (pbc, virtual_enabled) = header_flags(&self.header)?;
        let limits = ExyzParseLimits::default();
        if self.header.len() > limits.max_line_len
            || self.atoms.is_empty()
            || self.atoms.len() > limits.max_atoms
        {
            return Err(error(0, "header or atom count outside supported limits"));
        }
        if pbc != self.cell.is_some() {
            return Err(error(
                0,
                "%PBC must correspond to exactly one complete cell",
            ));
        }
        for atom in &self.atoms {
            if atom.position_angstrom.iter().any(|v| !v.is_finite()) {
                return Err(error(0, "non-finite atom coordinate"));
            }
            if atom.virtual_atom && !virtual_enabled {
                return Err(error(0, "VIRTUAL atom requires %VIRTUAL header"));
            }
        }
        if let Some(cell) = &self.cell {
            validate_cell(cell)?;
        }
        Ok(())
    }
}

fn validate_cell(cell: &ExyzCell) -> Result<(), ExyzError> {
    if cell.offset_angstrom.iter().any(|v| !v.is_finite()) {
        return Err(error(0, "non-finite cell offset"));
    }
    let mut rows = cell.lattice_angstrom;
    for row in &mut rows {
        let norm = row[0].hypot(row[1]).hypot(row[2]);
        if !norm.is_finite() || norm <= 0.0 {
            return Err(error(0, "cell vectors must have finite positive length"));
        }
        for v in row {
            *v /= norm;
        }
    }
    let [a, b, c] = rows;
    let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0]);
    if det.abs() <= 1e-6 {
        return Err(error(0, "singular or near-singular 3D lattice"));
    }
    Ok(())
}

/// Write current data with 17 decimal places in scientific notation.
/// Does not infer bonds, charge, spin, isotopes or implicit hydrogens.
pub fn write_exyz(geometry: &ExyzGeometry) -> Result<String, ExyzError> {
    geometry.validate()?;
    let mut out = format!("{}\n{}\n", geometry.atoms.len(), geometry.header);
    for atom in &geometry.atoms {
        let [x, y, z] = atom.position_angstrom;
        out.push_str(&format!(
            "{} {x:.17e} {y:.17e} {z:.17e}{}\n",
            atom.element.symbol(),
            if atom.virtual_atom { " VIRTUAL" } else { "" }
        ));
    }
    if let Some(cell) = &geometry.cell {
        out.push('\n');
        for (label, [x, y, z]) in ["Vector1", "Vector2", "Vector3", "Offset"].iter().zip(
            cell.lattice_angstrom
                .iter()
                .chain(std::iter::once(&cell.offset_angstrom)),
        ) {
            out.push_str(&format!("{label} {x:.17e} {y:.17e} {z:.17e}\n"));
        }
    }
    if out.len() > ExyzParseLimits::default().max_input_bytes {
        return Err(error(0, "serialized document exceeds byte limit"));
    }
    Ok(out)
}
