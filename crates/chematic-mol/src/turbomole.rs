//! Checked standalone TURBOMOLE coordinate documents (nonperiodic or 3D).
//!
//! Independently implemented from TURBOMOLE 8.0 and MCTC format documentation.
//! No topology, isotope, charge, spin or external-file inference is performed.
use chematic_core::Element;

mod reader;

pub use reader::{parse_turbomole, parse_turbomole_with_limits};

/// Bohr radius in angstroms, CODATA 2022 (central value).
pub const TURBOMOLE_BOHR_ANGSTROM: f64 = 0.529_177_210_544;

#[derive(Debug, Clone, Copy)]
pub struct TurbomoleParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
}
impl Default for TurbomoleParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 1024,
            max_atoms: 1_000_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TurbomoleError {
    /// One-based source line; zero denotes validation of a constructed geometry.
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for TurbomoleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "TURBOMOLE line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for TurbomoleError {}
fn error(line: usize, detail: impl Into<String>) -> TurbomoleError {
    TurbomoleError {
        line,
        detail: detail.into(),
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct TurbomoleAtom {
    pub element: Element,
    pub position_angstrom: [f64; 3],
    pub frozen: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TurbomoleGeometry {
    pub atoms: Vec<TurbomoleAtom>,
    /// Three Cartesian row vectors, preserving orientation and handedness.
    pub lattice_angstrom: Option<[[f64; 3]; 3]>,
    /// Original input; canonical writing does not reuse potentially stale text.
    pub source_text: String,
}
impl TurbomoleGeometry {
    pub fn periodicity(&self) -> u8 {
        if self.lattice_angstrom.is_some() {
            3
        } else {
            0
        }
    }
    /// Validate mutable public data before emitting a coordinate document.
    pub fn validate(&self) -> Result<(), TurbomoleError> {
        if self.atoms.is_empty() {
            return Err(error(0, "at least one atom is required"));
        }
        for atom in &self.atoms {
            if atom.position_angstrom.iter().any(|v| !v.is_finite()) {
                return Err(error(0, "non-finite coordinate"));
            }
        }
        if let Some(lattice) = &self.lattice_angstrom {
            validate_lattice(lattice)?;
        }
        Ok(())
    }
}

/// Require a finite nonsingular 3D cell without changing its basis.
/// The normalized determinant must exceed 1e-6 in magnitude (Gram > 1e-12).
fn validate_lattice(m: &[[f64; 3]; 3]) -> Result<(), TurbomoleError> {
    let mut normalized = [[0.0; 3]; 3];
    let mut volume_scale = 1.0;
    for (i, row) in m.iter().enumerate() {
        let norm = row[0].hypot(row[1]).hypot(row[2]);
        if !norm.is_finite() || norm <= 0.0 {
            return Err(error(0, "lattice vector must have finite positive length"));
        }
        volume_scale *= norm;
        for (j, value) in row.iter().enumerate() {
            normalized[i][j] = value / norm;
        }
    }
    let [a, b, c] = normalized;
    let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0]);
    let volume = det.abs() * volume_scale;
    if det.abs() <= 1e-6 || !volume.is_finite() || volume <= 0.0 {
        return Err(error(
            0,
            "singular, near-singular or non-finite lattice volume",
        ));
    }
    Ok(())
}

fn bohr(value: f64) -> Result<f64, TurbomoleError> {
    let result = value / TURBOMOLE_BOHR_ANGSTROM;
    if !result.is_finite() {
        return Err(error(
            0,
            "coordinate or lattice overflows on conversion to bohr",
        ));
    }
    Ok(result)
}

/// Canonical bohr coordinate document with 17 decimal places in scientific
/// notation. Preserves atom order, fixed-atom flags and all lattice components.
/// Coordinates are never wrapped or rotated; no bonds or implicit H are added.
pub fn write_turbomole(geometry: &TurbomoleGeometry) -> Result<String, TurbomoleError> {
    geometry.validate()?;
    let mut out = String::new();
    if let Some(lattice) = &geometry.lattice_angstrom {
        out.push_str("$periodic 3\n$lattice\n");
        for row in lattice {
            out.push_str(&format!(
                " {:.17e} {:.17e} {:.17e}\n",
                bohr(row[0])?,
                bohr(row[1])?,
                bohr(row[2])?
            ));
        }
    }
    out.push_str("$coord\n");
    for atom in &geometry.atoms {
        let p = atom.position_angstrom;
        out.push_str(&format!(
            " {:.17e} {:.17e} {:.17e} {}{}\n",
            bohr(p[0])?,
            bohr(p[1])?,
            bohr(p[2])?,
            atom.element.symbol().to_ascii_lowercase(),
            if atom.frozen { " f" } else { "" }
        ));
    }
    out.push_str("$end\n");
    Ok(out)
}
