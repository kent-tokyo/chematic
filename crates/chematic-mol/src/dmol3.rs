//! Checked standalone DMol3 .incoor geometry, independently implemented from
//! the public ASE format description. No Open Babel/ASE implementation copied.
use chematic_core::Element;
mod reader;
pub use reader::{parse_dmol3, parse_dmol3_with_limits};

/// Bohr radius in Angstrom, CODATA 2022 central value. Raw Bohr stays authoritative.
pub const DMOL3_BOHR_ANGSTROM: f64 = 0.529_177_210_544;

#[derive(Debug, Clone, Copy)]
pub struct Dmol3ParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
}
impl Default for Dmol3ParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 1024,
            max_atoms: 1_000_000,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Dmol3Error {
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for Dmol3Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "DMol3 line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for Dmol3Error {}
fn error(line: usize, detail: impl Into<String>) -> Dmol3Error {
    Dmol3Error {
        line,
        detail: detail.into(),
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Dmol3Atom {
    pub element: Element,
    pub position_bohr: [f64; 3],
}
#[derive(Debug, Clone, PartialEq)]
pub struct Dmol3Geometry {
    pub atoms: Vec<Dmol3Atom>,
    /// Three ordered Cartesian row vectors in Bohr; None means nonperiodic.
    pub cell_bohr: Option<[[f64; 3]; 3]>,
    pub source_text: String,
}
impl Dmol3Geometry {
    pub fn periodicity(&self) -> u8 {
        if self.cell_bohr.is_some() { 3 } else { 0 }
    }
    pub fn validate(&self) -> Result<(), Dmol3Error> {
        if self.atoms.is_empty() || self.atoms.len() > Dmol3ParseLimits::default().max_atoms {
            return Err(error(0, "atom count outside supported limits"));
        }
        if self
            .atoms
            .iter()
            .any(|a| a.position_bohr.iter().any(|v| !v.is_finite()))
        {
            return Err(error(0, "non-finite coordinate"));
        }
        if let Some(cell) = self.cell_bohr {
            validate_cell(&cell)?;
        }
        Ok(())
    }
    pub fn positions_angstrom(&self) -> Result<Vec<[f64; 3]>, Dmol3Error> {
        self.validate()?;
        let mut out = Vec::new();
        out.try_reserve_exact(self.atoms.len())
            .map_err(|e| error(0, format!("coordinate allocation failed: {e}")))?;
        out.extend(
            self.atoms
                .iter()
                .map(|a| a.position_bohr.map(|v| v * DMOL3_BOHR_ANGSTROM)),
        );
        Ok(out)
    }
    pub fn cell_angstrom(&self) -> Result<Option<[[f64; 3]; 3]>, Dmol3Error> {
        self.validate()?;
        Ok(self
            .cell_bohr
            .map(|m| m.map(|r| r.map(|v| v * DMOL3_BOHR_ANGSTROM))))
    }
}
fn validate_cell(m: &[[f64; 3]; 3]) -> Result<(), Dmol3Error> {
    let mut normalized = [[0.; 3]; 3];
    let mut scale = 1.;
    for (i, row) in m.iter().enumerate() {
        let norm = row[0].hypot(row[1]).hypot(row[2]);
        if !norm.is_finite() || norm <= 0. {
            return Err(error(0, "cell vector must have finite positive length"));
        }
        scale *= norm;
        for (j, v) in row.iter().enumerate() {
            normalized[i][j] = v / norm;
        }
    }
    let [a, b, c] = normalized;
    let det = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0]);
    let volume = det.abs() * scale;
    if det.abs() <= 1e-6 || !volume.is_finite() || volume <= 0. {
        return Err(error(
            0,
            "singular, near-singular or non-finite cell volume",
        ));
    }
    Ok(())
}

/// Serialize current validated Bohr values, never stale source_text. No wrapping,
/// rotation, cell normalization, bond perception or electronic-state inference.
pub fn write_dmol3(g: &Dmol3Geometry) -> Result<String, Dmol3Error> {
    g.validate()?;
    let mut out = String::new();
    if let Some(cell) = g.cell_bohr {
        out.push_str("$cell vectors\n");
        for row in cell {
            write_row(&mut out, row);
        }
    }
    out.push_str("$coordinates\n");
    for atom in &g.atoms {
        out.push_str(atom.element.symbol());
        write_row(&mut out, atom.position_bohr);
        if out.len() > Dmol3ParseLimits::default().max_input_bytes {
            return Err(error(0, "serialized byte limit exceeded"));
        }
    }
    out.push_str("$end\n");
    if out.len() > Dmol3ParseLimits::default().max_input_bytes {
        return Err(error(0, "serialized byte limit exceeded"));
    }
    Ok(out)
}
fn write_row(out: &mut String, p: [f64; 3]) {
    out.push_str(&format!(" {:.17e} {:.17e} {:.17e}\n", p[0], p[1], p[2]));
}
