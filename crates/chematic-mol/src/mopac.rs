//! Checked MOPAC Cartesian geometry from the public OpenMOPAC documentation.
//! Independently implemented; no MOPAC/Open Babel implementation code is used.
use chematic_core::Element;
mod reader;
mod writer;
pub use reader::{parse_mopac_cartesian, parse_mopac_cartesian_with_limits};
pub use writer::write_mopac_cartesian;

#[derive(Debug, Clone, Copy)]
pub struct MopacParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
}
impl Default for MopacParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 16_384,
            max_atoms: 1_000_000,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MopacError {
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for MopacError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "MOPAC line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for MopacError {}
fn error(line: usize, detail: impl Into<String>) -> MopacError {
    MopacError {
        line,
        detail: detail.into(),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MopacOptimizationFlag {
    MiniReaction,
    Reaction,
    Fixed,
    Optimize,
    MiniOptimize,
    TurningPoint,
}
impl MopacOptimizationFlag {
    pub fn code(self) -> &'static str {
        match self {
            Self::MiniReaction => "-2",
            Self::Reaction => "-1",
            Self::Fixed => "0",
            Self::Optimize => "1",
            Self::MiniOptimize => "2",
            Self::TurningPoint => "T",
        }
    }
    pub fn from_code(code: &str) -> Option<Self> {
        match code {
            "-2" => Some(Self::MiniReaction),
            "-1" => Some(Self::Reaction),
            "0" | "+0" => Some(Self::Fixed),
            "1" | "+1" => Some(Self::Optimize),
            "2" | "+2" => Some(Self::MiniOptimize),
            "T" => Some(Self::TurningPoint),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MopacCartesianRow {
    pub position_angstrom: [f64; 3],
    /// None retains an unflagged three-number row, without guessing defaults.
    pub optimization_flags: Option<[MopacOptimizationFlag; 3]>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct MopacAtom {
    pub element: Element,
    pub label: Option<String>,
    pub row: MopacCartesianRow,
}
#[derive(Debug, Clone, PartialEq)]
pub struct MopacCartesianGeometry {
    /// Raw keyword line; no charge, spin, Hamiltonian or optimizer is inferred.
    pub keywords: String,
    pub title: String,
    pub comment: String,
    pub atoms: Vec<MopacAtom>,
    /// Zero to three Cartesian translation vectors in source order.
    pub translation_vectors: Vec<MopacCartesianRow>,
    pub source_text: String,
}
impl MopacCartesianGeometry {
    pub fn periodicity(&self) -> usize {
        self.translation_vectors.len()
    }
    pub fn validate(&self) -> Result<(), MopacError> {
        for (i, text) in [&self.keywords, &self.title, &self.comment]
            .iter()
            .enumerate()
        {
            if text.len() > MopacParseLimits::default().max_line_len
                || text.chars().any(|c| c.is_control() && c != '\t')
            {
                return Err(error(i + 1, "header must fit a single supported line"));
            }
        }
        for token in self.keywords.split_whitespace() {
            let key = token.split('=').next().unwrap_or("").to_ascii_uppercase();
            if key == "A0" || matches!(key.as_str(), "+" | "++" | "&") {
                return Err(error(
                    1,
                    "atomic-unit input and keyword continuations are unsupported",
                ));
            }
        }
        if self.atoms.is_empty() || self.atoms.len() > MopacParseLimits::default().max_atoms {
            return Err(error(0, "atom count outside supported limits"));
        }
        for atom in &self.atoms {
            if matches!(atom.element.atomic_number(), 85 | 87 | 99) {
                return Err(error(
                    0,
                    "reserved MOPAC sparkle/dummy element encoding is unsupported",
                ));
            }
            if let Some(label) = &atom.label {
                validate_label(label)?;
            }
            validate_row(&atom.row)?;
        }
        for row in &self.translation_vectors {
            validate_row(row)?;
        }
        validate_vectors(&self.translation_vectors)
    }
}
fn validate_label(label: &str) -> Result<(), MopacError> {
    if !label.is_ascii()
        || label.len() > 38
        || label
            .chars()
            .any(|c| c.is_control() || c == '(' || c == ')')
    {
        return Err(error(
            0,
            "atom label must be <=38 ASCII characters without parentheses or controls",
        ));
    }
    Ok(())
}
fn validate_row(row: &MopacCartesianRow) -> Result<(), MopacError> {
    if row.position_angstrom.iter().any(|v| !v.is_finite()) {
        return Err(error(0, "non-finite Cartesian coordinate"));
    }
    Ok(())
}
fn validate_vectors(vectors: &[MopacCartesianRow]) -> Result<(), MopacError> {
    if vectors.len() > 3 {
        return Err(error(0, "at most three translation vectors are supported"));
    }
    let mut normalized = [[0.; 3]; 3];
    for (row, out) in vectors.iter().zip(&mut normalized) {
        let v = row.position_angstrom;
        let norm = v[0].hypot(v[1]).hypot(v[2]);
        if !norm.is_finite() || norm <= 0. {
            return Err(error(
                0,
                "translation vector must have finite positive length",
            ));
        }
        for (a, b) in out.iter_mut().zip(v) {
            *a = b / norm;
        }
    }
    let [a, b, c] = normalized;
    let cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let valid = match vectors.len() {
        2 => cross[0].hypot(cross[1]).hypot(cross[2]) > 1e-6,
        3 => (cross[0] * c[0] + cross[1] * c[1] + cross[2] * c[2]).abs() > 1e-6,
        _ => true,
    };
    if !valid {
        return Err(error(0, "dependent or near-dependent translation vectors"));
    }
    Ok(())
}
