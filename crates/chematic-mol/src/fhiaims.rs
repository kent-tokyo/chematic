//! Checked FHI-aims geometry.in, independently implemented from manual 261002.
//! No Open Babel or FHI-aims implementation source is incorporated.
use chematic_core::Element;
mod reader;
mod writer;
pub use reader::{parse_fhiaims, parse_fhiaims_with_limits};
pub use writer::write_fhiaims;

#[derive(Debug, Clone, Copy)]
pub struct FhiaimsParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
}
impl Default for FhiaimsParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 16_384,
            max_atoms: 1_000_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FhiaimsError {
    /// One-based source line; zero denotes whole-geometry validation.
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for FhiaimsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FHI-aims line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for FhiaimsError {}
fn error(line: usize, detail: impl Into<String>) -> FhiaimsError {
    FhiaimsError {
        line,
        detail: detail.into(),
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FhiaimsCoordinates {
    Cartesian([f64; 3]),
    Fractional([f64; 3]),
}
impl FhiaimsCoordinates {
    pub fn values(self) -> [f64; 3] {
        match self {
            Self::Cartesian(v) | Self::Fractional(v) => v,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct FhiaimsAtom {
    pub element: Element,
    pub coordinates: FhiaimsCoordinates,
    /// Fixed Cartesian components, in x/y/z order.
    pub constraints: [bool; 3],
    /// Starting-density charge in units of positive elementary charge.
    /// This is not a molecular formal charge or final atomic population.
    pub initial_charge: Option<f64>,
    /// Starting collinear N_up - N_down; not a spin multiplicity.
    pub initial_moment: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FhiaimsLatticeVector {
    pub vector_angstrom: [f64; 3],
    pub constraints: [bool; 3],
}

#[derive(Debug, Clone, PartialEq)]
pub struct FhiaimsGeometry {
    pub atoms: Vec<FhiaimsAtom>,
    /// Zero to three ordered independent Cartesian vectors; never padded.
    pub lattice_vectors: Vec<FhiaimsLatticeVector>,
    /// Input sidecar. Writing uses current structured data instead.
    pub source_text: String,
}
impl FhiaimsGeometry {
    pub fn periodicity(&self) -> usize {
        self.lattice_vectors.len()
    }

    pub fn validate(&self) -> Result<(), FhiaimsError> {
        if self.atoms.is_empty() || self.atoms.len() > FhiaimsParseLimits::default().max_atoms {
            return Err(error(0, "atom count outside supported limits"));
        }
        validate_lattice(&self.lattice_vectors)?;
        for atom in &self.atoms {
            if atom.coordinates.values().iter().any(|v| !v.is_finite())
                || atom.initial_charge.is_some_and(|v| !v.is_finite())
                || atom.initial_moment.is_some_and(|v| !v.is_finite())
            {
                return Err(error(
                    0,
                    "non-finite coordinate, initial charge or initial moment",
                ));
            }
            self.position_angstrom(atom)?;
        }
        Ok(())
    }

    fn position_angstrom(&self, atom: &FhiaimsAtom) -> Result<[f64; 3], FhiaimsError> {
        match atom.coordinates {
            FhiaimsCoordinates::Cartesian(v) => Ok(v),
            FhiaimsCoordinates::Fractional(v) => {
                if self.lattice_vectors.len() != 3 {
                    return Err(error(
                        0,
                        "fractional coordinates require exactly three lattice vectors",
                    ));
                }
                let mut out = [0.0; 3];
                for (coefficient, vector) in v.iter().zip(&self.lattice_vectors) {
                    for (component, basis) in out.iter_mut().zip(vector.vector_angstrom) {
                        *component += coefficient * basis;
                    }
                }
                if out.iter().any(|v| !v.is_finite()) {
                    return Err(error(0, "fractional-to-Cartesian conversion overflows"));
                }
                Ok(out)
            }
        }
    }

    /// Derived Cartesian angstrom positions, without wrapping or rotation.
    pub fn positions_angstrom(&self) -> Result<Vec<[f64; 3]>, FhiaimsError> {
        self.validate()?;
        self.atoms
            .iter()
            .map(|a| self.position_angstrom(a))
            .collect()
    }
}

fn validate_lattice(vectors: &[FhiaimsLatticeVector]) -> Result<(), FhiaimsError> {
    if vectors.len() > 3 {
        return Err(error(0, "at most three lattice vectors are supported"));
    }
    let mut normalized = [[0.0; 3]; 3];
    for (source, target) in vectors.iter().zip(&mut normalized) {
        let row = source.vector_angstrom;
        let norm = row[0].hypot(row[1]).hypot(row[2]);
        if !norm.is_finite() || norm <= 0.0 {
            return Err(error(0, "invalid lattice vector length"));
        }
        for (out, value) in target.iter_mut().zip(row) {
            *out = value / norm;
        }
    }
    let [a, b, c] = normalized;
    let cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let independent = match vectors.len() {
        2 => cross[0].hypot(cross[1]).hypot(cross[2]) > 1e-6,
        3 => (cross[0] * c[0] + cross[1] * c[1] + cross[2] * c[2]).abs() > 1e-6,
        _ => true,
    };
    if !independent {
        return Err(error(0, "dependent or near-dependent lattice vectors"));
    }
    Ok(())
}
