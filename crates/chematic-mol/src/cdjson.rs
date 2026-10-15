//! Checked molecular subset of ChemDoodle JSON, motivated by Open Babel #2763.
//! Coordinates retain drawing units; no scale or orientation is inferred.
use chematic_core::Molecule;

mod reader;
mod writer;
pub use reader::{parse_cdjson, parse_cdjson_with_limits};
pub use writer::write_cdjson_graph;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CdjsonParseLimits {
    pub max_input_bytes: usize,
    pub max_molecules: usize,
    pub max_atoms: usize,
    pub max_bonds: usize,
}

impl Default for CdjsonParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 16 * 1024 * 1024,
            max_molecules: 1024,
            max_atoms: 65_536,
            max_bonds: 131_072,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CdjsonError {
    Invalid(String),
    Unsupported(String),
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
}

impl std::fmt::Display for CdjsonError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Invalid(message) => write!(f, "CDJSON invalid: {message}"),
            Self::Unsupported(message) => write!(f, "CDJSON unsupported: {message}"),
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "CDJSON resource limit: {actual} {resource} exceeds {limit}"
            ),
        }
    }
}
impl std::error::Error for CdjsonError {}

/// Validated molecular projection. IDs remain separate from native atom maps.
#[derive(Clone)]
pub struct CdjsonMolecule {
    pub molecule: Molecule,
    pub coordinates: Vec<[f64; 3]>,
    pub id: Option<String>,
    pub atom_ids: Vec<Option<String>>,
    pub bond_ids: Vec<Option<String>>,
}

/// Immutable source record. Same-format preservation returns exact source text.
#[derive(Clone)]
pub struct CdjsonDocument {
    source: String,
    molecules: Vec<CdjsonMolecule>,
}
impl CdjsonDocument {
    pub fn source_text(&self) -> &str {
        &self.source
    }
    pub fn molecules(&self) -> &[CdjsonMolecule] {
        &self.molecules
    }
}

fn limit(
    resource: &'static str,
    actual: usize,
    requested: usize,
    hard: usize,
) -> Result<(), CdjsonError> {
    let limit = requested.min(hard);
    if actual > limit {
        return Err(CdjsonError::ResourceLimit {
            resource,
            actual,
            limit,
        });
    }
    Ok(())
}
