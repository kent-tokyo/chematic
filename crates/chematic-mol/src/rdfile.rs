//! Typed RDFile v1 records, independently implemented from CTfile 2026 Ch. 11.
//! Embedded V2000 RXN and V2000/V3000 MOL structures reuse existing parsers.
use serde::{Deserialize, Serialize};
mod reader;
mod writer;
pub use reader::{parse_rdfile, parse_rdfile_with_limits};
pub use writer::{write_rdfile, write_rdfile_with_options};

/// Select spec folding or the explicit unwrapped-data interchange profile.
#[derive(Debug, Clone, Copy)]
pub struct RdfileWriteOptions {
    pub wrap_data: bool,
}
impl Default for RdfileWriteOptions {
    fn default() -> Self {
        Self { wrap_data: true }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct RdfileParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_records: usize,
    pub max_fields_per_record: usize,
    pub max_reaction_components: usize,
}
impl Default for RdfileParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 1 << 20,
            max_records: 10_000,
            max_fields_per_record: 1_000,
            max_reaction_components: 1_000,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdfileError {
    /// One-based source line, or zero for document/write validation.
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for RdfileError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "RDFile line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for RdfileError {}
fn error(line: usize, detail: impl Into<String>) -> RdfileError {
    RdfileError {
        line,
        detail: detail.into(),
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RdfileKind {
    Molecule,
    Reaction,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RdfileData {
    pub name: String,
    /// Literal text, including significant spaces, blank lines and duplicates.
    /// Database field types are unavailable; numbers/references remain strings.
    pub value: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RdfileRecord {
    pub kind: RdfileKind,
    /// Decimal text, preserving leading zeros; no registry lookup is performed.
    pub internal_registry: Option<String>,
    pub external_registry: Option<String>,
    /// Exact embedded MOL/RXN bytes as UTF-8, or None for reference-only records.
    pub structure: Option<String>,
    /// Ordered fields; duplicate names are never collapsed into a map.
    pub data: Vec<RdfileData>,
}
impl RdfileRecord {
    /// Explicit graph projection. Returns None for unresolved registry references.
    /// Coordinates and fields remain in this record; they are not in SMILES.
    pub fn smiles(&self) -> Result<Option<String>, RdfileError> {
        let Some(block) = &self.structure else {
            return Ok(None);
        };
        let smiles = match self.kind {
            RdfileKind::Molecule => {
                let v3000 = block.lines().nth(3).is_some_and(|s| s.contains("V3000"));
                let parsed = if v3000 {
                    crate::mol3000::parse_mol_v3000(block)
                } else {
                    crate::mol2000::parse_mol(block)
                };
                let (mol, _) = parsed.map_err(|e| error(0, e.to_string()))?;
                chematic_smiles::canonical_smiles(&mol)
            }
            RdfileKind::Reaction => {
                let reaction = crate::parse_rxn_file(block).map_err(|e| error(0, e.to_string()))?;
                chematic_rxn::write_reaction(&reaction)
            }
        };
        Ok(Some(smiles))
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RdfileDocument {
    pub date_comment: String,
    pub records: Vec<RdfileRecord>,
    /// Original complete input. Canonical writing uses current fields instead.
    #[serde(default)]
    pub source_text: String,
}
