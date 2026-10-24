//! Editable nonperiodic BIOGRF 200 records with standard fixed-width fields.
//! Force-field types, valence capacity and lone pairs are source values;
//! none is inferred from coordinates or converted to formal charge.
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

mod reader;
mod writer;
pub use reader::{parse_bgf, parse_bgf_with_limits};
pub use writer::write_bgf;

pub const BGF_ATOM_FORMAT: &str = "(a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)";

#[derive(Debug, Clone, Copy)]
pub struct BgfLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
}
impl Default for BgfLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 16 << 20,
            max_line_len: 4096,
            max_atoms: 99_999,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BgfError {
    /// One-based source line; zero denotes record validation/global limits.
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for BgfError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "BGF line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for BgfError {}
fn error(line: usize, detail: impl Into<String>) -> BgfError {
    BgfError {
        line,
        detail: detail.into(),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BgfAtomKind {
    #[serde(rename = "ATOM")]
    Atom,
    #[serde(rename = "HETATM")]
    Hetatm,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BgfAtom {
    pub kind: BgfAtomKind,
    pub id: u32,
    pub name: String,
    pub residue_name: String,
    pub chain: String,
    /// A5 source identifier, not assumed to be an integer.
    pub residue_id: String,
    pub position_angstrom: [f64; 3],
    pub force_field_type: String,
    pub max_bonds: u16,
    pub lone_pairs: u8,
    pub partial_charge: f64,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BgfConnectivity {
    pub atom: u32,
    /// Source list order is retained; at most eleven neighbors per record.
    pub neighbors: Vec<u32>,
    /// None preserves implicit single bonds. Explicit 0 and 1 both mean single,
    /// but their source codes remain distinct. Only codes 0..=3 are supported.
    pub orders: Option<Vec<u8>>,
}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BgfRecord {
    /// This bounded profile accepts BIOGRF version 200 only.
    pub version: u16,
    pub description: Option<String>,
    pub remarks: Vec<String>,
    pub force_field: Option<String>,
    pub atoms: Vec<BgfAtom>,
    pub connectivity: Vec<BgfConnectivity>,
}
impl BgfRecord {
    pub fn validate(&self) -> Result<(), BgfError> {
        self.validate_with_limits(&BgfLimits::default())
    }
    pub(crate) fn validate_with_limits(&self, limits: &BgfLimits) -> Result<(), BgfError> {
        if self.version != 200 || self.atoms.is_empty() || self.atoms.len() > limits.max_atoms {
            return Err(error(
                0,
                "expected BIOGRF 200 and nonempty atom count within limits",
            ));
        }
        for value in self
            .description
            .iter()
            .chain(self.force_field.iter())
            .chain(&self.remarks)
        {
            if !value.is_ascii()
                || value.bytes().any(|b| b.is_ascii_control())
                || value.len() + 11 > limits.max_line_len
            {
                return Err(error(0, "metadata must be bounded printable ASCII"));
            }
        }
        let mut ids = BTreeSet::new();
        for atom in &self.atoms {
            if atom.id == 0 || atom.id > 99_999 || !ids.insert(atom.id) {
                return Err(error(0, "atom IDs must be unique in 1..99999"));
            }
            validate_atom(atom)?;
        }
        validate_connectivity(&self.connectivity, &ids)
    }
}
fn validate_atom(a: &BgfAtom) -> Result<(), BgfError> {
    for (s, width, required) in [
        (&a.name, 5, true),
        (&a.residue_name, 3, false),
        (&a.chain, 1, false),
        (&a.residue_id, 5, false),
        (&a.force_field_type, 5, true),
    ] {
        if s.len() > width || (required && s.is_empty()) || !s.bytes().all(|b| b.is_ascii_graphic())
        {
            return Err(error(
                0,
                "atom text field exceeds width or contains whitespace/non-ASCII",
            ));
        }
    }
    if !a.chain.bytes().all(|b| b.is_ascii_uppercase()) || a.max_bonds > 999 || a.lone_pairs > 99 {
        return Err(error(0, "invalid chain/max_bonds/lone_pairs field"));
    }
    for n in a.position_angstrom {
        fixed(n, 10)?;
    }
    fixed(a.partial_charge, 8)?;
    Ok(())
}
fn validate_connectivity(rows: &[BgfConnectivity], ids: &BTreeSet<u32>) -> Result<(), BgfError> {
    let mut centers = BTreeSet::new();
    let mut edges = BTreeMap::new();
    for row in rows {
        if !ids.contains(&row.atom) || !centers.insert(row.atom) || row.neighbors.len() > 11 {
            return Err(error(
                0,
                "unknown/duplicate CONECT center or too many neighbors",
            ));
        }
        if let Some(orders) = &row.orders
            && (orders.len() != row.neighbors.len() || orders.iter().any(|&n| n > 3))
        {
            return Err(error(0, "ORDER must match CONECT and use codes 0..3"));
        }
        let mut neighbors = BTreeSet::new();
        for (i, &n) in row.neighbors.iter().enumerate() {
            if n == row.atom || !ids.contains(&n) || !neighbors.insert(n) {
                return Err(error(0, "self/unknown/duplicate neighbor"));
            }
            let order = row.orders.as_ref().map_or(1, |v| v[i].max(1));
            if let Some(&reverse) = edges.get(&(n, row.atom))
                && reverse != order
            {
                return Err(error(0, "reciprocal bond orders disagree"));
            }
            edges.insert((row.atom, n), order);
        }
    }
    Ok(())
}
/// Refuse overflow and rounding rather than silently truncating edited values.
fn fixed(n: f64, width: usize) -> Result<String, BgfError> {
    let text = format!("{n:.5}");
    if !n.is_finite() || text.len() > width || text.parse::<f64>().ok() != Some(n) {
        return Err(error(
            0,
            "number is not exactly representable in fixed-width F*.5",
        ));
    }
    Ok(format!("{text:>width$}"))
}
pub fn bgf_record_from_json(json: &str) -> Result<BgfRecord, BgfError> {
    if json.len() > 64 << 20 {
        return Err(error(0, "record JSON exceeds byte limit"));
    }
    let record: BgfRecord = serde_json::from_str(json).map_err(|e| error(0, e.to_string()))?;
    record.validate()?;
    Ok(record)
}
