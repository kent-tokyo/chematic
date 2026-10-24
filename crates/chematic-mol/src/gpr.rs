//! Checked Ghemical gpr 100, one coordinate set, explicit S/D/T/C bonds.
//! Raw nm coordinate and partial-charge decimals remain authoritative strings.
//! No bond/charge perception, implicit H inference or molecular projection.
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
mod reader;
mod writer;
pub use reader::{parse_gpr, parse_gpr_with_limits};
pub use writer::{write_gpr, write_gpr_with_limits};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GprError {
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for GprError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "GPR line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for GprError {}
fn error(line: usize, detail: impl Into<String>) -> GprError {
    GprError {
        line,
        detail: detail.into(),
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GprParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
    pub max_bonds: usize,
}
impl Default for GprParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 4096,
            max_atoms: 1_000_000,
            max_bonds: 2_000_000,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GprAtom {
    pub id: u32,
    pub atomic_number: u8,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GprBondOrder {
    S,
    D,
    T,
    C,
}
impl GprBondOrder {
    pub fn code(self) -> &'static str {
        match self {
            Self::S => "S",
            Self::D => "D",
            Self::T => "T",
            Self::C => "C",
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GprBond {
    pub atoms: [u32; 2],
    pub order: GprBondOrder,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GprCoordinate {
    pub atom_id: u32,
    pub position_nm: [String; 3],
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GprPartialCharge {
    pub atom_id: u32,
    pub value: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GprRecord {
    /// Source atom-definition order; IDs must be consecutive 0..N-1.
    pub atoms: Vec<GprAtom>,
    /// Source bond order and endpoint orientation. C is the literal aromatic code.
    pub bonds: Vec<GprBond>,
    /// Exactly one row per atom, in consecutive atom-ID order.
    pub coordinates: Vec<GprCoordinate>,
    /// None retains absent !Charges; Some contains exactly one row per atom.
    /// These are partial charges, never formal charges, in consecutive atom-ID order.
    pub partial_charges: Option<Vec<GprPartialCharge>>,
}
impl GprRecord {
    pub fn validate(&self) -> Result<(), GprError> {
        self.validate_with_limits(&GprParseLimits::default())
    }
    pub fn validate_with_limits(&self, limits: &GprParseLimits) -> Result<(), GprError> {
        let n = self.atoms.len();
        if n == 0
            || n > limits.max_atoms
            || n > u32::MAX as usize
            || self.bonds.len() > limits.max_bonds
        {
            return Err(error(0, "atom or bond count outside limits"));
        }
        validate_ids(self.atoms.iter().map(|a| a.id), n)?;
        if self
            .atoms
            .iter()
            .any(|a| !(1..=118).contains(&a.atomic_number))
        {
            return Err(error(0, "atomic number must be 1..118"));
        }
        let mut pairs = HashSet::new();
        for bond in &self.bonds {
            let [a, b] = bond.atoms;
            if a as usize >= n || b as usize >= n || a == b || !pairs.insert((a.min(b), a.max(b))) {
                return Err(error(
                    0,
                    "unknown endpoint, self bond or duplicate unordered bond",
                ));
            }
        }
        validate_ids(self.coordinates.iter().map(|p| p.atom_id), n)?;
        for row in &self.coordinates {
            for value in &row.position_nm {
                real(value)?;
            }
        }
        if let Some(charges) = &self.partial_charges {
            validate_ids(charges.iter().map(|p| p.atom_id), n)?;
            for row in charges {
                real(&row.value)?;
            }
        }
        Ok(())
    }
}
fn validate_ids(ids: impl Iterator<Item = u32>, n: usize) -> Result<(), GprError> {
    let mut count = 0;
    for (index, id) in ids.enumerate() {
        if id as usize != index {
            return Err(error(
                0,
                "rows must use consecutive atom IDs 0..N-1 in order",
            ));
        }
        count += 1;
    }
    if count != n {
        return Err(error(0, "missing atom IDs"));
    }
    Ok(())
}
/// Strict editable-record JSON: unknown and recursively duplicate keys refuse.
pub fn parse_gpr_json(text: &str, limits: &GprParseLimits) -> Result<GprRecord, GprError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "JSON byte limit exceeded"));
    }
    let value: crate::cjson_record::UniqueValue =
        serde_json::from_str(text).map_err(|e| error(0, e.to_string()))?;
    let record: GprRecord = serde_json::from_value(value.0).map_err(|e| error(0, e.to_string()))?;
    record.validate_with_limits(limits)?;
    Ok(record)
}
/// Checked nominal value; decimal token stays authoritative. No unit conversion.
pub fn gpr_nominal_value(token: &str) -> Result<f64, GprError> {
    real(token)
}
fn real(token: &str) -> Result<f64, GprError> {
    if token.len() > 256 || !decimal(token) {
        return Err(error(0, "invalid decimal or token exceeds 256 bytes"));
    }
    let value: f64 = token.parse().map_err(|_| error(0, "invalid real value"))?;
    let mantissa = token.split(['e', 'E']).next().unwrap_or("");
    if !value.is_finite() || (value == 0. && mantissa.bytes().any(|b| matches!(b, b'1'..=b'9'))) {
        return Err(error(0, "non-finite or underflowed decimal"));
    }
    Ok(value)
}
fn decimal(s: &str) -> bool {
    let bytes = s.as_bytes();
    let mut i = usize::from(matches!(bytes.first(), Some(b'+' | b'-')));
    let mut digits = 0;
    while bytes.get(i).is_some_and(u8::is_ascii_digit) {
        i += 1;
        digits += 1;
    }
    if bytes.get(i) == Some(&b'.') {
        i += 1;
        while bytes.get(i).is_some_and(u8::is_ascii_digit) {
            i += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return false;
    }
    if matches!(bytes.get(i), Some(b'e' | b'E')) {
        i += 1;
        if matches!(bytes.get(i), Some(b'+' | b'-')) {
            i += 1;
        }
        let start = i;
        while bytes.get(i).is_some_and(u8::is_ascii_digit) {
            i += 1;
        }
        if i == start {
            return false;
        }
    }
    i == bytes.len()
}
