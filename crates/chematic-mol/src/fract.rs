//! Free Form Fractional (`fract`), a checked single-cell, element-site profile.
//!
//! Decimal tokens are authoritative: nominal precision, signed zero and opaque
//! parenthesized uncertainty digits survive edits and serialization. Fractional
//! positions are dimensionless; no wrapping, bond perception or cell defaults.
use serde::{Deserialize, Serialize};
mod reader;
pub use reader::{parse_fract, parse_fract_with_limits};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FractError {
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for FractError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fract line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for FractError {}
fn error(line: usize, detail: impl Into<String>) -> FractError {
    FractError {
        line,
        detail: detail.into(),
    }
}

#[derive(Debug, Clone, Copy)]
pub struct FractParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
}
impl Default for FractParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 4096,
            max_atoms: 1_000_000,
        }
    }
}

/// An ASCII decimal, optionally followed by `(digits)`. Uncertainty is retained
/// verbatim, not converted to an absolute standard deviation. At most 256 bytes.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(try_from = "String", into = "String")]
pub struct FractNumber(String);
impl FractNumber {
    pub fn new(token: impl Into<String>) -> Result<Self, FractError> {
        let token = token.into();
        checked_number(&token)?;
        Ok(Self(token))
    }
    pub fn token(&self) -> &str {
        &self.0
    }
    /// Finite nominal value; the original decimal token remains authoritative.
    pub fn value(&self) -> f64 {
        // Construction and deserialization validate this private field.
        nominal(&self.0).parse().expect("validated fract number")
    }
    pub fn uncertainty_digits(&self) -> Option<&str> {
        self.0
            .split_once('(')
            .map(|(_, tail)| &tail[..tail.len() - 1])
    }
}
impl TryFrom<String> for FractNumber {
    type Error = FractError;
    fn try_from(value: String) -> Result<Self, Self::Error> {
        Self::new(value)
    }
}
impl From<FractNumber> for String {
    fn from(value: FractNumber) -> Self {
        value.0
    }
}
fn nominal(token: &str) -> &str {
    token.split_once('(').map_or(token, |(v, _)| v)
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
fn checked_number(token: &str) -> Result<f64, FractError> {
    if token.len() > 256 {
        return Err(error(0, "numeric token exceeds 256 bytes"));
    }
    if let Some((_, tail)) = token.split_once('(') {
        let digits = tail.strip_suffix(')').unwrap_or("");
        if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
            return Err(error(
                0,
                "uncertainty must be trailing parenthesized ASCII digits",
            ));
        }
    }
    let raw = nominal(token);
    if !decimal(raw) {
        return Err(error(0, "invalid decimal token"));
    }
    let value: f64 = raw.parse().map_err(|_| error(0, "invalid nominal value"))?;
    let mantissa = raw.split(['e', 'E']).next().unwrap_or("");
    if !value.is_finite() || (value == 0. && mantissa.bytes().any(|b| matches!(b, b'1'..=b'9'))) {
        return Err(error(0, "non-finite or underflowed nominal value"));
    }
    Ok(value)
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FractAtom {
    /// Canonical chemical element symbol, not a site label or pseudoatom.
    pub element: String,
    pub fractional: [FractNumber; 3],
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FractRecord {
    pub title: String,
    /// a, b, c (Angstrom), alpha, beta, gamma (degrees), in source order.
    pub cell: [FractNumber; 6],
    pub atoms: Vec<FractAtom>,
}
impl FractRecord {
    pub fn validate(&self) -> Result<(), FractError> {
        self.validate_with_limits(&FractParseLimits::default())
    }
    pub fn validate_with_limits(&self, limits: &FractParseLimits) -> Result<(), FractError> {
        if self.title.len() > limits.max_line_len
            || self.title.chars().any(|c| c.is_control() && c != '\t')
        {
            return Err(error(
                1,
                "title contains a control character or exceeds line limit",
            ));
        }
        if self.atoms.is_empty() || self.atoms.len() > limits.max_atoms {
            return Err(error(0, "atom count outside supported limits"));
        }
        let [a, b, c, alpha, beta, gamma] = self.cell.each_ref().map(FractNumber::value);
        crate::cif::UnitCell {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        }
        .validate()
        .map_err(|e| error(2, e.to_string()))?;
        for (i, atom) in self.atoms.iter().enumerate() {
            validate_element(&atom.element, i + 3)?;
        }
        Ok(())
    }
}
/// Parse a strict editable record JSON, rejecting unknown fields and duplicate
/// keys at every depth (including escaped aliases). Numeric fields are strings.
pub fn parse_fract_json(text: &str, limits: &FractParseLimits) -> Result<FractRecord, FractError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "JSON byte limit exceeded"));
    }
    let value: crate::cjson_record::UniqueValue =
        serde_json::from_str(text).map_err(|e| error(0, format!("invalid record JSON: {e}")))?;
    let record: FractRecord = serde_json::from_value(value.0)
        .map_err(|e| error(0, format!("invalid record JSON: {e}")))?;
    record.validate_with_limits(limits)?;
    Ok(record)
}
/// Serialize current fields, retaining every decimal/uncertainty token. Delimiters
/// become spaces and line endings LF; original whitespace is not retained.
pub fn write_fract(record: &FractRecord) -> Result<String, FractError> {
    write_fract_with_limits(record, &FractParseLimits::default())
}
pub fn write_fract_with_limits(
    record: &FractRecord,
    limits: &FractParseLimits,
) -> Result<String, FractError> {
    record.validate_with_limits(limits)?;
    let mut out = String::new();
    append_line(&mut out, &record.title, limits)?;
    let cell = record.cell.each_ref().map(FractNumber::token).join(" ");
    append_line(&mut out, &cell, limits)?;
    for atom in &record.atoms {
        let row = format!(
            "{} {}",
            atom.element,
            atom.fractional.each_ref().map(FractNumber::token).join(" ")
        );
        append_line(&mut out, &row, limits)?;
    }
    Ok(out)
}
fn append_line(out: &mut String, line: &str, limits: &FractParseLimits) -> Result<(), FractError> {
    if line.len() > limits.max_line_len
        || out.len().saturating_add(line.len()).saturating_add(1) > limits.max_input_bytes
    {
        return Err(error(0, "serialized line or byte limit exceeded"));
    }
    out.push_str(line);
    out.push('\n');
    Ok(())
}

fn validate_element(symbol: &str, line: usize) -> Result<(), FractError> {
    if !chematic_core::Element::from_symbol(symbol).is_some_and(|e| e.symbol() == symbol) {
        return Err(error(line, "expected canonical element symbol"));
    }
    Ok(())
}
