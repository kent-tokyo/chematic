//! Checked single-record UniChem XYZ: title, atom count, atomic number + XYZ.
//! Decimal coordinates remain authoritative strings in Angstrom. No bond,
//! charge, isotope, H or stereo inference; trailing records/columns refuse.
use serde::{Deserialize, Serialize};
mod reader;
pub use reader::{parse_unixyz, parse_unixyz_with_limits};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnixyzError {
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for UnixyzError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "UniChem XYZ line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for UnixyzError {}
fn error(line: usize, detail: impl Into<String>) -> UnixyzError {
    UnixyzError {
        line,
        detail: detail.into(),
    }
}
#[derive(Debug, Clone, Copy)]
pub struct UnixyzParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
}
impl Default for UnixyzParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 4096,
            max_atoms: 1_000_000,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UnixyzAtom {
    pub atomic_number: u8,
    pub position_angstrom: [String; 3],
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UnixyzRecord {
    pub title: String,
    pub atoms: Vec<UnixyzAtom>,
}
impl UnixyzRecord {
    pub fn validate(&self) -> Result<(), UnixyzError> {
        self.validate_with_limits(&UnixyzParseLimits::default())
    }
    pub fn validate_with_limits(&self, limits: &UnixyzParseLimits) -> Result<(), UnixyzError> {
        check_line(&self.title, 1, limits)?;
        if self.atoms.is_empty() || self.atoms.len() > limits.max_atoms {
            return Err(error(0, "atom count outside limits"));
        }
        for (i, atom) in self.atoms.iter().enumerate() {
            if !(1..=118).contains(&atom.atomic_number) {
                return Err(error(i + 3, "atomic number must be 1..118"));
            }
            for v in &atom.position_angstrom {
                unixyz_nominal_value(v).map_err(|e| error(i + 3, e.detail))?;
            }
        }
        Ok(())
    }
}
/// Strict record JSON with recursively unique keys, including escaped aliases.
/// Coordinate entries must be decimal strings; unknown fields refuse.
pub fn parse_unixyz_json(
    text: &str,
    limits: &UnixyzParseLimits,
) -> Result<UnixyzRecord, UnixyzError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "JSON byte limit exceeded"));
    }
    let value: crate::cjson_record::UniqueValue =
        serde_json::from_str(text).map_err(|e| error(0, e.to_string()))?;
    let record: UnixyzRecord =
        serde_json::from_value(value.0).map_err(|e| error(0, e.to_string()))?;
    record.validate_with_limits(limits)?;
    Ok(record)
}
/// Checked binary64 nominal value. The unchanged decimal token stays authoritative.
pub fn unixyz_nominal_value(token: &str) -> Result<f64, UnixyzError> {
    if token.len() > 256 || !decimal(token) {
        return Err(error(0, "invalid decimal or token exceeds 256 bytes"));
    }
    let value: f64 = token
        .parse()
        .map_err(|_| error(0, "invalid nominal value"))?;
    let mantissa = token.split(['e', 'E']).next().unwrap_or("");
    if !value.is_finite() || (value == 0. && mantissa.bytes().any(|b| matches!(b, b'1'..=b'9'))) {
        return Err(error(0, "non-finite or underflowed coordinate"));
    }
    Ok(value)
}
fn decimal(s: &str) -> bool {
    let b = s.as_bytes();
    let mut i = usize::from(matches!(b.first(), Some(b'+' | b'-')));
    let mut digits = 0;
    while b.get(i).is_some_and(u8::is_ascii_digit) {
        i += 1;
        digits += 1;
    }
    if b.get(i) == Some(&b'.') {
        i += 1;
        while b.get(i).is_some_and(u8::is_ascii_digit) {
            i += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return false;
    }
    if matches!(b.get(i), Some(b'e' | b'E')) {
        i += 1;
        if matches!(b.get(i), Some(b'+' | b'-')) {
            i += 1;
        }
        let start = i;
        while b.get(i).is_some_and(u8::is_ascii_digit) {
            i += 1;
        }
        if i == start {
            return false;
        }
    }
    i == b.len()
}
fn check_line(line: &str, number: usize, limits: &UnixyzParseLimits) -> Result<(), UnixyzError> {
    if line.len() > limits.max_line_len || line.chars().any(|c| c.is_control() && c != '\t') {
        return Err(error(number, "line limit exceeded or control character"));
    }
    Ok(())
}
/// Write current fields. Count is regenerated, numeric tokens/title are unchanged;
/// atom-row spacing becomes single spaces, output uses LF and a final newline.
pub fn write_unixyz(record: &UnixyzRecord) -> Result<String, UnixyzError> {
    write_unixyz_with_limits(record, &UnixyzParseLimits::default())
}
pub fn write_unixyz_with_limits(
    record: &UnixyzRecord,
    limits: &UnixyzParseLimits,
) -> Result<String, UnixyzError> {
    record.validate_with_limits(limits)?;
    let mut out = String::new();
    append(&mut out, &record.title, limits)?;
    append(&mut out, &record.atoms.len().to_string(), limits)?;
    for atom in &record.atoms {
        append(
            &mut out,
            &format!(
                "{} {}",
                atom.atomic_number,
                atom.position_angstrom.join(" ")
            ),
            limits,
        )?;
    }
    Ok(out)
}
fn append(out: &mut String, line: &str, limits: &UnixyzParseLimits) -> Result<(), UnixyzError> {
    if line.len() > limits.max_line_len
        || out.len().saturating_add(line.len()).saturating_add(1) > limits.max_input_bytes
    {
        return Err(error(0, "serialized byte or line limit exceeded"));
    }
    out.push_str(line);
    out.push('\n');
    Ok(())
}
