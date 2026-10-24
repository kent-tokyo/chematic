//! Checked ViewMol geometry: $coord factor, optional $title, fixed-axis flags.
//! Raw coordinates and scale stay authoritative decimal strings. Unrepresented
//! vibrational, cell, orbital and trajectory groups refuse rather than disappear.
use serde::{Deserialize, Serialize};
mod reader;
mod writer;
pub use reader::{parse_viewmol, parse_viewmol_with_limits};
pub use writer::{write_viewmol, write_viewmol_with_limits};
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ViewmolError {
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for ViewmolError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ViewMol line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for ViewmolError {}
fn error(line: usize, detail: impl Into<String>) -> ViewmolError {
    ViewmolError {
        line,
        detail: detail.into(),
    }
}
#[derive(Debug, Clone, Copy)]
pub struct ViewmolParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    pub max_atoms: usize,
}
impl Default for ViewmolParseLimits {
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
pub struct ViewmolAtom {
    /// Canonical element symbol; no label or pseudoatom inference.
    pub element: String,
    /// Coordinates in source units; multiply by factor_to_angstrom to project.
    pub position: [String; 3],
    /// Optional nonempty combination of distinct x/y/z letters; order retained.
    pub fixed_axes: Option<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ViewmolRecord {
    /// None is absent; Some("") retains an explicitly empty title group.
    pub title: Option<String>,
    pub factor_to_angstrom: String,
    pub atoms: Vec<ViewmolAtom>,
}
impl ViewmolRecord {
    pub fn validate(&self) -> Result<(), ViewmolError> {
        self.validate_with_limits(&ViewmolParseLimits::default())
    }
    pub fn validate_with_limits(&self, limits: &ViewmolParseLimits) -> Result<(), ViewmolError> {
        if let Some(title) = &self.title {
            check_line(title, 0, limits)?;
            if title.trim_start().starts_with('$') {
                return Err(error(0, "title cannot begin with a data-group marker"));
            }
        }
        if number(&self.factor_to_angstrom)? <= 0. {
            return Err(error(0, "coordinate factor must be finite and positive"));
        }
        if self.atoms.is_empty() || self.atoms.len() > limits.max_atoms {
            return Err(error(0, "atom count outside limits"));
        }
        for atom in &self.atoms {
            element_number(&atom.element)?;
            for v in &atom.position {
                number(v)?;
            }
            validate_axes(atom.fixed_axes.as_deref())?;
        }
        Ok(())
    }
    /// Explicit projection only. Raw parse/write never changes the source scale.
    pub fn positions_angstrom(&self) -> Result<Vec<[f64; 3]>, ViewmolError> {
        self.positions_angstrom_with_limits(&ViewmolParseLimits::default())
    }
    pub fn positions_angstrom_with_limits(
        &self,
        limits: &ViewmolParseLimits,
    ) -> Result<Vec<[f64; 3]>, ViewmolError> {
        self.validate_with_limits(limits)?;
        let factor = number(&self.factor_to_angstrom)?;
        self.atoms
            .iter()
            .map(|atom| {
                let mut result = [0.; 3];
                for (i, token) in atom.position.iter().enumerate() {
                    let raw = number(token)?;
                    let scaled = raw * factor;
                    if !scaled.is_finite() || (raw != 0. && scaled == 0.) {
                        return Err(error(
                            0,
                            "Angstrom projection overflows or underflows to zero",
                        ));
                    }
                    result[i] = scaled;
                }
                Ok(result)
            })
            .collect()
    }
}
fn validate_axes(axes: Option<&str>) -> Result<(), ViewmolError> {
    if let Some(axes) = axes {
        let mut seen = 0u8;
        if axes.is_empty() {
            return Err(error(
                0,
                "fixed_axes must be absent or a nonempty x/y/z combination",
            ));
        }
        for c in axes.bytes() {
            let bit = match c {
                b'x' => 1,
                b'y' => 2,
                b'z' => 4,
                _ => return Err(error(0, "unsupported fixed-axis flag")),
            };
            if seen & bit != 0 {
                return Err(error(0, "duplicate fixed-axis flag"));
            }
            seen |= bit;
        }
    }
    Ok(())
}
fn element_number(symbol: &str) -> Result<u8, ViewmolError> {
    chematic_core::Element::from_symbol(symbol)
        .filter(|e| e.symbol() == symbol)
        .map(|e| e.atomic_number())
        .ok_or_else(|| error(0, "expected canonical element symbol"))
}
/// Strict editable-record JSON with recursively unique keys and known fields.
pub fn parse_viewmol_json(
    text: &str,
    limits: &ViewmolParseLimits,
) -> Result<ViewmolRecord, ViewmolError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "JSON byte limit exceeded"));
    }
    let value: crate::cjson_record::UniqueValue =
        serde_json::from_str(text).map_err(|e| error(0, e.to_string()))?;
    let record: ViewmolRecord =
        serde_json::from_value(value.0).map_err(|e| error(0, e.to_string()))?;
    record.validate_with_limits(limits)?;
    Ok(record)
}
fn number(token: &str) -> Result<f64, ViewmolError> {
    if token.len() > 256 || !decimal(token) {
        return Err(error(0, "invalid decimal or token exceeds 256 bytes"));
    }
    let v: f64 = token
        .parse()
        .map_err(|_| error(0, "invalid nominal value"))?;
    let mantissa = token.split(['e', 'E']).next().unwrap_or("");
    if !v.is_finite() || (v == 0. && mantissa.bytes().any(|b| matches!(b, b'1'..=b'9'))) {
        return Err(error(0, "non-finite or underflowed decimal"));
    }
    Ok(v)
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
fn check_line(line: &str, n: usize, limits: &ViewmolParseLimits) -> Result<(), ViewmolError> {
    if line.len() > limits.max_line_len || line.chars().any(|c| c.is_control() && c != '\t') {
        return Err(error(n, "line limit exceeded or control character"));
    }
    Ok(())
}
