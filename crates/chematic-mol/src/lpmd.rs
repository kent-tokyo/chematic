//! Bounded, editable LPMD 2.0 text trajectories. Decimal fields and column order
//! remain authoritative; Cartesian projection is an explicit, checked operation.
use serde::{Deserialize, Serialize};
mod numeric;
mod reader;
use numeric::{number, valid_cell};
pub use reader::{parse_lpmd, parse_lpmd_with_limits};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LpmdError {
    pub line: usize,
    pub detail: String,
}
impl core::fmt::Display for LpmdError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "LPMD line {}: {}", self.line, self.detail)
    }
}
impl std::error::Error for LpmdError {}
fn error(line: usize, detail: impl Into<String>) -> LpmdError {
    LpmdError {
        line,
        detail: detail.into(),
    }
}
#[derive(Debug, Clone, Copy)]
pub struct LpmdParseLimits {
    pub max_input_bytes: usize,
    pub max_line_len: usize,
    /// Total atom rows across all frames, rather than a per-frame allowance.
    pub max_atoms: usize,
    pub max_frames: usize,
}
impl Default for LpmdParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_line_len: 8192,
            max_atoms: 1_000_000,
            max_frames: 1024,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LpmdFrame {
    /// Three cell vectors in row order, in angstroms. No cell reconstruction.
    pub cell: [[String; 3]; 3],
    /// Each row follows the document's column order. Numeric values are strings.
    pub atoms: Vec<Vec<String>>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LpmdDocument {
    pub columns: Vec<String>,
    pub frames: Vec<LpmdFrame>,
}
impl LpmdDocument {
    pub fn validate(&self) -> Result<(), LpmdError> {
        self.validate_with_limits(&LpmdParseLimits::default())
    }
    pub fn validate_with_limits(&self, limits: &LpmdParseLimits) -> Result<(), LpmdError> {
        validate_columns(&self.columns)?;
        if self.frames.is_empty() || self.frames.len() > limits.max_frames {
            return Err(error(0, "frame count outside limits"));
        }
        let mut count = 0usize;
        for frame in &self.frames {
            count = count
                .checked_add(frame.atoms.len())
                .ok_or_else(|| error(0, "atom count overflow"))?;
            if frame.atoms.is_empty() || count > limits.max_atoms {
                return Err(error(0, "total atom count outside limits"));
            }
            valid_cell(&frame.cell)?;
            for row in &frame.atoms {
                validate_row(&self.columns, row)?;
            }
        }
        Ok(())
    }
    /// Project fractional XYZ using row cell vectors, without wrapping or bonds.
    /// No velocity, acceleration, force or charge conversion is implied.
    pub fn positions_angstrom(&self) -> Result<Vec<Vec<[f64; 3]>>, LpmdError> {
        self.positions_angstrom_with_limits(&LpmdParseLimits::default())
    }
    pub fn positions_angstrom_with_limits(
        &self,
        limits: &LpmdParseLimits,
    ) -> Result<Vec<Vec<[f64; 3]>>, LpmdError> {
        self.validate_with_limits(limits)?;
        let indices = ["X", "Y", "Z"].map(|c| {
            self.columns
                .iter()
                .position(|x| x == c)
                .expect("validated columns")
        });
        self.frames
            .iter()
            .map(|frame| {
                let mut cell = [[0.; 3]; 3];
                for (i, vector) in frame.cell.iter().enumerate() {
                    for (j, token) in vector.iter().enumerate() {
                        cell[i][j] = number(token)?;
                    }
                }
                frame
                    .atoms
                    .iter()
                    .map(|row| {
                        let mut position = [0.; 3];
                        for (j, result) in position.iter_mut().enumerate() {
                            let mut terms = [0.; 3];
                            for i in 0..3 {
                                let fraction = number(&row[indices[i]])?;
                                terms[i] = fraction * cell[i][j];
                                if !terms[i].is_finite()
                                    || (fraction != 0. && cell[i][j] != 0. && terms[i] == 0.)
                                {
                                    return Err(error(
                                        0,
                                        "Cartesian product overflows or underflows to zero",
                                    ));
                                }
                            }
                            *result = (terms[0] + terms[1]) + terms[2];
                            if !result.is_finite() {
                                return Err(error(0, "Cartesian sum overflows"));
                            }
                        }
                        Ok(position)
                    })
                    .collect()
            })
            .collect()
    }
}
fn validate_columns(columns: &[String]) -> Result<(), LpmdError> {
    for (i, column) in columns.iter().enumerate() {
        if !matches!(
            column.as_str(),
            "SYM"
                | "X"
                | "Y"
                | "Z"
                | "VX"
                | "VY"
                | "VZ"
                | "AX"
                | "AY"
                | "AZ"
                | "FX"
                | "FY"
                | "FZ"
                | "CHG"
                | "PHG"
                | "FHG"
                | "RGB"
                | "rgb"
                | "C"
        ) {
            return Err(error(0, "unsupported column"));
        }
        if columns[..i].iter().any(|c| c.eq_ignore_ascii_case(column)) {
            return Err(error(0, "duplicate column or RGB/rgb alias"));
        }
    }
    for required in ["SYM", "X", "Y", "Z"] {
        if !columns.iter().any(|c| c == required) {
            return Err(error(0, "missing SYM or complete XYZ columns"));
        }
    }
    for group in [["VX", "VY", "VZ"], ["AX", "AY", "AZ"], ["FX", "FY", "FZ"]] {
        let n = group
            .iter()
            .filter(|g| columns.iter().any(|c| c == **g))
            .count();
        if n != 0 && n != 3 {
            return Err(error(0, "incomplete vector columns"));
        }
    }
    Ok(())
}
fn validate_row(columns: &[String], row: &[String]) -> Result<(), LpmdError> {
    if columns.len() != row.len() {
        return Err(error(0, "atom row width differs from HDR"));
    }
    for (column, token) in columns.iter().zip(row) {
        match column.as_str() {
            "SYM" => {
                if !chematic_core::Element::from_symbol(token).is_some_and(|e| e.symbol() == token)
                {
                    return Err(error(0, "expected canonical element symbol"));
                }
            }
            "RGB" | "rgb" => {
                let inner = token
                    .strip_prefix('<')
                    .and_then(|s| s.strip_suffix('>'))
                    .ok_or_else(|| error(0, "expected <r,g,b>"))?;
                let values: Vec<_> = inner.split(',').collect();
                if values.len() != 3 {
                    return Err(error(0, "RGB requires three components"));
                }
                for value in values {
                    color(value)?;
                }
            }
            "C" => {
                color(token)?;
            }
            _ => {
                number(token)?;
            }
        }
    }
    Ok(())
}
fn color(token: &str) -> Result<(), LpmdError> {
    if !(0. ..=1.).contains(&number(token)?) {
        return Err(error(0, "color outside [0,1]"));
    }
    Ok(())
}
/// Editable JSON refuses duplicate keys recursively, unknown fields and numbers
/// in place of authoritative decimal strings.
pub fn parse_lpmd_json(text: &str, limits: &LpmdParseLimits) -> Result<LpmdDocument, LpmdError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "JSON byte limit exceeded"));
    }
    let value: crate::cjson_record::UniqueValue =
        serde_json::from_str(text).map_err(|e| error(0, e.to_string()))?;
    let document: LpmdDocument =
        serde_json::from_value(value.0).map_err(|e| error(0, e.to_string()))?;
    document.validate_with_limits(limits)?;
    Ok(document)
}
pub fn write_lpmd(document: &LpmdDocument) -> Result<String, LpmdError> {
    write_lpmd_with_limits(document, &LpmdParseLimits::default())
}
pub fn write_lpmd_with_limits(
    document: &LpmdDocument,
    limits: &LpmdParseLimits,
) -> Result<String, LpmdError> {
    document.validate_with_limits(limits)?;
    let mut output = String::new();
    let mut line = |text: &str| -> Result<(), LpmdError> {
        check_line(text, 0, limits)?;
        if output.len().saturating_add(text.len()).saturating_add(1) > limits.max_input_bytes {
            return Err(error(0, "text output byte limit exceeded"));
        }
        output.push_str(text);
        output.push('\n');
        Ok(())
    };
    line("LPMD 2.0 L")?;
    line(&format!("HDR {}", document.columns.join(" ")))?;
    for frame in &document.frames {
        line(&frame.atoms.len().to_string())?;
        let cell: Vec<_> = frame.cell.iter().flatten().cloned().collect();
        line(&cell.join(" "))?;
        for row in &frame.atoms {
            line(&row.join(" "))?;
        }
    }
    Ok(output)
}
fn check_line(line: &str, n: usize, limits: &LpmdParseLimits) -> Result<(), LpmdError> {
    if line.len() > limits.max_line_len || line.chars().any(|c| c.is_control() && c != '\t') {
        return Err(error(n, "line limit exceeded or control character"));
    }
    Ok(())
}
