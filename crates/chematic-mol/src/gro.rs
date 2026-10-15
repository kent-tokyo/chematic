//! Single-frame GROMACS GRO records, implemented from the GROMACS 2026.4 manual.
//! Atom names are not element symbols. No element or bond inference is performed.
use chematic_core::Element;

#[derive(Debug, Clone, Copy)]
pub struct GroParseLimits {
    pub max_input_bytes: usize,
    pub max_atoms: usize,
    pub max_line_len: usize,
}
impl Default for GroParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 << 20,
            max_atoms: 2_000_000,
            max_line_len: 1024,
        }
    }
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GroError {
    InputTooLarge { limit: usize },
    LineTooLong { line: usize, limit: usize },
    TooManyAtoms { limit: usize },
    InvalidField { line: usize, field: &'static str },
    MissingLine { line: usize, field: &'static str },
    TrailingData { line: usize },
    ElementCount { atoms: usize, elements: usize },
    AllocationFailed,
}
impl core::fmt::Display for GroError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InputTooLarge { limit } => write!(f, "GRO input exceeds {limit} bytes"),
            Self::LineTooLong { line, limit } => write!(f, "GRO line {line} exceeds {limit} bytes"),
            Self::TooManyAtoms { limit } => write!(f, "GRO atom count exceeds {limit}"),
            Self::InvalidField { line, field } => write!(f, "GRO line {line}: invalid {field}"),
            Self::MissingLine { line, field } => write!(f, "GRO line {line}: missing {field}"),
            Self::TrailingData { line } => write!(
                f,
                "GRO line {line}: additional frame or trailing data; one frame required"
            ),
            Self::ElementCount { atoms, elements } => write!(
                f,
                "GRO has {atoms} atoms but {elements} elements were supplied"
            ),
            Self::AllocationFailed => write!(f, "could not allocate GRO atom records"),
        }
    }
}
impl std::error::Error for GroError {}

#[derive(Debug, Clone, PartialEq)]
pub struct GroAtomRecord {
    pub residue_number: u32,
    pub residue_name: String,
    pub atom_name: String,
    pub atom_number: u32,
    pub position_nm: [f64; 3],
    pub velocity_nm_per_ps: Option<[f64; 3]>,
    /// Unknown unless explicitly supplied by the caller, never guessed from a name.
    pub element: Option<Element>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct GroFrame {
    pub title: String,
    pub atoms: Vec<GroAtomRecord>,
    /// Original three or nine values in GRO file order, in nm.
    pub box_values_nm: Vec<f64>,
    /// Row vectors v1, v2, v3, including zero off-diagonals omitted in a 3-value box.
    pub box_vectors_nm: [[f64; 3]; 3],
    /// Exact original input, including numeric spelling and line endings.
    pub source_text: String,
}
impl GroFrame {
    /// Resolve atom elements using authoritative caller data in source row order.
    /// Returns an independent copy. No name heuristic or topology is introduced.
    pub fn with_elements(&self, elements: &[Element]) -> Result<Self, GroError> {
        if elements.len() != self.atoms.len() {
            return Err(GroError::ElementCount {
                atoms: self.atoms.len(),
                elements: elements.len(),
            });
        }
        let mut result = self.clone();
        for (atom, &element) in result.atoms.iter_mut().zip(elements) {
            atom.element = Some(element);
        }
        Ok(result)
    }
}
fn invalid(line: usize, field: &'static str) -> GroError {
    GroError::InvalidField { line, field }
}
fn unsigned(value: &str, line: usize, field: &'static str) -> Result<usize, GroError> {
    let value = value.trim();
    if value.is_empty() || !value.bytes().all(|b| b.is_ascii_digit()) {
        return Err(invalid(line, field));
    }
    value.parse().map_err(|_| invalid(line, field))
}
fn name(value: &str, line: usize, field: &'static str) -> Result<String, GroError> {
    let value = value.trim();
    if value.is_empty() || !value.bytes().all(|b| b.is_ascii_graphic()) {
        return Err(invalid(line, field));
    }
    Ok(value.to_owned())
}
fn fixed_real(
    value: &str,
    decimals: usize,
    line: usize,
    field: &'static str,
) -> Result<f64, GroError> {
    let value = value.trim();
    let Some((whole, fraction)) = value.split_once('.') else {
        return Err(invalid(line, field));
    };
    let whole = whole
        .strip_prefix('-')
        .or_else(|| whole.strip_prefix('+'))
        .unwrap_or(whole);
    if whole.is_empty()
        || !whole.bytes().all(|b| b.is_ascii_digit())
        || fraction.len() != decimals
        || !fraction.bytes().all(|b| b.is_ascii_digit())
    {
        return Err(invalid(line, field));
    }
    let v: f64 = value.parse().map_err(|_| invalid(line, field))?;
    if !v.is_finite() {
        return Err(invalid(line, field));
    }
    Ok(v)
}
fn triple(
    value: &str,
    width: usize,
    decimals: usize,
    line: usize,
    field: &'static str,
) -> Result<[f64; 3], GroError> {
    Ok([
        fixed_real(&value[..width], decimals, line, field)?,
        fixed_real(&value[width..width * 2], decimals, line, field)?,
        fixed_real(&value[width * 2..width * 3], decimals, line, field)?,
    ])
}
fn atom_record(row: &str, line: usize) -> Result<GroAtomRecord, GroError> {
    let row = row.trim_end_matches([' ', '\t']);
    if !row.is_ascii() || row.len() < 20 {
        return Err(invalid(line, "fixed atom row"));
    }
    let numbers = &row[20..];
    let dots: Vec<_> = numbers
        .bytes()
        .enumerate()
        .filter_map(|(i, b)| (b == b'.').then_some(i))
        .take(3)
        .collect();
    if dots.len() != 3 {
        return Err(invalid(line, "coordinate precision"));
    }
    let width = dots[1] - dots[0];
    if width < 6 || dots[2] - dots[1] != width {
        return Err(invalid(line, "coordinate precision"));
    }
    if numbers.len() != 3 * width && numbers.len() != 6 * width {
        return Err(invalid(line, "coordinate/velocity field count"));
    }
    let position_nm = triple(&numbers[..width * 3], width, width - 5, line, "coordinates")?;
    let velocity_nm_per_ps = if numbers.len() == 6 * width {
        Some(triple(
            &numbers[width * 3..],
            width,
            width - 4,
            line,
            "velocities",
        )?)
    } else {
        None
    };
    Ok(GroAtomRecord {
        residue_number: unsigned(&row[..5], line, "residue number")? as u32,
        residue_name: name(&row[5..10], line, "residue name")?,
        atom_name: name(&row[10..15], line, "atom name")?,
        atom_number: unsigned(&row[15..20], line, "atom number")? as u32,
        position_nm,
        velocity_nm_per_ps,
        element: None,
    })
}
fn box_values(row: &str, line: usize) -> Result<(Vec<f64>, [[f64; 3]; 3]), GroError> {
    let values = row
        .split_whitespace()
        .map(|s| s.parse::<f64>().map_err(|_| invalid(line, "box vectors")))
        .collect::<Result<Vec<_>, _>>()?;
    if !matches!(values.len(), 3 | 9) || values.iter().any(|v| !v.is_finite()) {
        return Err(invalid(line, "box vectors"));
    }
    let mut m = [
        [values[0], 0.0, 0.0],
        [0.0, values[1], 0.0],
        [0.0, 0.0, values[2]],
    ];
    if values.len() == 9 {
        m[0][1] = values[3];
        m[0][2] = values[4];
        m[1][0] = values[5];
        m[1][2] = values[6];
        m[2][0] = values[7];
        m[2][1] = values[8];
    }
    Ok((values, m))
}
fn required_line<'a>(
    lines: &mut std::str::Lines<'a>,
    line: usize,
    field: &'static str,
    limits: &GroParseLimits,
) -> Result<&'a str, GroError> {
    let value = lines.next().ok_or(GroError::MissingLine { line, field })?;
    if value.len() > limits.max_line_len {
        return Err(GroError::LineTooLong {
            line,
            limit: limits.max_line_len,
        });
    }
    Ok(value)
}
/// Read exactly one GRO frame, keeping nm and nm/ps units and unresolved elements.
/// Supports decimal-inferred n+5 fields (n>=1), optional velocities and 3/9-value
/// box records. Additional nonblank frames/data are refused, never ignored.
pub fn parse_gro(text: &str) -> Result<GroFrame, GroError> {
    parse_gro_with_limits(text, &GroParseLimits::default())
}
pub fn parse_gro_with_limits(text: &str, limits: &GroParseLimits) -> Result<GroFrame, GroError> {
    if text.len() > limits.max_input_bytes {
        return Err(GroError::InputTooLarge {
            limit: limits.max_input_bytes,
        });
    }
    let mut lines = text.lines();
    let title = required_line(&mut lines, 1, "title", limits)?.to_owned();
    let count = unsigned(
        required_line(&mut lines, 2, "atom count", limits)?,
        2,
        "atom count",
    )?;
    if count > limits.max_atoms {
        return Err(GroError::TooManyAtoms {
            limit: limits.max_atoms,
        });
    }
    if lines.clone().count() <= count {
        return Err(GroError::MissingLine {
            line: count.saturating_add(3),
            field: "atom/box records",
        });
    }
    let mut atoms = Vec::new();
    atoms
        .try_reserve_exact(count)
        .map_err(|_| GroError::AllocationFailed)?;
    for i in 0..count {
        atoms.push(atom_record(
            required_line(&mut lines, i + 3, "atom record", limits)?,
            i + 3,
        )?);
    }
    let (box_values_nm, box_vectors_nm) = box_values(
        required_line(&mut lines, count + 3, "box record", limits)?,
        count + 3,
    )?;
    for (i, row) in lines.enumerate() {
        if row.len() > limits.max_line_len {
            return Err(GroError::LineTooLong {
                line: i + count + 4,
                limit: limits.max_line_len,
            });
        }
        if !row.trim().is_empty() {
            return Err(GroError::TrailingData {
                line: i + count + 4,
            });
        }
    }
    Ok(GroFrame {
        title,
        atoms,
        box_values_nm,
        box_vectors_nm,
        source_text: text.to_owned(),
    })
}
