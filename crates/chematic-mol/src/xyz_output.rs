//! Checked plain XYZ geometry output with caller-selected decimal formatting.
//! Independent output profile for Open Babel #1693; legacy writers are unchanged.
use crate::XyzAtom;

pub const MAX_XYZ_OUTPUT_ATOMS: usize = 100_000;
pub const MAX_XYZ_OUTPUT_BYTES: usize = 64 << 20;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XyzNotation {
    Fixed,
    Scientific,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct XyzWriteOptions {
    /// Number of digits after the decimal point, in either notation: 0..=17.
    pub decimal_places: usize,
    /// Minimum width per coordinate, 0..=64. Large values expand the field.
    pub min_width: usize,
    pub notation: XyzNotation,
    /// Aggregate document bytes, including header and padding; 1..=64 MiB.
    pub max_output_bytes: usize,
}
impl Default for XyzWriteOptions {
    fn default() -> Self {
        Self {
            decimal_places: 12,
            min_width: 20,
            notation: XyzNotation::Fixed,
            max_output_bytes: MAX_XYZ_OUTPUT_BYTES,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct XyzWriteError {
    /// Zero-based row, or None for header/options/allocation failures.
    pub atom: Option<usize>,
    pub detail: String,
}
impl core::fmt::Display for XyzWriteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "XYZ output {:?}: {}", self.atom, self.detail)
    }
}
impl std::error::Error for XyzWriteError {}
fn error(atom: Option<usize>, detail: impl Into<String>) -> XyzWriteError {
    XyzWriteError {
        atom,
        detail: detail.into(),
    }
}
impl XyzWriteOptions {
    /// Validate before collecting caller geometry. Only plain ASCII comments;
    /// '=' is refused to avoid advertising an unsupported extxyz schema/cell.
    pub fn validate(&self, comment: &str) -> Result<(), XyzWriteError> {
        if self.decimal_places > 17 || self.min_width > 64 {
            return Err(error(
                None,
                "decimal_places must be 0..=17 and min_width 0..=64",
            ));
        }
        if self.max_output_bytes == 0 || self.max_output_bytes > MAX_XYZ_OUTPUT_BYTES {
            return Err(error(None, "max_output_bytes must be 1..=64 MiB"));
        }
        if comment.len() > 4096
            || comment
                .bytes()
                .any(|b| !(b' '..=b'~').contains(&b) || b == b'=')
        {
            return Err(error(
                None,
                "comment must be <=4096 printable ASCII bytes without '='",
            ));
        }
        Ok(())
    }
}
fn append(out: &mut String, text: &str, limit: usize) -> Result<(), XyzWriteError> {
    if out.len().checked_add(text.len()).is_none_or(|n| n > limit) {
        return Err(error(None, "aggregate output byte limit exceeded"));
    }
    out.try_reserve(text.len())
        .map_err(|e| error(None, format!("output allocation failed: {e}")))?;
    out.push_str(text);
    Ok(())
}
fn coordinate(value: f64, options: &XyzWriteOptions, row: usize) -> Result<String, XyzWriteError> {
    let width = options.min_width;
    let precision = options.decimal_places;
    let text = match options.notation {
        XyzNotation::Fixed => format!("{value:>width$.precision$}"),
        XyzNotation::Scientific => format!("{value:>width$.precision$e}"),
    };
    if !text.trim().parse::<f64>().is_ok_and(f64::is_finite) {
        return Err(error(Some(row), "rounded coordinate is not finite"));
    }
    Ok(text)
}

/// Write exactly the supplied atom identities and Cartesian angstroms. Empty
/// geometry is allowed. No graph/H/bond/electronic-state inference or metadata
/// projection. Fixed precision is absolute decimal rounding; scientific with
/// at least 16 fractional digits round-trips finite f64 coordinates (including -0).
/// The minimum field width never truncates a value, and tokens always have an
/// additional separating space. Errors return no partial document.
pub fn write_xyz_checked(
    atoms: &[XyzAtom],
    comment: &str,
    options: &XyzWriteOptions,
) -> Result<String, XyzWriteError> {
    options.validate(comment)?;
    if atoms.len() > MAX_XYZ_OUTPUT_ATOMS {
        return Err(error(None, "atom count exceeds 100000"));
    }
    for (i, a) in atoms.iter().enumerate() {
        if ![a.x, a.y, a.z].iter().all(|v| v.is_finite()) {
            return Err(error(Some(i), "coordinates must be finite"));
        }
    }
    let mut out = String::new();
    append(
        &mut out,
        &format!("{}\n{comment}\n", atoms.len()),
        options.max_output_bytes,
    )?;
    for (i, atom) in atoms.iter().enumerate() {
        let line = format!(
            "{:<2} {} {} {}\n",
            atom.element.symbol(),
            coordinate(atom.x, options, i)?,
            coordinate(atom.y, options, i)?,
            coordinate(atom.z, options, i)?
        );
        append(&mut out, &line, options.max_output_bytes)?;
    }
    Ok(out)
}
