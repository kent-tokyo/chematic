//! Bounded, explicit Cartesian projection; no molecular graph is accepted.
use chematic_core::Element;

/// Units of the supplied coordinate numbers; values are never rescaled.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GjfCoordinateUnits {
    Angstrom,
    Bohr,
}

/// A refused Cartesian output. Existing Gaussian parser errors are unchanged.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GjfWriteError {
    pub detail: String,
}
impl core::fmt::Display for GjfWriteError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "GJF Cartesian output: {}", self.detail)
    }
}
impl std::error::Error for GjfWriteError {}
fn error(detail: impl Into<String>) -> GjfWriteError {
    GjfWriteError {
        detail: detail.into(),
    }
}

/// Explicit electronic state, coordinate units and single-line job headers.
/// The default route is a template; supply calculation keywords before execution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GjfWriteOptions {
    pub charge: i32,
    pub multiplicity: u32,
    pub units: GjfCoordinateUnits,
    pub title: String,
    /// Starts with #. Units keywords (including abbreviations) are refused;
    /// the writer appends one explicit Units=Ang or Units=AU keyword.
    pub route: String,
}
impl GjfWriteOptions {
    pub const MAX_ATOMS: usize = 100_000;
    pub fn new(charge: i32, multiplicity: u32, units: GjfCoordinateUnits) -> Self {
        Self {
            charge,
            multiplicity,
            units,
            title: "chematic Cartesian geometry".into(),
            route: "#".into(),
        }
    }
    fn validate(&self) -> Result<(), GjfWriteError> {
        if self.multiplicity == 0 {
            return Err(error("multiplicity must be positive"));
        }
        for field in [&self.title, &self.route] {
            if field.trim().is_empty()
                || field.len() > 1024
                || field
                    .chars()
                    .any(|c| c.is_control() || matches!(c, '\u{2028}' | '\u{2029}'))
            {
                return Err(error(
                    "headers must be nonempty single lines of at most 1024 bytes",
                ));
            }
        }
        if !self.route.starts_with('#')
            || !self.route.is_ascii()
            || self.route[1..].contains(['#', '!'])
        {
            return Err(error(
                "route must be ASCII, start with one #, and contain no ! comments",
            ));
        }
        if self
            .route
            .split(|c: char| !c.is_ascii_alphabetic())
            .filter(|s| !s.is_empty())
            .any(|word| "units".starts_with(&word.to_ascii_lowercase()))
        {
            return Err(error(
                "route must not override coordinate Units; use the units option",
            ));
        }
        Ok(())
    }
}

fn scalar(value: f64) -> Result<String, GjfWriteError> {
    if !value.is_finite() {
        return Err(error("coordinate values must be finite"));
    }
    let text = format!("{value:.16e}");
    if text.parse::<f64>().map(f64::to_bits).ok() != Some(value.to_bits()) {
        return Err(error(
            "coordinate cannot be serialized without changing its f64 value",
        ));
    }
    Ok(text)
}

/// Write explicit ordered elements, coordinate values and optional freeze codes.
///
/// All arrays must have identical lengths, with 1..=100000 atoms. Every finite
/// f64 is serialized with 17 significant digits and checked by bitwise readback,
/// including subnormals and signed zero. No missing values are synthesized.
/// Units=Ang/AU is emitted from `options.units`; coordinates are not converted.
/// This geometry-only API does not serialize a graph, isotopes or a complete job.
/// Header checks are lexical, not a Gaussian route grammar or execution check.
pub fn write_gjf_cartesian(
    elements: &[Element],
    coords: &[[f64; 3]],
    freeze_codes: &[Option<i32>],
    options: &GjfWriteOptions,
) -> Result<String, GjfWriteError> {
    options.validate()?;
    if elements.is_empty() || elements.len() > GjfWriteOptions::MAX_ATOMS {
        return Err(error("atom count must be 1..=100000"));
    }
    if coords.len() != elements.len() || freeze_codes.len() != elements.len() {
        return Err(error(
            "element, coordinate and freeze-code counts must match",
        ));
    }
    let units = match options.units {
        GjfCoordinateUnits::Angstrom => "Ang",
        GjfCoordinateUnits::Bohr => "AU",
    };
    let mut output = format!(
        "{} Units={}\n\n{}\n\n{} {}\n",
        options.route, units, options.title, options.charge, options.multiplicity
    );
    for ((element, position), code) in elements.iter().zip(coords).zip(freeze_codes) {
        if element.atomic_number() == 0 {
            return Err(error(
                "wildcard/dummy atoms are outside the elemental Cartesian profile",
            ));
        }
        output.push_str(element.symbol());
        if let Some(code) = code {
            output.push_str(&format!(" {code}"));
        }
        for &value in position {
            output.push(' ');
            output.push_str(&scalar(value)?);
        }
        output.push('\n');
    }
    output.push('\n');
    Ok(output)
}
