//! Immutable source retention for the single Cartesian job reader profile.
use super::{
    GaussianError, GaussianParseLimits, GjfGeometry, check_input_limits,
    parse_gjf_geometry_with_limits,
};

/// A validated Cartesian projection paired with its complete immutable source.
///
/// Source text is never reconstructed from the projection. Link0 commands,
/// route/title, original atom labels, numeric spelling, comments and following
/// job sections retain their original UTF-8 bytes. Opaque job sections and
/// calculation keywords are not interpreted or validated, and external files
/// are never opened. This is not a general Gaussian job validator.
pub struct GjfSourceRecord {
    source: String,
    geometry: GjfGeometry,
}
impl GjfSourceRecord {
    /// Exact original UTF-8 text, including BOM, line endings and final whitespace.
    pub fn source_text(&self) -> &str {
        &self.source
    }
    /// Immutable geometry projection; coordinates retain input numeric units.
    pub fn geometry(&self) -> &GjfGeometry {
        &self.geometry
    }
}

/// Retain one single-job Cartesian GJF source with the default Gaussian limits.
pub fn parse_gjf_source_record(input: &str) -> Result<GjfSourceRecord, GaussianError> {
    parse_gjf_source_record_with_limits(input, &GaussianParseLimits::default())
}

/// Retain a bounded Cartesian source without regenerating any job sections.
///
/// Accepts LF/CRLF and an optional leading UTF-8 BOM. Bare CR, NUL/other controls
/// (except TAB/LF/CRLF), and Unicode line/paragraph separators are refused.
/// Malformed/unsupported geometry and multiple Link1 jobs share the geometry
/// reader's refusals. Opaque sections are preserved without interpreting them.
pub fn parse_gjf_source_record_with_limits(
    input: &str,
    limits: &GaussianParseLimits,
) -> Result<GjfSourceRecord, GaussianError> {
    // Include the BOM in the byte/line limits before removing it for parsing.
    check_input_limits(input, limits)?;
    check_source_framing(input)?;
    let geometry =
        parse_gjf_geometry_with_limits(input.strip_prefix('\u{feff}').unwrap_or(input), limits)?;
    Ok(GjfSourceRecord {
        source: input.to_owned(),
        geometry,
    })
}

pub(super) fn check_source_framing(input: &str) -> Result<(), GaussianError> {
    for (offset, c) in input.char_indices() {
        if (c == '\r' && !input[offset..].starts_with("\r\n"))
            || (c.is_control() && !matches!(c, '\t' | '\n' | '\r'))
            || matches!(c, '\u{2028}' | '\u{2029}')
        {
            return Err(GaussianError::InvalidCoordinate(
                "unsupported GJF source text framing".into(),
            ));
        }
    }
    Ok(())
}

/// Return the exact immutable source; no geometry or job data is regenerated.
pub fn write_gjf_source_record(record: &GjfSourceRecord) -> String {
    record.source.clone()
}
