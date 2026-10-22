//! Explicit all-job source retention; single-job readers never select a job.
use super::{
    GaussianError, GaussianParseLimits, GjfGeometry, check_input_limits,
    parse_gjf_geometry_with_limits, source::check_source_framing,
};
use std::ops::Range;

/// Whole-document byte, line, section and atom budgets, plus a job count budget.
/// `gaussian.max_atoms` counts atoms across all jobs, not per job.
#[derive(Debug, Clone, Copy)]
pub struct GjfDocumentLimits {
    pub gaussian: GaussianParseLimits,
    pub max_jobs: usize,
}
impl Default for GjfDocumentLimits {
    fn default() -> Self {
        Self {
            gaussian: GaussianParseLimits::default(),
            max_jobs: 128,
        }
    }
}

/// One ordered Cartesian projection and its UTF-8 byte range in the document.
/// The range excludes the following Link1 delimiter but includes job whitespace.
pub struct GjfSourceJob {
    range: Range<usize>,
    geometry: GjfGeometry,
}
impl GjfSourceJob {
    pub fn source_range(&self) -> Range<usize> {
        self.range.clone()
    }
    pub fn geometry(&self) -> &GjfGeometry {
        &self.geometry
    }
}

/// Immutable complete source with every Link1 job checked, without job selection.
/// Every job must supply explicit Cartesian coordinates. Checkpoint-only geometry,
/// Z-matrices and empty/malformed jobs are refused. Opaque sections and route
/// keywords are not validated; references are never opened or bundled.
pub struct GjfSourceDocument {
    source: String,
    jobs: Vec<GjfSourceJob>,
}
impl GjfSourceDocument {
    pub fn source_text(&self) -> &str {
        &self.source
    }
    pub fn jobs(&self) -> &[GjfSourceJob] {
        &self.jobs
    }
    /// Ordered geometry inspection, not a Gaussian job validation report.
    /// Byte offsets index the original UTF-8 source, including any initial BOM.
    pub fn to_json(&self) -> serde_json::Value {
        serde_json::json!({
            "schema": "chematic.gjf-document.v1",
            "scope": "all_explicit_cartesian_jobs; opaque_sections_not_interpreted",
            "jobs": self.jobs.iter().enumerate().map(|(index, job)| {
                let geometry = job.geometry();
                serde_json::json!({
                    "index": index,
                    "source_start": job.range.start, "source_end": job.range.end,
                    "elements": geometry.mol.atoms().map(|(_, a)| a.element.symbol()).collect::<Vec<_>>(),
                    "coords": geometry.coords, "freeze_codes": geometry.freeze_codes,
                    "charge": geometry.charge, "multiplicity": geometry.multiplicity,
                })
            }).collect::<Vec<_>>()
        })
    }
}

pub fn parse_gjf_source_document(input: &str) -> Result<GjfSourceDocument, GaussianError> {
    parse_gjf_source_document_with_limits(input, &GjfDocumentLimits::default())
}

fn invalid(detail: &str) -> GaussianError {
    GaussianError::InvalidCoordinate(detail.into())
}

fn job_ranges(input: &str, limits: &GjfDocumentLimits) -> Result<Vec<Range<usize>>, GaussianError> {
    let mut ranges = Vec::new();
    let (mut start, mut offset, mut sections) = (0, 0, 0);
    let (mut previous_blank, mut in_section) = (false, false);
    for line in input.split_inclusive('\n') {
        let trimmed = line.trim();
        if trimmed.eq_ignore_ascii_case("--link1--") {
            if !previous_blank || input[start..offset].trim().is_empty() {
                return Err(invalid(
                    "Link1 requires a preceding nonempty job and blank line",
                ));
            }
            ranges.push(start..offset);
            start = offset + line.len();
            in_section = false;
            if ranges.len() >= limits.max_jobs {
                return Err(GaussianError::ResourceLimit {
                    resource: "jobs",
                    actual: ranges.len() + 1,
                    limit: limits.max_jobs,
                });
            }
        } else {
            if trimmed.to_ascii_lowercase().starts_with("--link1") {
                return Err(invalid("malformed Link1 delimiter"));
            }
            if !trimmed.is_empty() && !in_section {
                sections += 1;
                if sections > limits.gaussian.max_sections {
                    return Err(GaussianError::ResourceLimit {
                        resource: "document sections",
                        actual: sections,
                        limit: limits.gaussian.max_sections,
                    });
                }
            }
            in_section = !trimmed.is_empty();
        }
        previous_blank = trimmed.is_empty();
        offset += line.len();
    }
    if input[start..].trim().is_empty() {
        return Err(invalid("empty GJF document job"));
    }
    ranges.push(start..input.len());
    if ranges.len() > limits.max_jobs {
        return Err(GaussianError::ResourceLimit {
            resource: "jobs",
            actual: ranges.len(),
            limit: limits.max_jobs,
        });
    }
    Ok(ranges)
}

/// Validate every explicit Cartesian job before returning any document.
/// Link1 delimiters are standalone lines (case insensitive, surrounding
/// whitespace allowed), preceded by a blank line. Only the document may have
/// one leading BOM. All resource budgets are cumulative across jobs.
pub fn parse_gjf_source_document_with_limits(
    input: &str,
    limits: &GjfDocumentLimits,
) -> Result<GjfSourceDocument, GaussianError> {
    check_input_limits(input, &limits.gaussian)?;
    check_source_framing(input)?;
    let ranges = job_ranges(input, limits)?;
    let mut jobs = Vec::new();
    let mut remaining_atoms = limits.gaussian.max_atoms;
    for range in ranges {
        let text = &input[range.clone()];
        let text = if range.start == 0 {
            text.strip_prefix('\u{feff}').unwrap_or(text)
        } else {
            text
        };
        let first = text
            .lines()
            .find(|line| !line.trim().is_empty())
            .unwrap_or("")
            .trim();
        if !first.starts_with(['%', '#']) {
            return Err(invalid("GJF job must start with Link0 or route"));
        }
        let geometry = parse_gjf_geometry_with_limits(
            text,
            &GaussianParseLimits {
                max_atoms: remaining_atoms,
                ..limits.gaussian
            },
        )
        .map_err(|error| match error {
            GaussianError::ResourceLimit {
                resource: "atoms",
                actual,
                ..
            } => GaussianError::ResourceLimit {
                resource: "document atoms",
                actual: limits.gaussian.max_atoms - remaining_atoms + actual,
                limit: limits.gaussian.max_atoms,
            },
            other => other,
        })?;
        remaining_atoms -= geometry.mol.atom_count();
        jobs.push(GjfSourceJob { range, geometry });
    }
    Ok(GjfSourceDocument {
        source: input.to_owned(),
        jobs,
    })
}

/// Return the whole original source, including every delimiter and job.
pub fn write_gjf_source_document(document: &GjfSourceDocument) -> String {
    document.source.clone()
}
