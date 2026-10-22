use super::{MolParseError, SdfParseLimits};

/// Bounds for immutable SDF retention. Byte/line/record checks precede parsing;
/// atom counts are checked on each bounded native projection. Native CTAB caps
/// also apply, independently of these caller-selected limits.
#[derive(Debug, Clone, Copy)]
pub struct SdfSourceLimits {
    pub sdf: SdfParseLimits,
    pub max_atoms_per_record: usize,
    pub max_total_atoms: usize,
}
impl Default for SdfSourceLimits {
    fn default() -> Self {
        Self {
            sdf: SdfParseLimits {
                max_input_bytes: 16 * 1024 * 1024,
                max_record_bytes: 16 * 1024 * 1024,
                max_line_bytes: 1024 * 1024,
                max_records: 10_000,
            },
            max_atoms_per_record: 100_000,
            max_total_atoms: 100_000,
        }
    }
}

/// Counts from the bounded native projection, not a full chemistry validation.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SdfSourceSummary {
    pub records: usize,
    pub total_atoms: usize,
}

fn bound(resource: &'static str, actual: usize, limit: usize) -> Result<(), MolParseError> {
    if actual > limit {
        return Err(MolParseError::ResourceLimit {
            resource,
            actual,
            limit,
        });
    }
    Ok(())
}

fn ctab_end(block: &str) -> Result<usize, MolParseError> {
    let mut end = 0;
    for (i, line) in block.split_inclusive('\n').enumerate() {
        end += line.len();
        if i >= 4 && content(line) == "M  END" {
            return Ok(end);
        }
    }
    Err(MolParseError::UnexpectedEnd)
}

fn content(line: &str) -> &str {
    let line = line.strip_suffix('\n').unwrap_or(line);
    line.strip_suffix('\r').unwrap_or(line)
}

fn record(
    block: &str,
    limits: &SdfSourceLimits,
    summary: &mut SdfSourceSummary,
) -> Result<(), MolParseError> {
    if block.trim().is_empty() {
        return Err(MolParseError::InvalidHeader {
            line: 1,
            detail: "empty SDF record is outside the source-retention profile".into(),
        });
    }
    bound("record bytes", block.len(), limits.sdf.max_record_bytes)?;
    bound(
        "records",
        summary.records.saturating_add(1),
        limits.sdf.max_records,
    )?;
    ctab_end(block)?;
    let report = super::read_sdf_mol_with_diagnostics(block)?;
    let atoms = report.mol.atom_count();
    bound("atoms per record", atoms, limits.max_atoms_per_record)?;
    let total = summary.total_atoms.saturating_add(atoms);
    bound("total atoms", total, limits.max_total_atoms)?;
    summary.records += 1;
    summary.total_atoms = total;
    Ok(())
}

/// Validate a nonempty local-profile SDF document without retaining graphs.
/// Exact physical `$$$$` lines separate records. Empty records are refused;
/// whitespace after the final delimiter and a final unterminated record are
/// accepted. Each record needs an exact `M  END` and a native V2000/V3000 core
/// projection. SD fields/opaque chemistry/stereo are not verified or edited.
pub fn validate_sdf_source_with_limits(
    input: &str,
    limits: &SdfSourceLimits,
) -> Result<SdfSourceSummary, MolParseError> {
    bound("input bytes", input.len(), limits.sdf.max_input_bytes)?;
    let mut summary = SdfSourceSummary::default();
    let (mut start, mut offset) = (0, 0);
    for line in input.split_inclusive('\n') {
        bound("line bytes", line.len(), limits.sdf.max_line_bytes)?;
        if content(line) == "$$$$" {
            record(&input[start..offset], limits, &mut summary)?;
            start = offset + line.len();
        }
        offset += line.len();
    }
    if !input[start..].trim().is_empty() {
        record(&input[start..], limits, &mut summary)?;
    }
    if summary.records == 0 {
        return Err(MolParseError::UnexpectedEnd);
    }
    Ok(summary)
}

/// Retain the exact UTF-8 source of every validated record, including SD data,
/// coordinates, order, headers and line endings. No editable chemistry or
/// interpreted metadata/stereo fidelity is promised. Nothing is returned if
/// any record fails or exceeds bounds; this does not return a valid prefix.
pub fn roundtrip_sdf_source(input: &str) -> Result<String, MolParseError> {
    roundtrip_sdf_source_with_limits(input, &SdfSourceLimits::default())
}

pub fn roundtrip_sdf_source_with_limits(
    input: &str,
    limits: &SdfSourceLimits,
) -> Result<String, MolParseError> {
    validate_sdf_source_with_limits(input, limits)?;
    Ok(input.to_owned())
}

/// Extract one unchanged CTAB for a single-molecule bridge, refusing additional
/// SDF records instead of silently returning the first. SD fields and delimiter
/// are intentionally omitted; use source retention to retain them.
pub fn single_sdf_mol_block(input: &str) -> Result<&str, MolParseError> {
    single_sdf_mol_block_with_limits(input, &SdfSourceLimits::default())
}

pub fn single_sdf_mol_block_with_limits<'a>(
    input: &'a str,
    limits: &SdfSourceLimits,
) -> Result<&'a str, MolParseError> {
    let single = SdfSourceLimits {
        sdf: SdfParseLimits {
            max_records: limits.sdf.max_records.min(1),
            ..limits.sdf
        },
        ..*limits
    };
    validate_sdf_source_with_limits(input, &single)?;
    Ok(&input[..ctab_end(input)?])
}
