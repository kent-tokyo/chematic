//! Immutable single-block V3000 retention with a bounded graph projection.
use super::{collect_v30_lines, v3k_err, validate_v3000_physical_lines};
use crate::MolParseError;
use std::collections::HashSet;

/// Resource bounds for an immutable single V3000 MOL block.
/// The native reader's existing 100,000-atom cap also applies.
#[derive(Debug, Clone, Copy)]
pub struct V3000SourceLimits {
    pub max_input_bytes: usize,
    pub max_atoms: usize,
    pub max_lines: usize,
    pub max_line_bytes: usize,
}

impl Default for V3000SourceLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 16 * 1024 * 1024,
            max_atoms: 100_000,
            max_lines: 1_000_000,
            max_line_bytes: 1024 * 1024,
        }
    }
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

fn validate_ids(
    lines: &[super::LogicalLine<'_>],
    limits: &V3000SourceLimits,
) -> Result<(), MolParseError> {
    let mut atoms = HashSet::new();
    let mut bonds = HashSet::new();
    let mut section = None;
    for line in lines {
        let mut words = line.payload.split_whitespace();
        let (first, second) = (words.next(), words.next());
        match (first, second) {
            (Some("BEGIN"), Some("ATOM" | "BOND")) => {
                section = second;
                continue;
            }
            (Some("END"), Some("ATOM" | "BOND")) => {
                section = None;
                continue;
            }
            _ => {}
        }
        let ids = match section {
            Some("ATOM") => &mut atoms,
            Some("BOND") => &mut bonds,
            _ => continue,
        };
        let id = first.and_then(|s| s.parse::<u32>().ok());
        if !id.is_some_and(|id| id != 0 && ids.insert(id)) {
            return Err(v3k_err(
                line.line_num,
                "source atom/bond IDs must be positive and unique",
            ));
        }
        if section == Some("ATOM") {
            bound("V3000 atoms", atoms.len(), limits.max_atoms)?;
        }
    }
    Ok(())
}

/// Validate one local-profile V3000 block and return its exact UTF-8 source.
/// Coordinates, IDs, line endings and opaque fields survive without edits,
/// normalization or verification of opaque chemistry/stereo. Multiple CTABs,
/// trailing records and malformed core structures are rejected.
pub fn roundtrip_mol_v3000_source(input: &str) -> Result<String, MolParseError> {
    roundtrip_mol_v3000_source_with_limits(input, &V3000SourceLimits::default())
}

/// Source retention with caller-selected byte/line/atom bounds. Byte/line checks
/// precede graph parsing; atom IDs/counts are checked on source rows and the
/// native projection. Existing native graph/count caps also apply.
pub fn roundtrip_mol_v3000_source_with_limits(
    input: &str,
    limits: &V3000SourceLimits,
) -> Result<String, MolParseError> {
    bound("V3000 input bytes", input.len(), limits.max_input_bytes)?;
    for (index, line) in input.lines().enumerate() {
        bound("V3000 physical lines", index + 1, limits.max_lines)?;
        bound("V3000 line bytes", line.len(), limits.max_line_bytes)?;
    }
    let lines = validate_v3000_physical_lines(input)?;
    for pair in lines.windows(2).skip(4) {
        if pair[0].1.starts_with(super::V30_PREFIX)
            && pair[0].1.ends_with('-')
            && !pair[1].1.starts_with(super::V30_PREFIX)
        {
            return Err(v3k_err(
                pair[0].0,
                "source continuation requires a following V30 line",
            ));
        }
    }
    let logical = collect_v30_lines(&lines);
    validate_ids(&logical, limits)?;
    let markers = logical
        .iter()
        .map(|line| {
            let mut words = line.payload.split_whitespace();
            (words.next(), words.next())
        })
        .collect::<Vec<_>>();
    if markers.first() != Some(&(Some("BEGIN"), Some("CTAB")))
        || markers.last() != Some(&(Some("END"), Some("CTAB")))
        || markers
            .iter()
            .filter(|&&m| m == (Some("BEGIN"), Some("CTAB")))
            .count()
            != 1
        || markers
            .iter()
            .filter(|&&m| m == (Some("END"), Some("CTAB")))
            .count()
            != 1
    {
        return Err(v3k_err(
            5,
            "source retention requires exactly one complete CTAB",
        ));
    }
    let (mol, _) = super::parse_mol_v3000(input)?;
    bound("V3000 atoms", mol.atom_count(), limits.max_atoms)?;
    Ok(input.to_owned())
}
