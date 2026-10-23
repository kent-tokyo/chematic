use crate::{MolMetadata, MolParseError};
use chematic_core::Molecule;

fn coordinate(row: &str, atom_line: usize, axis: usize) -> Result<f64, MolParseError> {
    let invalid = || MolParseError::InvalidAtomLine {
        line: atom_line,
        detail: format!(
            "{} coordinate must be a present finite number",
            ["x", "y", "z"][axis]
        ),
    };
    let value = row
        .get(axis * 10..axis * 10 + 10)
        .and_then(|field| field.trim().parse::<f64>().ok())
        .ok_or_else(invalid)?;
    if !value.is_finite() {
        return Err(invalid());
    }
    Ok(value)
}

/// Parse V2000 with all finite source-order Cartesian coordinate triplets in Å.
///
/// Unlike the optional perception conformer, this always retains every source
/// z value, including tiny z, flat geometry and signed zero. All three 10-byte
/// fields must contain finite numbers; missing/malformed values are not filled.
/// The graph/stereo interpretation uses the existing reader unchanged. This
/// does not preserve source text, SD properties or opaque CTAB extensions.
/// Input is bounded to 16 MiB and native atom/bond limits also apply.
#[allow(clippy::type_complexity)]
pub fn parse_mol_with_coords_3d(
    input: &str,
) -> Result<(Molecule, MolMetadata, Vec<(f64, f64, f64)>), MolParseError> {
    const MAX_INPUT_BYTES: usize = 16 * 1024 * 1024;
    if input.len() > MAX_INPUT_BYTES {
        return Err(MolParseError::ResourceLimit {
            resource: "V2000 input bytes",
            actual: input.len(),
            limit: MAX_INPUT_BYTES,
        });
    }
    let mut lines = input.lines().enumerate();
    lines.nth(2).ok_or(MolParseError::UnexpectedEnd)?;
    let (atoms, _) = super::reader::parse_counts(&mut lines)?;
    let mut coords = Vec::with_capacity(atoms.min(input.len() / 34));
    for _ in 0..atoms {
        let (index, row) = lines.next().ok_or(MolParseError::UnexpectedEnd)?;
        coords.push((
            coordinate(row, index + 1, 0)?,
            coordinate(row, index + 1, 1)?,
            coordinate(row, index + 1, 2)?,
        ));
    }
    let report = super::read_mol_with_diagnostics(input)?;
    Ok((report.mol, report.metadata, coords))
}
