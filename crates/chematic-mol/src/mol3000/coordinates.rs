use crate::{MolMetadata, MolParseError};
use chematic_core::Molecule;

/// Parse V3000 with source-order x/y/z values, independently of the header's
/// dimension code and geometry-perception tolerance. Values are finite f64s in
/// Angstrom; signed zero is retained. The existing graph/stereo interpretation
/// is unchanged. This is not source-text or opaque-metadata preservation.
///
/// Input is bounded to 16 MiB; native atom/bond/count limits also apply.
#[allow(clippy::type_complexity)]
pub fn parse_mol_v3000_with_coords_3d(
    input: &str,
) -> Result<(Molecule, MolMetadata, Vec<(f64, f64, f64)>), MolParseError> {
    const MAX_INPUT_BYTES: usize = 16 * 1024 * 1024;
    if input.len() > MAX_INPUT_BYTES {
        return Err(MolParseError::ResourceLimit {
            resource: "V3000 input bytes",
            actual: input.len(),
            limit: MAX_INPUT_BYTES,
        });
    }
    let (report, z) = super::read_v3000_with_source_z(input)?;
    let coords = report
        .coords
        .into_iter()
        .zip(z)
        .map(|((x, y), z)| (x, y, z))
        .collect();
    Ok((report.mol, report.metadata, coords))
}
