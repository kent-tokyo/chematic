//! Reusable immutable setup for corrected-gradient single attempts.

use super::*;

/// Reuse distance bounds and torsion knowledge across independent corrected
/// embedding attempts. Each call accepts its explicit nonnegative seed.
pub struct PreparedCorrectedEmbedding {
    options: RdkitEmbedOptions,
    fragments: Option<Vec<PreparedFragment>>,
}

impl PreparedCorrectedEmbedding {
    /// Prepare immutable molecule data once. A smoothing failure is retained
    /// as an empty session, matching one-at-a-time ensemble behavior.
    pub fn new(molecule: &Molecule, options: &RdkitEmbedOptions) -> Result<Self, RdkitEmbedError> {
        validate_options(molecule, options)?;
        let fragments = match prepare_fragments(molecule, options) {
            Ok(fragments) => Some(fragments),
            Err(RdkitEmbedError::Failed) => None,
            Err(error) => return Err(error),
        };
        Ok(Self {
            options: options.clone(),
            fragments,
        })
    }

    /// Return one embedding or `None` when this attempt fails to embed.
    pub fn embed(&self, seed: i32) -> Result<Option<Vec<[f64; 3]>>, RdkitEmbedError> {
        if seed < 0 {
            return Err(RdkitEmbedError::Unsupported(
                "randomSeed must be >= 0 (RDKit's -1 is not reproducible)".into(),
            ));
        }
        let Some(fragments) = &self.fragments else {
            return Ok(None);
        };
        match embed_prepared_fragments_with(fragments, &self.options, seed, false, true) {
            Ok(coords) => Ok(Some(coords)),
            Err(RdkitEmbedError::Failed) => Ok(None),
            Err(error) => Err(error),
        }
    }
}
