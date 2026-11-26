//! Mixed checked conformers with an explicit MMFF-unavailable fallback.

use super::*;

/// An embedded conformer, including coordinates when MMFF setup failed.
#[derive(Clone, Debug)]
pub struct CheckedConformer {
    pub coords: Vec<[f64; 3]>,
    /// Zero means converged, one means unfinished, and minus one means that
    /// MMFF could not be prepared. The latter has no energy.
    pub status: i32,
    pub energy: Option<f64>,
    pub force_field_error: Option<String>,
}

/// Counts for a mixed checked ensemble. Detailed rejection counts are
/// available when MMFF setup fails and the embedded geometries are screened.
#[derive(Clone, Debug)]
pub struct CheckedEnsembleDiagnostics {
    pub conformers: Vec<CheckedConformer>,
    pub embedded: usize,
    pub force_field_error: Option<String>,
    pub stereo_violated: Option<usize>,
    pub stereo_unevaluable: Option<usize>,
    pub invalid_coordinates: Option<usize>,
}

pub(super) fn embed_mixed(
    mh: &Molecule,
    num_confs: usize,
    o: &RdkitEmbedOptions,
) -> Result<Vec<Vec<[f64; 3]>>, RdkitEmbedMmffError> {
    validate_options(mh, o)?;
    if num_confs == 0 {
        return Ok(Vec::new());
    }
    let prepared = match prepare_fragments(mh, o) {
        Ok(prepared) => prepared,
        Err(RdkitEmbedError::Failed) => return Ok(Vec::new()),
        Err(error) => return Err(error.into()),
    };
    let first_count = num_confs / 2;
    let mut embedded = embed_prepared_multiple_confs(&prepared, first_count, 0, o, false)?;
    embedded.extend(embed_prepared_multiple_confs(
        &prepared,
        num_confs - first_count,
        first_count,
        o,
        true,
    )?);
    Ok(embedded)
}

/// Generate a checked mixed ensemble while retaining stereo-checked embedded
/// coordinates if MMFF parameterization is unavailable.
pub fn etkdg_embed_multiple_confs_mixed_checked_with_fallback(
    mh: &Molecule,
    num_confs: usize,
    o: &RdkitEmbedOptions,
    max_iters: u32,
    retry_iters: u32,
) -> Result<Vec<CheckedConformer>, RdkitEmbedMmffError> {
    Ok(etkdg_embed_multiple_confs_mixed_checked_diagnostics(
        mh,
        num_confs,
        o,
        max_iters,
        retry_iters,
    )?
    .conformers)
}

/// Also report why embedded coordinates were rejected when MMFF setup fails.
pub fn etkdg_embed_multiple_confs_mixed_checked_diagnostics(
    mh: &Molecule,
    num_confs: usize,
    o: &RdkitEmbedOptions,
    max_iters: u32,
    retry_iters: u32,
) -> Result<CheckedEnsembleDiagnostics, RdkitEmbedMmffError> {
    let embedded = embed_mixed(mh, num_confs, o)?;
    optimize_embedded_checked_diagnostics_inner(mh, embedded, max_iters, retry_iters)
}

/// Optimize previously embedded coordinates with the checked MMFF fallback.
/// Rejects invalid dimensions or nonfinite coordinates before field creation.
pub fn optimize_embedded_checked_diagnostics(
    mh: &Molecule,
    embedded: Vec<Vec<[f64; 3]>>,
    max_iters: u32,
    retry_iters: u32,
) -> Result<CheckedEnsembleDiagnostics, RdkitEmbedMmffError> {
    if (!embedded.is_empty() && mh.atom_count() == 0)
        || embedded.iter().any(|coords| {
            coords.len() != mh.atom_count() || coords.iter().flatten().any(|v| !v.is_finite())
        })
    {
        return Err(RdkitEmbedError::Unsupported(
            "invalid embedded coordinate shape or nonfinite coordinate".into(),
        )
        .into());
    }
    optimize_embedded_checked_diagnostics_inner(mh, embedded, max_iters, retry_iters)
}

fn optimize_embedded_checked_diagnostics_inner(
    mh: &Molecule,
    embedded: Vec<Vec<[f64; 3]>>,
    max_iters: u32,
    retry_iters: u32,
) -> Result<CheckedEnsembleDiagnostics, RdkitEmbedMmffError> {
    use crate::coords::{Coords3D, Point3};
    use crate::stereo_constraints::verify_stereo;
    let count = embedded.len();
    let Some(first) = embedded.first() else {
        return Ok(CheckedEnsembleDiagnostics {
            conformers: Vec::new(),
            embedded: 0,
            force_field_error: None,
            stereo_violated: None,
            stereo_unevaluable: None,
            invalid_coordinates: None,
        });
    };
    let field = match chematic_ff::rdkit_mmff::RdkitMmffField::new(mh, first, 100.0, true) {
        Ok(field) => field,
        Err(error) => {
            let reason = error.to_string();
            let (mut violated, mut unevaluable, mut invalid) = (0, 0, 0);
            let mut conformers = Vec::new();
            for coords in embedded {
                if coords.len() != mh.atom_count()
                    || coords.iter().flatten().any(|v| !v.is_finite())
                {
                    invalid += 1;
                    continue;
                }
                let geometry = Coords3D {
                    points: coords
                        .iter()
                        .map(|&[x, y, z]| Point3::new(x, y, z))
                        .collect(),
                };
                let stereo = verify_stereo(mh, &geometry);
                if stereo.n_violations() != 0 {
                    violated += 1;
                } else if stereo.n_unevaluable() != 0 {
                    unevaluable += 1;
                } else {
                    conformers.push(CheckedConformer {
                        coords,
                        status: -1,
                        energy: None,
                        force_field_error: Some(reason.clone()),
                    });
                }
            }
            return Ok(CheckedEnsembleDiagnostics {
                conformers,
                embedded: count,
                force_field_error: Some(reason),
                stereo_violated: Some(violated),
                stereo_unevaluable: Some(unevaluable),
                invalid_coordinates: Some(invalid),
            });
        }
    };
    let mut conformers = optimize_checked_ensemble_with_field(mh, embedded, max_iters, &field)?;
    for conformer in &mut conformers {
        repair_linear_nitrogen(mh, conformer, retry_iters);
    }
    let conformers = conformers
        .into_iter()
        .map(|conformer| CheckedConformer {
            coords: conformer.coords,
            status: conformer.status,
            energy: Some(conformer.energy),
            force_field_error: None,
        })
        .collect();
    Ok(CheckedEnsembleDiagnostics {
        conformers,
        embedded: count,
        force_field_error: None,
        stereo_violated: None,
        stereo_unevaluable: None,
        invalid_coordinates: None,
    })
}
