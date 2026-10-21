//! Bounded MMFF94 retry for a quality ensemble with no converged conformer.

use chematic_core::{Coords3D, Molecule, Point3};

use crate::stereo_constraints::verify_stereo;

use super::{RdkitMmffConformer, optimize_mmff_conformer};

pub(super) fn rescue_quality_without_convergence(
    mh: &Molecule,
    mut candidates: Vec<RdkitMmffConformer>,
    max_iters: u32,
) -> Vec<RdkitMmffConformer> {
    let retry_iters = max_iters.saturating_mul(2);
    if retry_iters == 0 || candidates.iter().any(|candidate| candidate.status == 0) {
        return candidates;
    }
    let Some((best_index, best)) = candidates
        .iter()
        .enumerate()
        .filter(|(_, candidate)| candidate.energy.is_finite())
        .min_by(|(_, left), (_, right)| left.energy.total_cmp(&right.energy))
    else {
        return candidates;
    };
    let Ok(field) = chematic_ff::rdkit_mmff::RdkitMmffField::new(mh, &best.coords, 100.0, true)
    else {
        return candidates;
    };
    let retry = optimize_mmff_conformer(&field, best.coords.clone(), retry_iters);
    if retry.status == 0
        && retry.energy.is_finite()
        && retry.energy <= best.energy
        && retry.coords.iter().flatten().all(|value| value.is_finite())
    {
        let coords = Coords3D {
            points: retry
                .coords
                .iter()
                .map(|&[x, y, z]| Point3::new(x, y, z))
                .collect(),
        };
        let stereo = verify_stereo(mh, &coords);
        if stereo.n_violations() == 0 && stereo.n_unevaluable() == 0 {
            candidates[best_index] = retry;
        }
    }
    candidates
}
