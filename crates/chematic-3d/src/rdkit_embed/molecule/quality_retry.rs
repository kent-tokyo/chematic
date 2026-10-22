//! Bounded MMFF94 retry for the lowest-energy unfinished quality candidate.

use chematic_core::{Coords3D, Molecule, Point3};

use crate::stereo_constraints::verify_stereo;

use super::{RdkitMmffConformer, optimize_mmff_conformer};

pub(super) fn refine_unfinished_quality_minimum(
    mh: &Molecule,
    mut candidates: Vec<RdkitMmffConformer>,
    max_iters: u32,
) -> Vec<RdkitMmffConformer> {
    let retry_iters = max_iters.saturating_mul(2);
    if retry_iters == 0 {
        return candidates;
    }
    let Some(best_index) = lowest_unfinished_index(&candidates) else {
        return candidates;
    };
    let best = &candidates[best_index];
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

fn lowest_unfinished_index(candidates: &[RdkitMmffConformer]) -> Option<usize> {
    candidates
        .iter()
        .enumerate()
        .filter(|(_, candidate)| candidate.status != 0 && candidate.energy.is_finite())
        .min_by(|(_, left), (_, right)| left.energy.total_cmp(&right.energy))
        .map(|(index, _)| index)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn retry_selection_skips_converged_minimum() {
        let candidate = |energy, status| RdkitMmffConformer {
            coords: vec![],
            status,
            energy,
        };
        let candidates = [
            candidate(-10.0, 0),
            candidate(-7.0, 1),
            candidate(-9.0, 1),
            candidate(f64::NAN, 1),
        ];
        assert_eq!(lowest_unfinished_index(&candidates), Some(2));
        assert_eq!(lowest_unfinished_index(&candidates[..1]), None);
    }
}
