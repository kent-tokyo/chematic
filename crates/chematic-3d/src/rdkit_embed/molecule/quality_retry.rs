//! Bounded MMFF94 retry for the lowest-energy unfinished quality candidate.

use chematic_core::{Coords3D, Molecule, Point3};

use crate::stereo_constraints::verify_stereo;

use super::{RdkitMmffConformer, optimize_mmff_conformer};

/// Accept a refined geometry only when every declared stereochemical
/// constraint remains evaluable and satisfied. The adaptive retry APIs share
/// this gate with the quality profile.
pub(super) fn stereo_safe_retry(mh: &Molecule, coords: &[[f64; 3]]) -> bool {
    if coords.len() != mh.atom_count() || !coords.iter().flatten().all(|value| value.is_finite()) {
        return false;
    }
    let geometry = Coords3D {
        points: coords
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect(),
    };
    let stereo = verify_stereo(mh, &geometry);
    stereo.n_violations() == 0 && stereo.n_unevaluable() == 0
}

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
        && stereo_safe_retry(mh, &retry.coords)
    {
        candidates[best_index] = retry;
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
    fn adaptive_retry_rejects_mirrored_and_planar_stereo() {
        let mol = chematic_smiles::parse("F[C@](Cl)(Br)I").unwrap();
        let coords = [
            [1.0, 1.0, 1.0],
            [0.0, 0.0, 0.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ];
        let mirrored = coords.map(|point| [-point[0], point[1], point[2]]);
        assert_ne!(
            stereo_safe_retry(&mol, &coords),
            stereo_safe_retry(&mol, &mirrored)
        );
        assert!(!stereo_safe_retry(
            &mol,
            &coords.map(|point| [point[0], point[1], 0.0])
        ));
        assert!(!stereo_safe_retry(&mol, &coords[..4]));
    }

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
