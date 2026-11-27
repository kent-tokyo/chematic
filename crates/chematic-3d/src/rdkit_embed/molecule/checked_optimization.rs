//! Stereo-safe checked MMFF optimization and rejection diagnostics.

use super::*;

#[derive(Default)]
pub(super) struct RejectionCounts {
    pub stereo_violated: usize,
    pub stereo_unevaluable: usize,
    pub invalid_coordinates: usize,
}

pub(super) fn optimize_checked_ensemble_with_field(
    mh: &Molecule,
    embedded: Vec<Vec<[f64; 3]>>,
    max_iters: u32,
    field: &chematic_ff::rdkit_mmff::RdkitMmffField,
) -> Result<Vec<RdkitMmffConformer>, RdkitEmbedMmffError> {
    Ok(optimize_checked_ensemble_with_rejections(mh, embedded, max_iters, field)?.0)
}

pub(super) fn optimize_checked_ensemble_with_rejections(
    mh: &Molecule,
    embedded: Vec<Vec<[f64; 3]>>,
    max_iters: u32,
    field: &chematic_ff::rdkit_mmff::RdkitMmffField,
) -> Result<(Vec<RdkitMmffConformer>, RejectionCounts), RdkitEmbedMmffError> {
    let optimized = optimize_mmff_conformers(field, embedded.clone(), max_iters);
    let mut conformers = Vec::with_capacity(optimized.len());
    let mut rejected = RejectionCounts::default();
    for (mut candidate, original) in optimized.into_iter().zip(embedded) {
        if !stereo_safe_retry(mh, &candidate.coords) {
            if !stereo_safe_retry(mh, &original) {
                count_rejected_original(mh, &original, &mut rejected);
                continue;
            }
            let energy = field.energy(original.as_flattened());
            if !energy.is_finite() {
                continue;
            }
            candidate = RdkitMmffConformer {
                coords: original,
                status: 1,
                energy,
            };
        }
        conformers.push(candidate);
    }
    for conformer in &mut conformers {
        if conformer.status != 0 {
            continue;
        }
        if !conformer.energy.is_finite()
            || conformer
                .coords
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
        {
            conformer.status = 1;
            continue;
        }
        let field =
            chematic_ff::rdkit_mmff::RdkitMmffField::new(mh, &conformer.coords, 100.0, true)
                .map_err(|error| RdkitEmbedMmffError::ForceField(error.to_string()))?;
        let coordinates = conformer.coords.as_flattened();
        let mut gradient = vec![0.0; coordinates.len()];
        field.gradient(coordinates, &mut gradient);
        if gradient
            .iter()
            .any(|value| !value.is_finite() || value.abs() > 0.01)
        {
            conformer.status = 1;
        }
    }
    Ok((conformers, rejected))
}

fn count_rejected_original(mh: &Molecule, coords: &[[f64; 3]], counts: &mut RejectionCounts) {
    use crate::coords::{Coords3D, Point3};
    use crate::stereo_constraints::verify_stereo;
    if coords.len() != mh.atom_count() || coords.iter().flatten().any(|v| !v.is_finite()) {
        counts.invalid_coordinates += 1;
        return;
    }
    let geometry = Coords3D {
        points: coords
            .iter()
            .map(|&[x, y, z]| Point3::new(x, y, z))
            .collect(),
    };
    let stereo = verify_stereo(mh, &geometry);
    if stereo.n_violations() != 0 {
        counts.stereo_violated += 1;
    } else if stereo.n_unevaluable() != 0 {
        counts.stereo_unevaluable += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checked_mmff_reports_rejected_inverted_stereo_without_changing_output() {
        let mh = chematic_chem::add_hydrogens(&chematic_smiles::parse("C[C@H](O)F").unwrap());
        let options = RdkitEmbedOptions::default();
        let coords = etkdg_embed_multiple_confs_corrected(&mh, 1, &options)
            .unwrap()
            .remove(0);
        let mirror: Vec<_> = coords.iter().map(|&[x, y, z]| [-x, y, z]).collect();
        assert!(stereo_safe_retry(&mh, &coords));
        assert!(!stereo_safe_retry(&mh, &mirror));
        let field =
            chematic_ff::rdkit_mmff::RdkitMmffField::new(&mh, &coords, 100.0, true).unwrap();
        let input = vec![coords, mirror];
        let expected = optimize_checked_ensemble_with_field(&mh, input.clone(), 0, &field).unwrap();
        let report = optimize_embedded_checked_diagnostics(&mh, input, 0, 0).unwrap();
        assert_eq!(report.embedded, 2);
        assert_eq!(report.stereo_violated, Some(1));
        assert_eq!(report.stereo_unevaluable, Some(0));
        assert_eq!(report.invalid_coordinates, Some(0));
        assert_eq!(report.conformers.len(), 1);
        assert_eq!(report.conformers[0].coords, expected[0].coords);
        assert_eq!(report.conformers[0].status, expected[0].status);
        assert_eq!(report.conformers[0].energy, Some(expected[0].energy));
    }

    #[test]
    fn rejection_classifier_distinguishes_invalid_coordinates() {
        let mh = chematic_chem::add_hydrogens(&chematic_smiles::parse("CCO").unwrap());
        let mut counts = RejectionCounts::default();
        count_rejected_original(&mh, &[], &mut counts);
        count_rejected_original(
            &mh,
            &vec![[f64::NAN, 0.0, 0.0]; mh.atom_count()],
            &mut counts,
        );
        assert_eq!(counts.invalid_coordinates, 2);
        assert_eq!(counts.stereo_violated, 0);
        assert_eq!(counts.stereo_unevaluable, 0);
    }

    #[test]
    fn rejection_classifier_distinguishes_degenerate_stereo() {
        let mh = chematic_chem::add_hydrogens(&chematic_smiles::parse("C[C@H](O)F").unwrap());
        let mut counts = RejectionCounts::default();
        count_rejected_original(&mh, &vec![[0.0; 3]; mh.atom_count()], &mut counts);
        assert_eq!(counts.stereo_unevaluable, 1);
        assert_eq!(counts.stereo_violated, 0);
        assert_eq!(counts.invalid_coordinates, 0);
    }
}
