use super::identity_helpers::PyRdkitMmffConformer;
use super::*;

type FallbackConformer = (i32, Option<f64>, Vec<[f64; 3]>, Option<String>);
type FallbackDiagnostics = (
    Vec<FallbackConformer>,
    usize,
    Option<String>,
    Option<usize>,
    Option<usize>,
    Option<usize>,
);

#[pymethods]
impl Mol {
    /// Prepare one corrected-gradient embedding for the timed worker.
    #[pyo3(signature = (random_seed = 42, max_iterations = 0))]
    fn _timed_embed_corrected(
        &self,
        py: Python<'_>,
        random_seed: i32,
        max_iterations: u32,
    ) -> PyResult<Option<Vec<[f64; 3]>>> {
        use chematic_3d::rdkit_embed::{RdkitEmbedOptions, etkdg_embed_multiple_confs_corrected};
        let options = RdkitEmbedOptions {
            random_seed,
            max_iterations,
            ..Default::default()
        };
        let mol = self.inner.clone();
        py.detach(move || etkdg_embed_multiple_confs_corrected(&mol, 1, &options))
            .map(|mut values| values.pop())
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Check declared stereochemistry before publishing a timed preview.
    fn _timed_stereo_valid(&self, coords: Vec<[f64; 3]>) -> PyResult<bool> {
        use chematic_3d::coords::{Coords3D, Point3};
        use chematic_3d::stereo_constraints::verify_stereo;
        if coords.len() != self.inner.atom_count()
            || coords.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(PyValueError::new_err("invalid embedded coordinates"));
        }
        let geometry = Coords3D {
            points: coords
                .into_iter()
                .map(|[x, y, z]| Point3::new(x, y, z))
                .collect(),
        };
        let result = verify_stereo(&self.inner, &geometry);
        Ok(result.n_violations() == 0 && result.n_unevaluable() == 0)
    }

    /// Optimize a previously prepared embedding without re-embedding it.
    #[pyo3(signature = (coords, *, max_mmff_iters = 200, repair_mmff_iters = 0))]
    fn _timed_optimize_embedded(
        &self,
        py: Python<'_>,
        coords: Vec<[f64; 3]>,
        max_mmff_iters: u32,
        repair_mmff_iters: u32,
    ) -> PyResult<FallbackDiagnostics> {
        use chematic_3d::rdkit_embed::optimize_embedded_checked_diagnostics;
        let mol = self.inner.clone();
        py.detach(move || {
            optimize_embedded_checked_diagnostics(
                &mol,
                vec![coords],
                max_mmff_iters,
                repair_mmff_iters,
            )
        })
        .map(|result| {
            let rows = result
                .conformers
                .into_iter()
                .map(|c| (c.status, c.energy, c.coords, c.force_field_error))
                .collect();
            (
                rows,
                result.embedded,
                result.force_field_error,
                result.stereo_violated,
                result.stereo_unevaluable,
                result.invalid_coordinates,
            )
        })
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Return conformers, embedded count, MMFF setup error, and counts of
    /// stereo violations, unevaluable stereo and invalid coordinates.
    /// Rejection counts are available when MMFF setup fails.
    #[pyo3(signature = (num_confs, *, random_seed = 42, max_iterations = 0,
                        max_mmff_iters = 200, repair_mmff_iters = 0))]
    fn etkdg_embed_multiple_confs_mixed_checked_diagnostics(
        &self,
        py: Python<'_>,
        num_confs: usize,
        random_seed: i32,
        max_iterations: u32,
        max_mmff_iters: u32,
        repair_mmff_iters: u32,
    ) -> PyResult<FallbackDiagnostics> {
        use chematic_3d::rdkit_embed::{
            RdkitEmbedOptions, etkdg_embed_multiple_confs_mixed_checked_diagnostics,
        };
        let options = RdkitEmbedOptions {
            random_seed,
            max_iterations,
            ..Default::default()
        };
        let mol = self.inner.clone();
        py.detach(move || {
            etkdg_embed_multiple_confs_mixed_checked_diagnostics(
                &mol,
                num_confs,
                &options,
                max_mmff_iters,
                repair_mmff_iters,
            )
        })
        .map(|result| {
            let rows = result
                .conformers
                .into_iter()
                .map(|c| (c.status, c.energy, c.coords, c.force_field_error))
                .collect();
            (
                rows,
                result.embedded,
                result.force_field_error,
                result.stereo_violated,
                result.stereo_unevaluable,
                result.invalid_coordinates,
            )
        })
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }
    /// Return (status, optional energy, coordinates, optional MMFF error).
    /// Status -1 retains the embedded coordinates after MMFF setup fails.
    #[pyo3(signature = (num_confs, *, random_seed = 42, max_iterations = 0,
                        max_mmff_iters = 200, repair_mmff_iters = 0))]
    fn etkdg_embed_multiple_confs_mixed_checked_with_fallback(
        &self,
        py: Python<'_>,
        num_confs: usize,
        random_seed: i32,
        max_iterations: u32,
        max_mmff_iters: u32,
        repair_mmff_iters: u32,
    ) -> PyResult<Vec<FallbackConformer>> {
        use chematic_3d::rdkit_embed::{
            RdkitEmbedOptions, etkdg_embed_multiple_confs_mixed_checked_with_fallback,
        };
        let options = RdkitEmbedOptions {
            random_seed,
            max_iterations,
            ..Default::default()
        };
        let mol = self.inner.clone();
        py.detach(move || {
            etkdg_embed_multiple_confs_mixed_checked_with_fallback(
                &mol,
                num_confs,
                &options,
                max_mmff_iters,
                repair_mmff_iters,
            )
        })
        .map(|conformers| {
            conformers
                .into_iter()
                .map(|c| (c.status, c.energy, c.coords, c.force_field_error))
                .collect()
        })
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }
    /// Generate and optimize one ensemble using the corrected ETKDG torsion
    /// gradient. This uses the requested number of embedding attempts and
    /// MMFF94 minimizations without the extra batches in the quality profile.
    /// Status 0 additionally requires a finite final gradient whose largest
    /// absolute component is at most 0.01 kcal/mol/angstrom.
    /// If minimization violates declared stereo, the embedded coordinates are
    /// returned with status 1 when their stereo is evaluable and correct.
    /// Candidates whose embedded and optimized stereo both fail are omitted.
    /// `repair_mmff_iters` enables a bounded retry for near-linear nitrogen.
    #[pyo3(signature = (num_confs, *, random_seed = 42, max_iterations = 0,
                        max_mmff_iters = 200, repair_mmff_iters = 0))]
    fn etkdg_embed_multiple_confs_mmff_corrected(
        &self,
        py: Python<'_>,
        num_confs: usize,
        random_seed: i32,
        max_iterations: u32,
        max_mmff_iters: u32,
        repair_mmff_iters: u32,
    ) -> PyResult<Vec<PyRdkitMmffConformer>> {
        use chematic_3d::rdkit_embed::{
            RdkitEmbedOptions, etkdg_embed_multiple_confs_mmff_corrected_repair,
        };
        let options = RdkitEmbedOptions {
            random_seed,
            max_iterations,
            ..Default::default()
        };
        let mol = self.inner.clone();
        py.detach(move || {
            etkdg_embed_multiple_confs_mmff_corrected_repair(
                &mol,
                num_confs,
                &options,
                max_mmff_iters,
                repair_mmff_iters,
            )
        })
        .map(|conformers| {
            conformers
                .into_iter()
                .map(|conformer| (conformer.status, conformer.energy, conformer.coords))
                .collect()
        })
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Generate an RDKit-compatible-gradient ensemble with the same stereo,
    /// final-gradient and optional nitrogen checks as the corrected profile.
    #[pyo3(signature = (num_confs, *, random_seed = 42, max_iterations = 0,
                        max_mmff_iters = 200, repair_mmff_iters = 0))]
    fn etkdg_embed_multiple_confs_mmff_compatible_checked(
        &self,
        py: Python<'_>,
        num_confs: usize,
        random_seed: i32,
        max_iterations: u32,
        max_mmff_iters: u32,
        repair_mmff_iters: u32,
    ) -> PyResult<Vec<PyRdkitMmffConformer>> {
        use chematic_3d::rdkit_embed::{
            RdkitEmbedOptions, etkdg_embed_multiple_confs_mmff_compatible_checked,
        };
        let options = RdkitEmbedOptions {
            random_seed,
            max_iterations,
            ..Default::default()
        };
        let mol = self.inner.clone();
        py.detach(move || {
            etkdg_embed_multiple_confs_mmff_compatible_checked(
                &mol,
                num_confs,
                &options,
                max_mmff_iters,
                repair_mmff_iters,
            )
        })
        .map(|conformers| {
            conformers
                .into_iter()
                .map(|conformer| (conformer.status, conformer.energy, conformer.coords))
                .collect()
        })
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// Split one checked ensemble between compatibility and corrected gradients.
    #[pyo3(signature = (num_confs, *, random_seed = 42, max_iterations = 0,
                        max_mmff_iters = 200, repair_mmff_iters = 0))]
    fn etkdg_embed_multiple_confs_mmff_mixed_checked(
        &self,
        py: Python<'_>,
        num_confs: usize,
        random_seed: i32,
        max_iterations: u32,
        max_mmff_iters: u32,
        repair_mmff_iters: u32,
    ) -> PyResult<Vec<PyRdkitMmffConformer>> {
        use chematic_3d::rdkit_embed::{
            RdkitEmbedOptions, etkdg_embed_multiple_confs_mmff_mixed_checked,
        };
        let options = RdkitEmbedOptions {
            random_seed,
            max_iterations,
            ..Default::default()
        };
        let mol = self.inner.clone();
        py.detach(move || {
            etkdg_embed_multiple_confs_mmff_mixed_checked(
                &mol,
                num_confs,
                &options,
                max_mmff_iters,
                repair_mmff_iters,
            )
        })
        .map(|conformers| {
            conformers
                .into_iter()
                .map(|conformer| (conformer.status, conformer.energy, conformer.coords))
                .collect()
        })
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }
}
