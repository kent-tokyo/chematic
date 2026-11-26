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
