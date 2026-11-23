use super::identity_helpers::PyRdkitMmffConformer;
use super::*;

#[pymethods]
impl Mol {
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
