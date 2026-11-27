use super::*;

#[pymethods]
impl Mol {
    /// Embed and additionally screen nonplanarity and other declared stereo.
    /// The dimensionless planarity limit must be explicitly supplied.
    #[pyo3(signature = (*, max_unit_triple_product, random_seed = 42, max_iterations = 0,
                        trans_min_degrees = 150.0, cis_tolerance_degrees = 30.0))]
    fn rdkit_embed_square_planar_checked<'py>(
        &self,
        py: Python<'py>,
        max_unit_triple_product: f64,
        random_seed: i32,
        max_iterations: u32,
        trans_min_degrees: f64,
        cis_tolerance_degrees: f64,
    ) -> PyResult<Bound<'py, PyDict>> {
        use chematic_3d::rdkit_embed::RdkitEmbedOptions;
        use chematic_3d::square_planar_checked::{
            SquarePlanarCheckedConfig, SquarePlanarCheckedError, rdkit_embed_square_planar_checked,
        };
        use chematic_3d::square_planar_embedding::{
            SquarePlanarEmbedConfig, SquarePlanarEmbedError, SquarePlanarRetryError,
        };
        let mol = self.inner.clone();
        let options = RdkitEmbedOptions {
            random_seed,
            max_iterations,
            ..Default::default()
        };
        let config = SquarePlanarCheckedConfig {
            angles: SquarePlanarEmbedConfig {
                trans_min_degrees,
                cis_tolerance_degrees,
            },
            max_unit_triple_product,
        };
        let result = py
            .detach(move || rdkit_embed_square_planar_checked(&mol, &options, config))
            .map_err(|error| {
                let message = error.to_string();
                match error {
                    SquarePlanarCheckedError::InvalidPlanarityLimit
                    | SquarePlanarCheckedError::Embedding(SquarePlanarRetryError::Initial(
                        SquarePlanarEmbedError::InvalidThreshold
                        | SquarePlanarEmbedError::NoDeclaration
                        | SquarePlanarEmbedError::UnsupportedCenter { .. },
                    )) => PyValueError::new_err(message),
                    _ => pyo3::exceptions::PyRuntimeError::new_err(message),
                }
            })?;
        let item = PyDict::new(py);
        item.set_item("coordinates", result.embedding.coords)?;
        item.set_item("attempts", result.embedding.attempts)?;
        item.set_item(
            "used_tightened_bounds",
            result.embedding.used_tightened_bounds,
        )?;
        item.set_item("initial_rejection", result.embedding.initial_rejection)?;
        item.set_item("max_unit_triple_product", result.max_unit_triple_product)?;
        item.set_item("other_stereo_checked", result.other_stereo_checked)?;
        Ok(item)
    }

    /// Try square-planar embedding once, with one tighter-bounds retry after
    /// a final angle rejection. The original angle policy applies to both.
    #[pyo3(signature = (*, random_seed = 42, max_iterations = 0,
                        trans_min_degrees = 150.0, cis_tolerance_degrees = 30.0))]
    fn rdkit_embed_square_planar_with_retry<'py>(
        &self,
        py: Python<'py>,
        random_seed: i32,
        max_iterations: u32,
        trans_min_degrees: f64,
        cis_tolerance_degrees: f64,
    ) -> PyResult<Bound<'py, PyDict>> {
        use chematic_3d::rdkit_embed::RdkitEmbedOptions;
        use chematic_3d::square_planar_embedding::{
            SquarePlanarEmbedConfig, SquarePlanarEmbedError, SquarePlanarRetryError,
            rdkit_embed_square_planar_with_retry,
        };
        let mol = self.inner.clone();
        let options = RdkitEmbedOptions {
            random_seed,
            max_iterations,
            ..Default::default()
        };
        let config = SquarePlanarEmbedConfig {
            trans_min_degrees,
            cis_tolerance_degrees,
        };
        let outcome = py
            .detach(move || rdkit_embed_square_planar_with_retry(&mol, &options, config))
            .map_err(|error| {
                let message = error.to_string();
                match error {
                    SquarePlanarRetryError::Initial(
                        SquarePlanarEmbedError::InvalidThreshold
                        | SquarePlanarEmbedError::NoDeclaration
                        | SquarePlanarEmbedError::UnsupportedCenter { .. },
                    ) => PyValueError::new_err(message),
                    _ => pyo3::exceptions::PyRuntimeError::new_err(message),
                }
            })?;
        let item = PyDict::new(py);
        item.set_item("coordinates", outcome.coords)?;
        item.set_item("attempts", outcome.attempts)?;
        item.set_item("used_tightened_bounds", outcome.used_tightened_bounds)?;
        item.set_item("initial_rejection", outcome.initial_rejection)?;
        Ok(item)
    }

    /// Generate one experimental square-planar-constrained ETKDG conformer.
    /// Raises if bounds, embedding, or final angle validation fail.
    #[pyo3(signature = (*, random_seed = 42, max_iterations = 0,
                        trans_min_degrees = 150.0, cis_tolerance_degrees = 30.0))]
    fn rdkit_embed_square_planar(
        &self,
        py: Python<'_>,
        random_seed: i32,
        max_iterations: u32,
        trans_min_degrees: f64,
        cis_tolerance_degrees: f64,
    ) -> PyResult<Vec<[f64; 3]>> {
        use chematic_3d::rdkit_embed::RdkitEmbedOptions;
        use chematic_3d::square_planar_embedding::{
            SquarePlanarEmbedConfig, SquarePlanarEmbedError, rdkit_embed_square_planar,
        };
        let mol = self.inner.clone();
        let options = RdkitEmbedOptions {
            random_seed,
            max_iterations,
            ..Default::default()
        };
        let config = SquarePlanarEmbedConfig {
            trans_min_degrees,
            cis_tolerance_degrees,
        };
        py.detach(move || rdkit_embed_square_planar(&mol, &options, config))
            .map_err(|error| {
                let message = error.to_string();
                match error {
                    SquarePlanarEmbedError::InvalidThreshold
                    | SquarePlanarEmbedError::NoDeclaration
                    | SquarePlanarEmbedError::UnsupportedCenter { .. } => {
                        PyValueError::new_err(message)
                    }
                    _ => pyo3::exceptions::PyRuntimeError::new_err(message),
                }
            })
    }

    /// Measure declared square-planar ligand angles without judging validity.
    fn square_planar_geometry_diagnostics<'py>(
        &self,
        py: Python<'py>,
        coords: Vec<[f64; 3]>,
    ) -> PyResult<Vec<Bound<'py, PyDict>>> {
        use chematic_3d::square_planar_geometry::square_planar_angles;
        use chematic_core::SquarePlanarPermutation;
        if coords.len() != self.inner.atom_count()
            || coords.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(PyValueError::new_err(
                "coordinates must be finite and match the atom count",
            ));
        }
        let geometry = flat_to_coords3d(&coords);
        let reports = square_planar_angles(&self.inner, &geometry)
            .ok_or_else(|| PyValueError::new_err("invalid coordinates"))?;
        reports
            .into_iter()
            .map(|report| {
                let item = PyDict::new(py);
                item.set_item("center", report.center.0)?;
                let tag = report.permutation.map(|p| match p {
                    SquarePlanarPermutation::SP1 => "SP1",
                    SquarePlanarPermutation::SP2 => "SP2",
                    SquarePlanarPermutation::SP3 => "SP3",
                });
                item.set_item("permutation", tag)?;
                item.set_item("neighbors", report.neighbors)?;
                item.set_item("max_unit_triple_product", report.max_unit_triple_product)?;
                item.set_item(
                    "trans",
                    report
                        .trans
                        .into_iter()
                        .map(|a| (a.atoms[0], a.atoms[1], a.degrees))
                        .collect::<Vec<_>>(),
                )?;
                item.set_item(
                    "cis",
                    report
                        .cis
                        .into_iter()
                        .map(|a| (a.atoms[0], a.atoms[1], a.degrees))
                        .collect::<Vec<_>>(),
                )?;
                item.set_item("unavailable", report.unavailable)?;
                Ok(item)
            })
            .collect()
    }
}
