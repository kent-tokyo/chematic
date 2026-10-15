use super::*;

#[pymethods]
impl Mol {
    /// RDKit-compatible adjacency in sanitized RDKit atom order.
    /// Weighted dative entries are asymmetric. Only empty_value=0 is supported.
    /// Input caps: 1024 atoms / 16384 bonds, before model preparation.
    #[pyo3(signature = (use_bo = false, *, empty_value = 0))]
    fn rdkit_adjacency_matrix(&self, use_bo: bool, empty_value: i32) -> PyResult<Vec<Vec<f64>>> {
        chematic_smiles::rdkit_adjacency_matrix(&self.inner, use_bo, empty_value)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Unweighted shortest path, including endpoints in sanitized RDKit order.
    /// Disconnected endpoints return []; equal/invalid endpoints raise ValueError.
    /// Ties follow RDKit bond-index neighbor order. Input caps match adjacency.
    fn rdkit_shortest_path(&self, start: usize, end: usize) -> PyResult<Vec<usize>> {
        chematic_smiles::rdkit_shortest_path(&self.inner, start, end)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
}
