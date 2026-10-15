use super::*;
use ndarray::Array2;

pub(crate) fn raw_ogb_graph_dict<'py>(
    py: Python<'py>,
    graph: chematic_smiles::OgbGraphFeatures,
) -> PyResult<Bound<'py, PyDict>> {
    let mut sources = Vec::with_capacity(graph.edge_index.len());
    let mut destinations = Vec::with_capacity(graph.edge_index.len());
    for [source, destination] in &graph.edge_index {
        sources.push(*source);
        destinations.push(*destination);
    }
    let num_nodes = graph.node_feat.len();
    let result = PyDict::new(py);
    result.set_item("node_feat", graph.node_feat)?;
    result.set_item("edge_index", [sources, destinations])?;
    result.set_item("edge_feat", graph.edge_feat)?;
    result.set_item("num_nodes", num_nodes)?;
    Ok(result)
}

pub(crate) fn numpy_ogb_graph_dict<'py>(
    py: Python<'py>,
    graph: chematic_smiles::OgbGraphFeatures,
) -> PyResult<Bound<'py, PyDict>> {
    let mut nodes = Array2::<i64>::zeros((graph.node_feat.len(), 9));
    for (row, features) in graph.node_feat.iter().enumerate() {
        for (column, &value) in features.iter().enumerate() {
            nodes[(row, column)] = value;
        }
    }
    let mut edges = Array2::<i64>::zeros((2, graph.edge_index.len()));
    let mut edge_features = Array2::<i64>::zeros((graph.edge_feat.len(), 3));
    for (column, (&[begin, end], features)) in
        graph.edge_index.iter().zip(&graph.edge_feat).enumerate()
    {
        edges[(0, column)] = begin;
        edges[(1, column)] = end;
        for (feature, &value) in features.iter().enumerate() {
            edge_features[(column, feature)] = value;
        }
    }
    let result = PyDict::new(py);
    result.set_item("node_feat", nodes.into_pyarray(py))?;
    result.set_item("edge_index", edges.into_pyarray(py))?;
    result.set_item("edge_feat", edge_features.into_pyarray(py))?;
    result.set_item("num_nodes", graph.node_feat.len())?;
    Ok(result)
}

#[pymethods]
impl Mol {
    /// Benchmark hook for the same OGB graph before NumPy allocation.
    /// Its Python list conversion remains part of the measured cost.
    #[pyo3(signature = (remove_hs=true))]
    fn _ogb_graph_features_raw<'py>(
        &self,
        py: Python<'py>,
        remove_hs: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        let graph = chematic_smiles::ogb_graph_features_with_options(&self.inner, remove_hs)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        raw_ogb_graph_dict(py, graph)
    }

    /// OGB categorical molecular graph as NumPy arrays.
    ///
    /// The keys match ``ogb.utils.mol.smiles2graph``: ``node_feat`` (N, 9),
    /// ``edge_index`` (2, 2E), ``edge_feat`` (2E, 3), and ``num_nodes``.
    /// Arrays have dtype int64. Removable graph hydrogen atoms are dropped
    /// by default, as in OGB's ``smiles2graph``. Pass ``remove_hs=False`` to
    /// retain all input atoms and their indices.
    #[pyo3(signature = (remove_hs=true))]
    fn ogb_graph_features<'py>(
        &self,
        py: Python<'py>,
        remove_hs: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        let graph = chematic_smiles::ogb_graph_features_with_options(&self.inner, remove_hs)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        numpy_ogb_graph_dict(py, graph)
    }
}
