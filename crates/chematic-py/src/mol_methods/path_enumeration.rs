use super::*;

#[pymethods]
impl Mol {
    /// RDKit FindAllPathsOfLengthN with exact ordered bond-set deduplication.
    /// Indices refer to the sanitized RDKit model; limits fail without partial results.
    #[pyo3(signature=(length, *, use_bonds=true, use_hs=false, rooted_at_atom=None, only_shortest_paths=false, max_visits=1_000_000, max_frontier=100_000, max_results=100_000, max_path_indices=1_000_000))]
    #[allow(clippy::too_many_arguments)]
    fn rdkit_paths_of_length(
        &self,
        length: u32,
        use_bonds: bool,
        use_hs: bool,
        rooted_at_atom: Option<usize>,
        only_shortest_paths: bool,
        max_visits: usize,
        max_frontier: usize,
        max_results: usize,
        max_path_indices: usize,
    ) -> PyResult<Vec<Vec<usize>>> {
        let params = chematic_smiles::RdkitPathParams {
            use_bonds,
            use_hs,
            rooted_at_atom,
            only_shortest_paths,
        };
        let limits = chematic_smiles::RdkitPathLimits {
            max_visits,
            max_frontier,
            max_results,
            max_path_indices,
        };
        chematic_smiles::rdkit_paths_of_length_with_limits(&self.inner, length, &params, &limits)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
}
