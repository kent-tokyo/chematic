use super::*;

#[pymethods]
impl Mol {
    /// RDKit FindAtomEnvironmentOfRadiusN, returning ordered bond_indices and
    /// atom_distances in sanitized model numbering. Inputs: <=1024 atoms/16384
    /// bonds. Budgets fail without returning a truncated result or mutating Mol.
    #[pyo3(signature=(radius, rooted_at_atom, *, use_hs=false, enforce_size=true, max_visits=1_000_000, max_frontier=100_000))]
    #[allow(clippy::too_many_arguments)]
    fn rdkit_atom_environment<'py>(
        &self,
        py: Python<'py>,
        radius: u32,
        rooted_at_atom: usize,
        use_hs: bool,
        enforce_size: bool,
        max_visits: usize,
        max_frontier: usize,
    ) -> PyResult<Bound<'py, PyDict>> {
        let limits = chematic_smiles::RdkitAtomEnvironmentLimits {
            max_visits,
            max_frontier,
        };
        let result = chematic_smiles::rdkit_atom_environment_with_limits(
            &self.inner,
            radius,
            rooted_at_atom,
            use_hs,
            enforce_size,
            &limits,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let out = PyDict::new(py);
        out.set_item("bond_indices", result.bond_indices)?;
        out.set_item("atom_distances", result.atom_distances)?;
        Ok(out)
    }
}
