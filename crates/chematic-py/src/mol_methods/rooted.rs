use super::*;

#[pymethods]
impl Mol {
    /// Non-canonical rooted SMILES and original atom indices, with stereo checks.
    ///
    /// ``order[parsed_index]`` identifies the original graph atom. Supports
    /// <=256 atoms and recorded four-ligand tetrahedral frames; refuses missing
    /// or invalid frames, axial/non-tetrahedral/lone-pair stereo, enhanced groups
    /// and zero/query bonds. Never guesses a constructed graph's stereo frame.
    /// Raises ValueError on refusal. The molecule and compatibility APIs are unchanged.
    fn smiles_rooted_checked(&self, rooted_at_atom: usize) -> PyResult<(String, Vec<usize>)> {
        let (text, order) = chematic_smiles::write_rooted_checked(&self.inner, rooted_at_atom)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok((text, order.into_iter().map(|i| i.0 as usize).collect()))
    }
}
