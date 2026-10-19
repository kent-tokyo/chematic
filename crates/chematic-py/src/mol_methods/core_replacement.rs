use super::*;

type CoreOutput = (Mol, Vec<Option<u32>>, Vec<(u32, u32, u32, u32, u32, u32)>);

#[pymethods]
impl Mol {
    /// Internal bridge for explicit typed core removal.
    fn _replace_core_checked(&self, core: &Mol, matching: Vec<u32>) -> PyResult<CoreOutput> {
        let matching = matching
            .into_iter()
            .map(chematic_core::AtomIdx)
            .collect::<Vec<_>>();
        let result = chematic_chem::core_replacement::replace_core_checked(
            &self.inner,
            &core.inner,
            &matching,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let attachments = result
            .attachments
            .into_iter()
            .map(|a| {
                (
                    a.core_atom.0,
                    a.source_core_atom.0,
                    a.source_anchor.0,
                    a.source_bond.0,
                    a.anchor.0,
                    a.dummy.0,
                )
            })
            .collect();
        let remap = result
            .original_to_output
            .into_iter()
            .map(|a| a.map(|i| i.0))
            .collect();
        Ok((
            Mol {
                num_rings_cache: Default::default(),
                inner: Arc::new(result.molecule),
                props: self.props.clone(),
            },
            remap,
            attachments,
        ))
    }
}
