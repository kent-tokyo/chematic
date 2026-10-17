use super::*;

type FragmentOutput = (Mol, Vec<usize>, Vec<(u32, u32, u32, u16)>);

#[pymethods]
impl Mol {
    /// Internal bridge for the typed checked bond-fragmentation API.
    fn _fragment_single_bonds_checked(
        &self,
        bonds: Vec<u32>,
        labels: Option<Vec<(u16, u16)>>,
    ) -> PyResult<FragmentOutput> {
        let cuts = bonds
            .into_iter()
            .map(chematic_core::BondIdx)
            .collect::<Vec<_>>();
        let result = chematic_chem::bond_fragment::fragment_single_bonds_checked(
            &self.inner,
            &cuts,
            labels.as_deref(),
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let attachments = result
            .attachments
            .into_iter()
            .map(|a| (a.source_bond.0, a.anchor.0, a.dummy.0, a.isotope))
            .collect();
        let molecule = Mol {
            inner: Arc::new(result.molecule),
            props: self.props.clone(),
        };
        Ok((molecule, result.cuts_per_atom, attachments))
    }
}
