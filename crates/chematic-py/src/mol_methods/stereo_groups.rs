use super::*;
use chematic_core::{AtomIdx, BondIdx, StereoGroup, StereoGroupKind};

type PyGroupSpec = (String, u32, Vec<u32>, Vec<u32>);
type PyGroupInput<'py> = (
    String,
    Bound<'py, PyAny>,
    Vec<Bound<'py, PyAny>>,
    Vec<Bound<'py, PyAny>>,
);

fn group_integer(value: &Bound<'_, PyAny>) -> PyResult<u32> {
    if value.is_instance_of::<pyo3::types::PyBool>() {
        return Err(PyValueError::new_err(
            "stereo group numbers and indices must be integers, not bool",
        ));
    }
    value.extract::<u32>().map_err(|_| {
        PyValueError::new_err(
            "stereo group numbers and indices must be nonnegative 32-bit integers",
        )
    })
}
fn group_spec(input: PyGroupInput<'_>) -> PyResult<(StereoGroup, Vec<BondIdx>)> {
    let (kind, number, atoms, bonds) = input;
    let number = group_integer(&number)?;
    let kind = match kind.as_str() {
        "abs" if number == 0 => StereoGroupKind::Absolute,
        "abs" => {
            return Err(PyValueError::new_err(
                "absolute stereo group number must be 0",
            ));
        }
        "and" => StereoGroupKind::And(number),
        "or" => StereoGroupKind::Or(number),
        _ => {
            return Err(PyValueError::new_err(
                "stereo group kind must be 'abs', 'and' or 'or'",
            ));
        }
    };
    let atoms = atoms
        .iter()
        .map(|v| group_integer(v).map(AtomIdx))
        .collect::<PyResult<_>>()?;
    let bonds = bonds
        .iter()
        .map(|v| group_integer(v).map(BondIdx))
        .collect::<PyResult<_>>()?;
    Ok((StereoGroup::new(kind, atoms), bonds))
}

#[pymethods]
impl Mol {
    /// Enhanced groups as (kind, number, atom_indices, bond_indices) tuples.
    /// Indices address the stored graph; kind is abs/and/or and ABS number is 0.
    /// Returned lists are detached copies, including bond-only groups.
    #[getter]
    fn stereo_groups(&self) -> Vec<PyGroupSpec> {
        self.inner
            .stereo_groups()
            .iter()
            .enumerate()
            .map(|(i, group)| {
                let (kind, number) = match group.kind {
                    StereoGroupKind::Absolute => ("abs", 0),
                    StereoGroupKind::And(n) => ("and", n),
                    StereoGroupKind::Or(n) => ("or", n),
                };
                (
                    kind.into(),
                    number,
                    group.atom_indices.iter().map(|a| a.0).collect(),
                    self.inner
                        .stereo_group_bond_indices(i)
                        .unwrap_or_default()
                        .iter()
                        .map(|b| b.0)
                        .collect(),
                )
            })
            .collect()
    }

    /// Replace enhanced atom/bond stereo groups on a checked independent copy.
    /// No configuration is inverted. Existing stereo frames, graph metadata
    /// and typed molecule properties are retained. This validates grouping
    /// structure only; use V3000 to retain bond members on interchange.
    fn with_stereo_groups(&self, groups: Vec<PyGroupInput<'_>>) -> PyResult<Mol> {
        let groups = groups
            .into_iter()
            .map(group_spec)
            .collect::<PyResult<_>>()?;
        let inner = self
            .inner
            .with_checked_stereo_groups(groups)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Mol {
            inner: Arc::new(inner),
            props: self.props.clone(),
        })
    }
    /// Return the full alternative of an AND/OR tetrahedral/geometric group.
    /// Group membership and typed properties remain intact. Absolute groups,
    /// axes, incomplete configurations and unsupported directions are refused.
    fn invert_stereo_group(&self, index: Bound<'_, PyAny>) -> PyResult<Mol> {
        let index = group_integer(&index)? as usize;
        let inner = chematic_chem::invert_stereo_group(&self.inner, index)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Mol {
            inner: Arc::new(inner),
            props: self.props.clone(),
        })
    }
    /// Return (relative_group_indices, representatives), ordered by inversion mask.
    /// Each AND/OR group flips all tetrahedral/geometric members together.
    /// ABS, ungrouped and unspecified configurations stay fixed; metadata and
    /// typed properties are retained. No deduplication or compound counting.
    #[pyo3(signature=(max_alternatives=64))]
    fn enumerate_stereo_group_alternatives(
        &self,
        #[pyo3(from_py_with = group_integer)] max_alternatives: u32,
    ) -> PyResult<(Vec<usize>, Vec<Mol>)> {
        let max_alternatives = max_alternatives as usize;
        let result =
            chematic_chem::enumerate_stereo_group_alternatives(&self.inner, max_alternatives)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let molecules = result
            .molecules
            .into_iter()
            .map(|inner| Mol {
                inner: Arc::new(inner),
                props: self.props.clone(),
            })
            .collect();
        Ok((result.relative_group_indices, molecules))
    }
    /// Conservatively promote a redundant AND/OR group to ABS using exact
    /// stereo-preserving graph witnesses in every other-group assignment.
    /// Budget exhaustion returns an error; no partial reduction is returned.
    #[pyo3(signature=(index, *, max_visits=1000000))]
    fn reduce_meso_stereo_group<'py>(
        &self,
        py: Python<'py>,
        index: Bound<'py, PyAny>,
        #[pyo3(from_py_with = group_integer)] max_visits: u32,
    ) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        let index = group_integer(&index)? as usize;
        let report =
            chematic_chem::reduce_meso_stereo_group(&self.inner, index, u64::from(max_visits))
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let out = pyo3::types::PyDict::new(py);
        out.set_item(
            "molecule",
            Mol {
                inner: Arc::new(report.molecule),
                props: self.props.clone(),
            },
        )?;
        out.set_item("reduced", report.reduced)?;
        out.set_item("checked_pairs", report.checked_pairs)?;
        out.set_item("search_visits", report.search_visits)?;
        let witnesses = report
            .witnesses
            .into_iter()
            .map(|w| {
                (
                    w.mask,
                    w.atom_permutation
                        .into_iter()
                        .map(|a| a.0)
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<Vec<_>>();
        out.set_item("witnesses", witnesses)?;
        Ok(out)
    }
}
