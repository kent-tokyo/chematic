use super::*;
use pyo3::types::{PyBytes, PyList};

/// Read binary CDX or UTF-8 CDXML, retaining source bytes and CDX omissions.
#[pyfunction]
#[pyo3(signature = (data, *, format="AUTO", allow_metadata_loss=false))]
pub(super) fn from_chemdraw_block<'py>(
    py: Python<'py>,
    data: &[u8],
    format: &str,
    allow_metadata_loss: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let format = chematic_mol::ChemDrawFormat::from_name(format)
        .ok_or_else(|| PyValueError::new_err("format must be AUTO, CDX or CDXML"))?;
    let block = chematic_mol::parse_chemdraw_block(
        data,
        format,
        &chematic_mol::ChemDrawOptions {
            allow_metadata_loss,
        },
    )
    .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let out = PyDict::new(py);
    out.set_item("format", block.format.name())?;
    out.set_item("source_bytes", PyBytes::new(py, &block.source))?;
    out.set_item("projected_cdxml", block.projected_cdxml)?;
    let mut molecules = Vec::new();
    let mut coordinates = Vec::new();
    for (m, c) in block.molecules {
        molecules.push(Mol {
            inner: Arc::new(m),
            props: Default::default(),
        });
        coordinates.push(c.into_iter().map(|(x, y)| [x, y]).collect::<Vec<_>>());
    }
    let omissions = PyList::empty(py);
    for o in block.omissions {
        let row = PyDict::new(py);
        row.set_item("offset", o.offset)?;
        row.set_item("tag", o.tag)?;
        omissions.append(row)?;
    }
    out.set_item("molecules", molecules)?;
    out.set_item("coordinates", coordinates)?;
    out.set_item("omissions", omissions)?;
    Ok(out)
}
