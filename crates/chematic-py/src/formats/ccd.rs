use crate::Mol;
use pyo3::prelude::*;
use std::sync::Arc;

fn error(e: chematic_mol::CcdError) -> PyErr {
    match &e {
        chematic_mol::CcdError::ResourceLimit { .. } => {
            crate::errors::input_error("resource_limit", "ccd_limit", "ccd", e.to_string())
        }
        chematic_mol::CcdError::UnknownCharge { .. } => {
            crate::errors::input_error("unsupported", "ccd_unknown_charge", "ccd", e.to_string())
        }
        chematic_mol::CcdError::StereoLoss => {
            crate::errors::input_error("unsupported", "ccd_stereo_loss", "ccd", e.to_string())
        }
        _ => crate::errors::malformed("ccd", e.to_string()),
    }
}
/// Parse one CCD template to JSON; retains IDs, charges, bonds, both optional
/// coordinate families, stereo annotations and unmodeled atom/bond fields.
#[pyfunction]
#[pyo3(signature = (text, *, max_input_bytes=None, max_atoms=None, max_bonds=None))]
pub(super) fn parse_ccd_json(
    text: &str,
    max_input_bytes: Option<usize>,
    max_atoms: Option<usize>,
    max_bonds: Option<usize>,
) -> PyResult<String> {
    let mut limits = chematic_mol::CcdParseLimits::default();
    if let Some(n) = max_input_bytes {
        limits.max_input_bytes = n;
    }
    if let Some(n) = max_atoms {
        limits.max_atoms = n;
    }
    if let Some(n) = max_bonds {
        limits.max_bonds = n;
    }
    let record = chematic_mol::parse_ccd_component_with_limits(text, &limits).map_err(error)?;
    serde_json::to_string(&record)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
}
/// Explicit template graph. Refuses unassigned R/S or E/Z annotations unless
/// stereo loss is explicitly accepted. Does not infer implicit H or bonds.
#[pyfunction]
#[pyo3(signature = (text, *, allow_stereo_loss=false))]
pub(super) fn from_ccd(text: &str, allow_stereo_loss: bool) -> PyResult<Mol> {
    let record = chematic_mol::parse_ccd_component(text).map_err(error)?;
    let inner = record.to_molecule(allow_stereo_loss).map_err(error)?;
    Ok(Mol {
        inner: Arc::new(inner),
        props: Default::default(),
    })
}
