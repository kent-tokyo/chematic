use super::*;

/// Parse V2000 as (Mol, name, source-order [[x, y, z], ...]) in Å.
/// Every coordinate field must contain a finite number. Flat/tiny z and signed
/// zero survive independently of the dimensional header and perception cutoff.
/// Graph/stereo interpretation is unchanged; SD fields/opaque CTAB metadata
/// are not returned. Accepts str or UTF-8 bytes, bounded to 16 MiB.
#[pyfunction]
pub(super) fn from_mol_block_with_coords_3d(
    #[pyo3(from_py_with = serialized_text)] block: &str,
) -> PyResult<(Mol, String, Vec<Vec<f64>>)> {
    chematic_mol::parse_mol_with_coords_3d(block)
        .map(|(mol, metadata, coords)| {
            (
                Mol {
                    num_rings_cache: Default::default(),
                    inner: Arc::new(mol),
                    props: Default::default(),
                },
                metadata.name,
                coords.into_iter().map(|(x, y, z)| vec![x, y, z]).collect(),
            )
        })
        .map_err(|error| PyValueError::new_err(error.to_string()))
}
