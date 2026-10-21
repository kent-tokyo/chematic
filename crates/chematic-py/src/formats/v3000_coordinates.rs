use super::*;

/// Parse one V3000 block as (Mol, name, source-order [[x, y, z], ...]).
/// Finite source coordinates in Angstrom, including tiny z and signed zero,
/// survive independently of dimension headers and stereo-perception tolerance.
/// Graph/stereo interpretation is unchanged; opaque metadata is not returned.
/// Accepts str or UTF-8 bytes; input is bounded to 16 MiB.
#[pyfunction]
pub(super) fn from_mol_v3000_with_coords_3d(
    #[pyo3(from_py_with = serialized_text)] block: &str,
) -> PyResult<(Mol, String, Vec<Vec<f64>>)> {
    chematic_mol::parse_mol_v3000_with_coords_3d(block)
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
