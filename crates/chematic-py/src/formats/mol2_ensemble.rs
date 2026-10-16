//! Bounded MOL2 record and all-reference RMSD text APIs.
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

fn limits(
    input_bytes: usize,
    records: usize,
    record_bytes: usize,
    atoms: usize,
    bonds: usize,
) -> chematic_mol::Mol2StreamLimits {
    chematic_mol::Mol2StreamLimits {
        max_input_bytes: input_bytes,
        max_records: records,
        record: chematic_mol::Mol2ParseLimits {
            max_input_bytes: record_bytes,
            max_atoms: atoms,
            max_bonds: bonds,
            ..Default::default()
        },
    }
}

/// Inspect all concatenated MOL2 records, retaining exact source per record.
/// Input is text, not a path. Raises ValueError on any bad record or exceeded
/// limit; returns no partial result. This materializes the input and output.
#[pyfunction]
#[pyo3(signature = (text, *, max_input_bytes=67_108_864, max_records=10_000,
    max_record_bytes=16_777_216, max_atoms=100_000, max_bonds=200_000))]
pub(super) fn parse_mol2_records_json(
    text: &str,
    max_input_bytes: usize,
    max_records: usize,
    max_record_bytes: usize,
    max_atoms: usize,
    max_bonds: usize,
) -> PyResult<String> {
    let limits = limits(
        max_input_bytes,
        max_records,
        max_record_bytes,
        max_atoms,
        max_bonds,
    );
    let records = chematic_mol::Mol2SourceReader::with_limits(text, limits)
        .map(|record| record.map(|record| record.to_json()))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    Ok(serde_json::json!({"schema":"chematic.mol2-records.v1", "records":records}).to_string())
}

/// Heavy-atom RMSD for EVERY reference MOL2 record against ONE probe MOL2.
/// Input is text, not paths. Requires the same source atom IDs/names/types,
/// atom row order and indexed graph. No symmetry search or remapping. Raw RMSD
/// by default; align=True enables a proper rotation/translation (no reflection).
/// Any malformed or mismatched record raises ValueError without a partial result.
#[pyfunction]
#[pyo3(signature = (reference_text, probe_text, *, align=false, max_input_bytes=67_108_864,
    max_records=10_000, max_record_bytes=16_777_216, max_atoms=100_000, max_bonds=200_000))]
#[allow(clippy::too_many_arguments)]
pub(super) fn mol2_reference_rmsds_json(
    reference_text: &str,
    probe_text: &str,
    align: bool,
    max_input_bytes: usize,
    max_records: usize,
    max_record_bytes: usize,
    max_atoms: usize,
    max_bonds: usize,
) -> PyResult<String> {
    let limits = limits(
        max_input_bytes,
        max_records,
        max_record_bytes,
        max_atoms,
        max_bonds,
    );
    let results = chematic_mol::mol2_reference_rmsds(reference_text, probe_text, align, &limits)
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    Ok(serde_json::json!({
        "schema":"chematic.mol2-reference-rmsds.v1",
        "correspondence":"source_id_name_type_and_indexed_graph",
        "alignment":if align { "proper_rotation_translation" } else { "none" },
        "heavy_atoms_only":true, "results":results,
    })
    .to_string())
}
