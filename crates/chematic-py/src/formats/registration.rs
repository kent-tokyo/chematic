use super::*;

pub(super) fn register_quantum_formats(m: &Bound<'_, PyModule>) -> PyResult<()> {
    register_pyfunctions!(
        m,
        from_gjf,
        parse_gjf_geometry,
        roundtrip_gjf_record,
        roundtrip_gjf_document,
        parse_gjf_document_json,
        write_gjf_cartesian,
        parse_gaussian_log,
        write_gjf,
        parse_orca_input,
        parse_orca_cartesian,
        write_orca_input,
        parse_orca_output,
        parse_qcschema_molecule,
        write_qcschema_molecule,
        chematic_to_qc_molecule,
        qc_molecule_to_chematic,
        parse_atomic_input,
        write_atomic_input,
        parse_atomic_result,
        write_atomic_result,
    )
}
