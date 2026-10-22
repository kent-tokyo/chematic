use super::{
    ConversionError, ConversionLimits, ConversionObservation, ConversionPolicy, PreservationStatus,
};

pub(super) fn unemitted_charge_observation(
    record: &crate::Mol2SourceRecord,
    target: &str,
    emit_partial: bool,
    has_method: bool,
) -> Option<ConversionObservation> {
    let has_values = record
        .atoms()
        .iter()
        .any(|atom| atom.partial_charge.is_some());
    if emit_partial || !(has_values || (target == "cjson" && has_method)) {
        return None;
    }
    let unassigned_zeros = target == "cjson"
        && record.charge_type() == Some("NO_CHARGES")
        && has_values
        && record
            .atoms()
            .iter()
            .filter_map(|atom| atom.partial_charge)
            .all(|value| value == 0.0);
    Some(ConversionObservation {
        field: if unassigned_zeros {
            "unassigned_zero_columns"
        } else {
            "partial_charges"
        },
        status: if has_values {
            PreservationStatus::Dropped
        } else {
            PreservationStatus::Unverified
        },
        stage: "write",
        reason: if unassigned_zeros {
            "NO_CHARGES explicitly declares no assigned dataset; raw zero placeholders are omitted, not labeled as computed or assigned partial charges."
        } else if target == "cjson" {
            "Partial charges require a non-NO_CHARGES method label and one value for every atom; missing values are not filled and undeclared charge columns are not labeled as a computed dataset."
        } else {
            "Source charge-column values are not formal charges; the candidate writer does not emit them as partial charges."
        },
    })
}

/// Project one MOL2 record into CJSON, checking mapped graph, xyz, name and
/// complete declared partial charges on reparse. Missing/undeclared partial
/// values and unverified correspondence fail instead of silently dropping data.
/// NO_CHARGES zero placeholders are omitted as explicitly unassigned values.
///
/// This preserves selected fields only: MOL2 atom types/names/IDs, substructures,
/// comments, opaque sections, implicit hydrogens and stereo are not preserved
/// or verified. Use `convert_with_report` to inspect all losses, or immutable
/// same-format retention for source bytes. This is not the report's Strict policy.
pub fn project_mol2_to_cjson(text: &str) -> Result<String, ConversionError> {
    project_mol2_to_cjson_with_limits(text, &ConversionLimits::default())
}

/// Selected-field projection with shared input/atom/bond/output limits.
pub fn project_mol2_to_cjson_with_limits(
    text: &str,
    limits: &ConversionLimits,
) -> Result<String, ConversionError> {
    let report = super::convert_with_report_with_limits(
        text,
        "mol2",
        "cjson",
        ConversionPolicy::AllowLoss,
        limits,
    )?;
    let selected = [
        "coordinates",
        "molecule_name",
        "partial_charges",
        "charge_method",
        "projected_formal_charges",
    ];
    let coordinates_verified = report.observations.iter().any(|observation| {
        observation.field == "coordinates" && observation.status == PreservationStatus::Preserved
    });
    if !report.candidate_graph_verified || !coordinates_verified {
        return Err(ConversionError::Unsupported(
            "MOL2 to CJSON cannot verify mapped graph and coordinates; inspect convert_with_report with allow_loss".into(),
        ));
    }
    if let Some(observation) = report.observations.iter().find(|observation| {
        selected.contains(&observation.field) && observation.status != PreservationStatus::Preserved
    }) {
        return Err(ConversionError::Unsupported(format!(
            "MOL2 to CJSON cannot preserve {}; inspect convert_with_report with allow_loss",
            observation.field,
        )));
    }
    report.output.ok_or_else(|| {
        ConversionError::Unsupported("MOL2 to CJSON projection produced no output".into())
    })
}
