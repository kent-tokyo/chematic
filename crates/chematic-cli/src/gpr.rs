pub(super) fn convert(
    text: &str,
    inspect: bool,
    from_json: bool,
    max_bytes: usize,
) -> Result<String, String> {
    let limits = chematic_mol::GprParseLimits {
        max_input_bytes: max_bytes,
        ..Default::default()
    };
    let record = if from_json {
        chematic_mol::parse_gpr_json(text, &limits)
    } else {
        chematic_mol::parse_gpr_with_limits(text, &limits)
    }
    .map_err(|e| e.to_string())?;
    if inspect {
        Ok(format!(
            "{}\n",
            serde_json::to_string_pretty(&record).map_err(|e| e.to_string())?
        ))
    } else {
        chematic_mol::write_gpr_with_limits(&record, &limits).map_err(|e| e.to_string())
    }
}
