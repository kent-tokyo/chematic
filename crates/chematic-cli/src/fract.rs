pub(super) fn convert(
    text: &str,
    inspect: bool,
    from_json: bool,
    max_bytes: usize,
) -> Result<String, String> {
    let limits = chematic_mol::FractParseLimits {
        max_input_bytes: max_bytes,
        ..Default::default()
    };
    let record = if from_json {
        chematic_mol::parse_fract_json(text, &limits)
    } else {
        chematic_mol::parse_fract_with_limits(text, &limits)
    }
    .map_err(|e| e.to_string())?;
    if inspect {
        let json = serde_json::to_string_pretty(&record).map_err(|e| e.to_string())?;
        Ok(format!("{json}\n"))
    } else {
        chematic_mol::write_fract_with_limits(&record, &limits).map_err(|e| e.to_string())
    }
}
