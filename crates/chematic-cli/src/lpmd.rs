pub(super) fn convert(
    text: &str,
    inspect: bool,
    from_json: bool,
    angstrom: bool,
    max_bytes: usize,
) -> Result<String, String> {
    let l = chematic_mol::LpmdParseLimits {
        max_input_bytes: max_bytes,
        ..Default::default()
    };
    let r = if from_json {
        chematic_mol::parse_lpmd_json(text, &l)
    } else {
        chematic_mol::parse_lpmd_with_limits(text, &l)
    }
    .map_err(|e| e.to_string())?;
    if angstrom {
        let p = r
            .positions_angstrom_with_limits(&l)
            .map_err(|e| e.to_string())?;
        Ok(format!(
            "{}\n",
            serde_json::to_string_pretty(&p).map_err(|e| e.to_string())?
        ))
    } else if inspect {
        Ok(format!(
            "{}\n",
            serde_json::to_string_pretty(&r).map_err(|e| e.to_string())?
        ))
    } else {
        chematic_mol::write_lpmd_with_limits(&r, &l).map_err(|e| e.to_string())
    }
}
