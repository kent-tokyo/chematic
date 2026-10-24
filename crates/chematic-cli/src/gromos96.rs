pub(super) fn convert(
    text: &str,
    inspect: bool,
    from_json: bool,
    angstrom: bool,
    max_bytes: usize,
    length_unit: Option<&str>,
) -> Result<String, String> {
    let l = chematic_mol::Gromos96ParseLimits {
        max_input_bytes: max_bytes,
        ..Default::default()
    };
    let r = if from_json {
        chematic_mol::parse_gromos96_json(text, &l)
    } else {
        chematic_mol::parse_gromos96_with_limits(
            text,
            length_unit
                .unwrap_or("unspecified")
                .parse::<chematic_mol::Gromos96LengthUnit>()
                .map_err(|e| e.to_string())?,
            &l,
        )
    }
    .map_err(|e| e.to_string())?;
    if angstrom {
        let p = r
            .geometry_angstrom_with_limits(&l)
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
        chematic_mol::write_gromos96_with_limits(&r, &l).map_err(|e| e.to_string())
    }
}
