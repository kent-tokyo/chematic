use super::*;
pub fn write_viewmol(record: &ViewmolRecord) -> Result<String, ViewmolError> {
    write_viewmol_with_limits(record, &ViewmolParseLimits::default())
}
/// Current fields, title first if present, unchanged raw factor/coordinates/flags,
/// one space between fields, LF and a final $end line. No coordinate conversion.
pub fn write_viewmol_with_limits(
    record: &ViewmolRecord,
    limits: &ViewmolParseLimits,
) -> Result<String, ViewmolError> {
    record.validate_with_limits(limits)?;
    let mut out = String::new();
    if let Some(title) = &record.title {
        append(&mut out, "$title", limits)?;
        append(&mut out, title, limits)?;
    }
    append(
        &mut out,
        &format!("$coord {}", record.factor_to_angstrom),
        limits,
    )?;
    for atom in &record.atoms {
        let mut row = format!("{} {}", atom.position.join(" "), atom.element);
        if let Some(axes) = &atom.fixed_axes {
            row.push(' ');
            row.push_str(axes);
        }
        append(&mut out, &row, limits)?;
    }
    append(&mut out, "$end", limits)?;
    Ok(out)
}
fn append(out: &mut String, line: &str, limits: &ViewmolParseLimits) -> Result<(), ViewmolError> {
    if line.len() > limits.max_line_len
        || out.len().saturating_add(line.len()).saturating_add(1) > limits.max_input_bytes
    {
        return Err(error(0, "serialized byte or line limit exceeded"));
    }
    out.push_str(line);
    out.push('\n');
    Ok(())
}
