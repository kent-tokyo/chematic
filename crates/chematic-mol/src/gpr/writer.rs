use super::*;
pub fn write_gpr(record: &GprRecord) -> Result<String, GprError> {
    write_gpr_with_limits(record, &GprParseLimits::default())
}
/// Current validated fields, retaining IDs, row order, numeric tokens and absence
/// of partial charges. Output is gpr 100, !Info 1, canonical spaces and LF.
pub fn write_gpr_with_limits(
    record: &GprRecord,
    limits: &GprParseLimits,
) -> Result<String, GprError> {
    record.validate_with_limits(limits)?;
    let mut out = String::new();
    for line in [
        "!Header gpr 100".to_string(),
        "!Info 1".into(),
        format!("!Atoms {}", record.atoms.len()),
    ] {
        append(&mut out, &line, limits)?;
    }
    for a in &record.atoms {
        append(&mut out, &format!("{} {}", a.id, a.atomic_number), limits)?;
    }
    append(&mut out, &format!("!Bonds {}", record.bonds.len()), limits)?;
    for b in &record.bonds {
        append(
            &mut out,
            &format!("{} {} {}", b.atoms[0], b.atoms[1], b.order.code()),
            limits,
        )?;
    }
    append(&mut out, "!Coord", limits)?;
    for c in &record.coordinates {
        append(
            &mut out,
            &format!("{} {}", c.atom_id, c.position_nm.join(" ")),
            limits,
        )?;
    }
    if let Some(charges) = &record.partial_charges {
        append(&mut out, "!Charges", limits)?;
        for c in charges {
            append(&mut out, &format!("{} {}", c.atom_id, c.value), limits)?;
        }
    }
    append(&mut out, "!End", limits)?;
    Ok(out)
}
fn append(out: &mut String, line: &str, limits: &GprParseLimits) -> Result<(), GprError> {
    if line.len() > limits.max_line_len
        || out.len().saturating_add(line.len()).saturating_add(1) > limits.max_input_bytes
    {
        return Err(error(0, "serialized byte or line limit exceeded"));
    }
    out.push_str(line);
    out.push('\n');
    Ok(())
}
