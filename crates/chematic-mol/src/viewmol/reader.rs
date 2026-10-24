use super::*;
pub fn parse_viewmol(text: &str) -> Result<ViewmolRecord, ViewmolError> {
    parse_viewmol_with_limits(text, &ViewmolParseLimits::default())
}
pub fn parse_viewmol_with_limits(
    text: &str,
    limits: &ViewmolParseLimits,
) -> Result<ViewmolRecord, ViewmolError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input byte limit exceeded"));
    }
    let mut lines = text.lines().enumerate();
    let mut title = None;
    let mut factor = None;
    let mut atoms = Vec::new();
    let mut in_coord = false;
    let mut ended = false;
    while let Some((i, line)) = lines.next() {
        let n = i + 1;
        check_line(line, n, limits)?;
        let fields: Vec<_> = line.split_ascii_whitespace().collect();
        if fields.is_empty() {
            continue;
        }
        if ended {
            return Err(error(n, "trailing group or record after $end"));
        }
        if fields[0].starts_with('$') {
            in_coord = false;
            match fields.as_slice() {
                ["$title"] if title.is_none() => {
                    let (j, value) = lines.next().ok_or_else(|| error(n, "missing title line"))?;
                    check_line(value, j + 1, limits)?;
                    if value.trim_start().starts_with('$') {
                        return Err(error(
                            j + 1,
                            "expected a literal title line, not a data group",
                        ));
                    }
                    title = Some(value.to_string());
                }
                ["$coord", scale] if factor.is_none() => {
                    if number(scale).map_err(|e| error(n, e.detail))? <= 0. {
                        return Err(error(n, "coordinate factor must be positive"));
                    }
                    factor = Some(scale.to_string());
                    in_coord = true;
                }
                ["$end"] => ended = true,
                _ => return Err(error(n, "unsupported, malformed or duplicate data group")),
            }
        } else if in_coord {
            if atoms.len() >= limits.max_atoms {
                return Err(error(n, "atom limit exceeded"));
            }
            atoms.push(atom(&fields, n)?);
        } else {
            return Err(error(n, "unexpected row outside $coord"));
        }
    }
    if !ended {
        return Err(error(0, "missing $end terminator"));
    }
    let record = ViewmolRecord {
        title,
        factor_to_angstrom: factor.ok_or_else(|| error(0, "missing $coord group"))?,
        atoms,
    };
    record.validate_with_limits(limits)?;
    Ok(record)
}
fn atom(fields: &[&str], line: usize) -> Result<ViewmolAtom, ViewmolError> {
    if fields.len() != 4 && fields.len() != 5 {
        return Err(error(line, "expected x y z element [fixed_axes]"));
    }
    element_number(fields[3]).map_err(|e| error(line, e.detail))?;
    let mut position = [String::new(), String::new(), String::new()];
    for i in 0..3 {
        number(fields[i]).map_err(|e| error(line, e.detail))?;
        position[i] = fields[i].into();
    }
    let fixed_axes = fields.get(4).map(|s| s.to_string());
    validate_axes(fixed_axes.as_deref()).map_err(|e| error(line, e.detail))?;
    Ok(ViewmolAtom {
        element: fields[3].into(),
        position,
        fixed_axes,
    })
}
