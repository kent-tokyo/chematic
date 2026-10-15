use super::*;

/// Parse exactly one document; trailing records or unknown fields are errors.
pub fn parse_exyz(text: &str) -> Result<ExyzGeometry, ExyzError> {
    parse_exyz_with_limits(text, &ExyzParseLimits::default())
}

pub fn parse_exyz_with_limits(
    text: &str,
    limits: &ExyzParseLimits,
) -> Result<ExyzGeometry, ExyzError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input byte limit exceeded"));
    }
    for (i, line) in text.lines().enumerate() {
        if line.len() > limits.max_line_len {
            return Err(error(i + 1, "line length limit exceeded"));
        }
    }
    let mut lines = text.lines().enumerate();
    let first = next_line(&mut lines)?;
    let count: usize = first
        .1
        .trim()
        .parse()
        .map_err(|_| error(1, "invalid atom count"))?;
    if count == 0 || count > limits.max_atoms || count > ExyzParseLimits::default().max_atoms {
        return Err(error(1, "atom count outside supported limits"));
    }
    // A valid atom row needs at least eight bytes; reject huge declared counts
    // before any count-dependent allocation on tiny malformed inputs.
    if count > text.len() / 8 {
        return Err(error(1, "declared atom count exceeds input size"));
    }
    let header = next_line(&mut lines)?.1.to_string();
    let (pbc, virtual_enabled) = header_flags(&header)?;
    let mut atoms = Vec::new();
    for _ in 0..count {
        let (line, raw) = next_line(&mut lines)?;
        atoms.push(parse_atom(raw, line, virtual_enabled)?);
    }
    let cell = if pbc {
        Some(parse_cell(&mut lines)?)
    } else {
        None
    };
    for (i, line) in lines {
        if !line.trim().is_empty() {
            return Err(error(i + 1, "trailing data after single EXYZ document"));
        }
    }
    let geometry = ExyzGeometry {
        header,
        atoms,
        cell,
        source_text: text.to_string(),
    };
    geometry.validate()?;
    Ok(geometry)
}

fn next_line<'a>(
    lines: &mut impl Iterator<Item = (usize, &'a str)>,
) -> Result<(usize, &'a str), ExyzError> {
    lines
        .next()
        .map(|(i, s)| (i + 1, s))
        .ok_or_else(|| error(0, "unexpected end of document"))
}

fn triple(fields: &[&str], line: usize) -> Result<[f64; 3], ExyzError> {
    let mut result = [0.0; 3];
    for (out, field) in result.iter_mut().zip(fields) {
        *out = field
            .parse::<f64>()
            .map_err(|_| error(line, "invalid Cartesian number"))?;
        if !out.is_finite() {
            return Err(error(line, "non-finite Cartesian number"));
        }
    }
    Ok(result)
}

fn parse_atom(raw: &str, line: usize, virtual_enabled: bool) -> Result<ExyzAtom, ExyzError> {
    let fields: Vec<_> = raw.split_whitespace().collect();
    if fields.len() != 4 && !(fields.len() == 5 && fields[4] == "VIRTUAL" && virtual_enabled) {
        return Err(error(
            line,
            "expected symbol x y z and optional enabled VIRTUAL marker",
        ));
    }
    let element = Element::from_symbol(fields[0])
        .filter(|e| e.symbol() == fields[0])
        .ok_or_else(|| error(line, "expected a canonical element symbol"))?;
    Ok(ExyzAtom {
        element,
        position_angstrom: triple(&fields[1..4], line)?,
        virtual_atom: fields.len() == 5,
    })
}

fn parse_cell<'a>(
    lines: &mut impl Iterator<Item = (usize, &'a str)>,
) -> Result<ExyzCell, ExyzError> {
    let (line, separator) = next_line(lines)?;
    if !separator.trim().is_empty() {
        return Err(error(line, "expected blank separator before cell"));
    }
    let mut rows = [[0.0; 3]; 4];
    for (label, row) in ["Vector1", "Vector2", "Vector3", "Offset"]
        .iter()
        .zip(&mut rows)
    {
        let (line, raw) = next_line(lines)?;
        let fields: Vec<_> = raw.split_whitespace().collect();
        if fields.len() != 4 || fields[0] != *label {
            return Err(error(line, format!("expected {label} x y z")));
        }
        *row = triple(&fields[1..4], line)?;
    }
    Ok(ExyzCell {
        lattice_angstrom: [rows[0], rows[1], rows[2]],
        offset_angstrom: rows[3],
    })
}
