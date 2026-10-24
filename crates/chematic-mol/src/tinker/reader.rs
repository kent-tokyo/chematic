use super::*;

pub fn parse_tinker_xyz(text: &str) -> Result<TinkerXyz, TinkerError> {
    parse_tinker_xyz_with_limits(text, &TinkerParseLimits::default())
}

/// Read exactly one frame, refusing trailing records instead of dropping them.
pub fn parse_tinker_xyz_with_limits(
    text: &str,
    limits: &TinkerParseLimits,
) -> Result<TinkerXyz, TinkerError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input exceeds byte limit"));
    }
    for (i, line) in text.lines().enumerate() {
        if line.len() > limits.max_line_len || line.chars().any(|c| c.is_control() && c != '\t') {
            return Err(error(
                i + 1,
                "line exceeds limit or contains control characters",
            ));
        }
    }
    let mut lines = text.lines().enumerate();
    let (_, header) = lines.next().ok_or_else(|| error(1, "missing atom count"))?;
    let header = header.trim_start();
    let end = header.find(char::is_whitespace).unwrap_or(header.len());
    let count = header[..end]
        .parse::<usize>()
        .map_err(|_| error(1, "invalid atom count"))?;
    if count == 0 || count > limits.max_atoms {
        return Err(error(1, "atom count outside supported limits"));
    }
    let title = header[end..].trim_start().to_string();
    let (line_no, first) = lines.next().ok_or_else(|| error(2, "missing atom row"))?;
    let fields: Vec<_> = first.split_whitespace().collect();
    let cell_values: Option<Vec<f64>> = (fields.len() == 6)
        .then(|| {
            fields
                .iter()
                .map(|s| s.parse::<f64>())
                .collect::<Result<_, _>>()
                .ok()
        })
        .flatten();
    let mut atoms = Vec::with_capacity(count.min(4096));
    let cell = if let Some(v) = cell_values {
        Some(TinkerCell {
            lengths_angstrom: [v[0], v[1], v[2]],
            angles_degrees: [v[3], v[4], v[5]],
        })
    } else {
        atoms.push(parse_atom(first, line_no + 1, limits)?);
        None
    };
    while atoms.len() < count {
        let (i, line) = lines
            .next()
            .ok_or_else(|| error(0, "truncated atom table"))?;
        atoms.push(parse_atom(line, i + 1, limits)?);
    }
    if let Some((i, _)) = lines.find(|(_, line)| !line.trim().is_empty()) {
        return Err(error(
            i + 1,
            "trailing content: ARC archives and extra fields are unsupported",
        ));
    }
    let record = TinkerXyz { title, cell, atoms };
    record.validate_with_limits(limits)?;
    Ok(record)
}

fn parse_atom(
    line: &str,
    line_no: usize,
    limits: &TinkerParseLimits,
) -> Result<TinkerAtom, TinkerError> {
    let mut fields = line.split_whitespace();
    let mut next = || {
        fields
            .next()
            .ok_or_else(|| error(line_no, "atom row needs ID, name, x, y, z and type"))
    };
    let index = next()?
        .parse()
        .map_err(|_| error(line_no, "invalid atom ID"))?;
    let name = next()?.to_string();
    let mut position_angstrom = [0.0; 3];
    for v in &mut position_angstrom {
        *v = next()?
            .parse()
            .map_err(|_| error(line_no, "invalid coordinate"))?;
    }
    let atom_type = next()?
        .parse()
        .map_err(|_| error(line_no, "invalid nonnegative atom type"))?;
    let mut neighbors = Vec::new();
    for token in fields {
        if neighbors.len() >= limits.max_neighbors {
            return Err(error(line_no, "neighbor count exceeds limit"));
        }
        neighbors.push(
            token
                .parse()
                .map_err(|_| error(line_no, "invalid neighbor ID"))?,
        );
    }
    Ok(TinkerAtom {
        index,
        name,
        position_angstrom,
        atom_type,
        neighbors,
    })
}
