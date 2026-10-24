use super::*;

pub fn parse_dlpoly_config(text: &str) -> Result<DlpolyConfig, DlpolyError> {
    parse_dlpoly_config_with_limits(text, &DlpolyConfigLimits::default())
}

/// Exactly one CONFIG; legacy two-field headers and name-only atoms are kept.
pub fn parse_dlpoly_config_with_limits(
    text: &str,
    limits: &DlpolyConfigLimits,
) -> Result<DlpolyConfig, DlpolyError> {
    check_text(text, limits)?;
    let mut lines = text.lines().enumerate();
    let (_, title) = lines.next().ok_or_else(|| error(1, "missing title"))?;
    let (_, header) = lines
        .next()
        .ok_or_else(|| error(2, "missing control record"))?;
    let fields: Vec<_> = header.split_whitespace().collect();
    if !(2..=3).contains(&fields.len()) {
        return Err(error(
            2,
            "expected levcfg imcon [atom count]; extended headers are unsupported",
        ));
    }
    let levcfg = fields[0]
        .parse::<u8>()
        .map_err(|_| error(2, "invalid levcfg"))?;
    let imcon = fields[1]
        .parse::<u8>()
        .map_err(|_| error(2, "invalid imcon"))?;
    if levcfg > 2 || imcon > 3 {
        return Err(error(2, "supported levcfg is 0..2 and imcon is 0..3"));
    }
    let declared_atom_count = fields
        .get(2)
        .map(|s| {
            s.parse::<usize>()
                .map_err(|_| error(2, "invalid atom count"))
        })
        .transpose()?;
    if declared_atom_count.is_some_and(|n| n == 0 || n > limits.max_atoms) {
        return Err(error(2, "atom count outside limits"));
    }
    let cell_vectors = if imcon > 0 {
        Some([
            next_vector(&mut lines)?,
            next_vector(&mut lines)?,
            next_vector(&mut lines)?,
        ])
    } else {
        None
    };
    let mut atoms = Vec::new();
    while let Some((i, line)) = lines.next() {
        if atoms.len() >= limits.max_atoms || declared_atom_count.is_some_and(|n| atoms.len() >= n)
        {
            return Err(error(i + 1, "extra atom/content or atom limit exceeded"));
        }
        atoms.push(read_atom(line, i + 1, levcfg, &mut lines)?);
    }
    let record = DlpolyConfig {
        title: title.to_string(),
        levcfg,
        imcon,
        declared_atom_count,
        cell_vectors,
        atoms,
    };
    record.validate_with_limits(limits)?;
    Ok(record)
}

fn check_text(text: &str, limits: &DlpolyConfigLimits) -> Result<(), DlpolyError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input exceeds byte limit"));
    }
    for (i, line) in text.lines().enumerate() {
        if line.len() > 72 || !line.is_ascii() || line.chars().any(char::is_control) {
            return Err(error(
                i + 1,
                "record needs at most 72 ASCII bytes without controls",
            ));
        }
        if i > 0 && (line.trim().is_empty() || line.trim_start().starts_with(['#', '!'])) {
            return Err(error(i + 1, "blank/comment records are unsupported"));
        }
        if line.split_whitespace().any(|s| s.len() > 40) && i > 0 {
            return Err(error(i + 1, "token exceeds 40-byte format limit"));
        }
    }
    Ok(())
}

type Lines<'a> = std::iter::Enumerate<std::str::Lines<'a>>;

fn next_vector(lines: &mut Lines<'_>) -> Result<[f64; 3], DlpolyError> {
    let (i, line) = lines
        .next()
        .ok_or_else(|| error(0, "truncated vector block"))?;
    let fields: Vec<_> = line.split_whitespace().collect();
    if fields.len() != 3 {
        return Err(error(i + 1, "vector needs exactly three components"));
    }
    let mut v = [0.0; 3];
    for (value, token) in v.iter_mut().zip(fields) {
        *value = token
            .replace(['D', 'd'], "e")
            .parse()
            .map_err(|_| error(i + 1, "invalid vector number"))?;
    }
    finite_vector(&v).map_err(|e| error(i + 1, e.detail))?;
    Ok(v)
}

fn read_atom(
    line: &str,
    line_no: usize,
    level: u8,
    lines: &mut Lines<'_>,
) -> Result<DlpolyAtom, DlpolyError> {
    let fields: Vec<_> = line.split_whitespace().collect();
    if fields.is_empty() {
        return Err(error(line_no, "atom header needs a name"));
    }
    let index = fields
        .get(1)
        .map(|s| {
            s.parse::<usize>()
                .map_err(|_| error(line_no, "invalid atom index"))
        })
        .transpose()?;
    Ok(DlpolyAtom {
        name: fields[0].to_string(),
        header_extra: fields.iter().skip(2).map(|s| s.to_string()).collect(),
        index,
        position: next_vector(lines)?,
        velocity: if level >= 1 {
            Some(next_vector(lines)?)
        } else {
            None
        },
        force: if level == 2 {
            Some(next_vector(lines)?)
        } else {
            None
        },
    })
}
