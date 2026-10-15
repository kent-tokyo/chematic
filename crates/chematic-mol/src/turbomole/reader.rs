use super::*;
use crate::UnitCell;

type Rows<'a> = Vec<(usize, &'a str)>;
struct Group<'a> {
    name: &'a str,
    option: &'a str,
    line: usize,
    rows: Rows<'a>,
}

/// Parse one standalone coordinate document. Unknown groups, external files,
/// duplicate/conflicting cell definitions and unsupported periodicities fail.
pub fn parse_turbomole(text: &str) -> Result<TurbomoleGeometry, TurbomoleError> {
    parse_turbomole_with_limits(text, &TurbomoleParseLimits::default())
}

fn groups<'a>(
    text: &'a str,
    limits: &TurbomoleParseLimits,
) -> Result<Vec<Group<'a>>, TurbomoleError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input byte limit exceeded"));
    }
    let mut result: Vec<Group<'a>> = Vec::new();
    let mut ended = false;
    for (i, raw) in text.lines().enumerate() {
        let line = i + 1;
        if raw.len() > limits.max_line_len {
            return Err(error(line, "line length limit exceeded"));
        }
        let row = raw.split('#').next().unwrap_or("").trim();
        if row.is_empty() {
            continue;
        }
        if ended {
            return Err(error(line, "data after $end"));
        }
        if row.starts_with('$') {
            let mut parts = row.split_whitespace();
            let name = parts.next().unwrap();
            let option = parts.next().unwrap_or("");
            if parts.next().is_some() {
                return Err(error(line, "unsupported data-group options"));
            }
            if name == "$end" {
                if !option.is_empty() {
                    return Err(error(line, "$end takes no options"));
                }
                ended = true;
                continue;
            }
            if !matches!(name, "$coord" | "$periodic" | "$lattice" | "$cell") {
                return Err(error(line, format!("unsupported data group {name}")));
            }
            if result.iter().any(|g| g.name == name) {
                return Err(error(line, "duplicate data group"));
            }
            result.push(Group {
                name,
                option,
                line,
                rows: Vec::new(),
            });
        } else {
            let group = result
                .last_mut()
                .ok_or_else(|| error(line, "data outside a group"))?;
            if group.name == "$coord" && group.rows.len() >= limits.max_atoms {
                return Err(error(line, "atom count limit exceeded"));
            }
            if group.name != "$coord" && group.rows.len() >= 3 {
                return Err(error(line, "too many cell or periodicity rows"));
            }
            group
                .rows
                .try_reserve(1)
                .map_err(|_| error(line, "allocation failed"))?;
            group.rows.push((line, row));
        }
    }
    if !ended {
        return Err(error(0, "missing $end"));
    }
    Ok(result)
}

fn finite(s: &str, line: usize) -> Result<f64, TurbomoleError> {
    // Fortran D exponents occur in quantum-chemistry coordinate files.
    let v = s
        .replace(['D', 'd'], "e")
        .parse::<f64>()
        .map_err(|_| error(line, "invalid real number"))?;
    if !v.is_finite() {
        return Err(error(line, "non-finite real number"));
    }
    Ok(v)
}
fn scale(option: &str, line: usize) -> Result<f64, TurbomoleError> {
    match option {
        "" | "bohr" => Ok(TURBOMOLE_BOHR_ANGSTROM),
        "angs" => Ok(1.0),
        _ => Err(error(line, "unsupported unit or external-file option")),
    }
}
fn real_row(row: &str, line: usize, count: usize) -> Result<Vec<f64>, TurbomoleError> {
    let fields: Vec<_> = row.split_whitespace().collect();
    if fields.len() != count {
        return Err(error(line, "wrong numeric field count"));
    }
    fields.iter().map(|s| finite(s, line)).collect()
}
fn periodicity(group: Option<&Group<'_>>) -> Result<u8, TurbomoleError> {
    let Some(g) = group else {
        return Ok(0);
    };
    if !g.rows.is_empty() {
        return Err(error(g.line, "$periodic value must be on its header"));
    }
    match g.option {
        "0" => Ok(0),
        "3" => Ok(3),
        _ => Err(error(
            g.line,
            "only nonperiodic (0) and 3D (3) documents are supported",
        )),
    }
}
fn lattice_rows(g: &Group<'_>) -> Result<[[f64; 3]; 3], TurbomoleError> {
    let factor = scale(g.option, g.line)?;
    if g.rows.len() != 3 {
        return Err(error(g.line, "3D lattice requires three rows"));
    }
    let mut m = [[0.0; 3]; 3];
    for (i, &(line, row)) in g.rows.iter().enumerate() {
        let values = real_row(row, line, 3)?;
        for j in 0..3 {
            m[i][j] = values[j] * factor;
        }
    }
    Ok(m)
}
fn cell_rows(g: &Group<'_>) -> Result<[[f64; 3]; 3], TurbomoleError> {
    let factor = scale(g.option, g.line)?;
    if g.rows.len() != 1 {
        return Err(error(g.line, "3D cell requires one row of six values"));
    }
    let (line, row) = g.rows[0];
    let v = real_row(row, line, 6)?;
    let cell = UnitCell {
        a: v[0] * factor,
        b: v[1] * factor,
        c: v[2] * factor,
        alpha: v[3],
        beta: v[4],
        gamma: v[5],
    };
    cell.validate().map_err(|e| error(line, e.to_string()))?;
    Ok(
        [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)].map(|(x, y, z)| {
            let (x, y, z) = cell.frac_to_cart(x, y, z);
            [x, y, z]
        }),
    )
}
fn lattice(groups: &[Group<'_>]) -> Result<Option<[[f64; 3]; 3]>, TurbomoleError> {
    let periodic = periodicity(groups.iter().find(|g| g.name == "$periodic"))?;
    let l = groups.iter().find(|g| g.name == "$lattice");
    let c = groups.iter().find(|g| g.name == "$cell");
    if l.is_some() && c.is_some() {
        return Err(error(0, "$cell and $lattice conflict"));
    }
    if periodic == 0 && (l.is_some() || c.is_some()) {
        return Err(error(0, "cell requires explicit $periodic 3"));
    }
    if periodic == 0 {
        return Ok(None);
    }
    let m = if let Some(g) = l {
        lattice_rows(g)?
    } else if let Some(g) = c {
        cell_rows(g)?
    } else {
        return Err(error(0, "$periodic 3 requires $lattice or $cell"));
    };
    validate_lattice(&m)?;
    Ok(Some(m))
}
fn atom(row: &str, line: usize) -> Result<TurbomoleAtom, TurbomoleError> {
    let fields: Vec<_> = row.split_whitespace().collect();
    if fields.len() != 4 && fields.len() != 5 {
        return Err(error(line, "atom requires x y z element [f]"));
    }
    let frozen = fields.len() == 5;
    if frozen && fields[4] != "f" {
        return Err(error(line, "unsupported atom flag"));
    }
    let name = fields[3];
    if !name.is_ascii()
        || name.is_empty()
        || name.len() > 2
        || !name.bytes().all(|b| b.is_ascii_alphabetic())
    {
        return Err(error(
            line,
            "invalid element (isotopes/ghost atoms are unsupported)",
        ));
    }
    let canonical = format!(
        "{}{}",
        name[..1].to_ascii_uppercase(),
        name[1..].to_ascii_lowercase()
    );
    let element = Element::from_symbol(&canonical).ok_or_else(|| error(line, "unknown element"))?;
    Ok(TurbomoleAtom {
        element,
        position_angstrom: [
            finite(fields[0], line)?,
            finite(fields[1], line)?,
            finite(fields[2], line)?,
        ],
        frozen,
    })
}
fn converted_position(
    p: [f64; 3],
    g: &Group<'_>,
    lattice: &Option<[[f64; 3]; 3]>,
) -> Result<[f64; 3], TurbomoleError> {
    let result = if g.option == "frac" {
        let m = lattice
            .as_ref()
            .ok_or_else(|| error(g.line, "fractional coordinates require a 3D lattice"))?;
        [0, 1, 2].map(|i| p[0] * m[0][i] + p[1] * m[1][i] + p[2] * m[2][i])
    } else {
        let factor = scale(g.option, g.line)?;
        p.map(|v| v * factor)
    };
    if result.iter().any(|v| !v.is_finite()) {
        return Err(error(g.line, "coordinate conversion overflow"));
    }
    Ok(result)
}

pub fn parse_turbomole_with_limits(
    text: &str,
    limits: &TurbomoleParseLimits,
) -> Result<TurbomoleGeometry, TurbomoleError> {
    let groups = groups(text, limits)?;
    let lattice_angstrom = lattice(&groups)?;
    let g = groups
        .iter()
        .find(|g| g.name == "$coord")
        .ok_or_else(|| error(0, "missing $coord"))?;
    let mut atoms = Vec::new();
    atoms
        .try_reserve(g.rows.len())
        .map_err(|_| error(0, "allocation failed"))?;
    for &(line, row) in &g.rows {
        let mut a = atom(row, line)?;
        a.position_angstrom = converted_position(a.position_angstrom, g, &lattice_angstrom)?;
        atoms.push(a);
    }
    let geometry = TurbomoleGeometry {
        atoms,
        lattice_angstrom,
        source_text: text.to_owned(),
    };
    geometry.validate()?;
    Ok(geometry)
}
