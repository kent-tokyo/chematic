use super::*;

pub fn parse_dmol3(text: &str) -> Result<Dmol3Geometry, Dmol3Error> {
    parse_dmol3_with_limits(text, &Dmol3ParseLimits::default())
}
/// Strict standalone .incoor: optional three-row cell, coordinates, mandatory
/// $end. Preambles, output logs, repeated blocks/jobs and unknown fields fail.
pub fn parse_dmol3_with_limits(
    text: &str,
    limits: &Dmol3ParseLimits,
) -> Result<Dmol3Geometry, Dmol3Error> {
    check_input_limits(text, limits)?;
    let hard = Dmol3ParseLimits::default();
    let mut lines = text.lines().enumerate().filter_map(|(i, line)| {
        let line = line.trim();
        (!line.is_empty()).then_some((i + 1, line))
    });
    let (mut number, mut line) = lines
        .next()
        .ok_or_else(|| error(1, "expected $cell vectors or $coordinates"))?;
    let cell_bohr = if line == "$cell vectors" {
        let mut cell = [[0.; 3]; 3];
        for row in &mut cell {
            let (n, text) = lines
                .next()
                .ok_or_else(|| error(number, "incomplete three-vector cell"))?;
            *row = values(text, n)?;
        }
        (number, line) = lines
            .next()
            .ok_or_else(|| error(number, "expected $coordinates after cell"))?;
        Some(cell)
    } else {
        None
    };
    if line != "$coordinates" {
        return Err(error(number, "expected exact $coordinates header"));
    }
    let mut atoms = Vec::new();
    let mut ended = false;
    for (number, line) in lines {
        if ended {
            return Err(error(number, "trailing sections or jobs are unsupported"));
        }
        if line == "$end" {
            ended = true;
            continue;
        }
        if atoms.len() >= limits.max_atoms.min(hard.max_atoms) {
            return Err(error(number, "atom limit exceeded"));
        }
        let f: Vec<_> = line.split_whitespace().collect();
        if f.len() != 4 {
            return Err(error(
                number,
                "expected canonical element followed by three Bohr coordinates",
            ));
        }
        let element = Element::from_symbol(f[0])
            .filter(|e| e.symbol() == f[0])
            .ok_or_else(|| error(number, "invalid canonical element"))?;
        let position_bohr = [
            real(f[1], number)?,
            real(f[2], number)?,
            real(f[3], number)?,
        ];
        atoms
            .try_reserve(1)
            .map_err(|e| error(number, format!("atom allocation failed: {e}")))?;
        atoms.push(Dmol3Atom {
            element,
            position_bohr,
        });
    }
    if !ended {
        return Err(error(0, "missing $end terminator"));
    }
    let mut g = Dmol3Geometry {
        atoms,
        cell_bohr,
        source_text: String::new(),
    };
    g.validate()?;
    g.source_text = text.to_string();
    Ok(g)
}
fn check_input_limits(text: &str, limits: &Dmol3ParseLimits) -> Result<(), Dmol3Error> {
    let hard = Dmol3ParseLimits::default();
    if text.len() > limits.max_input_bytes.min(hard.max_input_bytes) {
        return Err(error(0, "input byte limit exceeded"));
    }
    for (i, line) in text.lines().enumerate() {
        if line.len() > limits.max_line_len.min(hard.max_line_len) {
            return Err(error(i + 1, "line length limit exceeded"));
        }
    }
    Ok(())
}
fn values(line: &str, number: usize) -> Result<[f64; 3], Dmol3Error> {
    let f: Vec<_> = line.split_whitespace().collect();
    if f.len() != 3 {
        return Err(error(number, "cell row requires exactly three Bohr values"));
    }
    Ok([
        real(f[0], number)?,
        real(f[1], number)?,
        real(f[2], number)?,
    ])
}
fn real(text: &str, line: usize) -> Result<f64, Dmol3Error> {
    let v = text
        .replace(['D', 'd'], "e")
        .parse::<f64>()
        .map_err(|_| error(line, "invalid real number"))?;
    if !v.is_finite() {
        return Err(error(line, "non-finite real number"));
    }
    Ok(v)
}
