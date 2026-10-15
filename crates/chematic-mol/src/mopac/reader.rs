use super::*;

pub fn parse_mopac_cartesian(text: &str) -> Result<MopacCartesianGeometry, MopacError> {
    parse_mopac_cartesian_with_limits(text, &MopacParseLimits::default())
}
pub fn parse_mopac_cartesian_with_limits(
    text: &str,
    limits: &MopacParseLimits,
) -> Result<MopacCartesianGeometry, MopacError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input byte limit exceeded"));
    }
    for (i, line) in text.lines().enumerate() {
        if line.len() > limits.max_line_len {
            return Err(error(i + 1, "line length limit exceeded"));
        }
    }
    let mut lines = text.lines();
    let mut header = Vec::new();
    for line in 1..=3 {
        header.push(
            lines
                .next()
                .ok_or_else(|| error(line, "three header lines are required"))?
                .to_string(),
        );
    }
    let mut g = MopacCartesianGeometry {
        keywords: header[0].clone(),
        title: header[1].clone(),
        comment: header[2].clone(),
        atoms: Vec::new(),
        translation_vectors: Vec::new(),
        source_text: String::new(),
    };
    let mut ended = false;
    for (i, line) in lines.enumerate() {
        let number = i + 4;
        if line.trim().is_empty() {
            ended = true;
            continue;
        }
        if ended {
            return Err(error(
                number,
                "additional sections after geometry are unsupported",
            ));
        }
        let (symbol, label, body) = identity(line, number)?;
        let row = parse_row(body, number)?;
        if symbol == "Tv" {
            if label.is_some() || g.translation_vectors.len() >= 3 {
                return Err(error(number, "expected up to three unlabeled Tv records"));
            }
            g.translation_vectors.push(row);
        } else {
            if !g.translation_vectors.is_empty() {
                return Err(error(number, "atoms cannot follow a Tv record"));
            }
            if g.atoms.len() >= limits.max_atoms.min(MopacParseLimits::default().max_atoms) {
                return Err(error(number, "atom limit exceeded"));
            }
            let element=element(symbol).ok_or_else(|| error(number,"expected canonical element/atomic number; isotopes and dummy/sparkle encodings unsupported"))?;
            if let Some(s) = &label {
                validate_label(s).map_err(|e| error(number, e.detail))?;
            }
            g.atoms.push(MopacAtom {
                element,
                label,
                row,
            });
        }
    }
    g.validate()?;
    g.source_text = text.to_string();
    Ok(g)
}

fn element(symbol: &str) -> Option<Element> {
    if symbol.bytes().all(|b| b.is_ascii_digit()) {
        return symbol
            .parse::<u8>()
            .ok()
            .and_then(Element::from_atomic_number);
    }
    Element::from_symbol(symbol).filter(|e| e.symbol() == symbol)
}

fn identity(raw: &str, line: usize) -> Result<(&str, Option<String>, &str), MopacError> {
    let raw = raw.trim_start();
    let boundary = raw
        .find(|c: char| c.is_whitespace() || c == '(')
        .unwrap_or(raw.len());
    let symbol = &raw[..boundary];
    let tail = &raw[boundary..];
    if let Some(label_tail) = tail.strip_prefix('(') {
        let end = label_tail
            .find(')')
            .ok_or_else(|| error(line, "unclosed atom label"))?;
        let label = &label_tail[..end];
        let body = &label_tail[end + 1..];
        if !body.starts_with(char::is_whitespace) {
            return Err(error(line, "whitespace required after atom label"));
        }
        Ok((symbol, Some(label.to_string()), body))
    } else {
        Ok((symbol, None, tail))
    }
}

fn number(text: &str, line: usize) -> Result<f64, MopacError> {
    let v = text
        .replace(['D', 'd'], "e")
        .parse::<f64>()
        .map_err(|_| error(line, "invalid real number"))?;
    if !v.is_finite() {
        return Err(error(line, "non-finite real number"));
    }
    Ok(v)
}
fn parse_row(body: &str, line: usize) -> Result<MopacCartesianRow, MopacError> {
    let f: Vec<_> = body.split_whitespace().collect();
    if f.len() == 3 {
        return Ok(MopacCartesianRow {
            position_angstrom: [
                number(f[0], line)?,
                number(f[1], line)?,
                number(f[2], line)?,
            ],
            optimization_flags: None,
        });
    }
    if f.len() != 6 {
        return Err(error(
            line,
            "expected x y z or x flag y flag z flag, without connectivity/extra fields",
        ));
    }
    let flag = |code| {
        MopacOptimizationFlag::from_code(code)
            .ok_or_else(|| error(line, "unsupported optimization flag"))
    };
    Ok(MopacCartesianRow {
        position_angstrom: [
            number(f[0], line)?,
            number(f[2], line)?,
            number(f[4], line)?,
        ],
        optimization_flags: Some([flag(f[1])?, flag(f[3])?, flag(f[5])?]),
    })
}
