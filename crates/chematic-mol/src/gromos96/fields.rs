use super::*;
pub(super) fn error(line: usize, detail: impl Into<String>) -> Gromos96Error {
    Gromos96Error {
        line,
        detail: detail.into(),
    }
}
pub(super) fn check_line(
    line: &str,
    n: usize,
    limits: &Gromos96ParseLimits,
) -> Result<(), Gromos96Error> {
    if line.len() > limits.max_line_len || line.chars().any(|c| c.is_control() && c != '\t') {
        return Err(error(n, "line limit exceeded or control character"));
    }
    Ok(())
}
pub(super) fn layout(line: &str) -> bool {
    line.trim().is_empty() || line.starts_with('#')
}
pub(super) fn integer(value: &str, signed: bool) -> Result<(), Gromos96Error> {
    let value = value.trim();
    let digits = if signed {
        value.strip_prefix(['-', '+']).unwrap_or(value)
    } else {
        value
    };
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return Err(error(0, "invalid fixed-width integer"));
    }
    Ok(())
}
pub(super) fn number(value: &str) -> Result<f64, Gromos96Error> {
    let token = value.trim();
    let bytes = token.as_bytes();
    let mut i = usize::from(matches!(bytes.first(), Some(b'+' | b'-')));
    let mut digits = 0;
    while bytes.get(i).is_some_and(u8::is_ascii_digit) {
        i += 1;
        digits += 1;
    }
    if bytes.get(i) == Some(&b'.') {
        i += 1;
        while bytes.get(i).is_some_and(u8::is_ascii_digit) {
            i += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return Err(error(0, "invalid decimal"));
    }
    if matches!(bytes.get(i), Some(b'e' | b'E')) {
        i += 1;
        if matches!(bytes.get(i), Some(b'+' | b'-')) {
            i += 1;
        }
        let start = i;
        while bytes.get(i).is_some_and(u8::is_ascii_digit) {
            i += 1;
        }
        if start == i {
            return Err(error(0, "invalid exponent"));
        }
    }
    if i != bytes.len() {
        return Err(error(0, "invalid decimal suffix"));
    }
    let result: f64 = token
        .parse()
        .map_err(|_| error(0, "invalid nominal decimal"))?;
    let nonzero = token
        .split(['e', 'E'])
        .next()
        .unwrap_or("")
        .bytes()
        .any(|b| matches!(b, b'1'..=b'9'));
    if !result.is_finite() || (nonzero && result == 0.) {
        return Err(error(0, "nonfinite or underflowed decimal"));
    }
    Ok(result)
}
pub(super) fn numbers(line: &str, count: usize) -> Result<Vec<f64>, Gromos96Error> {
    if !line.is_ascii() || line.len() != count * 15 {
        return Err(error(0, "expected complete 15-column numeric fields"));
    }
    (0..count)
        .map(|i| number(&line[15 * i..15 * (i + 1)]))
        .collect()
}
pub(super) fn identity(line: &str) -> Result<Gromos96Identity, Gromos96Error> {
    if !line.is_ascii()
        || line.len() < 24
        || line.as_bytes()[5] != b' '
        || line.as_bytes()[11] != b' '
    {
        return Err(error(0, "invalid 24-column atom identity prefix"));
    }
    integer(&line[..5], true)?;
    integer(&line[17..24], false)?;
    let name = |s: &str| -> Result<String, Gromos96Error> {
        let s = s.trim();
        if s.is_empty() || !s.bytes().all(|c| c.is_ascii_graphic()) {
            return Err(error(0, "invalid atom or residue name"));
        }
        Ok(s.into())
    };
    Ok(Gromos96Identity {
        residue_number: line[..5].trim().into(),
        residue_name: name(&line[6..11])?,
        atom_name: name(&line[12..17])?,
        atom_number: line[17..24].trim().into(),
    })
}
type AtomRows = (Vec<[f64; 3]>, Option<Vec<Gromos96Identity>>);

pub(super) fn rows(
    block: &Gromos96Block,
    limits: &Gromos96ParseLimits,
) -> Result<AtomRows, Gromos96Error> {
    let full = !block.kind.ends_with("RED");
    let mut positions = Vec::new();
    let mut identities = Vec::new();
    for line in &block.lines {
        if line.starts_with('#') {
            continue;
        }
        let shift = if full { 24 } else { 0 };
        if !line.is_ascii() || (line.len() != shift + 45 && line.len() != shift + 48) {
            return Err(error(0, "wrong atom row width"));
        }
        if positions.len() >= limits.max_atoms {
            return Err(error(0, "atom row limit exceeded"));
        }
        let values = if line.len() == shift + 45 {
            numbers(&line[shift..], 3)?
        } else {
            let mut values = Vec::new();
            for i in 0..3 {
                let start = shift + 16 * i;
                if line.as_bytes()[start] != b' ' {
                    return Err(error(0, "expected separator before 15-column field"));
                }
                values.push(number(&line[start + 1..start + 16])?);
            }
            values
        };
        positions.push([values[0], values[1], values[2]]);
        if full {
            identities.push(identity(line)?);
        }
    }
    if positions.is_empty() {
        return Err(error(0, "empty position or velocity block"));
    }
    Ok((positions, if full { Some(identities) } else { None }))
}
pub(super) fn data(block: &Gromos96Block) -> Vec<&str> {
    block
        .lines
        .iter()
        .filter(|line| !line.starts_with('#'))
        .map(String::as_str)
        .collect()
}
pub(super) fn cell(block: &Gromos96Block) -> Result<[[f64; 3]; 3], Gromos96Error> {
    let lines = data(block);
    if lines.len() != 1 {
        return Err(error(
            0,
            "BOX requires one row of three or nine fixed-width values",
        ));
    }
    let n = match lines[0].len() {
        45 => 3,
        135 => 9,
        _ => return Err(error(0, "unsupported BOX width")),
    };
    let v = numbers(lines[0], n)?;
    let mut cell = [[0.; 3]; 3];
    cell[0][0] = v[0];
    cell[1][1] = v[1];
    cell[2][2] = v[2];
    if n == 9 {
        cell[0][1] = v[3];
        cell[0][2] = v[4];
        cell[1][0] = v[5];
        cell[1][2] = v[6];
        cell[2][0] = v[7];
        cell[2][1] = v[8];
    }
    Ok(cell)
}
