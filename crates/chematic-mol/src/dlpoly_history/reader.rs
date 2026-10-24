use super::*;
use std::iter::{Enumerate, Peekable};
use std::str::Lines;

type Rows<'a> = Peekable<Enumerate<Lines<'a>>>;
pub fn parse_dlpoly_history(text: &str) -> Result<DlpolyHistory, DlpolyHistoryError> {
    parse_dlpoly_history_with_limits(text, &DlpolyHistoryLimits::default())
}
pub fn parse_dlpoly_history_with_limits(
    text: &str,
    limits: &DlpolyHistoryLimits,
) -> Result<DlpolyHistory, DlpolyHistoryError> {
    if text.len() > limits.max_input_bytes || !text.is_ascii() {
        return Err(error(0, "input byte limit or non-ASCII input"));
    }
    for (i, line) in text.lines().enumerate() {
        if line.len() > 4096 || line.bytes().any(|b| b.is_ascii_control() && b != b'\t') {
            return Err(error(i + 1, "line length/control character limit"));
        }
    }
    let mut rows = text.lines().enumerate().peekable();
    let (_, title) = next(&mut rows)?;
    let (line, raw) = next(&mut rows)?;
    let w = fields(raw, 3, &[10, 10, 10], line)?;
    let keytrj = integer(w[0], line)?;
    let imcon = integer(w[1], line)?;
    let atom_count = integer(w[2], line)?;
    let mut g = DlpolyHistory {
        title: title.to_owned(),
        keytrj,
        imcon,
        atom_count,
        frames: Vec::new(),
    };
    if keytrj > 2 || imcon > 3 || atom_count == 0 || atom_count > limits.max_atoms {
        return Err(error(line, "unsupported header key or atom count"));
    }
    while let Some((_, s)) = rows.peek() {
        if s.trim().is_empty() {
            if rows.any(|(_, s)| !s.trim().is_empty()) {
                return Err(error(0, "blank records within trajectory"));
            }
            break;
        }
        if g.frames.len() >= limits.max_frames
            || g.frames
                .len()
                .checked_add(1)
                .and_then(|n| n.checked_mul(atom_count))
                .is_none_or(|n| n > limits.max_atom_records)
        {
            return Err(error(0, "frame/total atom-record limit"));
        }
        g.frames.push(frame(&mut rows, &g)?);
    }
    g.validate_with_limits(limits)?;
    Ok(g)
}
fn next<'a>(rows: &mut Rows<'a>) -> Result<(usize, &'a str), DlpolyHistoryError> {
    rows.next()
        .map(|(i, s)| (i + 1, s))
        .ok_or_else(|| error(0, "truncated HISTORY record"))
}
fn frame(rows: &mut Rows<'_>, g: &DlpolyHistory) -> Result<DlpolyHistoryFrame, DlpolyHistoryError> {
    let (line, raw) = next(rows)?;
    let w = fields(raw, 6, &[8, 10, 10, 10, 10, 12], line)?;
    if !w[0].eq_ignore_ascii_case("timestep") {
        return Err(error(line, "expected timestep header"));
    }
    let step = integer(w[1], line)?;
    let n: usize = integer(w[2], line)?;
    let level: u8 = integer(w[3], line)?;
    let imcon: u8 = integer(w[4], line)?;
    if (n, level, imcon) != (g.atom_count, g.keytrj, g.imcon) {
        return Err(error(line, "frame/header count or keys disagree"));
    }
    let timestep = number(w[5], line)?;
    let cell_vectors = if imcon > 0 {
        Some([vector_row(rows)?, vector_row(rows)?, vector_row(rows)?])
    } else {
        None
    };
    let mut atoms = Vec::new();
    for _ in 0..n {
        atoms.push(atom(rows, level)?);
    }
    Ok(DlpolyHistoryFrame {
        step,
        timestep,
        cell_vectors,
        atoms,
    })
}
fn atom(rows: &mut Rows<'_>, level: u8) -> Result<DlpolyHistoryAtom, DlpolyHistoryError> {
    let (line, raw) = next(rows)?;
    let w = fields(raw, 4, &[8, 10, 12, 12], line)?;
    let name = w[0].to_owned();
    let index = integer(w[1], line)?;
    let mass = number(w[2], line)?;
    let charge = number(w[3], line)?;
    let position = vector_row(rows)?;
    let velocity = if level >= 1 {
        Some(vector_row(rows)?)
    } else {
        None
    };
    let force = if level == 2 {
        Some(vector_row(rows)?)
    } else {
        None
    };
    Ok(DlpolyHistoryAtom {
        name,
        index,
        mass,
        charge,
        position,
        velocity,
        force,
    })
}
fn vector_row(rows: &mut Rows<'_>) -> Result<[f64; 3], DlpolyHistoryError> {
    let (line, raw) = next(rows)?;
    let w = fields(raw, 3, &[12, 12, 12], line)?;
    Ok([
        number(w[0], line)?,
        number(w[1], line)?,
        number(w[2], line)?,
    ])
}
fn fields<'a>(
    raw: &'a str,
    count: usize,
    widths: &[usize],
    line: usize,
) -> Result<Vec<&'a str>, DlpolyHistoryError> {
    let words: Vec<_> = raw.split_whitespace().collect();
    if words.len() == count {
        return Ok(words);
    }
    let width: usize = widths.iter().sum();
    if raw.len() < width || !raw[width..].trim().is_empty() {
        return Err(error(
            line,
            "unexpected record field count or unsupported extension",
        ));
    }
    let mut start = 0;
    Ok(widths
        .iter()
        .map(|n| {
            let s = &raw[start..start + n];
            start += n;
            s.trim()
        })
        .collect())
}
fn integer<T: std::str::FromStr>(s: &str, line: usize) -> Result<T, DlpolyHistoryError> {
    s.parse().map_err(|_| error(line, "invalid integer field"))
}
fn number(s: &str, line: usize) -> Result<f64, DlpolyHistoryError> {
    if !s.contains('.') {
        return Err(error(
            line,
            "real fields need an explicit decimal point; implied decimals are not guessed",
        ));
    }
    let n: f64 = s
        .replace(['D', 'd'], "e")
        .parse()
        .map_err(|_| error(line, "invalid real field"))?;
    if !n.is_finite()
        || (n == 0.0
            && s.split(['e', 'E', 'd', 'D'])
                .next()
                .unwrap_or("")
                .bytes()
                .any(|b| matches!(b, b'1'..=b'9')))
    {
        return Err(error(line, "nonfinite or underflowing real field"));
    }
    Ok(n)
}
