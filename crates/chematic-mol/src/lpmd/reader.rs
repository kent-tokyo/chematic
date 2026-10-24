use super::*;
pub fn parse_lpmd(text: &str) -> Result<LpmdDocument, LpmdError> {
    parse_lpmd_with_limits(text, &LpmdParseLimits::default())
}
pub fn parse_lpmd_with_limits(
    text: &str,
    limits: &LpmdParseLimits,
) -> Result<LpmdDocument, LpmdError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input byte limit exceeded"));
    }
    // Physical LF records only: CRLF is accepted, embedded CR/control bytes refuse.
    let mut lines = text
        .split_terminator('\n')
        .enumerate()
        .map(|(i, l)| (i + 1, l.strip_suffix('\r').unwrap_or(l)))
        .peekable();
    let (n, header) = next(&mut lines, limits)?;
    if header != ["LPMD", "2.0", "L"] {
        return Err(error(n, "expected uncompressed LPMD 2.0 L"));
    }
    let (n, hdr) = next(&mut lines, limits)?;
    if hdr.first() != Some(&"HDR") {
        return Err(error(n, "expected HDR"));
    }
    let columns: Vec<_> = hdr[1..].iter().map(|s| s.to_string()).collect();
    validate_columns(&columns).map_err(|mut e| {
        e.line = n;
        e
    })?;
    let mut document = LpmdDocument {
        columns,
        frames: Vec::new(),
    };
    let mut total = 0usize;
    while lines.peek().is_some() {
        let (n, count) = next(&mut lines, limits)?;
        if count.is_empty() {
            // Blank records are allowed only as a trailing layout suffix.
            for (n, line) in lines {
                check_line(line, n, limits)?;
                if !line.trim().is_empty() {
                    return Err(error(n, "blank record before another frame"));
                }
            }
            break;
        }
        if document.frames.len() >= limits.max_frames {
            return Err(error(n, "frame limit exceeded"));
        }
        if count.len() != 1 || !count[0].bytes().all(|b| b.is_ascii_digit()) {
            return Err(error(n, "expected unsigned atom count"));
        }
        let count: usize = count[0]
            .parse()
            .map_err(|_| error(n, "atom count overflow"))?;
        total = total
            .checked_add(count)
            .ok_or_else(|| error(n, "atom count overflow"))?;
        if count == 0 || total > limits.max_atoms {
            return Err(error(n, "total atom count outside limits"));
        }
        document
            .frames
            .push(read_frame(&mut lines, &document.columns, count, limits)?);
    }
    document.validate_with_limits(limits)?;
    Ok(document)
}

fn next<'a>(
    lines: &mut impl Iterator<Item = (usize, &'a str)>,
    limits: &LpmdParseLimits,
) -> Result<(usize, Vec<&'a str>), LpmdError> {
    let (n, line) = lines.next().ok_or_else(|| error(0, "truncated document"))?;
    check_line(line, n, limits)?;
    Ok((n, line.split_ascii_whitespace().collect()))
}

fn read_frame<'a>(
    lines: &mut impl Iterator<Item = (usize, &'a str)>,
    columns: &[String],
    count: usize,
    limits: &LpmdParseLimits,
) -> Result<LpmdFrame, LpmdError> {
    let (n, values) = next(lines, limits)?;
    if values.len() != 9 {
        return Err(error(n, "expected nine explicit cell-vector components"));
    }
    let cell = std::array::from_fn(|i| std::array::from_fn(|j| values[3 * i + j].to_string()));
    valid_cell(&cell).map_err(|mut e| {
        e.line = n;
        e
    })?;
    let mut atoms = Vec::new();
    for _ in 0..count {
        let (n, values) = next(lines, limits)?;
        let row: Vec<_> = values.into_iter().map(str::to_string).collect();
        validate_row(columns, &row).map_err(|mut e| {
            e.line = n;
            e
        })?;
        atoms.push(row);
    }
    Ok(LpmdFrame { cell, atoms })
}
