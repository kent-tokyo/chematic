use super::*;
pub fn parse_unixyz(text: &str) -> Result<UnixyzRecord, UnixyzError> {
    parse_unixyz_with_limits(text, &UnixyzParseLimits::default())
}
pub fn parse_unixyz_with_limits(
    text: &str,
    limits: &UnixyzParseLimits,
) -> Result<UnixyzRecord, UnixyzError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input byte limit exceeded"));
    }
    let mut lines = text.lines().enumerate();
    let (_, title) = lines.next().ok_or_else(|| error(1, "missing title line"))?;
    check_line(title, 1, limits)?;
    let (_, count_line) = lines
        .next()
        .ok_or_else(|| error(2, "missing atom-count line"))?;
    check_line(count_line, 2, limits)?;
    let count = count_line.trim_matches([' ', '\t']);
    if count.is_empty() || !count.bytes().all(|b| b.is_ascii_digit()) {
        return Err(error(2, "expected positive decimal atom count"));
    }
    let n: usize = count
        .parse()
        .map_err(|_| error(2, "atom count outside usize range"))?;
    if n == 0 || n > limits.max_atoms {
        return Err(error(2, "atom count outside limits"));
    }
    let mut atoms = Vec::new();
    for index in 0..n {
        let number = index + 3;
        let (_, row) = lines
            .next()
            .ok_or_else(|| error(number, "truncated atom section"))?;
        check_line(row, number, limits)?;
        atoms.push(atom(row, number)?);
    }
    // Blank atom-section lines are not skipped; the declared count is physical.
    // Reject every trailing line, including whitespace, rather than guessing a
    // title/record boundary or silently selecting only the first molecule.
    if let Some((i, _)) = lines.next() {
        return Err(error(i + 1, "trailing record or extra line"));
    }
    let record = UnixyzRecord {
        title: title.into(),
        atoms,
    };
    record.validate_with_limits(limits)?;
    Ok(record)
}
fn atom(row: &str, line: usize) -> Result<UnixyzAtom, UnixyzError> {
    let f: Vec<_> = row.split_ascii_whitespace().collect();
    if f.len() != 4 {
        return Err(error(
            line,
            "expected atomic number and exactly three coordinates",
        ));
    }
    if !f[0].bytes().all(|b| b.is_ascii_digit()) {
        return Err(error(line, "expected decimal atomic number"));
    }
    let z: u8 = f[0]
        .parse()
        .map_err(|_| error(line, "invalid atomic number"))?;
    if !(1..=118).contains(&z) {
        return Err(error(line, "atomic number must be 1..118"));
    }
    let mut position_angstrom = [String::new(), String::new(), String::new()];
    for i in 0..3 {
        unixyz_nominal_value(f[i + 1]).map_err(|e| error(line, e.detail))?;
        position_angstrom[i] = f[i + 1].into();
    }
    Ok(UnixyzAtom {
        atomic_number: z,
        position_angstrom,
    })
}
