use super::*;

pub fn parse_fract(text: &str) -> Result<FractRecord, FractError> {
    parse_fract_with_limits(text, &FractParseLimits::default())
}
pub fn parse_fract_with_limits(
    text: &str,
    limits: &FractParseLimits,
) -> Result<FractRecord, FractError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input byte limit exceeded"));
    }
    let mut lines = text.lines().enumerate();
    let (_, title) = lines.next().ok_or_else(|| error(1, "missing title line"))?;
    let (_, cell_line) = lines.next().ok_or_else(|| error(2, "missing cell line"))?;
    check_line(title, 1, limits)?;
    check_line(cell_line, 2, limits)?;
    if cell_line.split(',').any(|part| part.trim().is_empty()) {
        return Err(error(2, "empty cell field"));
    }
    let cell_tokens: Vec<_> = cell_line
        .split([',', ' ', '\t'])
        .filter(|s| !s.is_empty())
        .collect();
    let cell = numbers::<6>(&cell_tokens, 2)?;
    let mut atoms = Vec::new();
    for (index, row) in lines {
        let line = index + 1;
        check_line(row, line, limits)?;
        if row.trim().is_empty() {
            continue;
        }
        if atoms.len() >= limits.max_atoms {
            return Err(error(line, "atom limit exceeded"));
        }
        let fields: Vec<_> = row.split_ascii_whitespace().collect();
        if fields.len() != 4 {
            return Err(error(
                line,
                "expected element and exactly three fractional coordinates",
            ));
        }
        validate_element(fields[0], line)?;
        atoms.push(FractAtom {
            element: fields[0].into(),
            fractional: numbers::<3>(&fields[1..], line)?,
        });
    }
    let record = FractRecord {
        title: title.into(),
        cell,
        atoms,
    };
    record.validate_with_limits(limits)?;
    Ok(record)
}
fn numbers<const N: usize>(tokens: &[&str], line: usize) -> Result<[FractNumber; N], FractError> {
    if tokens.len() != N {
        return Err(error(line, format!("expected {N} numeric fields")));
    }
    let values = tokens
        .iter()
        .map(|s| FractNumber::new(*s).map_err(|e| error(line, e.detail)))
        .collect::<Result<Vec<_>, _>>()?;
    values
        .try_into()
        .map_err(|_| error(line, "wrong numeric field count"))
}
fn check_line(line: &str, number: usize, limits: &FractParseLimits) -> Result<(), FractError> {
    if line.len() > limits.max_line_len || line.chars().any(|c| c.is_control() && c != '\t') {
        return Err(error(
            number,
            "control character or line length limit exceeded",
        ));
    }
    Ok(())
}
