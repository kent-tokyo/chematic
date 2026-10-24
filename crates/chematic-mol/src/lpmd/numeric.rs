use super::{LpmdError, error};
pub(super) fn number(token: &str) -> Result<f64, LpmdError> {
    if token.len() > 256 || !decimal(token) {
        return Err(error(0, "invalid decimal or token exceeds 256 bytes"));
    }
    let value: f64 = token.parse().map_err(|_| error(0, "invalid decimal"))?;
    let mantissa = token.split(['e', 'E']).next().unwrap_or("");
    if !value.is_finite() || (value == 0. && mantissa.bytes().any(|b| matches!(b, b'1'..=b'9'))) {
        return Err(error(0, "nonfinite or underflowed decimal"));
    }
    Ok(value)
}
fn decimal(s: &str) -> bool {
    let bytes = s.as_bytes();
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
        return false;
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
            return false;
        }
    }
    i == bytes.len()
}
pub(super) fn valid_cell(cell: &[[String; 3]; 3]) -> Result<(), LpmdError> {
    let mut scaled = [[0f64; 3]; 3];
    for (i, row) in cell.iter().enumerate() {
        for (j, token) in row.iter().enumerate() {
            scaled[i][j] = number(token)?;
        }
        let max = scaled[i].iter().map(|v| v.abs()).fold(0., f64::max);
        if max == 0. {
            return Err(error(0, "zero cell vector"));
        }
        for value in &mut scaled[i] {
            *value /= max;
        }
    }
    let [a, b, c] = scaled;
    let determinant = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0]);
    if determinant == 0. {
        return Err(error(0, "singular or numerically degenerate cell"));
    }
    Ok(())
}
