//! Checked atom fields which otherwise defaulted away when malformed.
use crate::MolParseError;

fn invalid(line: usize, field: &str, detail: &str) -> MolParseError {
    MolParseError::InvalidAtomLine {
        line,
        detail: format!("V3000 {field}: {detail}"),
    }
}

/// Identify attribute starts outside quoted strings and parenthesized values.
/// CTfile literal quotes inside a quoted value are doubled, not backslash escaped.
pub(super) fn scoped_tokens<'a>(
    tokens: &'a [&'a str],
) -> impl Iterator<Item = (bool, &'a str)> + 'a {
    let (mut quoted, mut depth) = (false, 0usize);
    tokens.iter().map(move |&token| {
        let top = !quoted && depth == 0;
        let bytes = token.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            if bytes[i] == b'"' {
                if quoted && bytes.get(i + 1) == Some(&b'"') {
                    i += 2;
                    continue;
                }
                quoted = !quoted;
            } else if !quoted {
                match bytes[i] {
                    b'(' => depth += 1,
                    b')' => depth = depth.saturating_sub(1),
                    _ => {}
                }
            }
            i += 1;
        }
        (top, token)
    })
}

fn scalar<T: std::str::FromStr>(
    tokens: &[&str],
    key: &str,
    line: usize,
) -> Result<Option<T>, MolParseError> {
    let mut result = None;
    for (top, token) in scoped_tokens(tokens) {
        if !top {
            continue;
        }
        if token == key {
            return Err(invalid(line, key, "missing '=value'"));
        }
        let Some(value) = token
            .strip_prefix(key)
            .and_then(|rest| rest.strip_prefix('='))
        else {
            continue;
        };
        if result.is_some() {
            return Err(invalid(line, key, "duplicate attribute"));
        }
        result = Some(value.parse().map_err(|_| {
            invalid(
                line,
                key,
                "expected an integer within the supported storage range",
            )
        })?);
    }
    Ok(result)
}

pub(super) fn charge_and_isotope(
    tokens: &[&str],
    line: usize,
) -> Result<(i8, Option<u16>), MolParseError> {
    let charge = scalar::<i8>(tokens, "CHG", line)?.unwrap_or(0);
    let isotope = scalar::<u16>(tokens, "MASS", line)?;
    if isotope == Some(0) {
        return Err(invalid(line, "MASS", "expected a positive isotope mass"));
    }
    Ok((charge, isotope))
}

pub(super) fn atom_map(text: &str, line: usize) -> Result<Option<u16>, MolParseError> {
    let number = text
        .parse::<u16>()
        .map_err(|_| invalid(line, "atom map", "expected an integer in 0..=65535"))?;
    Ok((number != 0).then_some(number))
}
