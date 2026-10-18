//! Physical keyword-line handling for the checked Cartesian '+' profile.
use super::*;

pub(super) const MAX_KEYWORD_LINES: usize = 64;
pub(super) const MAX_KEYWORD_BYTES: usize = 65_536;

/// Detect a standalone '+' outside double-quoted keyword values. Other
/// keyword values stay opaque; this does not interpret calculation settings.
pub(super) fn continues(line: &str, number: usize) -> Result<bool, MopacError> {
    if line.len() > MopacParseLimits::default().max_line_len {
        return Err(error(number, "keyword line length limit exceeded"));
    }
    let mut quoted = false;
    let mut start = None;
    let mut previous_plus = false;
    let mut continuation = false;
    for (index, ch) in line
        .char_indices()
        .chain(std::iter::once((line.len(), ' ')))
    {
        if ch.is_control() && ch != '\t' {
            return Err(error(number, "unsupported keyword line control"));
        }
        if ch == '"' {
            quoted = !quoted;
        }
        if ch == '+' && previous_plus && !quoted {
            return Err(error(
                number,
                "double-plus keyword continuation is unsupported",
            ));
        }
        previous_plus = ch == '+' && !quoted;
        if ch.is_whitespace() && !quoted {
            if let Some(begin) = start.take() {
                let token: &str = &line[begin..index];
                let key = token.split('=').next().unwrap_or("");
                if key.eq_ignore_ascii_case("A0") || key == "&" {
                    return Err(error(
                        number,
                        "atomic-unit input and ampersand continuation are unsupported",
                    ));
                }
                continuation |= token == "+";
            }
        } else if start.is_none() {
            start = Some(index);
        }
    }
    if quoted {
        return Err(error(number, "unclosed double-quoted keyword value"));
    }
    Ok(continuation)
}

/// Validate a caller-owned LF keyword block before serialization.
pub(super) fn validate(block: &str) -> Result<usize, MopacError> {
    if block.len() > MAX_KEYWORD_BYTES {
        return Err(error(0, "combined keyword byte limit exceeded"));
    }
    let mut lines = block.split('\n').peekable();
    let mut count = 0;
    while let Some(line) = lines.next() {
        count += 1;
        if count > MAX_KEYWORD_LINES {
            return Err(error(count, "keyword line count limit exceeded"));
        }
        let continuation = continues(line, count)?;
        if continuation != lines.peek().is_some() {
            return Err(error(
                count,
                "each nonfinal keyword line requires '+'; the final line must not continue",
            ));
        }
    }
    Ok(count)
}

/// Consume a variable-length keyword block and exactly two description lines.
pub(super) fn read<'a>(
    lines: &mut impl Iterator<Item = (usize, &'a str)>,
) -> Result<(String, String, String), MopacError> {
    let mut keywords = String::new();
    let mut keyword_lines = 0;
    loop {
        let (index, line) = lines
            .next()
            .ok_or_else(|| error(keyword_lines + 1, "missing keyword continuation line"))?;
        keyword_lines += 1;
        if keyword_lines > MAX_KEYWORD_LINES {
            return Err(error(index + 1, "keyword line count limit exceeded"));
        }
        let extra = usize::from(keyword_lines > 1);
        if keywords.len() + extra + line.len() > MAX_KEYWORD_BYTES {
            return Err(error(index + 1, "combined keyword byte limit exceeded"));
        }
        if extra != 0 {
            keywords.push('\n');
        }
        keywords.push_str(line);
        if !continues(line, index + 1)? {
            break;
        }
    }
    let title = lines
        .next()
        .ok_or_else(|| error(keyword_lines + 1, "two description lines are required"))?
        .1
        .to_owned();
    let comment = lines
        .next()
        .ok_or_else(|| error(keyword_lines + 2, "two description lines are required"))?
        .1
        .to_owned();
    Ok((keywords, title, comment))
}
