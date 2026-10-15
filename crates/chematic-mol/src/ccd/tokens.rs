//! Minimal checked CIF1 lexer. Quoted values retain their non-keyword status;
//! comments inside quoted/semicolon values are data, not stripped lines.
use super::{CcdError, CcdParseLimits, invalid};

#[derive(Clone, Copy)]
pub(super) struct Token<'a> {
    pub value: &'a str,
    pub bare: bool,
}
impl Token<'_> {
    pub fn keyword(&self, value: &str) -> bool {
        self.bare && self.value.eq_ignore_ascii_case(value)
    }
    pub fn tag(&self) -> bool {
        self.bare && self.value.starts_with('_')
    }
    pub fn control(&self) -> bool {
        self.tag()
            || self.keyword("loop_")
            || self.keyword("stop_")
            || self.keyword("global_")
            || (self.bare
                && ["data_", "save_"].iter().any(|prefix| {
                    self.value
                        .get(..prefix.len())
                        .is_some_and(|s| s.eq_ignore_ascii_case(prefix))
                }))
    }
}
fn quoted_end(bytes: &[u8], start: usize, quote: u8) -> Result<usize, CcdError> {
    let mut end = start;
    while end < bytes.len() {
        if bytes[end] == b'\n' || bytes[end] == b'\r' {
            return Err(invalid("multiline quoted value"));
        }
        if bytes[end] == quote && (end + 1 == bytes.len() || bytes[end + 1].is_ascii_whitespace()) {
            return Ok(end);
        }
        end += 1;
    }
    Err(invalid("unterminated quoted value"))
}
fn text_end(bytes: &[u8], start: usize) -> Result<usize, CcdError> {
    (start..bytes.len())
        .find(|&i| bytes[i] == b';' && i > 0 && bytes[i - 1] == b'\n')
        .ok_or_else(|| invalid("unterminated semicolon text"))
}
pub(super) fn tokenize<'a>(
    input: &'a str,
    limits: &CcdParseLimits,
) -> Result<Vec<Token<'a>>, CcdError> {
    let bytes = input.as_bytes();
    let mut result = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i].is_ascii_whitespace() {
            i += 1;
            continue;
        }
        if bytes[i] == b'#' {
            while i < bytes.len() && bytes[i] != b'\n' {
                i += 1;
            }
            continue;
        }
        if result.len() == limits.max_tokens {
            return Err(CcdError::ResourceLimit {
                resource: "tokens",
                limit: limits.max_tokens,
            });
        }
        let (start, end, bare) = if b"'\"".contains(&bytes[i]) {
            let start = i + 1;
            let end = quoted_end(bytes, start, bytes[i])?;
            i = end + 1;
            (start, end, false)
        } else if bytes[i] == b';' && (i == 0 || bytes[i - 1] == b'\n') {
            let start = i + 1;
            let end = text_end(bytes, start)?;
            i = end + 1;
            if i < bytes.len() && !bytes[i].is_ascii_whitespace() {
                return Err(invalid("text delimiter suffix"));
            }
            (start, end, false)
        } else {
            let start = i;
            while i < bytes.len() && !bytes[i].is_ascii_whitespace() {
                i += 1;
            }
            (start, i, true)
        };
        result.push(Token {
            value: &input[start..end],
            bare,
        });
    }
    Ok(result)
}
