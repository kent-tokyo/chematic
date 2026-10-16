use super::{ChemDrawError as Error, ChemDrawLimits};

/// Bound nesting before delegating to the existing molecular reader. Arbitrary
/// non-XML bytes must not be mistaken for a successfully read empty document.
pub(super) fn validate(input: &str, limits: &ChemDrawLimits) -> Result<(), Error> {
    let mut stack = Vec::<String>::new();
    let mut saw_root = false;
    let mut closed = false;
    let xml = input.strip_prefix('\u{feff}').unwrap_or(input);
    for (index, token) in crate::cdxml_document::logical_cdxml_lines(xml)
        .into_iter()
        .enumerate()
    {
        if index >= limits.max_items {
            return Err(Error::Limit {
                resource: "XML items",
                limit: limits.max_items,
            });
        }
        let t = token.trim();
        if t.starts_with("<?") || t.starts_with("<!") {
            continue;
        }
        let invalid = || Error::Invalid {
            offset: 0,
            detail: "expected one balanced CDXML document",
        };
        if !t.starts_with('<') {
            if stack.is_empty() && !t.is_empty() {
                return Err(invalid());
            } else {
                continue;
            }
        }
        if !t.ends_with('>') {
            return Err(invalid());
        }
        let end = t.starts_with("</");
        let name = t
            .trim_start_matches('<')
            .trim_start_matches('/')
            .split(|c: char| c.is_whitespace() || c == '/' || c == '>')
            .next()
            .unwrap_or("");
        if end {
            if stack.pop().as_deref() != Some(name) {
                return Err(invalid());
            }
            if stack.is_empty() {
                closed = true;
            }
        } else {
            if closed || (stack.is_empty() && (saw_root || name != "CDXML")) {
                return Err(invalid());
            }
            if stack.is_empty() {
                saw_root = true;
            }
            if stack.len() >= limits.max_depth {
                return Err(Error::Limit {
                    resource: "XML depth",
                    limit: limits.max_depth,
                });
            }
            if !t.ends_with("/>") {
                stack.push(name.to_owned());
            } else if stack.is_empty() {
                closed = true;
            }
        }
    }
    if !saw_root || !closed || !stack.is_empty() {
        return Err(Error::Invalid {
            offset: input.len(),
            detail: "expected one balanced CDXML document",
        });
    }
    Ok(())
}
