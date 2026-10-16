//! SD data fields following an SDF record's MOL CTAB.

/// Parse SD data fields from the section after `M  END`.
///
/// Each field starts with `> <FieldName>` and optional header information.
/// The value is
/// everything on subsequent lines until a blank line (or end of input).
pub(super) fn parse_sd_fields(data: &str) -> std::collections::HashMap<String, String> {
    let mut fields = std::collections::HashMap::new();
    let mut current_key: Option<String> = None;
    let mut current_value = String::new();

    for raw_line in data.lines() {
        let line = raw_line.trim_end_matches('\r');

        if let Some(key) = parse_sd_field_header(line) {
            // Flush previous field.
            if let Some(k) = current_key.take() {
                fields.insert(k, std::mem::take(&mut current_value));
            }
            current_key = Some(key);
        } else if line.is_empty() {
            // Blank line ends the current field's value.
            if let Some(k) = current_key.take() {
                fields.insert(k, std::mem::take(&mut current_value));
            }
        } else if current_key.is_some() {
            if !current_value.is_empty() {
                current_value.push('\n');
            }
            current_value.push_str(line);
        }
    }
    // Flush trailing field with no blank line.
    if let Some(k) = current_key {
        fields.insert(k, current_value);
    }

    fields
}

/// Parse `> <FieldName>` header lines, returning the field name or `None`.
fn parse_sd_field_header(line: &str) -> Option<String> {
    // SDF spec: field headers start with "> " and contain the name in `<...>`.
    // We accept "> <Name>" and also ">  <Name>" (extra spaces).
    let rest = line.strip_prefix('>')?;
    let rest = rest.trim();
    // RDKit's SDWriter appends a record identifier, e.g. `> <Name> (1)`.
    // Header information after the first closing bracket is not field data.
    let (inner, _) = rest.strip_prefix('<')?.split_once('>')?;
    Some(inner.to_string())
}

#[cfg(test)]
mod tests {
    use super::parse_sd_fields;

    #[test]
    fn optional_header_information_does_not_drop_or_rename_fields() {
        for suffix in ["", "  (1) ", " (999)", " registry metadata"] {
            for newline in ["\n", "\r\n"] {
                let input = format!(
                    ">  <Activity>{suffix}{newline}7.2{newline}{newline}> <名前>{suffix}{newline}日本語{newline}second line{newline}{newline}"
                );
                let fields = parse_sd_fields(&input);
                assert_eq!(fields.len(), 2);
                assert_eq!(fields["Activity"], "7.2");
                assert_eq!(fields["名前"], "日本語\nsecond line");
            }
        }
    }

    #[test]
    fn non_headers_and_field_like_value_content_keep_their_meaning() {
        let fields = parse_sd_fields(
            "> <unterminated\nignored\n\n> <value> (4)\ntext with > and <brackets>\n> not a header\n\n",
        );
        assert_eq!(fields.len(), 1);
        assert_eq!(
            fields["value"],
            "text with > and <brackets>\n> not a header"
        );
    }
}
