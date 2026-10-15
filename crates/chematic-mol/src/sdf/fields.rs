//! SD data fields following an SDF record's MOL CTAB.

/// Parse SD data fields from the section after `M  END`.
///
/// Each field starts with `> <FieldName>` on its own line. The value is
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
    let inner = rest.strip_prefix('<')?.strip_suffix('>')?;
    Some(inner.to_string())
}
