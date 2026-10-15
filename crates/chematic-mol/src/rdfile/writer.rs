use super::*;

fn scalar(s: &str, nonempty: bool) -> Result<(), RdfileError> {
    if !s.is_ascii() || s.contains(['\n', '\r', '\0']) || nonempty && s.is_empty() {
        return Err(error(
            0,
            "header fields must be single-line ASCII with a nonempty identifier",
        ));
    }
    Ok(())
}
fn header(out: &mut String, keyword: &str, value: &str) -> Result<(), RdfileError> {
    scalar(value, false)?;
    let line = if value.is_empty() {
        keyword.to_owned()
    } else {
        format!("{keyword} {value}")
    };
    if line.len() > 80 {
        return Err(error(0, "header exceeds supported 80-column layout"));
    }
    out.push_str(&line);
    out.push('\n');
    Ok(())
}
fn datum(out: &mut String, value: &str, wrap_data: bool) -> Result<(), RdfileError> {
    if !value.is_ascii() || value.contains(['\r', '\0']) {
        return Err(error(0, "data values must be ASCII text with LF newlines"));
    }
    for (i, line) in value.split('\n').enumerate() {
        if i > 0 && line.starts_with('$') {
            return Err(error(
                0,
                "data continuation cannot begin with a reserved '$' keyword",
            ));
        }
        let physical = if i == 0 {
            format!("$DATUM {line}")
        } else {
            line.to_owned()
        };
        let mut remaining = physical.as_str();
        while wrap_data && remaining.len() > 80 {
            out.push_str(&remaining[..80]);
            out.push_str("+\n");
            remaining = &remaining[80..];
            if remaining.starts_with('$') {
                return Err(error(
                    0,
                    "wrapped data would collide with an RDFile keyword",
                ));
            }
        }
        out.push_str(remaining);
        out.push('\n');
    }
    Ok(())
}
fn output_size(
    document: &RdfileDocument,
    limits: &RdfileParseLimits,
) -> Result<usize, RdfileError> {
    let mut estimate = 32usize.saturating_add(document.date_comment.len());
    for record in &document.records {
        if record.data.len() > limits.max_fields_per_record {
            return Err(error(0, "data-field count limit exceeded"));
        }
        for size in [
            record.structure.as_ref().map_or(0, String::len),
            record.internal_registry.as_ref().map_or(0, String::len),
            record.external_registry.as_ref().map_or(0, String::len),
            64,
        ] {
            estimate = estimate.saturating_add(size);
        }
        for field in &record.data {
            estimate = estimate
                .saturating_add(field.name.len())
                .saturating_add(field.value.len())
                .saturating_add(field.value.len() / 80 * 2)
                .saturating_add(24);
        }
        if estimate > limits.max_input_bytes {
            return Err(error(0, "output byte limit exceeded"));
        }
    }
    if estimate > limits.max_input_bytes {
        return Err(error(0, "output byte limit exceeded"));
    }
    Ok(estimate)
}

fn write_record(
    out: &mut String,
    record: &RdfileRecord,
    options: RdfileWriteOptions,
) -> Result<(), RdfileError> {
    if record.internal_registry.is_some() && record.external_registry.is_some() {
        return Err(error(
            0,
            "only one registry identifier per record is supported",
        ));
    }
    let (fmt, ireg, ereg) = match record.kind {
        RdfileKind::Molecule => ("$MFMT", "$MIREG", "$MEREG"),
        RdfileKind::Reaction => ("$RFMT", "$RIREG", "$REREG"),
    };
    let registry = if let Some(value) = &record.internal_registry {
        scalar(value, true)?;
        if !value.bytes().all(|b| b.is_ascii_digit()) {
            return Err(error(0, "internal registry number must be decimal"));
        }
        Some((ireg, value))
    } else if let Some(value) = &record.external_registry {
        scalar(value, true)?;
        Some((ereg, value))
    } else {
        None
    };
    if let Some(block) = &record.structure {
        let identifier = registry
            .map(|(k, v)| format!("{k} {v}"))
            .unwrap_or_default();
        header(out, fmt, &identifier)?;
        out.push_str(block);
        if !block.ends_with('\n') {
            out.push('\n');
        }
    } else if let Some((keyword, value)) = registry {
        header(out, keyword, value)?;
    } else {
        return Err(error(
            0,
            "reference-only record requires a registry identifier",
        ));
    }
    for field in &record.data {
        scalar(&field.name, true)?;
        header(out, "$DTYPE", &field.name)?;
        datum(out, &field.value, options.wrap_data)?;
    }
    Ok(())
}

/// Canonical LF RDFile serialization of current records, never stale source
/// replay. Embedded structure bytes remain intact. Refuses unsupported fields,
/// delimiter injection, invalid structures and output beyond default limits.
pub fn write_rdfile(document: &RdfileDocument) -> Result<String, RdfileError> {
    write_rdfile_with_options(document, RdfileWriteOptions::default())
}
/// Explicit unwrapped-data profile for consumers that do not implement '+'
/// folding. Round parsing still refuses ambiguous keyword/continuation data.
pub fn write_rdfile_with_options(
    document: &RdfileDocument,
    options: RdfileWriteOptions,
) -> Result<String, RdfileError> {
    let limits = RdfileParseLimits::default();
    if document.records.len() > limits.max_records {
        return Err(error(0, "record count limit exceeded"));
    }
    let estimate = output_size(document, &limits)?;
    let mut out = String::new();
    out.try_reserve(estimate)
        .map_err(|_| error(0, "output allocation failed"))?;
    out.push_str("$RDFILE 1\n");
    header(&mut out, "$DATM", &document.date_comment)?;
    for record in &document.records {
        write_record(&mut out, record, options)?;
    }
    // Round parsing catches injected record boundaries and verifies chemistry.
    let reparsed = parse_rdfile(&out)?;
    let mut expected = document.records.clone();
    for record in &mut expected {
        if let Some(block) = &mut record.structure
            && !block.ends_with('\n')
        {
            block.push('\n');
        }
    }
    if reparsed.date_comment != document.date_comment || reparsed.records != expected {
        return Err(error(
            0,
            "output would change record boundaries, fields or structure payloads",
        ));
    }
    Ok(out)
}
