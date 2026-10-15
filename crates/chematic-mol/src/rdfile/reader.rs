use super::*;

struct Cursor<'a> {
    text: &'a str,
    offset: usize,
    line: usize,
}
impl<'a> Cursor<'a> {
    fn peek(&self) -> Option<&'a str> {
        if self.offset == self.text.len() {
            return None;
        }
        Some(
            self.text[self.offset..]
                .split('\n')
                .next()
                .unwrap()
                .trim_end_matches('\r'),
        )
    }
    fn take(&mut self) -> Result<&'a str, RdfileError> {
        let row = self
            .peek()
            .ok_or_else(|| error(self.line, "unexpected end of input"))?;
        self.offset += self.text[self.offset..]
            .find('\n')
            .map_or(self.text.len() - self.offset, |i| i + 1);
        self.line += 1;
        Ok(row)
    }
}
fn argument<'a>(row: &'a str, keyword: &str) -> Option<&'a str> {
    let rest = row.strip_prefix(keyword)?;
    if rest.is_empty() {
        Some("")
    } else if rest.starts_with([' ', '\t']) {
        Some(rest.trim_start_matches([' ', '\t']))
    } else {
        None
    }
}
fn record_start(row: &str) -> bool {
    ["$MFMT", "$RFMT", "$MIREG", "$MEREG", "$RIREG", "$REREG"]
        .iter()
        .any(|k| argument(row, k).is_some())
}
fn short_header(row: &str, line: usize) -> Result<(), RdfileError> {
    if row.len() > 80 || !row.is_ascii() || row.contains(['\r', '\n', '\0']) {
        return Err(error(
            line,
            "identifier/date/field header must be ASCII within 80 columns",
        ));
    }
    Ok(())
}
fn identifier(row: &str, line: usize) -> Result<RdfileRecord, RdfileError> {
    short_header(row, line)?;
    for (keyword, kind, embedded, internal) in [
        ("$MFMT", RdfileKind::Molecule, true, false),
        ("$RFMT", RdfileKind::Reaction, true, false),
        ("$MIREG", RdfileKind::Molecule, false, true),
        ("$MEREG", RdfileKind::Molecule, false, false),
        ("$RIREG", RdfileKind::Reaction, false, true),
        ("$REREG", RdfileKind::Reaction, false, false),
    ] {
        let Some(mut value) = argument(row, keyword) else {
            continue;
        };
        let mut is_internal = internal;
        let mut registry_present = !embedded;
        if embedded && !value.is_empty() {
            registry_present = true;
            let (i, e) = match kind {
                RdfileKind::Molecule => ("$MIREG", "$MEREG"),
                RdfileKind::Reaction => ("$RIREG", "$REREG"),
            };
            if let Some(v) = argument(value, i) {
                is_internal = true;
                value = v;
            } else if let Some(v) = argument(value, e) {
                value = v;
            } else {
                return Err(error(line, "unsupported embedded registry identifier"));
            }
            if value.contains('$') {
                return Err(error(
                    line,
                    "only one registry identifier per record is supported",
                ));
            }
        }
        if registry_present
            && (value.is_empty() || is_internal && !value.bytes().all(|b| b.is_ascii_digit()))
        {
            return Err(error(
                line,
                "registry reference requires a nonempty external or decimal internal number",
            ));
        }
        let mut record = RdfileRecord {
            kind,
            internal_registry: None,
            external_registry: None,
            structure: None,
            data: Vec::new(),
        };
        if !value.is_empty() {
            if is_internal {
                record.internal_registry = Some(value.to_owned());
            } else {
                record.external_registry = Some(value.to_owned());
            }
        }
        return Ok(record);
    }
    Err(error(
        line,
        "expected a molecule or reaction record identifier",
    ))
}
fn mol_block<'a>(cursor: &mut Cursor<'a>) -> Result<&'a str, RdfileError> {
    let start = cursor.offset;
    for _ in 0..4 {
        cursor.take()?;
    }
    loop {
        let row = cursor.take()?;
        if row.trim_end() == "M  END" {
            break;
        }
        if row.starts_with('$') {
            return Err(error(
                cursor.line - 1,
                "missing MOL terminator before RDFile keyword",
            ));
        }
    }
    Ok(&cursor.text[start..cursor.offset])
}
fn rxn_block<'a>(
    cursor: &mut Cursor<'a>,
    limits: &RdfileParseLimits,
) -> Result<&'a str, RdfileError> {
    let start = cursor.offset;
    if cursor.take()? != "$RXN" {
        return Err(error(
            cursor.line - 1,
            "only embedded RXN V2000 is supported",
        ));
    }
    for _ in 0..3 {
        cursor.take()?;
    }
    let counts = cursor.take()?;
    let fields: Vec<_> = counts.split_whitespace().collect();
    if !(fields.len() == 2 || fields.len() == 3 && fields[2] == "0") {
        return Err(error(
            cursor.line - 1,
            "RXN requires exactly reactant and product counts; agents are unsupported",
        ));
    }
    let n = fields[..2]
        .iter()
        .try_fold(0usize, |n, s| {
            s.parse::<usize>().ok().and_then(|v| n.checked_add(v))
        })
        .ok_or_else(|| error(cursor.line - 1, "invalid RXN component count"))?;
    if n > limits.max_reaction_components {
        return Err(error(cursor.line - 1, "reaction component limit exceeded"));
    }
    for _ in 0..n {
        if cursor.take()? != "$MOL" {
            return Err(error(cursor.line - 1, "expected RXN $MOL component marker"));
        }
        let block = mol_block(cursor)?;
        if block.lines().nth(3).is_some_and(|s| s.contains("V3000")) {
            return Err(error(
                cursor.line - 1,
                "RXN V2000 components must be V2000 MOL blocks",
            ));
        }
    }
    Ok(&cursor.text[start..cursor.offset])
}
fn datum_line(row: &str, line: usize) -> Result<(&str, bool), RdfileError> {
    if !row.is_ascii() {
        return Err(error(line, "data rows must be ASCII"));
    }
    // Exact 81-column '+' rows use CTfile folding. Longer plain rows are a
    // common unwrapped-data extension (including Indigo output), kept verbatim.
    let join = row.len() == 81 && row.as_bytes()[80] == b'+';
    Ok((if join { &row[..80] } else { row }, join))
}
fn data(cursor: &mut Cursor<'_>) -> Result<RdfileData, RdfileError> {
    let line = cursor.line;
    let header = cursor.take()?;
    short_header(header, line)?;
    let name = argument(header, "$DTYPE")
        .filter(|s| !s.is_empty())
        .ok_or_else(|| error(line, "empty data field name"))?
        .to_owned();
    let first = cursor.take()?;
    let (first, mut join) = datum_line(first, cursor.line - 1)?;
    let value = argument(first, "$DATUM")
        .ok_or_else(|| error(cursor.line - 1, "$DTYPE must be followed by $DATUM"))?;
    if value.starts_with("$MFMT") {
        return Err(error(
            cursor.line - 1,
            "embedded molecules inside $DATUM are unsupported",
        ));
    }
    let mut value = value.to_owned();
    while let Some(row) = cursor.peek() {
        if row.starts_with('$') {
            if join {
                return Err(error(
                    cursor.line,
                    "missing data continuation before keyword",
                ));
            }
            if argument(row, "$DTYPE").is_none() && !record_start(row) {
                return Err(error(cursor.line, "unsupported keyword in data"));
            }
            break;
        }
        let (row, next_join) = datum_line(row, cursor.line)?;
        if !join {
            value.push('\n');
        }
        value.push_str(row);
        join = next_join;
        cursor.take()?;
    }
    if join {
        return Err(error(cursor.line, "missing final data continuation"));
    }
    Ok(RdfileData { name, value })
}

pub fn parse_rdfile(text: &str) -> Result<RdfileDocument, RdfileError> {
    parse_rdfile_with_limits(text, &RdfileParseLimits::default())
}
/// Checked in-memory RDFile v1. Source text and embedded blocks retain their
/// original bytes; ASCII data fields join documented column-81 continuations.
pub fn parse_rdfile_with_limits(
    text: &str,
    limits: &RdfileParseLimits,
) -> Result<RdfileDocument, RdfileError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input byte limit exceeded"));
    }
    if text.contains('\0') {
        return Err(error(0, "NUL bytes are unsupported"));
    }
    for (i, line) in text.lines().enumerate() {
        if line.len() > limits.max_line_len {
            return Err(error(i + 1, "physical line length limit exceeded"));
        }
    }
    let mut cursor = Cursor {
        text,
        offset: 0,
        line: 1,
    };
    if cursor.take()? != "$RDFILE 1" {
        return Err(error(1, "expected $RDFILE 1"));
    }
    let date = cursor.take()?;
    short_header(date, 2)?;
    let date_comment = argument(date, "$DATM")
        .ok_or_else(|| error(2, "missing $DATM"))?
        .to_owned();
    let mut records = Vec::new();
    while let Some(row) = cursor.peek() {
        if records.len() >= limits.max_records {
            return Err(error(cursor.line, "record count limit exceeded"));
        }
        let line = cursor.line;
        let embedded = argument(row, "$MFMT").is_some() || argument(row, "$RFMT").is_some();
        let mut record = identifier(cursor.take()?, line)?;
        if embedded {
            let block = match record.kind {
                RdfileKind::Molecule => mol_block(&mut cursor)?,
                RdfileKind::Reaction => rxn_block(&mut cursor, limits)?,
            };
            record.structure = Some(block.to_owned());
            record
                .smiles()
                .map_err(|e| error(line, format!("invalid embedded structure: {e}")))?;
        }
        while cursor
            .peek()
            .is_some_and(|s| argument(s, "$DTYPE").is_some())
        {
            if record.data.len() >= limits.max_fields_per_record {
                return Err(error(cursor.line, "data-field count limit exceeded"));
            }
            record
                .data
                .try_reserve(1)
                .map_err(|_| error(cursor.line, "allocation failed"))?;
            record.data.push(data(&mut cursor)?);
        }
        records
            .try_reserve(1)
            .map_err(|_| error(line, "allocation failed"))?;
        records.push(record);
    }
    Ok(RdfileDocument {
        date_comment,
        records,
        source_text: text.to_owned(),
    })
}
