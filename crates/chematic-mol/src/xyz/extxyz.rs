use super::*;

// ---------------------------------------------------------------------------
// Extended XYZ (extxyz)
// ---------------------------------------------------------------------------

/// Tokenize an extxyz comment line into ordered `(key, value)` pairs.
/// Quoted values (`key="a b c"`) may contain spaces; a backslash inside a
/// quoted value escapes the following character verbatim (so `\"` is a
/// literal quote, `\\` a literal backslash), matching ASE's own
/// `key_val_str_to_dict` char-by-char escape handling. Unquoted values run
/// to the next whitespace. A bare key with no `=value` at all defaults to
/// value `"T"` (ASE's own convention for an implicit boolean flag -- e.g.
/// `constrained` alone means `constrained=T`), rather than being an error.
/// `lineno` only annotates errors.
fn tokenize_info_line(line: &str, lineno: usize) -> Result<Vec<(String, String)>, XyzError> {
    let chars: Vec<char> = line.chars().collect();
    let n = chars.len();
    let mut i = 0;
    let mut pairs = Vec::new();

    while i < n {
        while i < n && chars[i].is_whitespace() {
            i += 1;
        }
        if i >= n {
            break;
        }
        let key_start = i;
        while i < n && chars[i] != '=' && !chars[i].is_whitespace() {
            i += 1;
        }
        let key: String = chars[key_start..i].iter().collect();
        if key.is_empty() {
            // Only reachable via a stray '=' with nothing before it (e.g.
            // "=foo" or "==bar") -- whitespace already advanced past above,
            // so a real key never starts empty.
            return Err(XyzError::InvalidInfoField {
                line: lineno,
                detail: "empty key before '='".to_string(),
            });
        }
        if i >= n || chars[i] != '=' {
            pairs.push((key, "T".to_string()));
            continue;
        }
        i += 1; // consume '='

        let value: String = if i < n && chars[i] == '"' {
            i += 1;
            let mut v = String::new();
            loop {
                if i >= n {
                    return Err(XyzError::UnterminatedQuote { line: lineno, key });
                }
                match chars[i] {
                    '"' => {
                        i += 1;
                        break;
                    }
                    '\\' if i + 1 < n => {
                        v.push(chars[i + 1]);
                        i += 2;
                    }
                    c => {
                        v.push(c);
                        i += 1;
                    }
                }
            }
            v
        } else {
            let val_start = i;
            while i < n && !chars[i].is_whitespace() {
                i += 1;
            }
            chars[val_start..i].iter().collect()
        };
        pairs.push((key, value));
    }
    Ok(pairs)
}

/// Parse a `Lattice=` value: exactly 9 whitespace-separated finite numbers.
fn parse_lattice(raw: &str, lineno: usize) -> Result<[f64; 9], XyzError> {
    let parts: Vec<&str> = raw.split_whitespace().collect();
    if parts.len() != 9 {
        return Err(XyzError::InvalidLattice {
            line: lineno,
            detail: format!("expected 9 numbers, found {}", parts.len()),
        });
    }
    let mut out = [0.0; 9];
    for (i, p) in parts.iter().enumerate() {
        let v: f64 = p.parse().map_err(|_| XyzError::InvalidLattice {
            line: lineno,
            detail: format!("cannot parse '{p}' as a number"),
        })?;
        if !v.is_finite() {
            return Err(XyzError::InvalidLattice {
                line: lineno,
                detail: format!("non-finite value '{p}'"),
            });
        }
        out[i] = v;
    }
    Ok(out)
}

/// One column parsed from a `Properties=` value: `(name, kind, count)`.
type PropertyColumn = (String, XyzPropertyKind, usize);

/// Parse a `Properties=` value into its column list. Requires the first two
/// columns to be exactly `species:S:1` and `pos:R:3` (universal in
/// ASE-written extxyz) -- any further columns become [`XyzFrame::properties`].
fn parse_properties(raw: &str, lineno: usize) -> Result<Vec<PropertyColumn>, XyzError> {
    let parts: Vec<&str> = raw.split(':').collect();
    if parts.is_empty() || !parts.len().is_multiple_of(3) {
        return Err(XyzError::InvalidProperties {
            line: lineno,
            detail: format!("'{raw}' is not a multiple-of-3 colon-separated list"),
        });
    }
    let mut columns = Vec::with_capacity(parts.len() / 3);
    for chunk in parts.chunks(3) {
        let name = chunk[0].to_string();
        let kind = match chunk[1] {
            "R" => XyzPropertyKind::Real,
            "I" => XyzPropertyKind::Integer,
            "S" => XyzPropertyKind::String,
            "L" => XyzPropertyKind::Logical,
            other => {
                return Err(XyzError::InvalidProperties {
                    line: lineno,
                    detail: format!("unknown column type '{other}' for '{name}'"),
                });
            }
        };
        let count: usize = chunk[2].parse().map_err(|_| XyzError::InvalidProperties {
            line: lineno,
            detail: format!("invalid column count '{}' for '{name}'", chunk[2]),
        })?;
        if count == 0 {
            return Err(XyzError::InvalidProperties {
                line: lineno,
                detail: format!("column '{name}' declares 0 components"),
            });
        }
        columns.push((name, kind, count));
    }
    let valid_prefix = columns.len() >= 2
        && columns[0].0 == "species"
        && columns[0].1 == XyzPropertyKind::String
        && columns[0].2 == 1
        && columns[1].0 == "pos"
        && columns[1].1 == XyzPropertyKind::Real
        && columns[1].2 == 3;
    if !valid_prefix {
        return Err(XyzError::InvalidProperties {
            line: lineno,
            detail: format!("'{raw}' must start with 'species:S:1:pos:R:3'"),
        });
    }
    Ok(columns)
}

type ExtxyzHeader = (Option<[f64; 9]>, Vec<PropertyColumn>, Vec<(String, String)>);

fn parse_extxyz_header(comment_line: &str) -> Result<ExtxyzHeader, XyzError> {
    let default_columns = || {
        vec![
            ("species".to_string(), XyzPropertyKind::String, 1),
            ("pos".to_string(), XyzPropertyKind::Real, 3),
        ]
    };
    if !comment_line.contains('=') {
        return Ok((None, default_columns(), Vec::new()));
    }

    let mut lattice = None;
    let mut columns = None;
    let mut info = Vec::new();
    let mut seen = HashSet::new();
    for (key, value) in tokenize_info_line(comment_line, 2)? {
        if !seen.insert(key.clone()) {
            return Err(XyzError::DuplicateInfoKey { line: 2, key });
        }
        match key.as_str() {
            "Lattice" => lattice = Some(parse_lattice(&value, 2)?),
            "Properties" => columns = Some(parse_properties(&value, 2)?),
            _ => info.push((key, value)),
        }
    }
    Ok((lattice, columns.unwrap_or_else(default_columns), info))
}

fn parse_property_value(
    raw: &str,
    kind: XyzPropertyKind,
    lineno: usize,
    property: &str,
) -> Result<XyzValue, XyzError> {
    let invalid = || XyzError::InvalidPropertyValue {
        line: lineno,
        property: property.to_string(),
        raw: raw.to_string(),
    };
    match kind {
        XyzPropertyKind::Real => {
            let value: f64 = raw.parse().map_err(|_| invalid())?;
            value
                .is_finite()
                .then_some(XyzValue::Real(value))
                .ok_or_else(invalid)
        }
        XyzPropertyKind::Integer => raw.parse().map(XyzValue::Integer).map_err(|_| invalid()),
        XyzPropertyKind::String => Ok(XyzValue::Str(raw.to_string())),
        XyzPropertyKind::Logical => match raw {
            "T" | "True" | "true" | "1" => Ok(XyzValue::Logical(true)),
            "F" | "False" | "false" | "0" => Ok(XyzValue::Logical(false)),
            _ => Err(invalid()),
        },
    }
}

fn parse_extxyz_atom_line(
    raw_line: &str,
    lineno: usize,
    columns: &[PropertyColumn],
) -> Result<(XyzAtom, Vec<Vec<XyzValue>>), XyzError> {
    let tokens: Vec<&str> = raw_line.split_whitespace().collect();
    let expected: usize = columns.iter().map(|(_, _, count)| count).sum();
    if tokens.len() != expected {
        return Err(XyzError::PropertyColumnMismatch {
            line: lineno,
            expected,
            found: tokens.len(),
        });
    }

    let element = element_from_xyz_symbol(tokens[0]).ok_or_else(|| XyzError::UnknownElement {
        symbol: tokens[0].to_string(),
        line: lineno,
    })?;
    let parse_coord = |raw: &str, axis| {
        let value: f64 = raw.parse().map_err(|_| XyzError::InvalidAtomLine {
            line: lineno,
            detail: format!("cannot parse {axis} coordinate '{raw}'"),
        })?;
        value
            .is_finite()
            .then_some(value)
            .ok_or_else(|| XyzError::NonFiniteCoordinate {
                line: lineno,
                axis,
                raw: raw.to_string(),
            })
    };
    let atom = XyzAtom {
        element,
        x: parse_coord(tokens[1], "x")?,
        y: parse_coord(tokens[2], "y")?,
        z: parse_coord(tokens[3], "z")?,
    };

    let mut token_index = 4;
    let mut extra_values = Vec::with_capacity(columns.len().saturating_sub(2));
    for (name, kind, count) in &columns[2..] {
        let mut values = Vec::with_capacity(*count);
        for _ in 0..*count {
            values.push(parse_property_value(
                tokens[token_index],
                *kind,
                lineno,
                name,
            )?);
            token_index += 1;
        }
        extra_values.push(values);
    }
    Ok((atom, extra_values))
}

/// Parse one extxyz frame starting at `text`, returning the frame and the
/// remaining unparsed input. Shared by [`parse_extxyz`] and [`ExtxyzReader`].
pub(super) fn parse_one_frame_ext(text: &str) -> Result<(XyzFrame, &str), XyzError> {
    let (count_line, rest) = next_line(text).ok_or(XyzError::InvalidCountLine {
        line: 1,
        raw: String::new(),
    })?;
    let declared: usize = count_line
        .trim()
        .parse()
        .map_err(|_| XyzError::InvalidCountLine {
            line: 1,
            raw: count_line.to_string(),
        })?;

    let (comment_line, mut rest) = next_line(rest).ok_or(XyzError::MissingCommentLine)?;
    let comment = comment_line.to_string();

    let (lattice, columns, info) = parse_extxyz_header(comment_line)?;
    let extra_columns = &columns[2..];

    let mut atoms = Vec::with_capacity(declared);
    let mut extra_values: Vec<Vec<Vec<XyzValue>>> = vec![Vec::new(); extra_columns.len()];

    for lineno in (3usize..).take(declared) {
        let Some((raw_line, next_rest)) = next_line(rest) else {
            break;
        };
        rest = next_rest;

        let (atom, atom_properties) = parse_extxyz_atom_line(raw_line, lineno, &columns)?;
        atoms.push(atom);
        for (column_values, values) in extra_values.iter_mut().zip(atom_properties) {
            column_values.push(values);
        }
    }

    if atoms.len() != declared {
        return Err(XyzError::AtomCountMismatch {
            declared,
            found: atoms.len(),
        });
    }

    let properties = extra_columns
        .iter()
        .zip(extra_values)
        .map(|((name, kind, count), values)| XyzProperty {
            name: name.clone(),
            kind: *kind,
            count: *count,
            values,
        })
        .collect();

    Ok((
        XyzFrame {
            atoms,
            comment,
            lattice,
            properties,
            info,
        },
        rest,
    ))
}

/// Parse a single extended-XYZ (extxyz) frame from the start of `text`.
/// Trailing content is ignored, same as [`parse_xyz`] -- use
/// [`ExtxyzReader`] to iterate every frame in a multi-frame file. See the
/// module-level "Extended XYZ" docs for what `Lattice=`/`Properties=`
/// support covers.
pub fn parse_extxyz(text: &str) -> Result<XyzFrame, XyzError> {
    parse_extxyz_with_limits(text, &XyzParseLimits::default())
}

/// Parse a single Extended XYZ frame while enforcing resource limits.
pub fn parse_extxyz_with_limits(text: &str, limits: &XyzParseLimits) -> Result<XyzFrame, XyzError> {
    validate_xyz_limits(text, limits)?;
    parse_one_frame_ext(text).map(|(frame, _rest)| frame)
}

/// Iterator over frames in a multi-frame extxyz (trajectory) string. Same
/// stop-on-first-error behavior as [`XyzReader`].
pub struct ExtxyzReader<'a> {
    remaining: &'a str,
    limits: XyzParseLimits,
    input_validated: bool,
    frames_read: usize,
}

impl<'a> ExtxyzReader<'a> {
    /// Create a new `ExtxyzReader` over the given multi-frame extxyz string.
    pub fn new(input: &'a str) -> Self {
        Self::with_limits(input, XyzParseLimits::default())
    }

    /// Create an in-memory Extended XYZ reader with explicit resource limits.
    pub fn with_limits(input: &'a str, limits: XyzParseLimits) -> Self {
        Self {
            remaining: input,
            limits,
            input_validated: false,
            frames_read: 0,
        }
    }
}

impl<'a> Iterator for ExtxyzReader<'a> {
    type Item = Result<XyzFrame, XyzError>;

    fn next(&mut self) -> Option<Self::Item> {
        if !self.input_validated {
            self.input_validated = true;
            if let Err(error) = validate_xyz_input_limits(self.remaining, &self.limits) {
                self.remaining = "";
                return Some(Err(error));
            }
        }
        if self.remaining.trim().is_empty() {
            return None;
        }
        let leading_blank =
            self.remaining.len() - self.remaining.trim_start_matches(['\n', '\r']).len();
        self.remaining = &self.remaining[leading_blank..];

        let count_line = next_line(self.remaining)
            .map(|(line, _)| line)
            .unwrap_or("");
        if let Ok(declared) = count_line.trim().parse::<usize>()
            && declared > self.limits.max_atoms_per_frame
        {
            self.remaining = "";
            return Some(Err(XyzError::ResourceLimit {
                resource: "atoms per frame",
                actual: declared,
                limit: self.limits.max_atoms_per_frame,
            }));
        }
        self.frames_read = self.frames_read.saturating_add(1);
        if self.frames_read > self.limits.max_frames {
            self.remaining = "";
            return Some(Err(XyzError::ResourceLimit {
                resource: "frames",
                actual: self.frames_read,
                limit: self.limits.max_frames,
            }));
        }

        match parse_one_frame_ext(self.remaining) {
            Ok((frame, rest)) => {
                self.remaining = rest;
                Some(Ok(frame))
            }
            Err(e) => {
                self.remaining = "";
                Some(Err(e))
            }
        }
    }
}

/// Parse every frame in a multi-frame extxyz string. Stops and returns an
/// error on the first parse failure.
pub fn parse_extxyz_all(input: &str) -> Result<Vec<XyzFrame>, XyzError> {
    ExtxyzReader::new(input).collect()
}

fn format_xyz_value(v: &XyzValue) -> String {
    match v {
        XyzValue::Real(x) => x.to_string(),
        XyzValue::Integer(i) => i.to_string(),
        XyzValue::Str(s) => s.clone(),
        XyzValue::Logical(b) => (if *b { "T" } else { "F" }).to_string(),
    }
}

/// Build the extxyz comment line for `frame`. When `lattice`, `properties`
/// and `info` are all empty, the original `comment` is used verbatim --
/// this is what keeps a plain XYZ frame's comment line byte-identical
/// through the extxyz writer (atom rows are always re-emitted through the
/// generic species/pos formatter, so *those* are structurally, not
/// byte-for-byte, unchanged -- same non-goal as [`write_xyz`] already has
/// relative to hand-written input). Otherwise the comment line is rebuilt
/// canonically as `Lattice=... Properties=... key=value ...` (in that
/// order), which is not guaranteed to match a third-party file's original
/// key ordering.
fn build_comment_line(frame: &XyzFrame) -> String {
    if frame.lattice.is_none() && frame.properties.is_empty() && frame.info.is_empty() {
        return frame.comment.clone();
    }
    let mut parts = Vec::new();
    if let Some(l) = frame.lattice {
        let nums = l
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(" ");
        parts.push(format!("Lattice=\"{nums}\""));
    }
    let mut props = String::from("species:S:1:pos:R:3");
    for p in &frame.properties {
        let kind_char = match p.kind {
            XyzPropertyKind::Real => 'R',
            XyzPropertyKind::Integer => 'I',
            XyzPropertyKind::String => 'S',
            XyzPropertyKind::Logical => 'L',
        };
        props.push_str(&format!(":{}:{}:{}", p.name, kind_char, p.count));
    }
    parts.push(format!("Properties={props}"));
    for (k, v) in &frame.info {
        if v.is_empty() || v.chars().any(char::is_whitespace) || v.contains('"') {
            // Escape '"' and '\' so the quoted value reparses back to
            // exactly `v` -- see `tokenize_info_line`'s matching escape
            // handling, both modeled on ASE's own `key_val_str_to_dict`.
            let escaped: String = v
                .chars()
                .flat_map(|c| match c {
                    '"' => vec!['\\', '"'],
                    '\\' => vec!['\\', '\\'],
                    other => vec![other],
                })
                .collect();
            parts.push(format!("{k}=\"{escaped}\""));
        } else {
            parts.push(format!("{k}={v}"));
        }
    }
    parts.join(" ")
}

/// Rejects an [`XyzFrame`] that [`write_extxyz`]/[`ExtxyzWriter::write_frame`]
/// cannot losslessly round-trip through extxyz `key=value` comment-line
/// syntax. Returns `Err(detail)` naming the first offending field. Checks:
/// - an [`XyzFrame::info`] key or a [`XyzProperty::name`] containing
///   whitespace, `:`, `=`, or `"` (all syntactically significant to the
///   tokenizer/column-spec grammar, and -- unlike values -- names/keys have
///   no quoting or escaping to fall back on);
/// - an info key literally `"Lattice"` or `"Properties"`, which
///   [`parse_one_frame_ext`] always special-cases regardless of whether it
///   came from [`XyzFrame::lattice`]/[`properties`](XyzFrame::properties) or
///   a generic info entry, so a generic entry with either name would
///   reparse as the *wrong* field (or fail as `InvalidLattice`/
///   `InvalidProperties` if its value isn't in that field's own syntax);
/// - an info value containing an embedded newline, which no amount of
///   quoting/escaping can represent inside a single comment line;
/// - a [`XyzProperty`] declaring `count == 0`, which `Properties=`' own
///   grammar (parsed by [`parse_properties`]) requires to be at least 1.
///
/// Info values containing `"` or `\` are otherwise fine: `build_comment_line`
/// escapes them the same way [`tokenize_info_line`] un-escapes them on read
/// (both modeled on ASE's own `key_val_str_to_dict`), so they are not
/// rejected here.
fn validate_extxyz_writable(frame: &XyzFrame) -> Result<(), String> {
    for (k, v) in &frame.info {
        if k.is_empty() || k.chars().any(|c| c.is_whitespace() || c == '=' || c == '"') {
            return Err(format!(
                "info key {k:?} cannot be written (must be non-empty and contain no whitespace, '=', or '\"')"
            ));
        }
        if k == "Lattice" || k == "Properties" {
            return Err(format!(
                "info key {k:?} is reserved by extxyz -- use `XyzFrame::lattice` or \
                 `XyzFrame::properties` instead of a generic info entry"
            ));
        }
        if v.contains('\n') || v.contains('\r') {
            return Err(format!(
                "info['{k}'] value {v:?} cannot be written (extxyz comment lines cannot \
                 contain an embedded newline)"
            ));
        }
    }
    for p in &frame.properties {
        if p.name.is_empty()
            || p.name
                .chars()
                .any(|c| c.is_whitespace() || c == ':' || c == '=' || c == '"')
        {
            return Err(format!(
                "property name {:?} cannot be written (must be non-empty and contain no whitespace, ':', '=', or '\"')",
                p.name
            ));
        }
        if p.count == 0 {
            return Err(format!(
                "property '{}' declares 0 components, which extxyz's Properties= grammar \
                 cannot represent (minimum is 1)",
                p.name
            ));
        }
    }
    Ok(())
}

/// Streaming writer for a (possibly multi-frame) extxyz file.
pub struct ExtxyzWriter<W: std::io::Write> {
    writer: W,
}

impl<W: std::io::Write> ExtxyzWriter<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Write one frame. Call repeatedly for a multi-frame trajectory.
    ///
    /// Returns an [`std::io::Error`] of kind [`std::io::ErrorKind::InvalidData`]
    /// (not a sink failure) if `frame.info`/`frame.properties` contains a
    /// key/value/name that cannot be represented in extxyz syntax -- see
    /// [`validate_extxyz_writable`]. [`write_extxyz`] surfaces the same
    /// check as a typed [`XyzError`] instead.
    ///
    /// # Panics
    /// Panics if `frame.properties` is inconsistent with `frame.atoms`
    /// (a `values` row missing for some atom, or a row whose length
    /// doesn't match its column's declared `count`) -- this is a caller
    /// contract on hand-built `XyzFrame` values, not something that can
    /// happen from any of this module's own parsers.
    pub fn write_frame(&mut self, frame: &XyzFrame) -> std::io::Result<()> {
        if let Err(detail) = validate_extxyz_writable(frame) {
            return Err(std::io::Error::new(std::io::ErrorKind::InvalidData, detail));
        }
        writeln!(self.writer, "{}", frame.atoms.len())?;
        writeln!(self.writer, "{}", build_comment_line(frame))?;
        for (i, atom) in frame.atoms.iter().enumerate() {
            write!(
                self.writer,
                "{} {} {} {}",
                atom.element.symbol(),
                atom.x,
                atom.y,
                atom.z
            )?;
            for prop in &frame.properties {
                let values = &prop.values[i];
                debug_assert_eq!(
                    values.len(),
                    prop.count,
                    "XyzProperty '{}' declares count {} but atom {i} has {} value(s)",
                    prop.name,
                    prop.count,
                    values.len()
                );
                for v in values {
                    write!(self.writer, " {}", format_xyz_value(v))?;
                }
            }
            writeln!(self.writer)?;
        }
        Ok(())
    }
}

/// Write a single extxyz frame to a string.
///
/// Returns [`XyzError::UnwritableMetadata`] if `frame.info`/`frame.properties`
/// contains a key/value/name that cannot be represented in extxyz syntax --
/// see [`validate_extxyz_writable`]. A frame produced by [`parse_extxyz`]
/// can never trigger this (the parser's own tokenizer already rejects those
/// characters on read); it is reachable only from a hand-built `XyzFrame`,
/// e.g. via the Python/WASM `to_extxyz` bindings.
pub fn write_extxyz(frame: &XyzFrame) -> Result<String, XyzError> {
    let mut buf = Vec::new();
    ExtxyzWriter::new(&mut buf)
        .write_frame(frame)
        .map_err(|e| match e.kind() {
            std::io::ErrorKind::InvalidData => XyzError::UnwritableMetadata {
                detail: e.to_string(),
            },
            _ => unreachable!("writing to a Vec<u8> only fails with InvalidData (see write_frame)"),
        })?;
    Ok(String::from_utf8(buf).expect("ExtxyzWriter only emits ASCII/UTF-8 text"))
}
