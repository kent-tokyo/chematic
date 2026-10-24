use super::*;

pub fn parse_bgf(text: &str) -> Result<BgfRecord, BgfError> {
    parse_bgf_with_limits(text, &BgfLimits::default())
}
pub fn parse_bgf_with_limits(text: &str, limits: &BgfLimits) -> Result<BgfRecord, BgfError> {
    if text.len() > limits.max_input_bytes || !text.is_ascii() {
        return Err(error(0, "input exceeds byte limit or is not ASCII"));
    }
    let mut record = BgfRecord {
        version: 200,
        description: None,
        remarks: Vec::new(),
        force_field: None,
        atoms: Vec::new(),
        connectivity: Vec::new(),
    };
    let mut started = false;
    let mut ended = false;
    let mut state = ParseState::default();
    for (i, raw) in text.lines().enumerate() {
        let line = i + 1;
        if raw.len() > limits.max_line_len || raw.bytes().any(|b| b.is_ascii_control()) {
            return Err(error(line, "line limit or control character"));
        }
        if raw.trim().is_empty() {
            continue;
        }
        if ended {
            return Err(error(
                line,
                "content after END; only one structure is supported",
            ));
        }
        if !started {
            let fields: Vec<_> = raw.split_whitespace().collect();
            if fields != ["BIOGRF", "200"] {
                return Err(error(line, "expected nonperiodic BIOGRF 200 header"));
            }
            started = true;
            continue;
        }
        read_record(&mut record, &mut state, &mut ended, raw, line, limits)?;
    }
    if !started || !ended {
        return Err(error(0, "missing header or END"));
    }
    record.validate_with_limits(limits)?;
    Ok(record)
}
fn metadata(record: &mut BgfRecord, raw: &str, line: usize) -> Result<(), BgfError> {
    if let Some(v) = raw.strip_prefix("REMARK ") {
        record.remarks.push(v.to_owned());
        return Ok(());
    }
    let (target, value) = if let Some(v) = raw.strip_prefix("DESCRP ") {
        (&mut record.description, v)
    } else if let Some(v) = raw.strip_prefix("FORCEFIELD ") {
        (&mut record.force_field, v)
    } else {
        return Err(error(line, "invalid metadata record"));
    };
    if target.replace(value.to_owned()).is_some() {
        return Err(error(line, "duplicate metadata"));
    }
    Ok(())
}
fn parse_atom(raw: &str, line: usize) -> Result<BgfAtom, BgfError> {
    if raw.len() < 80
        || !raw[80..].trim().is_empty()
        || [6, 12, 18, 22, 24, 60, 71]
            .iter()
            .any(|&i| raw.as_bytes()[i] != b' ')
    {
        return Err(error(
            line,
            "expected standard 80-column ATOM record without extensions",
        ));
    }
    let f = |a, b| raw[a..b].trim().to_owned();
    let integer = |a, b| {
        raw[a..b]
            .trim()
            .parse::<u32>()
            .map_err(|_| error(line, "invalid atom integer"))
    };
    Ok(BgfAtom {
        kind: if &raw[..6] == "ATOM  " {
            BgfAtomKind::Atom
        } else if &raw[..6] == "HETATM" {
            BgfAtomKind::Hetatm
        } else {
            return Err(error(line, "invalid atom keyword alignment"));
        },
        id: integer(7, 12)?,
        name: f(13, 18),
        residue_name: f(19, 22),
        chain: f(23, 24),
        residue_id: f(25, 30),
        position_angstrom: [
            decimal(&raw[30..40], line)?,
            decimal(&raw[40..50], line)?,
            decimal(&raw[50..60], line)?,
        ],
        force_field_type: f(61, 66),
        max_bonds: u16::try_from(integer(66, 69)?)
            .map_err(|_| error(line, "max_bonds overflow"))?,
        lone_pairs: u8::try_from(integer(69, 71)?)
            .map_err(|_| error(line, "lone_pairs overflow"))?,
        partial_charge: decimal(&raw[72..80], line)?,
    })
}
fn decimal(s: &str, line: usize) -> Result<f64, BgfError> {
    let s = s.trim();
    // No implied decimals or exponent conventions are guessed.
    if !s.contains('.') || !s.bytes().all(|b| b.is_ascii_digit() || b"+-.".contains(&b)) {
        return Err(error(line, "expected explicit decimal point in F field"));
    }
    let n: f64 = s.parse().map_err(|_| error(line, "invalid decimal"))?;
    fixed(n, s.len().max(10)).map_err(|e| error(line, e.detail))?;
    Ok(n)
}
fn integers(raw: &str, line: usize) -> Result<Vec<u32>, BgfError> {
    let raw = raw.trim_end();
    if raw.len() < 12 || raw.len() > 78 || !(raw.len() - 6).is_multiple_of(6) {
        return Err(error(line, "expected (a6,12i6) connection/order fields"));
    }
    raw.as_bytes()[6..]
        .chunks(6)
        .map(|s| {
            std::str::from_utf8(s)
                .unwrap()
                .trim()
                .parse()
                .map_err(|_| error(line, "invalid connection/order integer"))
        })
        .collect()
}

#[derive(Default)]
struct ParseState {
    bonds: bool,
    atom_format: bool,
    conect_format: bool,
}
fn read_record(
    record: &mut BgfRecord,
    state: &mut ParseState,
    ended: &mut bool,
    raw: &str,
    line: usize,
    limits: &BgfLimits,
) -> Result<(), BgfError> {
    let key = raw.get(..6).unwrap_or(raw).trim();
    match key {
        "ATOM" | "HETATM" => {
            if state.bonds || record.atoms.len() >= limits.max_atoms {
                return Err(error(line, "atom after connectivity or atom count limit"));
            }
            record.atoms.push(parse_atom(raw, line)?);
        }
        "CONECT" => {
            state.bonds = true;
            if record.connectivity.len() >= limits.max_atoms {
                return Err(error(line, "connectivity count limit"));
            }
            let values = integers(raw, line)?;
            record.connectivity.push(BgfConnectivity {
                atom: values[0],
                neighbors: values[1..].to_vec(),
                orders: None,
            });
        }
        "ORDER" => {
            let values = integers(raw, line)?;
            let row = record
                .connectivity
                .last_mut()
                .ok_or_else(|| error(line, "ORDER without preceding CONECT"))?;
            if row.atom != values[0] || row.orders.is_some() {
                return Err(error(line, "ORDER must follow its matching CONECT"));
            }
            row.orders = Some(
                values[1..]
                    .iter()
                    .map(|&n| u8::try_from(n).map_err(|_| error(line, "ORDER overflow")))
                    .collect::<Result<_, _>>()?,
            );
        }
        "FORMAT" => {
            let normalized = raw.to_ascii_lowercase().replace(' ', "");
            if normalized == format!("formatatom{BGF_ATOM_FORMAT}")
                && !state.atom_format
                && record.atoms.is_empty()
                && !state.bonds
            {
                state.atom_format = true;
            } else if normalized == "formatconect(a6,12i6)" && !state.conect_format && !state.bonds
            {
                state.conect_format = true;
            } else {
                return Err(error(line, "unsupported, duplicate or misplaced FORMAT"));
            }
        }
        "END" if raw.trim() == "END" => *ended = true,
        "DESCRP" | "REMARK" | "FORCEF" if record.atoms.is_empty() && !state.bonds => {
            metadata(record, raw, line)?
        }
        _ => {
            return Err(error(
                line,
                "unsupported or misplaced record; no fields discarded",
            ));
        }
    }
    Ok(())
}
