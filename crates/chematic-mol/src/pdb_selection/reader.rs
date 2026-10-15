use super::*;
use std::collections::HashSet;

struct State {
    model: bool,
    closed: bool,
    ended: bool,
    connect: bool,
    last: Option<usize>,
    anisou: HashSet<u32>,
}

pub(super) fn read<'a>(
    text: &'a str,
    allow_loss: bool,
    limits: &PdbSelectionLimits,
) -> Result<Document<'a>, PdbSelectionError> {
    check_limits(text, limits)?;
    let hard = PdbSelectionLimits::default();
    let atoms = chematic_3d::pdb::parse_pdb_atoms_strict(
        text,
        &chematic_3d::pdb::PdbParseLimits {
            max_input_bytes: limits.max_input_bytes.min(hard.max_input_bytes),
            max_line_bytes: limits.max_line_len.min(hard.max_line_len),
            max_atoms: limits.max_atoms.min(hard.max_atoms),
            max_models: 1,
        },
    )
    .map_err(|e| match e {
        chematic_3d::pdb::PdbStrictError::Parse(ref p) => error(p.line, e.to_string()),
        _ => error(0, e.to_string()),
    })?;
    let mut doc = Document {
        rows: Vec::new(),
        atom_lines: Vec::new(),
        atoms,
        by_serial: HashMap::new(),
        losses: Vec::new(),
        model_num: None,
    };
    let mut state = State {
        model: false,
        closed: false,
        ended: false,
        connect: false,
        last: None,
        anisou: HashSet::new(),
    };
    for (i, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let number = i + 1;
        if state.ended {
            return Err(error(number, "records after END are unsupported"));
        }
        let name = line.get(..6).unwrap_or(line).trim_end();
        let kind = match name {
            "ATOM" | "HETATM" => Kind::Atom(atom(&mut doc, &mut state, line, number)?),
            "ANISOU" => Kind::Anisou(anisou(&doc, &mut state, line, number)?),
            "TER" => Kind::Ter(ter(&doc, &mut state, line, number)?),
            "CONECT" => {
                if state.model {
                    return Err(error(number, "CONECT must follow ENDMDL"));
                }
                state.connect = true;
                Kind::Conect(conect(line, number)?)
            }
            "MODEL" => {
                model(&mut doc, &mut state, line, number)?;
                Kind::Model
            }
            "ENDMDL" => {
                if !state.model {
                    return Err(error(number, "unmatched ENDMDL"));
                }
                state.model = false;
                state.closed = true;
                state.last = None;
                Kind::EndModel
            }
            "END" => {
                if state.model {
                    return Err(error(number, "MODEL lacks ENDMDL"));
                }
                state.ended = true;
                Kind::End
            }
            _ => {
                if !allow_loss {
                    return Err(error(
                        number,
                        format!(
                            "record {name:?} outside coordinate selection contract; metadata-loss consent required"
                        ),
                    ));
                }
                reserve(&mut doc.losses, 1)?;
                doc.losses.push(PdbSelectionLoss {
                    line: number,
                    record_name: name.to_string(),
                });
                continue;
            }
        };
        reserve(&mut doc.rows, 1)?;
        doc.rows.push(Row {
            line: number,
            text: line,
            kind,
        });
    }
    if state.model {
        return Err(error(0, "MODEL lacks ENDMDL"));
    }
    check_references(&doc)?;
    Ok(doc)
}

fn check_limits(text: &str, limits: &PdbSelectionLimits) -> Result<(), PdbSelectionError> {
    let hard = PdbSelectionLimits::default();
    if text.len() > limits.max_input_bytes.min(hard.max_input_bytes) {
        return Err(error(0, "input byte limit exceeded"));
    }
    for (i, line) in text.lines().enumerate() {
        if i >= limits.max_records.min(hard.max_records) {
            return Err(error(i + 1, "record count limit exceeded"));
        }
        if line.len() > limits.max_line_len.min(hard.max_line_len) {
            return Err(error(i + 1, "line length limit exceeded"));
        }
        if line.bytes().any(|b| !(b' '..=b'~').contains(&b)) {
            return Err(error(
                i + 1,
                "PDB records must use printable ASCII and spaces",
            ));
        }
    }
    Ok(())
}
fn serial(line: &str, number: usize) -> Result<u32, PdbSelectionError> {
    line.get(6..11)
        .and_then(|s| s.trim().parse::<u32>().ok())
        .filter(|v| *v > 0 && *v <= 99_999)
        .ok_or_else(|| error(number, "invalid positive decimal serial"))
}
fn coordinate_phase(state: &State, number: usize) -> Result<(), PdbSelectionError> {
    if state.closed || state.connect {
        return Err(error(
            number,
            "coordinate record outside single-model coordinate block",
        ));
    }
    Ok(())
}
fn atom<'a>(
    doc: &mut Document<'a>,
    state: &mut State,
    line: &'a str,
    number: usize,
) -> Result<usize, PdbSelectionError> {
    coordinate_phase(state, number)?;
    if line.len() > 80 && !line[80..].trim().is_empty() {
        return Err(error(number, "extra atom fields after column 80"));
    }
    if line.get(12..16).is_none_or(|s| s.trim().is_empty())
        || line.get(17..20).is_none_or(|s| s.trim().is_empty())
    {
        return Err(error(number, "atom and residue names are required"));
    }
    for (start, end) in [(54, 60), (60, 66)] {
        if let Some(s) = line.get(start..end) {
            if !s.trim().is_empty() && !s.trim().parse::<f64>().is_ok_and(f64::is_finite) {
                return Err(error(number, "invalid occupancy or temperature factor"));
            }
        } else if line.get(start..).is_some_and(|s| !s.trim().is_empty()) {
            return Err(error(
                number,
                "incomplete occupancy or temperature factor field",
            ));
        }
    }
    let serial = serial(line, number)?;
    let index = doc.atom_lines.len();
    doc.by_serial
        .try_reserve(1)
        .map_err(|e| error(number, format!("serial allocation failed: {e}")))?;
    if doc.by_serial.insert(serial, index).is_some() {
        return Err(error(number, "duplicate atom serial"));
    }
    reserve(&mut doc.atom_lines, 1)?;
    doc.atom_lines.push(line);
    state.last = Some(index);
    Ok(index)
}
fn anisou(
    doc: &Document<'_>,
    state: &mut State,
    line: &str,
    number: usize,
) -> Result<usize, PdbSelectionError> {
    coordinate_phase(state, number)?;
    let serial = serial(line, number)?;
    let index = *doc
        .by_serial
        .get(&serial)
        .ok_or_else(|| error(number, "ANISOU must follow its atom record"))?;
    if line.get(6..27) != doc.atom_lines[index].get(6..27) {
        return Err(error(number, "ANISOU identity disagrees with atom"));
    }
    for start in (28..70).step_by(7) {
        if line
            .get(start..start + 7)
            .and_then(|s| s.trim().parse::<i32>().ok())
            .is_none()
        {
            return Err(error(number, "invalid ANISOU tensor field"));
        }
    }
    if line.len() > 80 && !line[80..].trim().is_empty() {
        return Err(error(number, "extra ANISOU fields"));
    }
    state
        .anisou
        .try_reserve(1)
        .map_err(|e| error(number, format!("ANISOU allocation failed: {e}")))?;
    if !state.anisou.insert(serial) {
        return Err(error(number, "duplicate ANISOU record"));
    }
    Ok(index)
}
fn ter(
    doc: &Document<'_>,
    state: &mut State,
    line: &str,
    number: usize,
) -> Result<u32, PdbSelectionError> {
    coordinate_phase(state, number)?;
    let serial = serial(line, number)?;
    let index = state
        .last
        .ok_or_else(|| error(number, "TER has no preceding segment"))?;
    if line.get(17..27).is_none() || line.get(17..27) != doc.atom_lines[index].get(17..27) {
        return Err(error(number, "TER residue disagrees with preceding atom"));
    }
    if doc.atoms[index].serial.checked_add(1) != Some(serial) {
        return Err(error(
            number,
            "TER serial must follow the terminal atom serial",
        ));
    }
    state.last = None;
    Ok(serial)
}
fn model(
    doc: &mut Document<'_>,
    state: &mut State,
    line: &str,
    number: usize,
) -> Result<(), PdbSelectionError> {
    if doc.model_num.is_some() || !doc.atom_lines.is_empty() || state.connect {
        return Err(error(
            number,
            "only one complete explicit model or one implicit model is supported",
        ));
    }
    let id = line
        .get(10..14)
        .and_then(|s| s.trim().parse::<u32>().ok())
        .filter(|v| *v > 0)
        .ok_or_else(|| error(number, "invalid MODEL number"))?;
    doc.model_num = Some(id);
    state.model = true;
    Ok(())
}
fn conect(line: &str, number: usize) -> Result<Vec<u32>, PdbSelectionError> {
    if line.len() > 31 && !line[31..].trim().is_empty() {
        return Err(error(
            number,
            "legacy CONECT fields beyond four covalent neighbors are unsupported",
        ));
    }
    let mut values = Vec::new();
    for start in (6..31).step_by(5) {
        let field = line.get(start..line.len().min(start + 5)).unwrap_or("");
        if field.trim().is_empty() {
            continue;
        }
        let value = field
            .trim()
            .parse::<u32>()
            .ok()
            .filter(|v| *v > 0 && *v <= 99_999)
            .ok_or_else(|| error(number, "invalid CONECT serial"))?;
        reserve(&mut values, 1)?;
        values.push(value);
    }
    if values.is_empty() {
        return Err(error(number, "CONECT source is required"));
    }
    // Source must occupy its fixed field even if later fields are populated.
    if serial(line, number)? != values[0] {
        return Err(error(number, "invalid CONECT source"));
    }
    Ok(values)
}
fn check_references(doc: &Document<'_>) -> Result<(), PdbSelectionError> {
    let mut ter_serials = HashSet::new();
    for row in &doc.rows {
        if let Kind::Ter(serial) = row.kind {
            ter_serials
                .try_reserve(1)
                .map_err(|e| error(row.line, format!("TER allocation failed: {e}")))?;
            if doc.by_serial.contains_key(&serial) || !ter_serials.insert(serial) {
                return Err(error(
                    row.line,
                    "TER serial collides with another coordinate record",
                ));
            }
        }
        if let Kind::Conect(values) = &row.kind {
            for &serial in values {
                if !doc.by_serial.contains_key(&serial) {
                    return Err(error(row.line, "CONECT refers to an unknown atom serial"));
                }
            }
            if values[1..].contains(&values[0]) {
                return Err(error(row.line, "self CONECT reference"));
            }
        }
    }
    Ok(())
}
