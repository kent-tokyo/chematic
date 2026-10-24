use super::*;
use std::collections::BTreeMap;
type Rows<'a> = Vec<(usize, &'a str)>;
struct Section<'a> {
    argument: &'a str,
    rows: Rows<'a>,
}

pub fn parse_molden(text: &str) -> Result<MoldenRecord, MoldenError> {
    parse_molden_with_limits(text, &MoldenLimits::default())
}

pub fn parse_molden_with_limits(
    text: &str,
    limits: &MoldenLimits,
) -> Result<MoldenRecord, MoldenError> {
    let sections = sections(text, limits)?;
    let atoms = sections
        .get("ATOMS")
        .map(|s| parse_atoms(s, limits))
        .transpose()?;
    let fr_coords = sections
        .get("FR-COORD")
        .map(|s| parse_fr_coords(s, limits))
        .transpose()?;
    let frequencies = sections
        .get("FREQ")
        .map(|s| parse_frequencies(s, limits))
        .transpose()?;
    let normal_modes = sections
        .get("FR-NORM-COORD")
        .map(|s| parse_modes(s, fr_coords.as_ref().map_or(0, Vec::len), limits))
        .transpose()?;
    let intensities = sections
        .get("INT")
        .map(|s| parse_intensities(s, limits))
        .transpose()?;
    let record = MoldenRecord {
        atoms,
        fr_coords,
        frequencies,
        normal_modes,
        intensities,
    };
    record.validate_with_limits(limits)?;
    Ok(record)
}

fn sections<'a>(
    text: &'a str,
    limits: &MoldenLimits,
) -> Result<BTreeMap<String, Section<'a>>, MoldenError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input exceeds byte limit"));
    }
    let max_rows = limits
        .max_displacements
        .saturating_add(limits.max_atoms.saturating_mul(2))
        .saturating_add(limits.max_modes.saturating_mul(3))
        .saturating_add(32);
    let mut result = BTreeMap::<String, Section>::new();
    let mut current = None;
    let mut magic = false;
    let mut row_count = 0usize;
    for (i, raw) in text.lines().enumerate() {
        if raw.len() > limits.max_line_len || raw.chars().any(|c| c.is_control() && c != '\t') {
            return Err(error(
                i + 1,
                "line exceeds limit or contains control characters",
            ));
        }
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if !magic {
            if !line.eq_ignore_ascii_case("[Molden Format]") {
                return Err(error(i + 1, "expected [Molden Format] first"));
            }
            magic = true;
            continue;
        }
        if line.starts_with('[') {
            let (key, arg) = section_header(line, i + 1)?;
            if result
                .insert(
                    key.clone(),
                    Section {
                        argument: arg,
                        rows: Vec::new(),
                    },
                )
                .is_some()
            {
                return Err(error(i + 1, "duplicate section"));
            }
            current = Some(key);
        } else {
            row_count += 1;
            if row_count > max_rows {
                return Err(error(i + 1, "record count exceeds limits"));
            }
            let key = current
                .as_ref()
                .ok_or_else(|| error(i + 1, "data outside a supported section"))?;
            result
                .get_mut(key)
                .expect("current section")
                .rows
                .push((i + 1, line));
        }
    }
    if !magic {
        return Err(error(1, "missing [Molden Format]"));
    }
    if result.values().any(|s| s.rows.is_empty()) {
        return Err(error(0, "empty sections are unsupported"));
    }
    Ok(result)
}
fn section_header(line: &str, line_no: usize) -> Result<(String, &str), MoldenError> {
    let end = line
        .find(']')
        .ok_or_else(|| error(line_no, "unclosed section header"))?;
    let key = line[1..end].to_ascii_uppercase();
    let arg = line[end + 1..].trim();
    if !matches!(
        key.as_str(),
        "ATOMS" | "FREQ" | "FR-COORD" | "FR-NORM-COORD" | "INT"
    ) {
        return Err(error(
            line_no,
            "unsupported Molden section; electronic, convergence and trajectory data cannot be discarded",
        ));
    }
    if key != "ATOMS" && !arg.is_empty() {
        return Err(error(line_no, "unexpected section argument"));
    }
    Ok((key, arg))
}
fn number(token: &str, line: usize) -> Result<f64, MoldenError> {
    let v: f64 = token
        .replace(['D', 'd'], "e")
        .parse()
        .map_err(|_| error(line, "invalid number"))?;
    if !v.is_finite() {
        return Err(error(line, "non-finite number"));
    }
    Ok(v)
}
fn vector(fields: &[&str], line: usize) -> Result<[f64; 3], MoldenError> {
    if fields.len() != 3 {
        return Err(error(line, "expected exactly three vector components"));
    }
    Ok([
        number(fields[0], line)?,
        number(fields[1], line)?,
        number(fields[2], line)?,
    ])
}
fn parse_atoms(s: &Section<'_>, limits: &MoldenLimits) -> Result<MoldenAtoms, MoldenError> {
    let unit = if s.argument.eq_ignore_ascii_case("Angs") {
        MoldenLengthUnit::Angstrom
    } else if s.argument.eq_ignore_ascii_case("AU") {
        MoldenLengthUnit::Bohr
    } else {
        return Err(error(0, "Atoms requires explicit Angs or AU unit"));
    };
    if s.rows.len() > limits.max_atoms {
        return Err(error(0, "Atoms count exceeds limit"));
    }
    let mut atoms = Vec::new();
    for &(line, row) in &s.rows {
        let f: Vec<_> = row.split_whitespace().collect();
        if f.len() != 6 {
            return Err(error(
                line,
                "Atoms row needs symbol, index, atomic number, x, y, z",
            ));
        }
        atoms.push(MoldenAtom {
            symbol: f[0].into(),
            index: f[1]
                .parse()
                .map_err(|_| error(line, "invalid atom index"))?,
            atomic_number: f[2]
                .parse()
                .map_err(|_| error(line, "invalid atomic number"))?,
            position: vector(&f[3..], line)?,
        });
    }
    Ok(MoldenAtoms { unit, atoms })
}
fn parse_fr_coords(
    s: &Section<'_>,
    limits: &MoldenLimits,
) -> Result<Vec<MoldenVibrationAtom>, MoldenError> {
    if s.rows.len() > limits.max_atoms {
        return Err(error(0, "FR-COORD atom count exceeds limit"));
    }
    let mut atoms = Vec::new();
    for &(line, row) in &s.rows {
        let f: Vec<_> = row.split_whitespace().collect();
        if f.len() != 4 {
            return Err(error(line, "FR-COORD row needs symbol, x, y, z"));
        }
        atoms.push(MoldenVibrationAtom {
            symbol: f[0].into(),
            position_bohr: vector(&f[1..], line)?,
        });
    }
    Ok(atoms)
}
fn parse_frequencies(s: &Section<'_>, limits: &MoldenLimits) -> Result<Vec<f64>, MoldenError> {
    if s.rows.len() > limits.max_modes {
        return Err(error(0, "frequency count exceeds limit"));
    }
    s.rows
        .iter()
        .map(|&(line, row)| {
            if row.split_whitespace().count() != 1 {
                return Err(error(line, "frequency row needs exactly one value"));
            }
            number(row, line)
        })
        .collect()
}
fn parse_intensities(
    s: &Section<'_>,
    limits: &MoldenLimits,
) -> Result<Vec<MoldenIntensity>, MoldenError> {
    if s.rows.len() > limits.max_modes {
        return Err(error(0, "intensity count exceeds limit"));
    }
    s.rows
        .iter()
        .map(|&(line, row)| {
            let f: Vec<_> = row.split_whitespace().collect();
            if !(1..=2).contains(&f.len()) {
                return Err(error(line, "INT row needs IR and optional Raman value"));
            }
            Ok(MoldenIntensity {
                ir: number(f[0], line)?,
                raman: f.get(1).map(|v| number(v, line)).transpose()?,
            })
        })
        .collect()
}
fn parse_modes(
    s: &Section<'_>,
    atoms: usize,
    limits: &MoldenLimits,
) -> Result<Vec<MoldenNormalMode>, MoldenError> {
    if atoms == 0 {
        return Err(error(0, "normal modes require FR-COORD"));
    }
    let mut modes = Vec::new();
    let mut rows = s.rows.iter();
    let mut total = 0usize;
    while let Some(&(line, row)) = rows.next() {
        let f: Vec<_> = row.split_whitespace().collect();
        if f.len() != 2 || !f[0].eq_ignore_ascii_case("vibration") {
            return Err(error(line, "expected vibration and mode index"));
        }
        total = total
            .checked_add(atoms)
            .ok_or_else(|| error(line, "displacement count overflow"))?;
        if modes.len() >= limits.max_modes || total > limits.max_displacements {
            return Err(error(line, "mode/displacement count exceeds limits"));
        }
        let index = f[1]
            .parse()
            .map_err(|_| error(line, "invalid mode index"))?;
        let mut displacements_bohr = Vec::new();
        for _ in 0..atoms {
            let &(n, row) = rows
                .next()
                .ok_or_else(|| error(line, "truncated normal-mode block"))?;
            displacements_bohr.push(vector(&row.split_whitespace().collect::<Vec<_>>(), n)?);
        }
        modes.push(MoldenNormalMode {
            index,
            displacements_bohr,
        });
    }
    Ok(modes)
}
