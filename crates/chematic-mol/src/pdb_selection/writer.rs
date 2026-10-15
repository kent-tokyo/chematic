use super::*;
use std::collections::HashSet;

pub(super) fn select(
    text: &str,
    doc: &Document<'_>,
    keep: &[bool],
) -> Result<PdbSelectionReport, PdbSelectionError> {
    let mut report = PdbSelectionReport {
        pdb: String::new(),
        source_text: String::new(),
        source_atom_count: keep.len(),
        selected_source_indices: Vec::new(),
        old_to_new: Vec::new(),
        atoms: Vec::new(),
        losses: doc.losses.clone(),
        retargeted_ter_lines: Vec::new(),
        filtered_conect_lines: Vec::new(),
        model_num: doc.model_num,
    };
    reserve(&mut report.old_to_new, keep.len())?;
    report.old_to_new.resize(keep.len(), None);
    for (i, &selected) in keep.iter().enumerate() {
        if selected {
            report.old_to_new[i] = Some(report.atoms.len());
            reserve(&mut report.atoms, 1)?;
            report.atoms.push(identity(doc, i));
            reserve(&mut report.selected_source_indices, 1)?;
            report.selected_source_indices.push(i);
        }
    }
    let mut last = None;
    let mut ter_serials = HashSet::new();
    for row in &doc.rows {
        match &row.kind {
            Kind::Atom(i) => {
                if keep[*i] {
                    append(&mut report.pdb, row.text)?;
                    last = Some(*i);
                }
            }
            Kind::Anisou(i) => {
                if keep[*i] {
                    append(&mut report.pdb, row.text)?;
                }
            }
            Kind::Ter(_) => {
                if let Some(i) = last {
                    write_ter(&mut report, row, doc, i, &mut ter_serials)?;
                }
                last = None;
            }
            Kind::Conect(values) => write_conect(&mut report, row, values, doc, keep)?,
            Kind::Model | Kind::EndModel => {
                if !report.atoms.is_empty() {
                    append(&mut report.pdb, row.text)?;
                }
                last = None;
            }
            Kind::End => {}
        }
    }
    append(&mut report.pdb, "END")?;
    report.source_text = text.to_string();
    Ok(report)
}
fn append(out: &mut String, line: &str) -> Result<(), PdbSelectionError> {
    let n = out
        .len()
        .checked_add(line.len())
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| error(0, "output size overflow"))?;
    if n > PdbSelectionLimits::default().max_input_bytes {
        return Err(error(0, "output byte limit exceeded"));
    }
    out.try_reserve(line.len() + 1)
        .map_err(|e| error(0, format!("output allocation failed: {e}")))?;
    out.push_str(line);
    out.push('\n');
    Ok(())
}
fn identity(doc: &Document<'_>, i: usize) -> PdbSelectedAtom {
    let a = &doc.atoms[i];
    let line = doc.atom_lines[i];
    PdbSelectedAtom {
        source_index: i,
        serial: a.serial,
        record_type: line[..6].trim_end().into(),
        atom_name: a.name.clone(),
        alt_loc: line.as_bytes()[16] as char,
        residue_name: a.res_name.clone(),
        chain_id: a.chain_id,
        residue_number: a.res_seq,
        insertion_code: line.as_bytes()[26] as char,
        position_angstrom: [a.x, a.y, a.z],
    }
}
fn write_ter(
    report: &mut PdbSelectionReport,
    row: &Row<'_>,
    doc: &Document<'_>,
    index: usize,
    used: &mut HashSet<u32>,
) -> Result<(), PdbSelectionError> {
    let atom = doc.atom_lines[index];
    let serial = doc.atoms[index]
        .serial
        .checked_add(1)
        .filter(|v| *v <= 99_999)
        .ok_or_else(|| error(row.line, "selected TER serial overflows PDB width"))?;
    if doc
        .by_serial
        .get(&serial)
        .is_some_and(|&i| report.old_to_new[i].is_some())
    {
        return Err(error(
            row.line,
            "selected TER serial collides with a retained atom",
        ));
    }
    used.try_reserve(1)
        .map_err(|e| error(row.line, format!("TER allocation failed: {e}")))?;
    if !used.insert(serial) {
        return Err(error(row.line, "selected TER serial collision"));
    }
    let original = match row.kind {
        Kind::Ter(v) => v,
        _ => return Err(error(row.line, "expected TER record")),
    };
    if original == serial && row.text.get(17..27) == atom.get(17..27) {
        return append(&mut report.pdb, row.text);
    }
    let line = format!(
        "TER   {serial:5}{}{}{}",
        &row.text[11..17],
        &atom[17..27],
        &row.text[27..]
    );
    reserve(&mut report.retargeted_ter_lines, 1)?;
    report.retargeted_ter_lines.push(row.line);
    append(&mut report.pdb, &line)
}
fn write_conect(
    report: &mut PdbSelectionReport,
    row: &Row<'_>,
    values: &[u32],
    doc: &Document<'_>,
    keep: &[bool],
) -> Result<(), PdbSelectionError> {
    // Reader verified every serial exists; use checked lookup to keep failure explicit.
    let selected = |serial: &u32| doc.by_serial.get(serial).is_some_and(|&i| keep[i]);
    let targets: Vec<_> = values[1..]
        .iter()
        .filter(|s| selected(s))
        .copied()
        .collect();
    if selected(&values[0]) && targets.len() == values.len() - 1 {
        return append(&mut report.pdb, row.text);
    }
    reserve(&mut report.filtered_conect_lines, 1)?;
    report.filtered_conect_lines.push(row.line);
    if !selected(&values[0]) || targets.is_empty() {
        return Ok(());
    }
    let mut line = format!("CONECT{:5}", values[0]);
    for target in targets {
        line.push_str(&format!("{target:5}"));
    }
    append(&mut report.pdb, &line)
}
