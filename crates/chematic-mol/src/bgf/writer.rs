use super::*;
use std::fmt::Write;

pub fn write_bgf(record: &BgfRecord) -> Result<String, BgfError> {
    record.validate()?;
    let mut out = String::from("BIOGRF 200\n");
    if let Some(v) = &record.description {
        writeln!(out, "DESCRP {v}").unwrap();
    }
    for v in &record.remarks {
        writeln!(out, "REMARK {v}").unwrap();
    }
    if let Some(v) = &record.force_field {
        writeln!(out, "FORCEFIELD {v}").unwrap();
    }
    writeln!(out, "FORMAT ATOM   {BGF_ATOM_FORMAT}").unwrap();
    for a in &record.atoms {
        let kind = match a.kind {
            BgfAtomKind::Atom => "ATOM  ",
            BgfAtomKind::Hetatm => "HETATM",
        };
        writeln!(
            out,
            "{kind} {:>5} {:<5} {:>3} {:1} {:>5}{}{}{} {:<5}{:>3}{:>2} {}",
            a.id,
            a.name,
            a.residue_name,
            a.chain,
            a.residue_id,
            fixed(a.position_angstrom[0], 10)?,
            fixed(a.position_angstrom[1], 10)?,
            fixed(a.position_angstrom[2], 10)?,
            a.force_field_type,
            a.max_bonds,
            a.lone_pairs,
            fixed(a.partial_charge, 8)?
        )
        .unwrap();
    }
    out.push_str("FORMAT CONECT (a6,12i6)\n");
    for row in &record.connectivity {
        write!(out, "CONECT{:>6}", row.atom).unwrap();
        for n in &row.neighbors {
            write!(out, "{n:>6}").unwrap();
        }
        out.push('\n');
        if let Some(orders) = &row.orders {
            write!(out, "ORDER {:>6}", row.atom).unwrap();
            for n in orders {
                write!(out, "{n:>6}").unwrap();
            }
            out.push('\n');
        }
    }
    out.push_str("END\n");
    if out.len() > BgfLimits::default().max_input_bytes {
        return Err(error(0, "output exceeds byte limit"));
    }
    Ok(out)
}
