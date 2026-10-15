use super::*;

/// Write current Cartesian data, preserving keyword/title/comment text, labels,
/// flags (including absence), and ordered translation vectors. No calculation
/// settings, formal charge, spin state, isotope or topology are inferred.
pub fn write_mopac_cartesian(g: &MopacCartesianGeometry) -> Result<String, MopacError> {
    g.validate()?;
    let mut out = format!("{}\n{}\n{}\n", g.keywords, g.title, g.comment);
    for atom in &g.atoms {
        out.push_str(atom.element.symbol());
        if let Some(label) = &atom.label {
            out.push('(');
            out.push_str(label);
            out.push(')');
        }
        write_row(&mut out, &atom.row);
        if out.len() > MopacParseLimits::default().max_input_bytes {
            return Err(error(0, "serialized byte limit exceeded"));
        }
    }
    for row in &g.translation_vectors {
        out.push_str("Tv");
        write_row(&mut out, row);
    }
    out.push('\n');
    if out.len() > MopacParseLimits::default().max_input_bytes {
        return Err(error(0, "serialized byte limit exceeded"));
    }
    Ok(out)
}
fn write_row(out: &mut String, row: &MopacCartesianRow) {
    for (i, value) in row.position_angstrom.iter().enumerate() {
        out.push_str(&format!(" {value:.17e}"));
        if let Some(flags) = row.optimization_flags {
            out.push(' ');
            out.push_str(flags[i].code());
        }
    }
    out.push('\n');
}
