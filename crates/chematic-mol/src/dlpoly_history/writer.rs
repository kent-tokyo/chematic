use super::*;
use std::fmt::Write;

pub fn write_dlpoly_history(g: &DlpolyHistory) -> Result<String, DlpolyHistoryError> {
    g.validate()?;
    let mut out = format!(
        "{}\n{:>10}{:>10}{:>10}\n",
        g.title, g.keytrj, g.imcon, g.atom_count
    );
    for f in &g.frames {
        writeln!(
            out,
            "timestep{:>10}{:>10}{:>10}{:>10}{}",
            f.step,
            g.atom_count,
            g.keytrj,
            g.imcon,
            real(f.timestep)?
        )
        .unwrap();
        if let Some(cell) = &f.cell_vectors {
            for v in cell {
                vector_row(&mut out, v)?;
            }
        }
        for a in &f.atoms {
            writeln!(
                out,
                "{:<8}{:>10}{}{}",
                a.name,
                a.index,
                real(a.mass)?,
                real(a.charge)?
            )
            .unwrap();
            vector_row(&mut out, &a.position)?;
            if let Some(v) = &a.velocity {
                vector_row(&mut out, v)?;
            }
            if let Some(v) = &a.force {
                vector_row(&mut out, v)?;
            }
            if out.len() > DlpolyHistoryLimits::default().max_input_bytes {
                return Err(error(0, "output byte limit"));
            }
        }
    }
    Ok(out)
}
fn vector_row(out: &mut String, v: &[f64; 3]) -> Result<(), DlpolyHistoryError> {
    writeln!(out, "{}{}{}", real(v[0])?, real(v[1])?, real(v[2])?).unwrap();
    Ok(())
}
