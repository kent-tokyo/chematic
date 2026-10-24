use super::*;
use std::fmt::Write;

/// Serialize current data. Source numbers/units/signs/IDs remain unchanged;
/// layout and section order are canonicalized, never modes or atom rows.
pub fn write_molden(record: &MoldenRecord) -> Result<String, MoldenError> {
    record.validate()?;
    let mut out = String::from("[Molden Format]\n");
    if let Some(atoms) = &record.atoms {
        let unit = match atoms.unit {
            MoldenLengthUnit::Angstrom => "Angs",
            MoldenLengthUnit::Bohr => "AU",
        };
        writeln!(out, "[Atoms] {unit}").expect("writing a String");
        for atom in &atoms.atoms {
            write!(
                out,
                "{} {} {} ",
                atom.symbol, atom.index, atom.atomic_number
            )
            .expect("writing a String");
            write_vector(&mut out, &atom.position);
        }
    }
    if let Some(frequencies) = &record.frequencies {
        out.push_str("[FREQ]\n");
        for n in frequencies {
            writeln!(out, "{n:.17e}").expect("writing a String");
        }
    }
    if let Some(coords) = &record.fr_coords {
        out.push_str("[FR-COORD]\n");
        for atom in coords {
            write!(out, "{} ", atom.symbol).expect("writing a String");
            write_vector(&mut out, &atom.position_bohr);
        }
    }
    if let Some(modes) = &record.normal_modes {
        out.push_str("[FR-NORM-COORD]\n");
        for mode in modes {
            writeln!(out, "vibration {}", mode.index).expect("writing a String");
            for v in &mode.displacements_bohr {
                write_vector(&mut out, v);
            }
            if out.len() > MoldenLimits::default().max_input_bytes {
                return Err(error(0, "output exceeds byte limit"));
            }
        }
    }
    if let Some(intensities) = &record.intensities {
        out.push_str("[INT]\n");
        for v in intensities {
            write!(out, "{:.17e}", v.ir).expect("writing a String");
            if let Some(r) = v.raman {
                write!(out, " {r:.17e}").expect("writing a String");
            }
            out.push('\n');
        }
    }
    if out.len() > MoldenLimits::default().max_input_bytes {
        return Err(error(0, "output exceeds byte limit"));
    }
    Ok(out)
}
fn write_vector(out: &mut String, v: &[f64; 3]) {
    writeln!(out, "{:.17e} {:.17e} {:.17e}", v[0], v[1], v[2]).expect("writing a String");
}
