use super::*;
use std::fmt::Write;

/// Canonical block, retaining f64 values and explicit-vs-omitted options.
/// Comments move to an ordered leading block. No calculation is executed.
pub fn write_nwchem_geometry(g: &NwchemGeometry) -> Result<String, NwchemError> {
    g.validate()?;
    let mut out = String::new();
    for comment in &g.comments {
        writeln!(out, "#{comment}").unwrap();
        check_size(&out)?;
    }
    write_header(&mut out, g);
    for c in &g.centers {
        write!(out, "  {}", c.tag).unwrap();
        vector(&mut out, c.position);
        if let Some(v) = c.velocity_au {
            vector(&mut out, v);
        }
        if let Some(v) = c.nuclear_charge {
            write!(out, " charge {v:.17e}").unwrap();
        }
        if let Some(v) = c.mass_amu {
            write!(out, " mass {v:.17e}").unwrap();
        }
        if let Some(v) = c.nucleus {
            write!(out, " nucleus {}", model(v)).unwrap();
        }
        out.push('\n');
        check_size(&out)?;
    }
    out.push_str("end\n");
    check_size(&out)?;
    Ok(out)
}
fn write_header(out: &mut String, g: &NwchemGeometry) {
    out.push_str("geometry");
    if let Some(v) = &g.name {
        write!(out, " {v}").unwrap();
    }
    if let Some(v) = &g.units {
        write!(out, " units {v}").unwrap();
    }
    if let Some(v) = g.angstrom_to_au {
        write!(out, " angstrom_to_au {v:.17e}").unwrap();
    }
    for (value, yes, no) in [
        (g.center, "center", "nocenter"),
        (g.autosym, "autosym", "noautosym"),
    ] {
        if let Some(v) = value {
            write!(out, " {}", if v { yes } else { no }).unwrap();
        }
    }
    if let Some(v) = g.symmetry_tolerance {
        write!(out, " {v:.17e}").unwrap();
    }
    if let Some(v) = g.autoz {
        write!(out, " {}", if v { "autoz" } else { "noautoz" }).unwrap();
    }
    if let Some(v) = g.print {
        out.push_str(match v {
            NwchemPrint::Normal => " print",
            NwchemPrint::Xyz => " print xyz",
            NwchemPrint::None => " noprint",
        });
    }
    if g.bqbq {
        out.push_str(" bqbq");
    }
    if let Some(v) = g.nucleus {
        write!(out, " nucleus {}", model(v)).unwrap();
    }
    out.push('\n');
}
fn model(m: NwchemNucleus) -> &'static str {
    match m {
        NwchemNucleus::Point => "point",
        NwchemNucleus::Finite => "finite",
    }
}
fn vector(out: &mut String, v: [f64; 3]) {
    for n in v {
        write!(out, " {n:.17e}").unwrap();
    }
}
fn check_size(out: &str) -> Result<(), NwchemError> {
    if out.len() > NwchemLimits::default().max_input_bytes {
        return Err(error(0, "output byte limit"));
    }
    Ok(())
}
