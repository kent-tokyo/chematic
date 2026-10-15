use super::*;

pub fn parse_fhiaims(text: &str) -> Result<FhiaimsGeometry, FhiaimsError> {
    parse_fhiaims_with_limits(text, &FhiaimsParseLimits::default())
}

pub fn parse_fhiaims_with_limits(
    text: &str,
    limits: &FhiaimsParseLimits,
) -> Result<FhiaimsGeometry, FhiaimsError> {
    if text.len() > limits.max_input_bytes {
        return Err(error(0, "input byte limit exceeded"));
    }
    let mut geometry = FhiaimsGeometry {
        atoms: Vec::new(),
        lattice_vectors: Vec::new(),
        source_text: String::new(),
    };
    let mut target = None;
    for (index, raw) in text.lines().enumerate() {
        let line = index + 1;
        if raw.len() > limits.max_line_len {
            return Err(error(line, "line length limit exceeded"));
        }
        let fields: Vec<_> = raw
            .split('#')
            .next()
            .unwrap_or("")
            .split_whitespace()
            .collect();
        if fields.is_empty() {
            continue;
        }
        read_record(&fields, line, limits, &mut geometry, &mut target)?;
    }
    geometry.validate()?;
    geometry.source_text = text.to_string();
    Ok(geometry)
}

#[derive(Clone, Copy)]
enum Target {
    Atom(usize),
    Lattice(usize),
}

fn read_record(
    fields: &[&str],
    line: usize,
    limits: &FhiaimsParseLimits,
    g: &mut FhiaimsGeometry,
    target: &mut Option<Target>,
) -> Result<(), FhiaimsError> {
    match fields[0] {
        "atom" | "atom_frac" => {
            if fields.len() != 5 {
                return Err(error(
                    line,
                    "expected atom/atom_frac x y z canonical-element",
                ));
            }
            if g.atoms.len()
                >= limits
                    .max_atoms
                    .min(FhiaimsParseLimits::default().max_atoms)
            {
                return Err(error(line, "atom count limit exceeded"));
            }
            let element = Element::from_symbol(fields[4]).filter(|e| e.symbol() == fields[4])
                .ok_or_else(|| error(line, "custom species, isotope aliases and dummy sites require external definitions; canonical elements only"))?;
            let values = triple(&fields[1..4], line)?;
            *target = Some(Target::Atom(g.atoms.len()));
            g.atoms.push(FhiaimsAtom {
                element,
                coordinates: if fields[0] == "atom" {
                    FhiaimsCoordinates::Cartesian(values)
                } else {
                    FhiaimsCoordinates::Fractional(values)
                },
                constraints: [false; 3],
                initial_charge: None,
                initial_moment: None,
            });
        }
        "lattice_vector" => {
            if fields.len() != 4 || g.lattice_vectors.len() >= 3 {
                return Err(error(
                    line,
                    "expected at most three lattice_vector x y z records",
                ));
            }
            *target = Some(Target::Lattice(g.lattice_vectors.len()));
            g.lattice_vectors.push(FhiaimsLatticeVector {
                vector_angstrom: triple(&fields[1..], line)?,
                constraints: [false; 3],
            });
        }
        "constrain_relaxation" => {
            if fields.len() != 2 {
                return Err(error(line, "one constraint per line is required"));
            }
            let mask = match target {
                Some(Target::Atom(i)) => &mut g.atoms[*i].constraints,
                Some(Target::Lattice(i)) => &mut g.lattice_vectors[*i].constraints,
                None => {
                    return Err(error(
                        line,
                        "constraint has no preceding atom or lattice vector",
                    ));
                }
            };
            constrain(mask, fields[1], line)?;
        }
        "initial_charge" | "initial_moment" => initial(fields, line, g, *target)?,
        _ => {
            return Err(error(
                line,
                format!("unsupported geometry directive '{}'", fields[0]),
            ));
        }
    }
    Ok(())
}

fn initial(
    fields: &[&str],
    line: usize,
    g: &mut FhiaimsGeometry,
    target: Option<Target>,
) -> Result<(), FhiaimsError> {
    if fields.len() != 2 {
        return Err(error(line, "initial charge/moment requires one scalar"));
    }
    let Some(Target::Atom(index)) = target else {
        return Err(error(
            line,
            "initial charge/moment requires a preceding atom",
        ));
    };
    let atom = &mut g.atoms[index];
    let out = if fields[0] == "initial_charge" {
        &mut atom.initial_charge
    } else {
        &mut atom.initial_moment
    };
    if out.is_some() {
        return Err(error(line, "duplicate initial charge/moment"));
    }
    *out = Some(number(fields[1], line)?);
    Ok(())
}

fn constrain(mask: &mut [bool; 3], value: &str, line: usize) -> Result<(), FhiaimsError> {
    match value {
        ".true." => *mask = [true; 3],
        ".false." if mask.iter().any(|v| *v) => {
            return Err(error(
                line,
                "conflicting .false. after an active constraint",
            ));
        }
        ".false." => {}
        "x" => mask[0] = true,
        "y" => mask[1] = true,
        "z" => mask[2] = true,
        _ => {
            return Err(error(
                line,
                "unsupported constraint; expected .true., .false., x, y or z",
            ));
        }
    }
    Ok(())
}

fn number(text: &str, line: usize) -> Result<f64, FhiaimsError> {
    let value = text
        .replace(['D', 'd'], "e")
        .parse::<f64>()
        .map_err(|_| error(line, "invalid real number"))?;
    if !value.is_finite() {
        return Err(error(line, "non-finite real number"));
    }
    Ok(value)
}
fn triple(fields: &[&str], line: usize) -> Result<[f64; 3], FhiaimsError> {
    Ok([
        number(fields[0], line)?,
        number(fields[1], line)?,
        number(fields[2], line)?,
    ])
}
