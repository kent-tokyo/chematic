use super::*;

/// Canonical LF output retaining coordinate kind, 0D–3D basis, constraints,
/// and optional initial density values. Numbers use 17 decimal places in
/// scientific notation. Source text is never replayed after structured edits.
pub fn write_fhiaims(g: &FhiaimsGeometry) -> Result<String, FhiaimsError> {
    g.validate()?;
    let mut out = String::new();
    for vector in &g.lattice_vectors {
        let [x, y, z] = vector.vector_angstrom;
        out.push_str(&format!("lattice_vector {x:.17e} {y:.17e} {z:.17e}\n"));
        write_constraints(&mut out, vector.constraints);
    }
    for atom in &g.atoms {
        let keyword = match atom.coordinates {
            FhiaimsCoordinates::Cartesian(_) => "atom",
            FhiaimsCoordinates::Fractional(_) => "atom_frac",
        };
        let [x, y, z] = atom.coordinates.values();
        out.push_str(&format!(
            "{keyword} {x:.17e} {y:.17e} {z:.17e} {}\n",
            atom.element.symbol()
        ));
        write_constraints(&mut out, atom.constraints);
        if let Some(value) = atom.initial_charge {
            out.push_str(&format!("initial_charge {value:.17e}\n"));
        }
        if let Some(value) = atom.initial_moment {
            out.push_str(&format!("initial_moment {value:.17e}\n"));
        }
        if out.len() > FhiaimsParseLimits::default().max_input_bytes {
            return Err(error(0, "serialized byte limit exceeded"));
        }
    }
    Ok(out)
}

fn write_constraints(out: &mut String, mask: [bool; 3]) {
    if mask == [true; 3] {
        out.push_str("constrain_relaxation .true.\n");
    } else {
        for (axis, constrained) in ["x", "y", "z"].iter().zip(mask) {
            if constrained {
                out.push_str(&format!("constrain_relaxation {axis}\n"));
            }
        }
    }
}
