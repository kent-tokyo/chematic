//! Checked PQR to rigid PDBQT adaptation without bond or charge perception.
//!
//! Independently implemented from APBS PQR and AutoDock 4.2 field layouts.
//! Callers supply validated AutoDock types in source atom order. PQR has no
//! bond table, so this adapter cannot assign chemically reliable docking types.
use crate::pqr::PqrResult;
use chematic_core::Element;

/// Conversion failures use zero-based source atom indices.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PqrPdbqtError {
    EmptyInput,
    AtomTypeCount {
        atoms: usize,
        types: usize,
    },
    /// PDBQT has no per-atom radius field; explicit consent is required.
    RadiusLossRequiresConsent,
    InvalidField {
        atom: usize,
        field: &'static str,
    },
    InvalidAtomType {
        atom: usize,
        atom_type: String,
    },
    ElementMismatch {
        atom: usize,
        atom_type: String,
    },
}

impl core::fmt::Display for PqrPdbqtError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "PQR to PDBQT requires at least one atom"),
            Self::AtomTypeCount { atoms, types } => write!(
                f,
                "PQR has {atoms} atoms but {types} AutoDock types were supplied"
            ),
            Self::RadiusLossRequiresConsent => write!(
                f,
                "PDBQT cannot store PQR radii; set allow_radius_loss=true to retain radii in the sidecar report"
            ),
            Self::InvalidField { atom, field } => write!(
                f,
                "PQR atom {atom}: {field} is invalid or cannot fit its PDBQT field"
            ),
            Self::InvalidAtomType { atom, atom_type } => write!(
                f,
                "PQR atom {atom}: unsupported AutoDock type '{atom_type}'"
            ),
            Self::ElementMismatch { atom, atom_type } => write!(
                f,
                "PQR atom {atom}: AutoDock type '{atom_type}' does not match the source element"
            ),
        }
    }
}
impl std::error::Error for PqrPdbqtError {}

/// Source-order sidecar. PDBQT retains charges to three decimal places and
/// coordinates to three decimal places in Å. Original values remain here.
#[derive(Debug, Clone, PartialEq)]
pub struct PqrPdbqtReport {
    pub pdbqt: String,
    pub source_charges: Vec<f64>,
    pub source_radii: Vec<f64>,
    pub source_coords: Vec<(f64, f64, f64)>,
    /// Always true: the PDBQT payload has no radius field.
    pub radii_omitted: bool,
    pub max_charge_rounding_error: f64,
    pub max_coordinate_rounding_error: f64,
}

fn type_element(atom_type: &str) -> Option<Element> {
    let symbol = match atom_type {
        "A" | "C" => "C",
        "N" | "NA" => "N",
        "O" | "OA" | "OS" => "O",
        "S" | "SA" => "S",
        "H" | "HD" => "H",
        "P" | "F" | "Cl" | "Br" | "I" | "Mg" | "Mn" | "Zn" | "Ca" | "Fe" => atom_type,
        _ => return None,
    };
    Element::from_symbol(symbol)
}

fn text_fits(s: &str, width: usize) -> bool {
    !s.is_empty() && s.len() <= width && s.bytes().all(|b| b.is_ascii_graphic())
}

fn rounded_field(
    value: f64,
    width: usize,
    atom: usize,
    field: &'static str,
) -> Result<(String, f64), PqrPdbqtError> {
    let invalid = || PqrPdbqtError::InvalidField { atom, field };
    if !value.is_finite() {
        return Err(invalid());
    }
    let text = format!("{value:>width$.3}");
    if text.len() != width {
        return Err(invalid());
    }
    let rounded: f64 = text.trim().parse().map_err(|_| invalid())?;
    Ok((text, (rounded - value).abs()))
}

/// Convert source records to rigid receptor PDBQT, preserving source order,
/// serials, names, residue identity and partial charges without recalculation.
///
/// Explicit `atom_types` avoid guessing chemistry from a bondless PQR graph.
/// Their element compatibility is checked; acceptor/aromaticity assignments
/// remain the caller's responsibility. No hydrogens are merged or removed.
/// `allow_radius_loss` acknowledges that radii survive only in the report.
/// Unrepresentable metadata/numeric fields always fail, even with consent.
/// This does not prepare a flexible ligand or calculate a torsion tree.
pub fn pqr_to_pdbqt(
    source: &PqrResult,
    atom_types: &[&str],
    allow_radius_loss: bool,
) -> Result<PqrPdbqtReport, PqrPdbqtError> {
    if source.atoms.is_empty() {
        return Err(PqrPdbqtError::EmptyInput);
    }
    if source.atoms.len() != atom_types.len() {
        return Err(PqrPdbqtError::AtomTypeCount {
            atoms: source.atoms.len(),
            types: atom_types.len(),
        });
    }
    if !allow_radius_loss {
        return Err(PqrPdbqtError::RadiusLossRequiresConsent);
    }
    // PDB serials have five columns; keep the output bounded even for records
    // constructed directly rather than obtained from the limited PQR parser.
    if source.atoms.len() > 99_999 {
        return Err(PqrPdbqtError::InvalidField {
            atom: 99_999,
            field: "atom count",
        });
    }
    let mut result = PqrPdbqtReport {
        pdbqt: String::from("REMARK  PQR charges retained; radii omitted (see sidecar)\n"),
        source_charges: Vec::with_capacity(source.atoms.len()),
        source_radii: Vec::with_capacity(source.atoms.len()),
        source_coords: Vec::with_capacity(source.atoms.len()),
        radii_omitted: true,
        max_charge_rounding_error: 0.0,
        max_coordinate_rounding_error: 0.0,
    };
    let mut serials = std::collections::HashSet::new();
    for (i, (a, atom_type)) in source.atoms.iter().zip(atom_types).enumerate() {
        validate_identity(a, i, &mut serials)?;
        let element = type_element(atom_type).ok_or_else(|| PqrPdbqtError::InvalidAtomType {
            atom: i,
            atom_type: (*atom_type).into(),
        })?;
        if element != a.element {
            return Err(PqrPdbqtError::ElementMismatch {
                atom: i,
                atom_type: (*atom_type).into(),
            });
        }
        if !a.radius.is_finite() || a.radius < 0.0 {
            return Err(PqrPdbqtError::InvalidField {
                atom: i,
                field: "radius",
            });
        }
        let (x, dx) = rounded_field(a.x, 8, i, "x")?;
        let (y, dy) = rounded_field(a.y, 8, i, "y")?;
        let (z, dz) = rounded_field(a.z, 8, i, "z")?;
        let (q, dq) = rounded_field(a.charge, 6, i, "charge")?;
        result.max_charge_rounding_error = result.max_charge_rounding_error.max(dq);
        result.max_coordinate_rounding_error =
            result.max_coordinate_rounding_error.max(dx.max(dy).max(dz));
        let name = if a.atom_name.len() < 4 && a.element.symbol().len() == 1 {
            format!(" {:<3}", a.atom_name)
        } else {
            format!("{:<4}", a.atom_name)
        };
        result.pdbqt.push_str(&format!(
            "{:<6}{:>5} {} {:>3} {}{:>4}{}   {}{}{}{:>6.2}{:>6.2}    {} {:<2}\n",
            a.group_pdb,
            a.serial,
            name,
            a.res_name,
            a.chain_id.as_deref().unwrap_or(" "),
            a.res_seq,
            a.icode.unwrap_or(' '),
            x,
            y,
            z,
            1.0,
            0.0,
            q,
            atom_type,
        ));
        result.source_charges.push(a.charge);
        result.source_radii.push(a.radius);
        result.source_coords.push((a.x, a.y, a.z));
    }
    result.pdbqt.push_str("END\n");
    Ok(result)
}

fn validate_identity(
    a: &crate::pqr::PqrAtomRecord,
    i: usize,
    serials: &mut std::collections::HashSet<i64>,
) -> Result<(), PqrPdbqtError> {
    for (valid, field) in [
        (matches!(a.group_pdb.as_str(), "ATOM" | "HETATM"), "record"),
        (
            (1..=99_999).contains(&a.serial) && serials.insert(a.serial),
            "serial",
        ),
        (text_fits(&a.atom_name, 4), "atom name"),
        (text_fits(&a.res_name, 3), "residue name"),
        (
            a.chain_id.as_deref().is_none_or(|s| text_fits(s, 1)),
            "chain",
        ),
        ((-999..=9999).contains(&a.res_seq), "residue number"),
        (
            a.icode.is_none_or(|c| c.is_ascii_alphabetic()),
            "insertion code",
        ),
    ] {
        if !valid {
            return Err(PqrPdbqtError::InvalidField { atom: i, field });
        }
    }
    Ok(())
}
