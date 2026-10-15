use super::{CdjsonError, CdjsonParseLimits, limit};
use chematic_core::{BondOrder, Chirality, Molecule, implicit_hcount};
use serde_json::{Value, json};

fn atom_values(mol: &Molecule, coordinates: &[[f64; 3]]) -> Result<Vec<Value>, CdjsonError> {
    let mut atoms = Vec::new();
    for ((idx, atom), coordinate) in mol.atoms().zip(coordinates) {
        if atom.element.atomic_number() == 0
            || atom.aromatic
            || atom.wildcard
            || atom.chirality != Chirality::None
            || atom.atom_map.is_some()
        {
            return Err(CdjsonError::Unsupported(
                "aromatic/query/chiral/mapped atom in graph writer".into(),
            ));
        }
        if coordinate.iter().any(|v| !v.is_finite()) {
            return Err(CdjsonError::Invalid("nonfinite coordinate".into()));
        }
        if atom.isotope == Some(0) {
            return Err(CdjsonError::Invalid("zero isotope mass".into()));
        }
        let mut value = json!({"x":coordinate[0],"y":coordinate[1],"z":coordinate[2],
            "l":atom.element.symbol(),"c":atom.charge,"h":implicit_hcount(mol, idx)});
        if let Some(mass) = atom.isotope {
            value["m"] = json!(mass);
        }
        atoms.push(value);
    }
    Ok(atoms)
}

/// Write one molecular graph with explicit inferred/stored H overrides.
///
/// Coordinates retain caller drawing units. The writer creates no persistent
/// IDs and does not export tags or unrelated source metadata. Use the immutable
/// parsed document's `source_text()` for exact same-format preservation.
pub fn write_cdjson_graph(mol: &Molecule, coordinates: &[[f64; 3]]) -> Result<String, CdjsonError> {
    let hard = CdjsonParseLimits::default();
    limit("atoms", mol.atom_count(), hard.max_atoms, hard.max_atoms)?;
    limit("bonds", mol.bond_count(), hard.max_bonds, hard.max_bonds)?;
    if coordinates.len() != mol.atom_count() {
        return Err(CdjsonError::Invalid(
            "coordinate count must equal atom count".into(),
        ));
    }
    if !mol.stereo_groups().is_empty() {
        return Err(CdjsonError::Unsupported("enhanced stereo groups".into()));
    }
    let atoms = atom_values(mol, coordinates)?;
    let mut bonds = Vec::new();
    for (idx, bond) in mol.bonds() {
        if bond.atom1 == bond.atom2 {
            return Err(CdjsonError::Invalid("self bond".into()));
        }
        if mol.bond_direction(idx).is_some() {
            return Err(CdjsonError::Unsupported("directional bond".into()));
        }
        let order = match bond.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            _ => return Err(CdjsonError::Unsupported("bond order outside 1/2/3".into())),
        };
        bonds.push(json!({"b":bond.atom1.0,"e":bond.atom2.0,"o":order}));
    }
    let output = serde_json::to_string(&json!({"m":[{"a":atoms,"b":bonds}]}))
        .map_err(|e| CdjsonError::Invalid(e.to_string()))?;
    limit(
        "output bytes",
        output.len(),
        hard.max_input_bytes,
        hard.max_input_bytes,
    )?;
    Ok(output)
}
