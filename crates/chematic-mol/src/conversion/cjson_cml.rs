//! Strict selected-field CJSON -> CML projection, without opaque field loss.
use super::*;

fn unsupported(message: impl Into<String>) -> ConversionError {
    ConversionError::Unsupported(format!("CJSON to CML: {}", message.into()))
}

fn xml_10_character(character: char) -> bool {
    matches!(
        character as u32,
        0x9 | 0xA | 0xD | 0x20..=0xD7FF | 0xE000..=0xFFFD | 0x10000..=0x10FFFF
    )
}

fn validate_source(record: &crate::CjsonRecord) -> Result<(), ConversionError> {
    if record.molecule().atom_count() == 0 {
        return Err(unsupported("requires a nonempty molecule"));
    }
    if let Some(field) = record.opaque_fields().first() {
        return Err(unsupported(format!("cannot preserve opaque field {field}")));
    }
    if !record.partial_charge_datasets().is_empty() {
        return Err(unsupported("cannot preserve partialCharges datasets"));
    }
    if record.has_conflicting_names() {
        return Err(unsupported("conflicting source names"));
    }
    if record
        .name()
        .is_some_and(|name| !name.chars().all(xml_10_character))
    {
        return Err(unsupported("name contains characters outside XML 1.0"));
    }
    Ok(())
}

fn verify(
    record: &crate::CjsonRecord,
    output: &str,
    limits: &ConversionLimits,
) -> Result<(), ConversionError> {
    let (parsed, coords) = crate::parse_cml_with_3d_with_limits(
        output,
        &crate::CmlParseLimits {
            max_input_bytes: limits.max_output_bytes,
            max_atoms: limits.mol2.max_atoms,
            max_bonds: limits.mol2.max_bonds,
            ..Default::default()
        },
    )
    .map_err(|error| unsupported(format!("cannot reparse candidate: {error}")))?;
    if !graph_matches(record.molecule(), &parsed) {
        return Err(unsupported("cannot verify indexed graph"));
    }
    if coords.is_some() != record.coordinates_present() {
        return Err(unsupported("cannot verify coordinate presence"));
    }
    if let Some(coords) = coords {
        let expected = record.coordinates().iter().flat_map(|&(x, y, z)| [x, y, z]);
        if coords.len() != record.coordinates().len()
            || !coords
                .into_iter()
                .flatten()
                .zip(expected)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        {
            return Err(unsupported("cannot verify coordinate bits"));
        }
    }
    let xml = roxmltree::Document::parse(output)
        .map_err(|error| unsupported(format!("invalid candidate XML: {error}")))?;
    let name = xml
        .root_element()
        .children()
        .find(|n| n.has_tag_name(("http://www.xml-cml.org/schema", "name")));
    if name.map(|n| n.text().unwrap_or("")) != record.name() {
        return Err(unsupported("cannot verify name"));
    }
    Ok(())
}

// CJSON carries bond orders, not atom aromaticity flags. CML derives those
// flags from A bonds; compare only the explicitly projected source fields.
fn graph_matches(source: &Molecule, parsed: &Molecule) -> bool {
    source.atom_count() == parsed.atom_count()
        && source.bond_count() == parsed.bond_count()
        && source.atoms().zip(parsed.atoms()).all(|((a, x), (b, y))| {
            a == b && x.element == y.element && x.charge == y.charge && x.isotope == y.isotope
        })
        && source
            .bonds()
            .zip(parsed.bonds())
            .all(|((_, x), (_, y))| x == y)
}

/// Project the local CJSON graph/Cartesian/name profile to CML, or refuse.
///
/// Preserves source atom indices as generated a1..aN IDs, indexed elements,
/// formal charges, isotope masses, bond orders, optional Cartesian coordinate
/// bits and a consistent root/properties name. Absent 3D remains absent.
/// Rejects unknown/opaque fields (including stereo and unit cells), partial
/// charges, conflicting names, empty graphs and non-XML-1.0 name characters.
/// The existing CML writer adds inferred total H counts; chemical validity,
/// geometry-derived stereo and original JSON spelling are not verified.
/// CML derives atom aromaticity flags from aromatic bonds; CJSON does not carry
/// those flags. They are not part of the selected-field equivalence check.
/// Candidate graph, XYZ and name are checked after reparse before returning.
pub fn project_cjson_to_cml(text: &str) -> Result<String, ConversionError> {
    project_cjson_to_cml_with_limits(text, &ConversionLimits::default())
}

/// Checked CJSON -> CML projection with shared input/atom/bond/output bounds.
pub fn project_cjson_to_cml_with_limits(
    text: &str,
    limits: &ConversionLimits,
) -> Result<String, ConversionError> {
    let record = crate::parse_cjson_record_with_limits(
        text,
        &crate::CjsonParseLimits {
            max_input_bytes: limits.mol2.max_input_bytes,
            max_atoms: limits.mol2.max_atoms,
            max_bonds: limits.mol2.max_bonds,
            ..Default::default()
        },
    )
    .map_err(ConversionError::CjsonParse)?;
    validate_source(&record)?;
    let output = crate::cml::write_cml_geometry(
        record.molecule(),
        None,
        record.coordinates_present().then_some(record.coordinates()),
        record.name(),
    );
    if output.len() > limits.max_output_bytes {
        return Err(ConversionError::ResourceLimit("output bytes"));
    }
    verify(&record, &output, limits)?;
    Ok(output)
}
