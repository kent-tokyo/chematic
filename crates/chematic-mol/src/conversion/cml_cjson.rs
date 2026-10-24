//! Refuse unrepresented CML payload before producing selected-field CJSON.
use super::*;
use roxmltree::{Document, Node};

const CML_NS: &str = "http://www.xml-cml.org/schema";

fn unsupported(detail: impl Into<String>) -> ConversionError {
    ConversionError::Unsupported(format!("CML to CJSON: {}", detail.into()))
}

fn attributes(node: Node<'_, '_>, allowed: &[&str]) -> Result<(), ConversionError> {
    for attr in node.attributes() {
        if attr.namespace().is_some() || !allowed.contains(&attr.name()) {
            return Err(unsupported(format!(
                "cannot preserve {} attribute {}",
                node.tag_name().name(),
                attr.name()
            )));
        }
    }
    Ok(())
}

fn element_profile(
    node: Node<'_, '_>,
    mol: Node<'_, '_>,
    root: Node<'_, '_>,
) -> Result<(), ConversionError> {
    let parent = node.parent();
    let tag = node.tag_name().name();
    let allowed: &[&str] = match tag {
        "cml" if node == root => &[],
        "molecule" if node == mol && (node == root || parent == Some(root)) => &["id", "title"],
        "name" if parent == Some(mol) => &[],
        "atomArray" | "bondArray" if parent == Some(mol) => &[],
        "atom"
            if parent == Some(mol)
                || parent.is_some_and(|p| {
                    p.tag_name().namespace() == mol.tag_name().namespace()
                        && p.tag_name().name() == "atomArray"
                        && p.parent() == Some(mol)
                }) =>
        {
            &[
                "id",
                "elementType",
                "formalCharge",
                "isotope",
                "isotopeNumber",
                "hydrogenCount",
                "x3",
                "y3",
                "z3",
            ]
        }
        "bond"
            if parent == Some(mol)
                || parent.is_some_and(|p| {
                    p.tag_name().namespace() == mol.tag_name().namespace()
                        && p.tag_name().name() == "bondArray"
                        && p.parent() == Some(mol)
                }) =>
        {
            &["id", "atomRefs2", "order"]
        }
        _ => return Err(unsupported(format!("cannot preserve element {tag}"))),
    };
    if node != root && node != mol && node.tag_name().namespace() != mol.tag_name().namespace() {
        return Err(unsupported("cannot preserve mixed namespace payload"));
    }
    attributes(node, allowed)?;
    Ok(())
}

fn profile_name(text: &str) -> Result<Option<String>, ConversionError> {
    // The bounded geometry reader has already checked XML nodes/depth/DTD.
    let doc = Document::parse(text).map_err(|e| unsupported(e.to_string()))?;
    let root = doc.root_element();
    let mol = if root.tag_name().name() == "cml" {
        root.children()
            .find(Node::is_element)
            .ok_or_else(|| unsupported("missing molecule"))?
    } else {
        root
    };
    let mut name = mol.attribute("title").map(str::to_owned);
    for node in doc.descendants() {
        if node.is_text() {
            if node.parent().is_some_and(|p| {
                p.tag_name().namespace() == mol.tag_name().namespace()
                    && p.tag_name().name() == "name"
            }) {
                continue;
            }
            if !node
                .text()
                .unwrap_or("")
                .chars()
                .all(|c| matches!(c, ' ' | '\t' | '\r' | '\n'))
            {
                return Err(unsupported("cannot preserve text outside name"));
            }
        }
        if !node.is_element() {
            continue;
        }
        if !matches!(node.tag_name().namespace(), None | Some(CML_NS)) {
            return Err(unsupported("cannot preserve foreign namespace payload"));
        }
        element_profile(node, mol, root)?;
        let tag = node.tag_name().name();
        if tag == "name" {
            if node.children().any(|n| n.is_element()) {
                return Err(unsupported("name must contain text only"));
            }
            let value: String = node
                .children()
                .filter(Node::is_text)
                .filter_map(|n| n.text())
                .collect();
            if name.as_ref().is_some_and(|n| *n != value) {
                return Err(unsupported("conflicting source names"));
            }
            name = Some(value);
        }
        if tag == "bond"
            && !matches!(
                node.attribute("order").map(str::trim),
                Some("1" | "S" | "single" | "2" | "D" | "double" | "3" | "T" | "triple")
            )
        {
            return Err(unsupported(
                "bond order is outside single/double/triple profile",
            ));
        }
        if tag == "atom"
            && ["isotope", "isotopeNumber"].iter().any(|key| {
                node.attribute(*key)
                    .is_some_and(|v| v.trim().parse::<u16>() == Ok(0))
            })
        {
            return Err(unsupported("zero isotope mass is not a label"));
        }
    }
    Ok(name)
}

fn verify_graph(source: &Molecule, parsed: &Molecule) -> Result<(), ConversionError> {
    if source.atom_count() != parsed.atom_count()
        || source.bond_count() != parsed.bond_count()
        || !source
            .bonds()
            .zip(parsed.bonds())
            .all(|((_, a), (_, b))| a == b)
        || !source.atoms().zip(parsed.atoms()).all(|((i, a), (j, b))| {
            i == j && a.element == b.element && a.charge == b.charge && a.isotope == b.isotope
        })
    {
        return Err(unsupported("cannot verify indexed graph"));
    }
    for (idx, _) in source.atoms() {
        if chematic_core::implicit_hcount(source, idx)
            != chematic_core::implicit_hcount(parsed, idx)
        {
            return Err(unsupported(format!(
                "cannot preserve hydrogen count at atom {}",
                idx.0
            )));
        }
    }
    Ok(())
}

/// Project a bounded CML graph/Cartesian/name profile, or refuse before output.
///
/// Preserves source atom order, elements, charges, isotope masses, single/double/
/// triple bonds, optional finite binary64 XYZ (including signed zero), and
/// consistent title/name text. IDs resolve graph references but their spelling
/// is not retained. Explicit hydrogenCount is accepted only if the local CJSON
/// reader infers the same non-node H count. This does not certify chemistry or
/// stereo inferred from geometry. Rejects 2D, aromatic/query/zero bonds, stereo,
/// every other unrepresented attribute/element. Molecule/atom/bond ID spellings
/// are document-local references and are not retained. Comments, PIs,
/// whitespace, XML prefixes and lexical spelling are representation details.
/// Legacy graph readers remain available; source retention preserves opaque XML.
pub fn project_cml_to_cjson(text: &str) -> Result<String, ConversionError> {
    project_cml_to_cjson_with_limits(text, &ConversionLimits::default())
}

/// Checked CML -> CJSON with shared byte/line/atom/bond/output bounds.
pub fn project_cml_to_cjson_with_limits(
    text: &str,
    limits: &ConversionLimits,
) -> Result<String, ConversionError> {
    let (mol, xyz) = crate::parse_cml_with_3d_with_limits(
        text,
        &crate::CmlParseLimits {
            max_input_bytes: limits.mol2.max_input_bytes,
            max_lines: limits.mol2.max_lines,
            max_line_bytes: limits.mol2.max_line_bytes,
            max_atoms: limits.mol2.max_atoms,
            max_bonds: limits.mol2.max_bonds,
            ..Default::default()
        },
    )
    .map_err(|e| unsupported(e.to_string()))?;
    super::cml_source::validate_prolog(text)?;
    let name = profile_name(text)?;
    let points = xyz
        .as_ref()
        .map(|p| p.iter().map(|&[x, y, z]| (x, y, z)).collect::<Vec<_>>())
        .unwrap_or_default();
    let raw = crate::write_cjson(&mol, &points);
    let mut value: serde_json::Value =
        serde_json::from_str(&raw).map_err(|e| unsupported(e.to_string()))?;
    if let Some(name) = &name {
        value["name"] = serde_json::Value::String(name.clone());
    }
    let output = value.to_string();
    if output.len() > limits.max_output_bytes {
        return Err(ConversionError::ResourceLimit("output bytes"));
    }
    let record = crate::parse_cjson_record_with_limits(
        &output,
        &crate::CjsonParseLimits {
            max_input_bytes: limits.max_output_bytes,
            max_atoms: limits.mol2.max_atoms,
            max_bonds: limits.mol2.max_bonds,
            ..Default::default()
        },
    )
    .map_err(ConversionError::CjsonParse)?;
    verify_graph(&mol, record.molecule())?;
    if record.coordinates_present() != xyz.is_some()
        || record.name() != name.as_deref()
        || record
            .coordinates()
            .iter()
            .flat_map(|&(x, y, z)| [x, y, z])
            .zip(points.iter().flat_map(|&(x, y, z)| [x, y, z]))
            .any(|(a, b)| a.to_bits() != b.to_bits())
    {
        return Err(unsupported("cannot verify coordinates/name"));
    }
    Ok(output)
}
