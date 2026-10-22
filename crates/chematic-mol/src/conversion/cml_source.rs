//! One well-formed CML molecule document, retained without graph projection.
use super::*;

const CML_NS: &str = "http://www.xml-cml.org/schema";
const MAX_XML_NODES: u32 = 200_000;
const MAX_XML_DEPTH: usize = 64;

fn malformed(message: impl Into<String>) -> ConversionError {
    ConversionError::Unsupported(format!("invalid CML source: {}", message.into()))
}

// Inspect only the prolog, skipping comments/PIs instead of matching XML-like
// text inside opaque payload. roxmltree subsequently checks the whole document.
fn validate_prolog(text: &str) -> Result<(), ConversionError> {
    let ascii = text.is_ascii();
    let mut rest = text.strip_prefix('\u{feff}').unwrap_or(text);
    loop {
        rest = rest.trim_start();
        let (closing, declaration) = if rest.starts_with("<!--") {
            ("-->", false)
        } else if rest.starts_with("<?") {
            (
                "?>",
                rest.strip_prefix("<?xml")
                    .is_some_and(|s| s.starts_with(char::is_whitespace)),
            )
        } else if rest.starts_with("<!DOCTYPE") {
            return Err(malformed(
                "CML source retention does not support DTD declarations",
            ));
        } else {
            break;
        };
        let Some(end) = rest.find(closing).map(|i| i + closing.len()) else {
            break; // Whole-document parsing reports the unterminated construct.
        };
        if declaration {
            if rest[..end].contains('&') {
                return Err(malformed(
                    "XML declarations cannot contain entity references",
                ));
            }
            let attrs = crate::cml::parse_xml_attrs(&rest[..end]);
            if attrs.get("version").map(String::as_str) != Some("1.0") {
                return Err(malformed("CML source retention requires XML 1.0"));
            }
            if let Some(encoding) = attrs.get("encoding") {
                match encoding.to_ascii_lowercase().as_str() {
                    "utf-8" => {}
                    "us-ascii" if ascii => {}
                    _ => {
                        return Err(malformed(
                            "CML source retention requires UTF-8 or ASCII text",
                        ));
                    }
                }
            }
        }
        rest = &rest[end..];
    }
    Ok(())
}

fn validate(text: &str, limits: &ConversionLimits) -> Result<(), ConversionError> {
    if text.len() > limits.mol2.max_input_bytes {
        return Err(ConversionError::ResourceLimit("input bytes"));
    }
    if text.lines().count() > limits.mol2.max_lines {
        return Err(ConversionError::ResourceLimit("lines"));
    }
    if text
        .lines()
        .any(|line| line.len() > limits.mol2.max_line_bytes)
    {
        return Err(ConversionError::ResourceLimit("line bytes"));
    }
    validate_prolog(text)?;
    let doc = roxmltree::Document::parse_with_options(
        text,
        roxmltree::ParsingOptions {
            allow_dtd: false,
            nodes_limit: MAX_XML_NODES,
        },
    )
    .map_err(|error| match error {
        roxmltree::Error::NodesLimitReached => ConversionError::ResourceLimit("CML XML nodes"),
        _ => malformed(error.to_string()),
    })?;
    let root = doc.root_element();
    if !root.has_tag_name((CML_NS, "molecule")) {
        return Err(malformed(
            "CML source retention requires one namespaced molecule root",
        ));
    }
    let (mut atoms, mut bonds) = (0usize, 0usize);
    for node in root.descendants().filter(|n| n.is_element()) {
        if node
            .ancestors()
            .filter(|n| n.is_element())
            .take(MAX_XML_DEPTH + 1)
            .count()
            > MAX_XML_DEPTH
        {
            return Err(ConversionError::ResourceLimit("CML XML depth"));
        }
        if node != root && node.has_tag_name((CML_NS, "molecule")) {
            return Err(malformed(
                "nested CML molecule hierarchies are outside source-retention profile",
            ));
        }
        atoms += usize::from(node.has_tag_name((CML_NS, "atom")));
        bonds += usize::from(node.has_tag_name((CML_NS, "bond")));
        if atoms > limits.mol2.max_atoms || bonds > limits.mol2.max_bonds {
            return Err(ConversionError::ResourceLimit("CML atom/bond elements"));
        }
    }
    Ok(())
}

pub(super) fn convert(
    text: &str,
    policy: ConversionPolicy,
    limits: &ConversionLimits,
) -> Result<ConversionReport, ConversionError> {
    validate(text, limits)?;
    if text.len() > limits.max_output_bytes {
        return Err(ConversionError::ResourceLimit("output bytes"));
    }
    Ok(ConversionReport {
        schema: "chematic.conversion-report.v1",
        input_format: "cml",
        output_format: "cml",
        policy,
        status: "completed",
        source_record_index: 0,
        verification_scope: "exact_source_bytes",
        candidate_graph_verified: false,
        output: Some(text.to_owned()),
        observations: vec![
            ConversionObservation {
                field: "source_text",
                status: PreservationStatus::Preserved,
                stage: "source_retention",
                reason: "Exact UTF-8 source bytes retained after bounded XML parsing and single-molecule framing checks; no edits or chemical interpretation promised.",
            },
            ConversionObservation {
                field: "source_semantics",
                status: PreservationStatus::OpaqueRetained,
                stage: "source_retention",
                reason: "Coordinate spellings/dimensionality, IDs, properties, namespaces, comments, stereo and extensions remain in their source representation; chemistry and references are not interpreted or verified.",
            },
        ],
        candidate_atom_mapping: Vec::new(),
        candidate_bond_mapping: Vec::new(),
    })
}
