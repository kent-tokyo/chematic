//! Namespace-aware, bounded XML selection for the single-molecule projection.
use super::{CmlError, CmlParseLimits};
use roxmltree::{Document, Node};
use std::collections::HashMap;

const CML_NS: &str = "http://www.xml-cml.org/schema";
const MAX_DEPTH: usize = 64;
const MAX_NODES: usize = 4_000_000;
type Attributes = HashMap<String, String>;
type Records = (Vec<Attributes>, Vec<Attributes>);

fn malformed(detail: &str) -> CmlError {
    CmlError::MalformedXml(detail.into())
}
fn limit(resource: &'static str, actual: usize, limit: usize) -> CmlError {
    CmlError::ResourceLimit {
        resource,
        actual,
        limit,
    }
}
fn cml_namespace(node: Node<'_, '_>) -> bool {
    matches!(node.tag_name().namespace(), None | Some(CML_NS))
}

fn in_scope(node: Node<'_, '_>, mol: Node<'_, '_>) -> bool {
    node.ancestors()
        .take_while(|n| *n != mol)
        .filter(Node::is_element)
        .all(|n| n.tag_name().namespace() == mol.tag_name().namespace())
}

fn validate_tree(doc: &Document<'_>, limits: &CmlParseLimits) -> Result<(), CmlError> {
    let mut count = 0;
    for node in doc.descendants().filter(Node::is_element) {
        count += 1;
        if count > limits.max_xml_elements {
            return Err(limit("XML elements", count, limits.max_xml_elements));
        }
        if node
            .ancestors()
            .filter(Node::is_element)
            .take(MAX_DEPTH + 1)
            .count()
            > MAX_DEPTH
        {
            return Err(limit("XML depth", MAX_DEPTH + 1, MAX_DEPTH));
        }
    }
    Ok(())
}

fn molecule<'a, 'input>(doc: &'a Document<'input>) -> Result<Node<'a, 'input>, CmlError> {
    let root = doc.root_element();
    if !cml_namespace(root) || !matches!(root.tag_name().name(), "molecule" | "cml") {
        return Err(malformed(
            "CML projection requires a CML molecule or cml root",
        ));
    }
    let selected = if root.tag_name().name() == "molecule" {
        root
    } else {
        let mut molecules = root.children().filter(|n| {
            n.is_element()
                && n.tag_name().name() == "molecule"
                && cml_namespace(*n)
                && (root.tag_name().namespace().is_none()
                    || n.tag_name().namespace() == root.tag_name().namespace())
        });
        let first = molecules.next().ok_or(CmlError::MissingMolecule)?;
        if molecules.next().is_some() {
            return Err(malformed(
                "CML 3D reader requires exactly one direct molecule",
            ));
        }
        first
    };
    if selected
        .descendants()
        .skip(1)
        .any(|n| n.is_element() && n.tag_name().name() == "molecule" && in_scope(n, selected))
    {
        return Err(malformed(
            "CML 3D reader requires exactly one direct molecule",
        ));
    }
    Ok(selected)
}

fn graph_child(node: Node<'_, '_>, mol: Node<'_, '_>, array: &str) -> bool {
    let Some(parent) = node.parent() else {
        return false;
    };
    // Direct atom/bond children retain the legacy unwrapped-record extension.
    parent == mol || (parent.tag_name().name() == array && parent.parent() == Some(mol))
}

pub(super) fn records(input: &str, limits: &CmlParseLimits) -> Result<Records, CmlError> {
    let nodes_limit = limits
        .max_xml_elements
        .saturating_mul(4)
        .saturating_add(8)
        .min(MAX_NODES);
    let doc = Document::parse_with_options(
        input,
        roxmltree::ParsingOptions {
            allow_dtd: false,
            nodes_limit: nodes_limit as u32,
        },
    )
    .map_err(|error| match error {
        roxmltree::Error::NodesLimitReached => limit("XML nodes", nodes_limit + 1, nodes_limit),
        _ => CmlError::MalformedXml(error.to_string()),
    })?;
    validate_tree(&doc, limits)?;
    let mol = molecule(&doc)?;
    let (mut atoms, mut bonds) = (Vec::new(), Vec::new());
    for node in mol.descendants().filter(Node::is_element) {
        if !in_scope(node, mol) {
            continue;
        }
        let (target, cap, resource) = match node.tag_name().name() {
            "atom" if graph_child(node, mol, "atomArray") => {
                (&mut atoms, limits.max_atoms, "atom elements")
            }
            "bond" if graph_child(node, mol, "bondArray") => {
                (&mut bonds, limits.max_bonds, "bond elements")
            }
            _ => continue,
        };
        if target.len() >= cap {
            return Err(limit(resource, target.len() + 1, cap));
        }
        target.push(
            node.attributes()
                .filter(|a| a.namespace().is_none())
                .map(|a| (a.name().to_owned(), a.value().to_owned()))
                .collect(),
        );
    }
    Ok((atoms, bonds))
}
