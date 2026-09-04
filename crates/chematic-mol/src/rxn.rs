//! MDL RXN V2000 file format parser and writer.
//!
//! The MDL RXN format stores reactions as a sequence of MOL blocks.
//!
//! File structure:
//! ```text
//! $RXN
//! <blank line>
//! <program/date line>
//! <comment line>
//! nreactants nproducts
//! $MOL
//! <MOL block for reactant 1>
//! $MOL
//! <MOL block for reactant 2>
//! …
//! $MOL
//! <MOL block for product 1>
//! …
//! ```

use chematic_rxn::{
    ComponentRole, Reaction, ReactionDocument, ReactionDocumentError, ReactionLoss,
};

use crate::error::MolParseError;
use crate::mol2000::parse_mol;

/// Resource limits for parsing an MDL RXN V2000 file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RxnFileParseLimits {
    /// Maximum UTF-8 input size in bytes.
    pub max_input_bytes: usize,
    /// Maximum number of reactants declared by the RXN header.
    pub max_reactants: usize,
    /// Maximum number of products declared by the RXN header.
    pub max_products: usize,
    /// Maximum number of `$MOL` blocks retained from the file.
    pub max_molecules: usize,
}

impl Default for RxnFileParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 64 * 1024 * 1024,
            max_reactants: 10_000,
            max_products: 10_000,
            max_molecules: 20_000,
        }
    }
}

/// Error produced by [`parse_rxn_file`].
#[derive(Debug)]
pub enum RxnParseError {
    /// A configured input or reaction-component limit was exceeded.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
    /// The file does not start with `$RXN`.
    MissingHeader,
    /// The reactant/product count line could not be parsed.
    BadCountLine,
    /// A MOL block inside the RXN file failed to parse.
    MolParse(MolParseError),
}

/// Error returned by the loss-aware RXN/document adapter.
#[derive(Debug)]
pub enum RxnDocumentError {
    /// The upstream-backed RXN V2000 parser rejected the input.
    Rxn(RxnParseError),
    /// The typed document was invalid or could not be represented losslessly.
    Document(ReactionDocumentError),
}

impl core::fmt::Display for RxnDocumentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Rxn(error) => write!(f, "RXN document parse error: {error}"),
            Self::Document(error) => write!(f, "RXN document conversion error: {error}"),
        }
    }
}

impl std::error::Error for RxnDocumentError {}

impl From<RxnParseError> for RxnDocumentError {
    fn from(error: RxnParseError) -> Self {
        Self::Rxn(error)
    }
}

impl From<ReactionDocumentError> for RxnDocumentError {
    fn from(error: ReactionDocumentError) -> Self {
        Self::Document(error)
    }
}

impl core::fmt::Display for RxnParseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(f, "RXN {resource} exceeds limit {limit} (got {actual})"),
            Self::MissingHeader => write!(f, "RXN file must start with $RXN"),
            Self::BadCountLine => write!(f, "cannot parse reactant/product count line"),
            Self::MolParse(e) => write!(f, "MOL parse error in RXN: {e}"),
        }
    }
}

impl std::error::Error for RxnParseError {}

impl From<MolParseError> for RxnParseError {
    fn from(e: MolParseError) -> Self {
        Self::MolParse(e)
    }
}

/// Parse an MDL RXN V2000 string into a [`Reaction`].
pub fn parse_rxn_file(text: &str) -> Result<Reaction, RxnParseError> {
    parse_rxn_file_with_limits(text, RxnFileParseLimits::default())
}

/// Parse an upstream-backed MDL RXN V2000 file into a typed reaction document.
///
/// The classic RXN dialect has reactant and product MOL blocks but no agents,
/// conditions, provenance, or stoichiometric coefficient channel. Components
/// are therefore marked as derived and receive deterministic IDs.
pub fn parse_rxn_document(text: &str) -> Result<ReactionDocument, RxnDocumentError> {
    let reaction = parse_rxn_file(text)?;
    Ok(ReactionDocument::from_reaction(&reaction))
}

/// Write a typed reaction document as MDL RXN V2000 without silent data loss.
///
/// Agents and rich document metadata are rejected with a typed `Losses` error
/// because RXN V2000 has no representational channel for them.
pub fn write_rxn_document(document: &ReactionDocument) -> Result<String, RxnDocumentError> {
    document.validate()?;
    let mut losses = Vec::new();
    if document.steps.len() != 1 {
        losses.push(ReactionLoss {
            field: "steps".to_string(),
            detail: "RXN V2000 has one reaction boundary".to_string(),
        });
    }
    if !document.provenance.is_empty() {
        losses.push(ReactionLoss {
            field: "provenance".to_string(),
            detail: "RXN V2000 has no document provenance channel".to_string(),
        });
    }
    let step = &document.steps[0];
    if !step.conditions.is_empty() {
        losses.push(ReactionLoss {
            field: "conditions".to_string(),
            detail: "RXN V2000 has no reaction conditions channel".to_string(),
        });
    }
    if !step.provenance.is_empty() {
        losses.push(ReactionLoss {
            field: "step.provenance".to_string(),
            detail: "RXN V2000 has no step provenance channel".to_string(),
        });
    }
    for component in &step.components {
        if component.coefficient != 1 {
            losses.push(ReactionLoss {
                field: format!("{}.coefficient", component.id),
                detail: "RXN V2000 has no stoichiometric coefficient channel".to_string(),
            });
        }
        if component.role == ComponentRole::Agent {
            losses.push(ReactionLoss {
                field: format!("{}.role", component.id),
                detail: "RXN V2000 has no agent channel".to_string(),
            });
        }
    }
    if !losses.is_empty() {
        return Err(RxnDocumentError::Document(ReactionDocumentError::Losses(
            losses,
        )));
    }
    let reaction_smiles = document.to_reaction_smiles()?;
    let reaction = chematic_rxn::parse_reaction(&reaction_smiles)
        .map_err(|error| ReactionDocumentError::parse_message(error.to_string()))?;
    Ok(write_rxn_file(&reaction))
}

/// Parse an MDL RXN V2000 string with explicit resource limits.
pub fn parse_rxn_file_with_limits(
    text: &str,
    limits: RxnFileParseLimits,
) -> Result<Reaction, RxnParseError> {
    if text.len() > limits.max_input_bytes {
        return Err(RxnParseError::ResourceLimit {
            resource: "input bytes",
            actual: text.len(),
            limit: limits.max_input_bytes,
        });
    }

    let mut lines = text.lines();

    // Line 1: $RXN
    match lines.next() {
        Some(l) if l.trim() == "$RXN" => {}
        _ => return Err(RxnParseError::MissingHeader),
    }

    // Lines 2-4: blank, program, comment (skip)
    for _ in 0..3 {
        lines.next();
    }

    // Line 5: "  nreactants  nproducts  …"
    let count_line = lines.next().unwrap_or("");
    let counts: Vec<i64> = count_line
        .split_whitespace()
        .filter_map(|s| s.parse().ok())
        .collect();
    if counts.len() < 2 {
        return Err(RxnParseError::BadCountLine);
    }
    let n_reactants = counts[0].max(0) as usize;
    let n_products = counts[1].max(0) as usize;
    if n_reactants > limits.max_reactants {
        return Err(RxnParseError::ResourceLimit {
            resource: "reactants",
            actual: n_reactants,
            limit: limits.max_reactants,
        });
    }
    if n_products > limits.max_products {
        return Err(RxnParseError::ResourceLimit {
            resource: "products",
            actual: n_products,
            limit: limits.max_products,
        });
    }
    let declared_molecules = n_reactants.saturating_add(n_products);
    if declared_molecules > limits.max_molecules {
        return Err(RxnParseError::ResourceLimit {
            resource: "molecules",
            actual: declared_molecules,
            limit: limits.max_molecules,
        });
    }

    // Work directly on the remaining text to find "$MOL" blocks.
    // Find the position of the first "$MOL" in the original text.
    let first_mol_pos = match text.find("$MOL") {
        Some(p) => p,
        None => {
            return Ok(Reaction {
                reactants: vec![],
                agents: vec![],
                products: vec![],
            });
        }
    };
    let mol_section = &text[first_mol_pos..];

    // Split on "$MOL\n" to get individual MOL blocks.
    let mol_blocks = mol_section.split("$MOL\n").skip(1);

    let mut reactants = Vec::with_capacity(n_reactants);
    let mut products = Vec::with_capacity(n_products);

    for (i, block) in mol_blocks.enumerate() {
        if i >= limits.max_molecules {
            return Err(RxnParseError::ResourceLimit {
                resource: "molecules",
                actual: i.saturating_add(1),
                limit: limits.max_molecules,
            });
        }
        // Each block is already a valid MOL V2000 block (3 header lines + data).
        let (mol, _meta) = parse_mol(block)?;
        if i < n_reactants {
            reactants.push(mol);
        } else if i < declared_molecules {
            products.push(mol);
        }
    }

    Ok(Reaction {
        reactants,
        agents: vec![],
        products,
    })
}

/// Write a [`Reaction`] as an MDL RXN V2000 string.
pub fn write_rxn_file(rxn: &Reaction) -> String {
    use crate::mol2000::{MolMetadata, write_mol};

    let mut out = String::new();
    out.push_str("$RXN\n");
    out.push('\n'); // program line (blank)
    out.push_str("     chematic\n"); // program/date
    out.push('\n'); // comment (blank)
    out.push_str(&format!(
        "{:3}{:3}\n",
        rxn.reactants.len(),
        rxn.products.len()
    ));

    let meta = MolMetadata::default();
    for mol in rxn.reactants.iter().chain(rxn.products.iter()) {
        out.push_str("$MOL\n");
        out.push_str(&write_mol(mol, &meta));
    }
    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn minimal_rxn_block() -> String {
        // Build ethane→ethanol reaction by writing real MOL blocks via write_mol.
        use crate::mol2000::{MolMetadata, write_mol};
        use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};

        let mut b = MoleculeBuilder::new();
        let c1 = b.add_atom(Atom::new(Element::C));
        let c2 = b.add_atom(Atom::new(Element::C));
        b.add_bond(c1, c2, BondOrder::Single).unwrap();
        let ethane = b.build();

        let mut b2 = MoleculeBuilder::new();
        let c1 = b2.add_atom(Atom::new(Element::C));
        let c2 = b2.add_atom(Atom::new(Element::C));
        let o = b2.add_atom(Atom::new(Element::O));
        b2.add_bond(c1, c2, BondOrder::Single).unwrap();
        b2.add_bond(c2, o, BondOrder::Single).unwrap();
        let ethanol = b2.build();

        let meta = MolMetadata::default();
        format!(
            "$RXN\n\n     test\n\n  1  1\n$MOL\n{}$MOL\n{}",
            write_mol(&ethane, &meta),
            write_mol(&ethanol, &meta),
        )
    }

    #[test]
    fn test_parse_rxn_file_counts() {
        let rxn = parse_rxn_file(&minimal_rxn_block()).unwrap();
        assert_eq!(rxn.reactants.len(), 1);
        assert_eq!(rxn.products.len(), 1);
        assert_eq!(rxn.reactants[0].atom_count(), 2); // ethane
        assert_eq!(rxn.products[0].atom_count(), 3); // ethanol
    }

    #[test]
    fn test_parse_rxn_missing_header() {
        let err = parse_rxn_file("not a rxn file\n");
        assert!(matches!(err, Err(RxnParseError::MissingHeader)));
    }

    #[test]
    fn test_parse_rxn_resource_limits() {
        let text = minimal_rxn_block();
        let err = parse_rxn_file_with_limits(
            &text,
            RxnFileParseLimits {
                max_input_bytes: text.len() - 1,
                ..Default::default()
            },
        );
        assert!(matches!(
            err,
            Err(RxnParseError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));

        let err = parse_rxn_file_with_limits(
            &text,
            RxnFileParseLimits {
                max_reactants: 0,
                ..Default::default()
            },
        );
        assert!(matches!(
            err,
            Err(RxnParseError::ResourceLimit {
                resource: "reactants",
                ..
            })
        ));
    }

    #[test]
    fn test_write_rxn_file_roundtrip() {
        let rxn = parse_rxn_file(&minimal_rxn_block()).unwrap();
        let written = write_rxn_file(&rxn);
        assert!(written.starts_with("$RXN"));
        // Round-trip: re-parse and check counts.
        let rxn2 = parse_rxn_file(&written).unwrap();
        assert_eq!(rxn2.reactants.len(), 1);
        assert_eq!(rxn2.products.len(), 1);
    }

    #[test]
    fn rxn_document_adapter_round_trips_through_upstream_parser() {
        let source = minimal_rxn_block();
        let document = parse_rxn_document(&source).unwrap();
        assert_eq!(document.steps.len(), 1);
        assert_eq!(document.steps[0].components.len(), 2);
        let written = write_rxn_document(&document).unwrap();
        let decoded = parse_rxn_document(&written).unwrap();
        assert_eq!(decoded.steps[0].components.len(), 2);
    }

    #[test]
    fn rxn_document_rejects_agents_instead_of_dropping_them() {
        let mut document = ReactionDocument::from_reaction_smiles("CC>O>CC").unwrap();
        document.steps[0].components[0].role = ComponentRole::Agent;
        let error = write_rxn_document(&document).unwrap_err();
        assert!(matches!(
            error,
            RxnDocumentError::Document(ReactionDocumentError::Losses(_))
        ));
    }
}
