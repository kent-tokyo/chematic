use chematic_core::Molecule;
use chematic_smiles::{parse as parse_smiles, write as write_smiles};

/// A chemical reaction with reactants, agents, and products.
pub struct Reaction {
    pub reactants: Vec<Molecule>,
    pub agents: Vec<Molecule>,
    pub products: Vec<Molecule>,
}

/// Error type for reaction SMILES parsing.
#[derive(Debug)]
pub enum RxnError {
    /// The string does not contain two `>` delimiters (reaction arrow).
    MissingArrow,
    /// A SMILES component failed to parse.
    SmilesParse { part: String, source: String },
}

impl core::fmt::Display for RxnError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            RxnError::MissingArrow => write!(f, "reaction SMILES must contain '>>'"),
            RxnError::SmilesParse { part, source } => {
                write!(f, "failed to parse SMILES '{part}': {source}")
            }
        }
    }
}

/// Parse a reaction SMILES string of the form `"reactants>agents>products"`.
///
/// - `>` splits the string into 3 parts (reactants, agents, products).
/// - Each part is a dot-separated list of SMILES; empty parts yield empty `Vec`s.
/// - `"R>>P"` is the standard form with empty agents section.
/// - Returns `Err(RxnError::MissingArrow)` if fewer than two `>` delimiters are found.
pub fn parse_reaction(s: &str) -> Result<Reaction, RxnError> {
    // Split on ">" into at most 3 parts: reactants, agents, products
    let parts: Vec<&str> = s.splitn(3, '>').collect();
    if parts.len() < 3 {
        return Err(RxnError::MissingArrow);
    }
    // parts[0] = reactants, parts[1] = agents, parts[2] = products
    // "CC>>CC"   → ["CC",  "",  "CC"]  ✓
    // "CC>O>CC=O"→ ["CC", "O", "CC=O"] ✓
    // "CC>CC"    → ["CC", "CC"]         → MissingArrow ✓
    // ">>"       → ["",   "",  ""]      ✓

    let parse_part = |s: &str| -> Result<Vec<Molecule>, RxnError> {
        if s.is_empty() {
            return Ok(Vec::new());
        }
        s.split('.')
            .filter(|p| !p.is_empty())
            .map(|p| {
                parse_smiles(p).map_err(|e| RxnError::SmilesParse {
                    part: p.to_string(),
                    source: e.to_string(),
                })
            })
            .collect()
    };

    Ok(Reaction {
        reactants: parse_part(parts[0])?,
        agents:    parse_part(parts[1])?,
        products:  parse_part(parts[2])?,
    })
}

/// Serialize a `Reaction` back to a reaction SMILES string.
pub fn write_reaction(rxn: &Reaction) -> String {
    let join = |mols: &[Molecule]| -> String {
        mols.iter().map(|m| write_smiles(m)).collect::<Vec<_>>().join(".")
    };
    format!("{}>{}>{}",
        join(&rxn.reactants),
        join(&rxn.agents),
        join(&rxn.products),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_core::AtomIdx;

    #[test]
    fn test_simple_reaction() {
        let rxn = parse_reaction("[Na+].[Cl-]>>[Na+].[Cl-]").unwrap();
        assert_eq!(rxn.reactants.len(), 2);
        assert_eq!(rxn.agents.len(), 0);
        assert_eq!(rxn.products.len(), 2);
    }

    #[test]
    fn test_one_to_one() {
        let rxn = parse_reaction("CC>>CC").unwrap();
        assert_eq!(rxn.reactants.len(), 1);
        assert_eq!(rxn.agents.len(), 0);
        assert_eq!(rxn.products.len(), 1);
    }

    #[test]
    fn test_with_agent() {
        let rxn = parse_reaction("CC>O>CC=O").unwrap();
        assert_eq!(rxn.reactants.len(), 1);
        assert_eq!(rxn.agents.len(), 1);
        assert_eq!(rxn.products.len(), 1);
        // agent is water (O) — one heavy atom
        assert_eq!(rxn.agents[0].atom_count(), 1);
    }

    #[test]
    fn test_two_reactants() {
        let rxn = parse_reaction("C.C>>CC").unwrap();
        assert_eq!(rxn.reactants.len(), 2);
        assert_eq!(rxn.agents.len(), 0);
        assert_eq!(rxn.products.len(), 1);
    }

    #[test]
    fn test_empty_reaction() {
        let rxn = parse_reaction(">>").unwrap();
        assert_eq!(rxn.reactants.len(), 0);
        assert_eq!(rxn.agents.len(), 0);
        assert_eq!(rxn.products.len(), 0);
    }

    #[test]
    fn test_missing_arrow_error() {
        // "CC>CC" has only one ">" → splitn(3, '>') gives 2 parts → MissingArrow
        let err = parse_reaction("CC>CC");
        assert!(matches!(err, Err(RxnError::MissingArrow)));
    }

    #[test]
    fn test_invalid_smiles_error() {
        // "[X]" has an unrecognised element symbol → SmilesParse error
        let err = parse_reaction("[X]>>C");
        assert!(matches!(err, Err(RxnError::SmilesParse { .. })));
    }

    #[test]
    fn test_atom_map_preserved() {
        // [CH3:1] has atom_map = 1
        let rxn = parse_reaction("[CH3:1]>>[CH3:1]").unwrap();
        assert_eq!(rxn.reactants[0].atom(AtomIdx(0)).atom_map, Some(1));
        assert_eq!(rxn.products[0].atom(AtomIdx(0)).atom_map, Some(1));
    }

    #[test]
    fn test_product_atom_count() {
        let rxn = parse_reaction("CC>>CCO").unwrap();
        assert_eq!(rxn.products[0].atom_count(), 3); // C, C, O
    }

    #[test]
    fn test_empty_agents() {
        let rxn = parse_reaction("C>>C").unwrap();
        assert_eq!(rxn.agents.len(), 0);
    }

    #[test]
    fn test_complex_reaction() {
        // aspirin synthesis-like
        let rxn = parse_reaction("OC(=O)c1ccccc1.CC(=O)O>>CC(=O)Oc1ccccc1C(=O)O");
        assert!(rxn.is_ok());
        let rxn = rxn.unwrap();
        assert_eq!(rxn.reactants.len(), 2);
        assert_eq!(rxn.products.len(), 1);
    }

    #[test]
    fn test_reactant_only() {
        let rxn = parse_reaction("C>>").unwrap();
        assert_eq!(rxn.reactants.len(), 1);
        assert_eq!(rxn.products.len(), 0);
    }

    #[test]
    fn test_multi_product() {
        let rxn = parse_reaction(">>C.C").unwrap();
        assert_eq!(rxn.products.len(), 2);
    }

    #[test]
    fn test_roundtrip() {
        let s = "CC>>CC=O";
        let rxn = parse_reaction(s).unwrap();
        let written = write_reaction(&rxn);
        let rxn2 = parse_reaction(&written).unwrap();
        assert_eq!(rxn.reactants.len(), rxn2.reactants.len());
        assert_eq!(rxn.products.len(), rxn2.products.len());
        assert_eq!(rxn.reactants[0].atom_count(), rxn2.reactants[0].atom_count());
    }

    #[test]
    fn test_write_reaction_format() {
        let rxn = parse_reaction("CC>O>CC=O").unwrap();
        let s = write_reaction(&rxn);
        // Standard reaction SMILES format: "reactants>agents>products"
        assert!(s.contains('>'), "written reaction should contain '>'");
        // Should have exactly 2 '>' characters
        assert_eq!(s.chars().filter(|&c| c == '>').count(), 2);
    }
}
