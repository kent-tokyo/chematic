use std::collections::HashMap;

use chematic_core::{AtomIdx, Molecule};
use chematic_smiles::{parse as parse_smiles, write as write_smiles};

/// A chemical reaction with reactants, agents, and products.
pub struct Reaction {
    pub reactants: Vec<Molecule>,
    pub agents: Vec<Molecule>,
    pub products: Vec<Molecule>,
}

/// Resource limits for parsing untrusted reaction SMILES.
#[derive(Debug, Clone, Copy)]
pub struct ReactionParseLimits {
    /// Maximum UTF-8 input size in bytes.
    pub max_input_bytes: usize,
    /// Maximum number of non-empty dot-separated components per side.
    pub max_components_per_side: usize,
    /// Maximum atoms in one component.
    pub max_atoms_per_molecule: usize,
    /// Maximum bonds in one component.
    pub max_bonds_per_molecule: usize,
}

impl Default for ReactionParseLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 16 * 1024 * 1024,
            max_components_per_side: 10_000,
            max_atoms_per_molecule: 100_000,
            max_bonds_per_molecule: 200_000,
        }
    }
}

/// Error type for reaction SMILES parsing.
#[derive(Debug)]
pub enum RxnError {
    /// The string does not contain two `>` delimiters (reaction arrow).
    MissingArrow,
    /// The reaction exceeded a configured resource limit.
    ResourceLimit {
        resource: &'static str,
        actual: usize,
        limit: usize,
    },
    /// A SMILES component failed to parse.
    SmilesParse { part: String, source: String },
}

impl core::fmt::Display for RxnError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MissingArrow => write!(f, "reaction SMILES must contain '>>'"),
            Self::ResourceLimit {
                resource,
                actual,
                limit,
            } => write!(
                f,
                "reaction {resource} exceeds limit {limit} (got {actual})"
            ),
            Self::SmilesParse { part, source } => {
                write!(f, "failed to parse SMILES '{part}': {source}")
            }
        }
    }
}

impl std::error::Error for RxnError {}

/// Parse a reaction SMILES string of the form `"reactants>agents>products"`.
///
/// - `>` splits the string into 3 parts (reactants, agents, products).
/// - Each part is a dot-separated list of SMILES; empty parts yield empty `Vec`s.
/// - `"R>>P"` is the standard form with empty agents section.
/// - Returns `Err(RxnError::MissingArrow)` if fewer than two `>` delimiters are found.
pub fn parse_reaction(s: &str) -> Result<Reaction, RxnError> {
    parse_reaction_with_limits(s, &ReactionParseLimits::default())
}

/// Parse a reaction SMILES while enforcing input and per-component limits.
pub fn parse_reaction_with_limits(
    s: &str,
    limits: &ReactionParseLimits,
) -> Result<Reaction, RxnError> {
    if s.len() > limits.max_input_bytes {
        return Err(RxnError::ResourceLimit {
            resource: "input bytes",
            actual: s.len(),
            limit: limits.max_input_bytes,
        });
    }
    // splitn(3, '>') yields 3 parts when at least two `>` are present;
    // fewer parts means the reaction arrow is missing.
    let parts: Vec<&str> = s.splitn(3, '>').collect();
    if parts.len() < 3 {
        return Err(RxnError::MissingArrow);
    }

    let parse_part = |side: &'static str, s: &str| -> Result<Vec<Molecule>, RxnError> {
        if s.is_empty() {
            return Ok(Vec::new());
        }
        let components: Vec<&str> = s.split('.').filter(|p| !p.is_empty()).collect();
        if components.len() > limits.max_components_per_side {
            return Err(RxnError::ResourceLimit {
                resource: side,
                actual: components.len(),
                limit: limits.max_components_per_side,
            });
        }
        components
            .into_iter()
            .filter(|p| !p.is_empty())
            .map(|p| {
                let molecule = parse_smiles(p).map_err(|e| RxnError::SmilesParse {
                    part: p.to_string(),
                    source: e.to_string(),
                })?;
                if molecule.atom_count() > limits.max_atoms_per_molecule {
                    return Err(RxnError::ResourceLimit {
                        resource: "atoms per molecule",
                        actual: molecule.atom_count(),
                        limit: limits.max_atoms_per_molecule,
                    });
                }
                if molecule.bond_count() > limits.max_bonds_per_molecule {
                    return Err(RxnError::ResourceLimit {
                        resource: "bonds per molecule",
                        actual: molecule.bond_count(),
                        limit: limits.max_bonds_per_molecule,
                    });
                }
                Ok(molecule)
            })
            .collect()
    };

    Ok(Reaction {
        reactants: parse_part("reactants", parts[0])?,
        agents: parse_part("agents", parts[1])?,
        products: parse_part("products", parts[2])?,
    })
}

/// Serialize a `Reaction` back to a reaction SMILES string.
pub fn write_reaction(rxn: &Reaction) -> String {
    let join = |mols: &[Molecule]| -> String {
        mols.iter().map(write_smiles).collect::<Vec<_>>().join(".")
    };
    format!(
        "{}>{}>{}",
        join(&rxn.reactants),
        join(&rxn.agents),
        join(&rxn.products),
    )
}

/// The center of a chemical reaction: changed atoms and bonds.
#[derive(Debug, Clone)]
pub struct ReactionCenter {
    /// Bonds present in reactants but not in products (broken bonds).
    pub broken_bonds: Vec<(AtomIdx, AtomIdx)>,
    /// Bonds present in products but not in reactants (formed bonds).
    pub formed_bonds: Vec<(AtomIdx, AtomIdx)>,
    /// Atoms whose element, charge, or aromaticity changed between reactants and products.
    /// Uses reactant-side indexing.
    pub changed_atoms: Vec<AtomIdx>,
}

/// Identify the reaction center: bonds broken/formed and atoms changed.
///
/// Uses atom_map numbers (if present) to match reactant atoms to product atoms.
/// For each mapped atom pair:
/// - Compares bond connectivity to identify broken/formed bonds.
/// - Compares element, charge, aromaticity to identify changed atoms.
///
/// Returns empty vecs if atoms lack atom_map annotations.
pub fn find_reaction_center(rxn: &Reaction) -> ReactionCenter {
    let mut broken_bonds = Vec::new();
    let mut formed_bonds = Vec::new();
    let mut changed_atoms = Vec::new();

    // Build atom_map -> (mol_idx, atom_idx) for reactants and products
    let mut reactant_map: HashMap<u16, (usize, AtomIdx)> = HashMap::new();
    for (mol_idx, mol) in rxn.reactants.iter().enumerate() {
        for (atom_idx, atom) in mol.atoms() {
            if let Some(map_num) = atom.atom_map {
                reactant_map.insert(map_num, (mol_idx, atom_idx));
            }
        }
    }

    let mut product_map: HashMap<u16, (usize, AtomIdx)> = HashMap::new();
    for (mol_idx, mol) in rxn.products.iter().enumerate() {
        for (atom_idx, atom) in mol.atoms() {
            if let Some(map_num) = atom.atom_map {
                product_map.insert(map_num, (mol_idx, atom_idx));
            }
        }
    }

    // If no atom_map, return empty
    if reactant_map.is_empty() {
        return ReactionCenter {
            broken_bonds,
            formed_bonds,
            changed_atoms,
        };
    }

    // Identify broken bonds (edges in reactants not in products)
    for map_num in reactant_map.keys() {
        let (r_mol_idx, r_atom_idx) = reactant_map[map_num];
        let r_mol = &rxn.reactants[r_mol_idx];
        for (r_neighbor, _) in r_mol.neighbors(r_atom_idx) {
            if let Some(neighbor_map) = r_mol.atom(r_neighbor).atom_map
                && neighbor_map > *map_num
            {
                // Check if this bond exists in products
                if let Some((p_mol_idx, p_atom_idx)) = product_map.get(map_num) {
                    let p_mol = &rxn.products[*p_mol_idx];
                    if let Some((_, p_neighbor)) = product_map.get(&neighbor_map) {
                        let bond_exists = p_mol.bond_between(*p_atom_idx, *p_neighbor).is_some();
                        if !bond_exists {
                            broken_bonds.push((r_atom_idx, r_neighbor));
                        }
                    } else {
                        broken_bonds.push((r_atom_idx, r_neighbor));
                    }
                } else {
                    broken_bonds.push((r_atom_idx, r_neighbor));
                }
            }
        }
    }

    // Identify formed bonds (edges in products not in reactants)
    for map_num in product_map.keys() {
        let (p_mol_idx, p_atom_idx) = product_map[map_num];
        let p_mol = &rxn.products[p_mol_idx];
        for (p_neighbor, _) in p_mol.neighbors(p_atom_idx) {
            if let Some(neighbor_map) = p_mol.atom(p_neighbor).atom_map
                && neighbor_map > *map_num
            {
                // Check if this bond exists in reactants
                if let Some((r_mol_idx, r_atom_idx)) = reactant_map.get(map_num) {
                    let r_mol = &rxn.reactants[*r_mol_idx];
                    if let Some((_, r_neighbor)) = reactant_map.get(&neighbor_map) {
                        let bond_exists = r_mol.bond_between(*r_atom_idx, *r_neighbor).is_some();
                        if !bond_exists {
                            formed_bonds.push((p_atom_idx, p_neighbor));
                        }
                    } else {
                        formed_bonds.push((p_atom_idx, p_neighbor));
                    }
                } else {
                    formed_bonds.push((p_atom_idx, p_neighbor));
                }
            }
        }
    }

    // Identify changed atoms (element/charge/aromaticity changes)
    for map_num in reactant_map.keys() {
        let (r_mol_idx, r_atom_idx) = reactant_map[map_num];
        let r_atom = rxn.reactants[r_mol_idx].atom(r_atom_idx);

        if let Some((p_mol_idx, p_atom_idx)) = product_map.get(map_num) {
            let p_atom = rxn.products[*p_mol_idx].atom(*p_atom_idx);

            if r_atom.element != p_atom.element
                || r_atom.charge != p_atom.charge
                || r_atom.aromatic != p_atom.aromatic
            {
                changed_atoms.push(r_atom_idx);
            }
        }
    }

    ReactionCenter {
        broken_bonds,
        formed_bonds,
        changed_atoms,
    }
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
        assert_eq!(
            rxn.reactants[0].atom_count(),
            rxn2.reactants[0].atom_count()
        );
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

    #[test]
    fn test_parse_reaction_limits_reject_input_and_components() {
        let limits = ReactionParseLimits {
            max_input_bytes: 3,
            ..ReactionParseLimits::default()
        };
        assert!(matches!(
            parse_reaction_with_limits("C>>C", &limits),
            Err(RxnError::ResourceLimit {
                resource: "input bytes",
                ..
            })
        ));

        let limits = ReactionParseLimits {
            max_components_per_side: 1,
            ..ReactionParseLimits::default()
        };
        assert!(matches!(
            parse_reaction_with_limits("C.C>>C", &limits),
            Err(RxnError::ResourceLimit {
                resource: "reactants",
                ..
            })
        ));
    }

    #[test]
    fn test_parse_reaction_limits_reject_component_size() {
        let limits = ReactionParseLimits {
            max_atoms_per_molecule: 1,
            ..ReactionParseLimits::default()
        };
        assert!(matches!(
            parse_reaction_with_limits("CC>>C", &limits),
            Err(RxnError::ResourceLimit {
                resource: "atoms per molecule",
                ..
            })
        ));
    }
}
