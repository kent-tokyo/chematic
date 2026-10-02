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
    /// An atomic-number SMARTS primitive could not be expanded safely.
    UnsupportedAtomicNumberPrimitive { primitive: String },
    /// Expanding aromatic/aliphatic alternatives would exceed the safety cap.
    AtomicNumberExpansionLimit { actual: usize, limit: usize },
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
            Self::UnsupportedAtomicNumberPrimitive { primitive } => write!(
                f,
                "unsupported atomic-number SMARTS primitive '{primitive}'; expected [#N] or [#N:map]"
            ),
            Self::AtomicNumberExpansionLimit { actual, limit } => write!(
                f,
                "atomic-number SMARTS expansion exceeds limit {limit} (got {actual})"
            ),
        }
    }
}

/// Expand simple atomic-number SMARTS atoms into SMILES-compatible variants.
///
/// Reaction application uses the SMILES parser for template construction, while
/// the SMARTS parser already accepts query atoms such as `[#7:2]`. This helper
/// is the explicit compatibility boundary between those representations. The
/// atom map is retained verbatim and the returned order is deterministic.
///
/// The supported forms are `[#number]`, `[#number:map]`, and the common
/// explicit-single-hydrogen forms `[#number;H1]` / `[#number;H1:map]`.
/// More complex primitives must remain on the SMARTS query path and return an
/// explicit error rather than silently changing semantics.
///
/// Reactant-side atoms (every part before the last `>`) are queries: each
/// primitive expands to its aliphatic form and, where chemically valid, its
/// aromatic form (`[#7:2]` becomes `[N:2]` and `[n:2]`), and the variants
/// enumerate every combination. A string without `>` is treated as a
/// reactant-side pattern.
///
/// Product-side atoms (the last part) are literals, not queries (issue #679):
///
/// * an unmapped primitive becomes one aliphatic atom with implicit
///   hydrogens — organic-subset symbols are written bare (`[#6](=[#8])[#6]`
///   becomes `C(=O)C`), others in brackets (`[#14]` becomes `[Si]`), and
///   `;H1` keeps its explicit hydrogen (`[#6;H1]` becomes `[CH]`);
/// * a mapped primitive keeps the aromaticity of the matched reactant atom:
///   when the reactant side spells the same map as an atomic-number primitive
///   the two share one choice per variant (`[#7:1]>>[#7:1]C` becomes
///   `[N:1]>>[N:1]C` and `[n:1]>>[n:1]C`); otherwise it follows the case of
///   the reactant-side atom with that map, or is aliphatic when there is none.
pub fn expand_atomic_number_primitives(s: &str) -> Result<Vec<String>, RxnError> {
    const MAX_VARIANTS: usize = 256;
    const ORGANIC_SUBSET: [&str; 10] = ["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"];

    // Split into `>`-separated parts outside brackets; the last part is the
    // product side when there is at least one `>`.
    let bytes = s.as_bytes();
    let mut part_starts = vec![0usize];
    let mut depth = 0usize;
    for (i, &b) in bytes.iter().enumerate() {
        match b {
            b'[' => depth += 1,
            b']' => depth = depth.saturating_sub(1),
            b'>' if depth == 0 => part_starts.push(i + 1),
            _ => {}
        }
    }
    let has_products = part_starts.len() > 1;
    let product_start = *part_starts.last().unwrap_or(&0);
    let reactant_text = if has_products {
        &s[..product_start]
    } else {
        ""
    };

    // Aromaticity spelled for a mapped reactant-side bracket atom (`[n:1]` →
    // aromatic, `[N:1]` → aliphatic); `None` when the map is absent or the
    // atom is itself an atomic-number primitive.
    let reactant_case = |map: &str| -> Option<bool> {
        let mut rest = reactant_text;
        while let Some(open) = rest.find('[') {
            let Some(close_rel) = rest[open + 1..].find(']') else {
                break;
            };
            let inner = &rest[open + 1..open + 1 + close_rel];
            rest = &rest[open + 1 + close_rel + 1..];
            if !inner.ends_with(map) || !inner[..inner.len() - map.len()].ends_with(':') {
                continue;
            }
            let body = inner.trim_start_matches(|c: char| c.is_ascii_digit());
            match body.bytes().next() {
                Some(b'#') => return None,
                Some(c) if c.is_ascii_lowercase() => return Some(true),
                Some(c) if c.is_ascii_uppercase() => return Some(false),
                _ => return None,
            }
        }
        None
    };

    // Text segments and expansion slots, in order. A slot holds its
    // alternatives and the index of the choice group it belongs to (mapped
    // primitives with the same map share a group).
    enum Piece {
        Text(String),
        Slot { options: Vec<String>, group: usize },
    }
    let mut pieces: Vec<Piece> = Vec::new();
    let mut groups: Vec<usize> = Vec::new(); // option count per group
    let mut group_of_map: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    let mut text = String::new();

    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] != b'[' {
            text.push(bytes[i] as char);
            i += 1;
            continue;
        }
        let start = i;
        let Some(end_rel) = s[i + 1..].find(']') else {
            text.push_str(&s[i..]);
            break;
        };
        let end = i + 1 + end_rel;
        let primitive = &s[start..=end];
        let inner = &s[start + 1..end];
        if !inner.starts_with('#') {
            text.push_str(primitive);
            i = end + 1;
            continue;
        }
        let number_end = inner[1..]
            .bytes()
            .position(|b| !b.is_ascii_digit())
            .map(|p| p + 1)
            .unwrap_or(inner.len());
        let number_text = &inner[1..number_end];
        let suffix = &inner[number_end..];
        let (hydrogen_count, map_suffix) = if suffix.is_empty() {
            (None, "")
        } else if suffix.starts_with(':')
            && suffix.len() > 1
            && suffix[1..].bytes().all(|b| b.is_ascii_digit())
        {
            (None, suffix)
        } else if let Some(map) = suffix.strip_prefix(";H1:")
            && !map.is_empty()
            && map.bytes().all(|b| b.is_ascii_digit())
        {
            (Some(1_u8), &suffix[3..])
        } else if suffix == ";H1" {
            (Some(1_u8), "")
        } else {
            (None, "")
        };
        if number_text.is_empty()
            || (hydrogen_count.is_none() && map_suffix.is_empty() && !suffix.is_empty())
        {
            return Err(RxnError::UnsupportedAtomicNumberPrimitive {
                primitive: primitive.to_string(),
            });
        }
        let atomic_number = number_text.parse::<u8>().ok();
        let Some(atomic_number) = atomic_number.filter(|n| (1..=118).contains(n)) else {
            return Err(RxnError::UnsupportedAtomicNumberPrimitive {
                primitive: primitive.to_string(),
            });
        };
        let Some(element) = chematic_core::Element::from_atomic_number(atomic_number) else {
            return Err(RxnError::UnsupportedAtomicNumberPrimitive {
                primitive: primitive.to_string(),
            });
        };
        let symbol = element.symbol();
        let hydrogen = hydrogen_count.map_or("", |_| "H");
        let can_be_aromatic = matches!(atomic_number, 5 | 6 | 7 | 8 | 15 | 16);
        let aliphatic = format!("[{symbol}{hydrogen}{map_suffix}]");
        let aromatic = format!(
            "[{}{}{}]",
            symbol.to_ascii_lowercase(),
            hydrogen,
            map_suffix
        );
        let on_product_side = has_products && start >= product_start;
        let map = map_suffix.strip_prefix(':');

        let options: Vec<String> = if !on_product_side {
            if can_be_aromatic {
                vec![aliphatic, aromatic]
            } else {
                vec![aliphatic]
            }
        } else if let Some(map) = map {
            if group_of_map.contains_key(map) {
                // Shares the reactant-side choice (same option order).
                if can_be_aromatic {
                    vec![aliphatic, aromatic]
                } else {
                    vec![aliphatic]
                }
            } else {
                match reactant_case(map) {
                    Some(true) if can_be_aromatic => vec![aromatic],
                    _ => vec![aliphatic],
                }
            }
        } else if hydrogen_count.is_none() && ORGANIC_SUBSET.contains(&symbol) {
            vec![symbol.to_string()]
        } else {
            vec![aliphatic]
        };

        if options.len() == 1 {
            text.push_str(&options[0]);
            i = end + 1;
            continue;
        }
        let group = match map {
            Some(map) if !on_product_side || group_of_map.contains_key(map) => {
                *group_of_map.entry(map.to_string()).or_insert_with(|| {
                    groups.push(options.len());
                    groups.len() - 1
                })
            }
            _ => {
                groups.push(options.len());
                groups.len() - 1
            }
        };
        if !text.is_empty() {
            pieces.push(Piece::Text(std::mem::take(&mut text)));
        }
        pieces.push(Piece::Slot { options, group });
        i = end + 1;
    }
    if !text.is_empty() {
        pieces.push(Piece::Text(text));
    }

    let total: usize = groups
        .iter()
        .try_fold(1usize, |acc, &n| {
            let next = acc.saturating_mul(n);
            (next <= MAX_VARIANTS).then_some(next)
        })
        .unwrap_or(usize::MAX);
    if total > MAX_VARIANTS {
        return Err(RxnError::AtomicNumberExpansionLimit {
            actual: total,
            limit: MAX_VARIANTS,
        });
    }

    // Enumerate choice vectors in lexicographic order (first group slowest).
    let mut choice = vec![0usize; groups.len()];
    let mut variants = Vec::with_capacity(total);
    loop {
        let mut out = String::with_capacity(s.len());
        for piece in &pieces {
            match piece {
                Piece::Text(t) => out.push_str(t),
                Piece::Slot { options, group } => out.push_str(&options[choice[*group]]),
            }
        }
        variants.push(out);
        let mut k = groups.len();
        loop {
            if k == 0 {
                return Ok(variants);
            }
            k -= 1;
            choice[k] += 1;
            if choice[k] < groups[k] {
                break;
            }
            choice[k] = 0;
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

    #[test]
    fn expands_atomic_number_primitives_and_preserves_maps() {
        let variants = expand_atomic_number_primitives("[#7:2][#6]").unwrap();
        assert_eq!(
            variants,
            vec![
                "[N:2][C]".to_string(),
                "[N:2][c]".to_string(),
                "[n:2][C]".to_string(),
                "[n:2][c]".to_string(),
            ]
        );
    }

    #[test]
    fn rejects_compound_atomic_number_primitives() {
        assert!(matches!(
            expand_atomic_number_primitives("[#7;H2]>>[#7;H2]"),
            Err(RxnError::UnsupportedAtomicNumberPrimitive { .. })
        ));
    }

    #[test]
    fn expands_single_hydrogen_constraint_and_preserves_maps() {
        let variants = expand_atomic_number_primitives("[#7;H1:2][#6;H1]").unwrap();
        assert_eq!(
            variants,
            vec![
                "[NH:2][CH]".to_string(),
                "[NH:2][cH]".to_string(),
                "[nH:2][CH]".to_string(),
                "[nH:2][cH]".to_string(),
            ]
        );
    }
}
