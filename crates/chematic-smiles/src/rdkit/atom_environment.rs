//! Atom environments following RDKit Release_2026_09_1 Subgraphs.cpp.
//! Copyright (C) 2003-2026 Greg Landrum and other RDKit contributors.
//! BSD-3-Clause; see THIRD_PARTY_NOTICES.md. The ordered bond frontier, including
//! duplicate pending bonds, is necessary for upstream enforceSize behavior.
use super::{RdkitSmilesError, graph_queries, mol::Mol};
use chematic_core::Molecule;
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy)]
pub struct RdkitAtomEnvironmentLimits {
    /// Maximum number of popped frontier entries, including duplicates.
    pub max_visits: usize,
    /// Maximum pending entries in either one-layer frontier.
    pub max_frontier: usize,
}
impl Default for RdkitAtomEnvironmentLimits {
    fn default() -> Self {
        Self {
            max_visits: 1_000_000,
            max_frontier: 100_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RdkitAtomEnvironment {
    /// Ordered model bond indices, as in FindAtomEnvironmentOfRadiusN.
    pub bond_indices: Vec<usize>,
    /// Model atom index -> minimum radius; returned independently of caller maps.
    pub atom_distances: BTreeMap<usize, u32>,
}

fn unsupported(detail: impl Into<String>) -> RdkitSmilesError {
    RdkitSmilesError::Unsupported(detail.into())
}

/// FindAtomEnvironmentOfRadiusN on the existing sanitized RDKit model.
/// Graph H removal, atom/bond order and input caps match rdkit_shortest_path.
/// With use_hs=false only *neighbor* hydrogens are skipped; an H root can still
/// reach a heavy atom. Radius zero returns no bonds and {root: 0}.
pub fn rdkit_atom_environment(
    mol: &Molecule,
    radius: u32,
    rooted_at_atom: usize,
    use_hs: bool,
    enforce_size: bool,
) -> Result<RdkitAtomEnvironment, RdkitSmilesError> {
    rdkit_atom_environment_with_limits(
        mol,
        radius,
        rooted_at_atom,
        use_hs,
        enforce_size,
        &RdkitAtomEnvironmentLimits::default(),
    )
}

/// Bounded ordered bond expansion. Exhausting a requested radius returns empty
/// bonds and map when enforce_size=true, or the accumulated environment otherwise.
/// Pending duplicate bonds can satisfy a final layer without adding new bonds,
/// as in upstream cyclic environments. Limits return errors, never partial data.
pub fn rdkit_atom_environment_with_limits(
    mol: &Molecule,
    radius: u32,
    rooted_at_atom: usize,
    use_hs: bool,
    enforce_size: bool,
    limits: &RdkitAtomEnvironmentLimits,
) -> Result<RdkitAtomEnvironment, RdkitSmilesError> {
    let model = graph_queries::prepare(mol)?;
    if rooted_at_atom >= model.atoms.len() {
        return Err(unsupported("atom environment root is out of range"));
    }
    let mut result = RdkitAtomEnvironment {
        bond_indices: Vec::new(),
        atom_distances: BTreeMap::from([(rooted_at_atom, 0)]),
    };
    if radius == 0 {
        return Ok(result);
    }
    let mut frontier = Vec::new();
    let mut next = Vec::new();
    let mut seen = Vec::new();
    seen.try_reserve_exact(model.bonds.len())
        .map_err(|e| unsupported(format!("environment allocation failed: {e}")))?;
    seen.resize(model.bonds.len(), false);
    result
        .bond_indices
        .try_reserve_exact(model.bonds.len())
        .map_err(|e| unsupported(format!("environment allocation failed: {e}")))?;
    append_neighbors(&model, rooted_at_atom, use_hs, &seen, &mut frontier, limits)?;
    let mut visits = 0;
    let mut processed = 0;
    for layer in 0..radius {
        if frontier.is_empty() {
            break;
        }
        next.clear();
        for &(start, bond_index) in &frontier {
            if visits >= limits.max_visits {
                return Err(unsupported(
                    "atom environment visit resource limit exceeded",
                ));
            }
            visits += 1;
            if seen[bond_index] {
                continue;
            }
            seen[bond_index] = true;
            result.bond_indices.push(bond_index);
            let other = model.bonds[bond_index].other(start);
            result
                .atom_distances
                .entry(other)
                .and_modify(|v| *v = (*v).min(layer + 1))
                .or_insert(layer + 1);
            if layer < radius - 1 {
                append_neighbors(&model, other, use_hs, &seen, &mut next, limits)?;
            }
        }
        processed += 1;
        std::mem::swap(&mut frontier, &mut next);
    }
    if processed != radius && enforce_size {
        result.bond_indices.clear();
        result.atom_distances.clear();
    }
    Ok(result)
}

fn append_neighbors(
    model: &Mol,
    atom: usize,
    use_hs: bool,
    seen: &[bool],
    frontier: &mut Vec<(usize, usize)>,
    limits: &RdkitAtomEnvironmentLimits,
) -> Result<(), RdkitSmilesError> {
    for &bond in &model.atom_bonds[atom] {
        if seen[bond] || (!use_hs && model.atoms[model.bonds[bond].other(atom)].anum == 1) {
            continue;
        }
        if frontier.len() >= limits.max_frontier {
            return Err(unsupported(
                "atom environment frontier resource limit exceeded",
            ));
        }
        frontier
            .try_reserve(1)
            .map_err(|e| unsupported(format!("environment allocation failed: {e}")))?;
        frontier.push((atom, bond));
    }
    Ok(())
}
