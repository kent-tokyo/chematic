//! Bounded, opt-in MOL2 bond-token conversion for Open Babel #2668.
use std::collections::{HashMap, HashSet};

use chematic_core::{AtomIdx, BondIdx, BondOrder, KekuleError, Molecule, implicit_hcount};

use crate::{Mol2Error, Mol2ParseLimits, parse_mol2_source_record_with_limits};

#[derive(Debug)]
#[non_exhaustive]
pub enum Mol2KekuleError {
    Parse(Mol2Error),
    Unsupported(String),
    Kekulization(KekuleError),
}
impl core::fmt::Display for Mol2KekuleError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Parse(e) => e.fmt(f),
            Self::Unsupported(s) => write!(f, "MOL2 Kekule profile: {s}"),
            Self::Kekulization(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for Mol2KekuleError {}

fn aromatic_cycle(mol: &Molecule, excluded: BondIdx) -> bool {
    let bond = mol.bond(excluded);
    let mut seen = HashSet::from([bond.atom1]);
    let mut todo = vec![bond.atom1];
    while let Some(atom) = todo.pop() {
        for (next, edge) in mol.neighbors(atom) {
            if edge != excluded && mol.bond(edge).order == BondOrder::Aromatic && seen.insert(next)
            {
                if next == bond.atom2 {
                    return true;
                }
                todo.push(next);
            }
        }
    }
    false
}

fn profile(mol: &Molecule) -> Result<HashSet<AtomIdx>, Mol2KekuleError> {
    let mut atoms = HashSet::new();
    for (idx, bond) in mol.bonds() {
        if !matches!(
            bond.order,
            BondOrder::Single
                | BondOrder::Double
                | BondOrder::Triple
                | BondOrder::Aromatic
                | BondOrder::Up
                | BondOrder::Down
        ) {
            return Err(Mol2KekuleError::Unsupported("query/zero bonds".into()));
        }
        if bond.order == BondOrder::Aromatic {
            if !aromatic_cycle(mol, idx) {
                return Err(Mol2KekuleError::Unsupported(
                    "ar bond outside an aromatic cycle".into(),
                ));
            }
            atoms.extend([bond.atom1, bond.atom2]);
        }
    }
    for &idx in &atoms {
        let atom = mol.atom(idx);
        if atom.charge != 0
            || !matches!(atom.element.atomic_number(), 6 | 7 | 8 | 16)
            || !(2..=3).contains(&mol.degree(idx))
        {
            return Err(Mol2KekuleError::Unsupported(format!(
                "aromatic atom {} requires neutral C/N/O/S with degree 2 or 3",
                idx.0
            )));
        }
    }
    Ok(atoms)
}

/// Convert `ar` bond tokens to an explicit single/double assignment.
///
/// Only those tokens change: whitespace, line endings, atom types, coordinates,
/// names, IDs, charges and other source fields remain byte-for-byte identical.
/// Aromatic atom types deliberately remain (e.g. C.ar); this is a bond encoding
/// choice, not atom retyping or loss of aromatic chemical identity.
///
/// One record, <=8 MiB / 512 atoms / 1024 bonds. The conservative profile accepts
/// neutral aromatic C/N/O/S cycles with degree 2 or 3, ordinary covalent bonds,
/// unchanged per-atom implicit H counts, and re-perceived aromaticity of every
/// converted bond after conversion. Pyrrole-like H
/// must be explicit. Failed matching/unsupported inputs return no partial text.
/// Inputs without `ar` bonds are validated bounded exact no-ops. The ordinary
/// source round-trip and graph-only writers are unchanged.
pub fn kekulize_mol2_record(text: &str) -> Result<String, Mol2KekuleError> {
    let record = parse_mol2_source_record_with_limits(
        text,
        &Mol2ParseLimits {
            max_input_bytes: 8 * 1024 * 1024,
            max_atoms: 512,
            max_bonds: 1024,
            ..Default::default()
        },
    )
    .map_err(Mol2KekuleError::Parse)?;
    let mol = record.molecule();
    if !mol.bonds().any(|(_, b)| b.order == BondOrder::Aromatic) {
        return Ok(text.to_owned());
    }
    let atoms = profile(mol)?;
    let orders = chematic_core::kekulize(mol).map_err(Mol2KekuleError::Kekulization)?;
    let mut explicit = chematic_core::apply_kekule(mol, &orders);
    for &idx in &atoms {
        explicit.set_atom_aromatic(idx, true);
        let valence: u32 = explicit
            .neighbors(idx)
            .map(|(_, bond)| u32::from(explicit.bond(bond).order.order_int()))
            .sum();
        let maximum = match explicit.atom(idx).element.atomic_number() {
            6 => 4,
            7 => 3,
            _ => 2,
        };
        if valence > maximum {
            return Err(Mol2KekuleError::Unsupported(format!(
                "explicit valence exceeds neutral profile at atom {}",
                idx.0
            )));
        }
        if implicit_hcount(mol, idx) != implicit_hcount(&explicit, idx) {
            return Err(Mol2KekuleError::Unsupported(format!(
                "implicit hydrogen count changes at atom {}",
                idx.0
            )));
        }
    }
    let aromaticity = chematic_perception::assign_aromaticity(&explicit);
    if orders.keys().any(|&idx| !aromaticity.is_bond_aromatic(idx)) {
        return Err(Mol2KekuleError::Unsupported(
            "explicit assignment does not preserve every declared aromatic bond".into(),
        ));
    }
    let replacements: HashMap<_, _> = orders
        .iter()
        .map(|(idx, order)| {
            (
                record.bonds()[idx.0 as usize].source_id,
                if *order == BondOrder::Double {
                    "2"
                } else {
                    "1"
                },
            )
        })
        .collect();
    rewrite_bonds(text, &replacements)
}

fn rewrite_bonds(text: &str, replacements: &HashMap<u32, &str>) -> Result<String, Mol2KekuleError> {
    let mut output = String::with_capacity(text.len());
    let mut in_bonds = false;
    let mut changed = 0;
    for line in text.split_inclusive('\n') {
        let trimmed = line.trim();
        if trimmed
            .get(..9)
            .is_some_and(|s| s.eq_ignore_ascii_case("@<TRIPOS>"))
        {
            in_bonds = trimmed.eq_ignore_ascii_case("@<TRIPOS>BOND");
        }
        let mut tokens = line.split_whitespace();
        let replacement = in_bonds
            .then(|| tokens.next()?.parse().ok())
            .flatten()
            .and_then(|id| replacements.get(&id));
        if let Some(order) = replacement {
            let token = line
                .split_whitespace()
                .nth(3)
                .ok_or_else(|| Mol2KekuleError::Unsupported("missing bond token".into()))?;
            let start = token.as_ptr() as usize - line.as_ptr() as usize;
            output.push_str(&line[..start]);
            output.push_str(order);
            output.push_str(&line[start + token.len()..]);
            changed += 1;
        } else {
            output.push_str(line);
        }
    }
    if changed != replacements.len() {
        return Err(Mol2KekuleError::Unsupported(
            "bond source mapping mismatch".into(),
        ));
    }
    Ok(output)
}
