//! Bounded local cation resonance, motivated by RDKit #9458.
//!
//! Enumerates `A-C=B+ <-> A+=C-B` for non-aromatic N/O/S termini and a
//! neutral carbon. Neutral terminal valences are N=3, O=2, S=2; the positive
//! terminal has one additional bond-order unit. Hydrogen placement never
//! changes. This is a local electron-pair shift, not tautomer enumeration.
//!
//! The result is exhaustive only for this rule and the supplied graph when
//! its status is `CompleteWithinScope`. General aromatic/Kekulé resonance,
//! anions, charge separation, stereochemistry and energy ranking are outside
//! the scope. Stereo-bearing inputs are rejected instead of silently altered.

use std::collections::HashSet;

use chematic_core::{AtomIdx, BondIdx, BondOrder, Chirality, Element, Molecule, implicit_hcount};

/// Resource bounds, validated at execution even for struct-literal callers.
#[derive(Debug, Clone, Copy)]
pub struct LocalResonanceLimits {
    /// Includes the starting form. Accepted range: 1..=1024; default 32.
    pub max_forms: usize,
    /// Number of applied shifts, including duplicate states. 1..=100000; default 1024.
    pub max_transforms: usize,
    /// Input atom count cap. Accepted range: 1..=1024; default 256.
    pub max_atoms: usize,
}

impl Default for LocalResonanceLimits {
    fn default() -> Self {
        Self {
            max_forms: 32,
            max_transforms: 1024,
            max_atoms: 256,
        }
    }
}

/// Whether the local-rule search exhausted its reachable indexed states.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocalResonanceStatus {
    CompleteWithinScope,
    FormLimitReached,
    TransformLimitReached,
}

/// Search result. Forms keep the input atom/bond indices, isotope labels and
/// atom maps. The starting form is first; subsequent forms follow BFS order.
///
/// Deduplication uses indexed charge/bond-order states. Symmetry-equivalent
/// charge placements on different atoms are retained, as they describe
/// distinct placements in the supplied graph. No canonical SMILES search is
/// needed by this enumerator.
#[derive(Clone)]
pub struct LocalResonanceReport {
    pub forms: Vec<Molecule>,
    pub status: LocalResonanceStatus,
    pub transforms_attempted: usize,
}

/// Invalid configuration or an input outside the supported representation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LocalResonanceError {
    InvalidLimit {
        field: &'static str,
        value: usize,
        maximum: usize,
    },
    EmptyMolecule,
    AtomLimitExceeded {
        actual: usize,
        limit: usize,
    },
    UnsupportedStereo,
    QueryAtom {
        atom: AtomIdx,
    },
    QueryBond {
        bond: BondIdx,
    },
}

impl std::fmt::Display for LocalResonanceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidLimit {
                field,
                value,
                maximum,
            } => write!(f, "{field} must be in 1..={maximum}, got {value}"),
            Self::EmptyMolecule => write!(f, "local resonance requires a non-empty molecule"),
            Self::AtomLimitExceeded { actual, limit } => write!(
                f,
                "local resonance atom count {actual} exceeds limit {limit}"
            ),
            Self::UnsupportedStereo => write!(
                f,
                "local cation resonance does not support stereochemical annotations"
            ),
            Self::QueryAtom { atom } => {
                write!(f, "local resonance does not support query atom {}", atom.0)
            }
            Self::QueryBond { bond } => {
                write!(f, "local resonance does not support query bond {}", bond.0)
            }
        }
    }
}

impl std::error::Error for LocalResonanceError {}

fn validate(mol: &Molecule, limits: &LocalResonanceLimits) -> Result<(), LocalResonanceError> {
    for (field, value, maximum) in [
        ("max_forms", limits.max_forms, 1024),
        ("max_transforms", limits.max_transforms, 100_000),
        ("max_atoms", limits.max_atoms, 1024),
    ] {
        if value == 0 || value > maximum {
            return Err(LocalResonanceError::InvalidLimit {
                field,
                value,
                maximum,
            });
        }
    }
    if mol.atom_count() == 0 {
        return Err(LocalResonanceError::EmptyMolecule);
    }
    if mol.atom_count() > limits.max_atoms {
        return Err(LocalResonanceError::AtomLimitExceeded {
            actual: mol.atom_count(),
            limit: limits.max_atoms,
        });
    }
    if mol.has_bond_directions()
        || !mol.stereo_groups().is_empty()
        || mol
            .atoms()
            .any(|(_, a)| a.chirality != Chirality::None || a.cip_code.is_some())
        || mol
            .bonds()
            .any(|(_, b)| matches!(b.order, BondOrder::Up | BondOrder::Down))
    {
        return Err(LocalResonanceError::UnsupportedStereo);
    }
    if let Some((atom, _)) = mol.atoms().find(|(_, a)| a.wildcard) {
        return Err(LocalResonanceError::QueryAtom { atom });
    }
    if let Some((bond, _)) = mol.bonds().find(|(_, b)| {
        matches!(
            b.order,
            BondOrder::QueryAny
                | BondOrder::QuerySingleOrDouble
                | BondOrder::QuerySingleOrAromatic
                | BondOrder::QueryDoubleOrAromatic
        )
    }) {
        return Err(LocalResonanceError::QueryBond { bond });
    }
    Ok(())
}

fn valence(mol: &Molecule, atom: AtomIdx) -> Option<u16> {
    mol.neighbors(atom)
        .try_fold(u16::from(implicit_hcount(mol, atom)), |sum, (_, bond)| {
            let order = match mol.bond(bond).order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                _ => return None,
            };
            sum.checked_add(order)
        })
}

fn terminal(mol: &Molecule, idx: AtomIdx, charge: i8) -> bool {
    let atom = mol.atom(idx);
    let neutral_valence = match atom.element {
        Element::N => 3,
        Element::O | Element::S => 2,
        _ => return false,
    };
    !atom.aromatic
        && atom.charge == charge
        && valence(mol, idx) == Some(neutral_valence + charge as u16)
}

#[derive(Clone, Copy)]
struct Shift {
    donor: AtomIdx,
    acceptor: AtomIdx,
    single: BondIdx,
    double: BondIdx,
}

fn shifts(mol: &Molecule) -> Vec<Shift> {
    let mut shifts = Vec::new();
    for (center, atom) in mol.atoms() {
        if atom.element != Element::C
            || atom.charge != 0
            || atom.aromatic
            || valence(mol, center) != Some(4)
        {
            continue;
        }
        for (acceptor, double) in mol.neighbors(center) {
            if mol.bond(double).order != BondOrder::Double || !terminal(mol, acceptor, 1) {
                continue;
            }
            for (donor, single) in mol.neighbors(center) {
                if mol.bond(single).order == BondOrder::Single && terminal(mol, donor, 0) {
                    shifts.push(Shift {
                        donor,
                        acceptor,
                        single,
                        double,
                    });
                }
            }
        }
    }
    shifts
}

fn apply(mol: &Molecule, shift: Shift) -> Molecule {
    let mut next = mol.clone();
    // Freeze the original implicit counts before changing bond order/charge:
    // a resonance shift moves electrons, never a proton or an H graph node.
    for atom in [shift.donor, shift.acceptor] {
        next.set_hydrogen_count(atom, Some(implicit_hcount(mol, atom)));
    }
    next.set_bond_order(shift.single, BondOrder::Double);
    next.set_bond_order(shift.double, BondOrder::Single);
    next.set_charge(shift.donor, 1);
    next.set_charge(shift.acceptor, 0);
    next
}

type StateKey = (Vec<i8>, Vec<BondOrder>);

fn key(mol: &Molecule) -> StateKey {
    (
        mol.atoms().map(|(_, a)| a.charge).collect(),
        mol.bonds().map(|(_, b)| b.order).collect(),
    )
}

/// Enumerate indexed forms reachable by local N/O/S–C–N/O/S cation shifts.
///
/// A limit status is returned only when an additional transition/new form
/// would exceed the budget, so a search that finishes exactly at the form
/// cap can still report `CompleteWithinScope`. Unsupported chemistry may have
/// no eligible moves; this status does not certify general resonance coverage
/// or validate unrelated atoms' valences.
pub fn enumerate_local_cation_resonance(
    mol: &Molecule,
    limits: &LocalResonanceLimits,
) -> Result<LocalResonanceReport, LocalResonanceError> {
    validate(mol, limits)?;
    let mut seen = HashSet::from([key(mol)]);
    let mut report = LocalResonanceReport {
        forms: vec![mol.clone()],
        status: LocalResonanceStatus::CompleteWithinScope,
        transforms_attempted: 0,
    };
    let mut cursor = 0;
    while cursor < report.forms.len() {
        for shift in shifts(&report.forms[cursor]) {
            if report.transforms_attempted == limits.max_transforms {
                report.status = LocalResonanceStatus::TransformLimitReached;
                return Ok(report);
            }
            report.transforms_attempted += 1;
            let candidate = apply(&report.forms[cursor], shift);
            let candidate_key = key(&candidate);
            if seen.contains(&candidate_key) {
                continue;
            }
            if report.forms.len() == limits.max_forms {
                report.status = LocalResonanceStatus::FormLimitReached;
                return Ok(report);
            }
            seen.insert(candidate_key);
            report.forms.push(candidate);
        }
        cursor += 1;
    }
    Ok(report)
}
