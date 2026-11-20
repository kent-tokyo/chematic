//! Accurate CIP priority at the two ends of an ordinary double bond.

use super::{
    AtomIdx, BondIdx, BondOrder, CipCode, Element, LazyEzRanker, Molecule,
    double_bond_in_small_ring, end_has_marked_substituent, implicit_hcount, substituent_is_up,
};

/// Absolute or pseudoasymmetric CIP descriptor for a configured double bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DoubleBondCip {
    E,
    Z,
    SeqCis,
    SeqTrans,
}

impl DoubleBondCip {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::E => "E",
            Self::Z => "Z",
            Self::SeqCis => "seqCis",
            Self::SeqTrans => "seqTrans",
        }
    }

    /// The atom-keyed CIP API uses the older absolute E/Z representation.
    pub const fn absolute(self) -> Option<CipCode> {
        match self {
            Self::E => Some(CipCode::E),
            Self::Z => Some(CipCode::Z),
            Self::SeqCis | Self::SeqTrans => None,
        }
    }
}

/// Accurate CIP assignments keyed by the configured double bond.
pub fn assign_accurate_double_bond_cip(mol: &Molecule) -> Vec<(BondIdx, DoubleBondCip)> {
    let ranker = LazyEzRanker::new(mol, chematic_cip::CipBudget::default_budget());
    mol.bonds()
        .filter_map(|(bond, _)| {
            assign_configured_double_bond(mol, bond, &ranker).map(|code| (bond, code))
        })
        .collect()
}

pub(super) fn absolute_bond_labels(
    labels: Vec<(BondIdx, DoubleBondCip)>,
) -> Vec<(BondIdx, CipCode)> {
    labels
        .into_iter()
        .filter_map(|(bond, code)| Some((bond, code.absolute()?)))
        .collect()
}

pub(super) fn assign_configured_double_bond(
    mol: &Molecule,
    bond_idx: BondIdx,
    ranker: &LazyEzRanker<'_>,
) -> Option<DoubleBondCip> {
    let bond = mol.bond(bond_idx);
    // A ring of fewer than eight atoms holds its double bond cis: no E/Z
    // (RDKit drops such stereo too).
    if bond.order != BondOrder::Double
        || !end_has_marked_substituent(mol, bond.atom1, bond.atom2)
        || !end_has_marked_substituent(mol, bond.atom2, bond.atom1)
        || double_bond_in_small_ring(mol, bond_idx)
    {
        return None;
    }
    let a1 = bond.atom1;
    let a2 = bond.atom2;
    let subs = |end: AtomIdx, other: AtomIdx| -> Vec<AtomIdx> {
        mol.neighbors(end)
            .filter(|&(nb, bidx)| nb != other && mol.bond(bidx).order != BondOrder::Double)
            .map(|(nb, _)| nb)
            .collect()
    };
    let subs_a1 = subs(a1, a2);
    let subs_a2 = subs(a2, a1);
    if subs_a1.is_empty() || subs_a2.is_empty() {
        return None;
    }
    let (up_a1, pseudo_a1) = highest_stereo_sub_accurate(mol, a1, &subs_a1, ranker, true, |sub| {
        substituent_is_up(mol, a1, sub)
    })?;
    let (up_a2, pseudo_a2) = highest_stereo_sub_accurate(mol, a2, &subs_a2, ranker, true, |sub| {
        substituent_is_up(mol, a2, sub)
    })?;
    Some(match (up_a1 == up_a2, pseudo_a1 != pseudo_a2) {
        (true, false) => DoubleBondCip::Z,
        (false, false) => DoubleBondCip::E,
        (true, true) => DoubleBondCip::SeqCis,
        (false, true) => DoubleBondCip::SeqTrans,
    })
}

/// Direction of the higher-priority ligand and whether priority is pseudoasymmetric.
pub(super) fn highest_stereo_sub_accurate(
    mol: &Molecule,
    alkene_end: AtomIdx,
    subs: &[AtomIdx],
    ranker: &LazyEzRanker<'_>,
    configured: bool,
    side: impl Fn(AtomIdx) -> Option<bool>,
) -> Option<(bool, bool)> {
    // No marker on either substituent: nothing to read (see
    // `highest_stereo_sub`).
    if subs.iter().all(|&sub| side(sub).is_none()) {
        return None;
    }
    let (top, other, pseudo) = match subs {
        [only] => {
            let up = side(*only)?;
            if mol.atom(*only).element != Element::H || implicit_hcount(mol, alkene_end) == 0 {
                return Some((up, false));
            }
            // The natural-abundance implicit H can outrank explicit protium.
            return match ranker
                .get()
                .compare_to_implicit_hydrogen(mol, alkene_end, *only)
                .ok()??
            {
                std::cmp::Ordering::Greater => Some((up, false)),
                std::cmp::Ordering::Less => Some((!up, false)),
                std::cmp::Ordering::Equal => None,
            };
        }
        [a, b] => {
            let (order, pseudo) = if configured {
                ranker
                    .get()
                    .compare_configured(mol, alkene_end, *a, *b)
                    .ok()??
            } else {
                (ranker.get().compare(mol, alkene_end, *a, *b).ok()??, false)
            };
            match order {
                std::cmp::Ordering::Greater => (*a, *b, pseudo),
                std::cmp::Ordering::Less => (*b, *a, pseudo),
                std::cmp::Ordering::Equal => return None,
            }
        }
        _ => return None,
    };
    let up = side(top).or_else(|| side(other).map(|up| !up))?;
    Some((up, pseudo))
}
