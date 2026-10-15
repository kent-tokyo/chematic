//! Checked full alternative of a relative tetrahedral/geometric stereo group.
//! Uses the issue #9705 proposal; no canonicalization or axial assignment.
use chematic_core::ez_markers::{read_ez_facts, write_ez_facts};
use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Chirality, Molecule, STEREO_H_SENTINEL, StereoGroupError,
    StereoGroupKind,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StereoGroupInversionError {
    GroupIndex { index: usize, groups: usize },
    AbsoluteGroup,
    InvalidGroups(StereoGroupError),
    UnsupportedAtom(AtomIdx),
    IncompleteFrame(AtomIdx),
    UnsupportedBond(BondIdx),
    UnspecifiedBond(BondIdx),
    SmallRingBond(BondIdx),
    ConflictingDirections,
    UnrepresentableDirections,
}
impl core::fmt::Display for StereoGroupInversionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::GroupIndex { index, groups } => {
                write!(f, "stereo group index {index} is outside {groups} groups")
            }
            Self::AbsoluteGroup => {
                write!(f, "an absolute stereo group has no relative alternative")
            }
            Self::InvalidGroups(e) => e.fmt(f),
            Self::UnsupportedAtom(a) => write!(
                f,
                "group atom {} requires specified tetrahedral chirality",
                a.0
            ),
            Self::IncompleteFrame(a) => write!(
                f,
                "group atom {} requires a complete four-ligand stereo frame",
                a.0
            ),
            Self::UnsupportedBond(b) => write!(
                f,
                "group bond {} requires a supported non-aromatic double bond; axes are unsupported",
                b.0
            ),
            Self::UnspecifiedBond(b) => write!(
                f,
                "group bond {} has no fully specified geometric configuration",
                b.0
            ),
            Self::SmallRingBond(b) => write!(
                f,
                "group bond {} belongs to a ring smaller than eight atoms",
                b.0
            ),
            Self::ConflictingDirections => write!(
                f,
                "conflicting geometric bond directions cannot be inverted safely"
            ),
            Self::UnrepresentableDirections => write!(
                f,
                "geometric directions cannot be represented without losing configuration"
            ),
        }
    }
}
impl std::error::Error for StereoGroupInversionError {}

fn tetrahedral_frame(mol: &Molecule, idx: AtomIdx) -> Result<(), StereoGroupInversionError> {
    if !matches!(
        mol.atom(idx).chirality,
        Chirality::Clockwise | Chirality::CounterClockwise
    ) {
        return Err(StereoGroupInversionError::UnsupportedAtom(idx));
    }
    let mut expected: Vec<u32> = mol.neighbors(idx).map(|(a, _)| a.0).collect();
    expected.extend(std::iter::repeat_n(
        STEREO_H_SENTINEL,
        chematic_core::implicit_hcount(mol, idx) as usize,
    ));
    let mut frame = mol.stereo_neighbor_order(idx).unwrap_or_default().to_vec();
    expected.sort_unstable();
    frame.sort_unstable();
    if expected.len() != 4 || frame != expected {
        return Err(StereoGroupInversionError::IncompleteFrame(idx));
    }
    Ok(())
}
fn small_ring(mol: &Molecule, idx: BondIdx) -> bool {
    let bond = mol.bond(idx);
    let mut seen = vec![false; mol.atom_count()];
    seen[bond.atom1.0 as usize] = true;
    let mut queue = std::collections::VecDeque::from([(bond.atom1, 0)]);
    while let Some((atom, distance)) = queue.pop_front() {
        if distance == 6 {
            continue;
        }
        for (next, b) in mol.neighbors(atom) {
            if b == idx {
                continue;
            }
            if next == bond.atom2 {
                return true;
            }
            if !seen[next.0 as usize] {
                seen[next.0 as usize] = true;
                queue.push_back((next, distance + 1));
            }
        }
    }
    false
}
fn invert_bonds(
    mol: &Molecule,
    result: &mut Molecule,
    bonds: &[BondIdx],
) -> Result<(), StereoGroupInversionError> {
    if bonds.is_empty() {
        return Ok(());
    }
    let Some((mut facts, false)) = read_ez_facts(mol) else {
        return Err(StereoGroupInversionError::ConflictingDirections);
    };
    for &idx in bonds {
        let bond = mol.bond(idx);
        if bond.order != BondOrder::Double
            || mol.atom(bond.atom1).aromatic && mol.atom(bond.atom2).aromatic
        {
            return Err(StereoGroupInversionError::UnsupportedBond(idx));
        }
        if small_ring(mol, idx) {
            return Err(StereoGroupInversionError::SmallRingBond(idx));
        }
        let fact = facts
            .iter_mut()
            .find(|fact| fact.double == idx)
            .ok_or(StereoGroupInversionError::UnspecifiedBond(idx))?;
        fact.ends[0].up = !fact.ends[0].up;
    }
    // Preserve every other geometric fact, including bonds sharing a marker.
    // The existing parity solver performs a joint assignment of all markers.
    if !write_ez_facts(result, &facts) {
        return Err(StereoGroupInversionError::UnrepresentableDirections);
    }
    Ok(())
}

/// Return the full drawn alternative of one AND/OR group: flip every declared
/// tetrahedral member and geometric double-bond member together. Preserve group
/// membership and all unselected configurations. The indexed graph and stored
/// frames remain intact; CIP annotations are cleared for recomputation.
///
/// Reject ABS, incomplete/unspecified frames, axes/non-tetrahedral atoms, small
/// ring double bonds and unsupported direction representations. No coordinates
/// are generated, and neither meso reduction nor canonicalization is performed.
pub fn invert_stereo_group(
    mol: &Molecule,
    index: usize,
) -> Result<Molecule, StereoGroupInversionError> {
    let group = mol
        .stereo_groups()
        .get(index)
        .ok_or(StereoGroupInversionError::GroupIndex {
            index,
            groups: mol.stereo_groups().len(),
        })?;
    if group.kind == StereoGroupKind::Absolute {
        return Err(StereoGroupInversionError::AbsoluteGroup);
    }
    let groups = mol
        .stereo_groups()
        .iter()
        .enumerate()
        .map(|(i, g)| {
            (
                g.clone(),
                mol.stereo_group_bond_indices(i)
                    .unwrap_or_default()
                    .to_vec(),
            )
        })
        .collect();
    let mut result = mol
        .with_checked_stereo_groups(groups)
        .map_err(StereoGroupInversionError::InvalidGroups)?;
    // Validation normalizes constructor inputs; preserve stored metadata verbatim.
    result.copy_stereo_groups_from(mol);
    for &idx in &group.atom_indices {
        tetrahedral_frame(mol, idx)?;
    }
    let bonds = mol.stereo_group_bond_indices(index).unwrap_or_default();
    invert_bonds(mol, &mut result, bonds)?;
    for &idx in &group.atom_indices {
        let inverted = match mol.atom(idx).chirality {
            Chirality::Clockwise => Chirality::CounterClockwise,
            _ => Chirality::Clockwise,
        };
        result.set_chirality(idx, inverted);
    }
    // A changed neighbouring assignment can affect CIP Rule 5 priorities.
    for (idx, atom) in mol.atoms() {
        if atom.cip_code.is_some() {
            result.set_cip_code(idx, None);
        }
    }
    Ok(result)
}
