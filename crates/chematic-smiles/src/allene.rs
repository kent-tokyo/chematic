use chematic_core::{AtomIdx, BondOrder, Molecule};

/// One side of an allene/cumulene axis.
#[derive(Debug, Clone, Copy)]
pub(crate) struct AlleneArm {
    pub terminal: AtomIdx,
    pub axis_neighbor: AtomIdx,
}

/// Follow both double-bond arms from `center` to their terminal atoms.
///
/// `axis_order` is the order in which the two bonds at the tagged center were
/// encountered or will be written.  Equal arm lengths are required by the
/// OpenSMILES allene class.  Cycles and non-double axis bonds are rejected.
pub(crate) fn allene_arms(
    mol: &Molecule,
    center: AtomIdx,
    axis_order: &[u32],
) -> Option<[AlleneArm; 2]> {
    let mut first_steps = Vec::with_capacity(2);
    for &raw in axis_order {
        let atom = AtomIdx(raw);
        if mol
            .neighbors(center)
            .any(|(neighbor, bond)| neighbor == atom && mol.bond(bond).order == BondOrder::Double)
            && !first_steps.contains(&atom)
        {
            first_steps.push(atom);
        }
    }
    if first_steps.len() != 2 || mol.degree(center) != 2 {
        return None;
    }

    let mut lengths = [0usize; 2];
    let mut arms = [
        AlleneArm {
            terminal: center,
            axis_neighbor: center,
        },
        AlleneArm {
            terminal: center,
            axis_neighbor: center,
        },
    ];
    for (slot, &first) in first_steps.iter().enumerate() {
        let mut previous = center;
        let mut current = first;
        let mut length = 1usize;
        let mut seen = vec![false; mol.atom_count()];
        seen[center.0 as usize] = true;
        while mol.degree(current) == 2 {
            let mut next = None;
            for (neighbor, bond) in mol.neighbors(current) {
                if neighbor != previous && mol.bond(bond).order == BondOrder::Double {
                    next = Some(neighbor);
                    break;
                }
            }
            let Some(next_atom) = next else { break };
            if seen[next_atom.0 as usize] {
                return None;
            }
            seen[current.0 as usize] = true;
            previous = current;
            current = next_atom;
            length += 1;
        }
        if !(2..=3).contains(&mol.degree(current)) {
            return None;
        }
        for (neighbor, bond) in mol.neighbors(current) {
            if neighbor != previous && mol.bond(bond).order != BondOrder::Single {
                return None;
            }
        }
        lengths[slot] = length;
        arms[slot] = AlleneArm {
            terminal: current,
            axis_neighbor: previous,
        };
    }
    (lengths[0] == lengths[1]).then_some(arms)
}

/// Whether two encoded allene terminal orders differ by one terminal-ligand
/// swap. Each order is `[terminal, ligand_a, ligand_b]` repeated twice.
/// Swapping the two terminal blocks preserves axial configuration; swapping
/// the ligand pair at exactly one end inverts it.
pub(crate) fn allene_order_is_inverted(original: &[u32], written: &[u32]) -> Option<bool> {
    let original: &[u32; 6] = original.try_into().ok()?;
    let written: &[u32; 6] = written.try_into().ok()?;
    let mut inverted = false;
    for written_offset in [0, 3] {
        let original_offset = if original[0] == written[written_offset] {
            0
        } else if original[3] == written[written_offset] {
            3
        } else {
            return None;
        };
        let original_ligands = &original[original_offset + 1..original_offset + 3];
        let written_ligands = &written[written_offset + 1..written_offset + 3];
        if original_ligands == written_ligands {
            continue;
        }
        if original_ligands[0] == written_ligands[1] && original_ligands[1] == written_ligands[0] {
            inverted = !inverted;
        } else {
            return None;
        }
    }
    Some(inverted)
}
