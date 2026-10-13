use chematic_core::{AtomIdx, BondOrder, Chirality, Molecule, STEREO_H_SENTINEL, implicit_hcount};

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

/// Normalize parsed central allene tags into their four virtual terminal
/// ligands. Component starts are recovered from the completed graph: SMILES
/// parsing appends each dot-separated component, so its lowest atom index is
/// the atom that had no preceding neighbour in parse order.
pub(crate) fn normalize_allene_stereo(molecule: &mut Molecule, tagged_atoms: &[AtomIdx]) {
    // Ordinary tetrahedral stereocentres are common; actual allene/cumulene
    // centres are not. Select the only possible central geometry before doing
    // the component walk needed to reconstruct terminal parse order. This
    // keeps allene support off the hot path for ordinary chiral SMILES.
    let candidates: Vec<_> = tagged_atoms
        .iter()
        .copied()
        .filter(|&atom| {
            molecule.atom(atom).chirality.is_tetrahedral()
                && molecule.degree(atom) == 2
                && molecule
                    .neighbors(atom)
                    .filter(|(_, bond)| molecule.bond(*bond).order == BondOrder::Double)
                    .count()
                    == 2
        })
        .collect();
    if candidates.is_empty() {
        return;
    }
    let fragment_starts = fragment_start_flags(molecule);

    for center in candidates {
        let stored = molecule.atom(center).chirality;
        let direct = molecule
            .stereo_neighbor_order(center)
            .map(<[u32]>::to_vec)
            .unwrap_or_else(|| molecule.neighbors(center).map(|(n, _)| n.0).collect());
        let arms = allene_arms(molecule, center, &direct);
        let Some(arms) = arms else {
            // Valid class on incompatible geometry is cleaned; malformed class
            // numbers are rejected earlier by the token parser.
            molecule.set_chirality(center, Chirality::None);
            continue;
        };

        let mut virtual_order = Vec::with_capacity(6);
        for arm in arms {
            let mut incident: Vec<_> = molecule.neighbors(arm.terminal).collect();
            incident.sort_by_key(|(_, bond)| bond.0);
            let parent = (!fragment_starts[arm.terminal.0 as usize])
                .then(|| incident.first().map(|(atom, _)| *atom))
                .flatten();
            let mut local = Vec::with_capacity(3);
            if let Some(atom) = parent {
                local.push(atom.0);
            }
            if implicit_hcount(molecule, arm.terminal) > 0 {
                local.push(STEREO_H_SENTINEL);
            }
            local.extend(
                incident
                    .into_iter()
                    .filter(|(neighbor, _)| Some(*neighbor) != parent)
                    .map(|(neighbor, _)| neighbor.0),
            );
            local.retain(|&atom| atom != arm.axis_neighbor.0);
            if local.len() != 2 {
                virtual_order.clear();
                break;
            }
            virtual_order.push(arm.terminal.0);
            virtual_order.extend(local);
        }
        if virtual_order.len() == 6 {
            molecule.set_chirality(center, stored);
            molecule.set_stereo_neighbor_order(center, virtual_order);
        } else {
            molecule.set_chirality(center, Chirality::None);
        }
    }
}

fn fragment_start_flags(molecule: &Molecule) -> Vec<bool> {
    let mut visited = vec![false; molecule.atom_count()];
    let mut starts = vec![false; molecule.atom_count()];
    for raw in 0..molecule.atom_count() {
        if visited[raw] {
            continue;
        }
        starts[raw] = true;
        visited[raw] = true;
        let mut stack = vec![AtomIdx(raw as u32)];
        while let Some(atom) = stack.pop() {
            for (neighbor, _) in molecule.neighbors(atom) {
                let slot = neighbor.0 as usize;
                if !visited[slot] {
                    visited[slot] = true;
                    stack.push(neighbor);
                }
            }
        }
    }
    starts
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
