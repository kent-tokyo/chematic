//! Additional geometry-derived ring centers shared with MOL2.

use super::*;

fn same_ring_component(mol: &Molecule, ring_bonds: &[bool], start: AtomIdx, goal: AtomIdx) -> bool {
    let mut seen = std::collections::HashSet::from([start]);
    let mut frontier = vec![start];
    while let Some(atom) = frontier.pop() {
        if atom == goal {
            return true;
        }
        for (neighbor, bond) in mol.neighbors(atom) {
            if ring_bonds[bond.0 as usize] && seen.insert(neighbor) {
                frontier.push(neighbor);
            }
        }
    }
    false
}

fn ring_component_labels(mol: &Molecule, ring_bonds: &[bool]) -> Vec<usize> {
    let mut labels = vec![usize::MAX; mol.atom_count()];
    let mut next_label = 0;
    for (start, _) in mol.atoms() {
        if labels[start.0 as usize] != usize::MAX {
            continue;
        }
        let mut frontier = vec![start];
        labels[start.0 as usize] = next_label;
        while let Some(atom) = frontier.pop() {
            for (neighbor, bond) in mol.neighbors(atom) {
                if ring_bonds[bond.0 as usize] && labels[neighbor.0 as usize] == usize::MAX {
                    labels[neighbor.0 as usize] = next_label;
                    frontier.push(neighbor);
                }
            }
        }
        next_label += 1;
    }
    labels
}

/// A single substituted adamantane bridgehead distinguishes a rigid set of
/// three hydrogen-bearing bridgeheads. Require the exact subdivided K4 ring
/// graph: four degree-three cage carbons and six CH2 links, one for every
/// bridgehead pair. The topology gate avoids assigning arbitrary geometry
/// tags to other multi-bridgehead cages.
fn is_monosubstituted_adamantane(
    mol: &Molecule,
    ring_bonds: &[bool],
    hydrogen_bridgeheads: &[AtomIdx],
    substituted_bridgeheads: &[AtomIdx],
) -> bool {
    if hydrogen_bridgeheads.len() != 3 || substituted_bridgeheads.len() != 1 {
        return false;
    }
    let centers: std::collections::HashSet<_> = hydrogen_bridgeheads
        .iter()
        .chain(substituted_bridgeheads)
        .copied()
        .collect();
    let mut links = std::collections::HashSet::new();
    let mut pairs = std::collections::HashSet::new();
    for &center in &centers {
        let ring_neighbors: Vec<_> = mol
            .neighbors(center)
            .filter(|&(_, bond)| ring_bonds[bond.0 as usize])
            .collect();
        if ring_neighbors.len() != 3 {
            return false;
        }
        for (link, bond) in ring_neighbors {
            if mol.bond(bond).order != BondOrder::Single || centers.contains(&link) {
                return false;
            }
            if !links.insert(link) {
                continue;
            }
            let atom = mol.atom(link);
            if atom.element != Element::C || atom.aromatic {
                return false;
            }
            let (ring, outside): (Vec<_>, Vec<_>) = mol
                .neighbors(link)
                .partition(|&(_, bond)| ring_bonds[bond.0 as usize]);
            if ring.len() != 2
                || ring.iter().any(|&(neighbor, bond)| {
                    !centers.contains(&neighbor) || mol.bond(bond).order != BondOrder::Single
                })
                || outside.iter().any(|&(neighbor, bond)| {
                    mol.atom(neighbor).element != Element::H
                        || mol.bond(bond).order != BondOrder::Single
                })
                || outside.len() + usize::from(chematic_core::implicit_hcount(mol, link)) != 2
            {
                return false;
            }
            let a = ring[0].0.0;
            let b = ring[1].0.0;
            if !pairs.insert((a.min(b), a.max(b))) {
                return false;
            }
        }
    }
    links.len() == 6 && pairs.len() == 6
}

fn substituted_adamantane_centers(
    mol: &Molecule,
    ring_bonds: &[bool],
    hydrogen_bridgeheads: &[AtomIdx],
) -> Vec<AtomIdx> {
    if hydrogen_bridgeheads.len() != 3 {
        return Vec::new();
    }
    let substituted_bridgeheads: Vec<_> = mol
        .atoms()
        .filter(|(idx, atom)| {
            if atom.element != Element::C || atom.aromatic || mol.degree(*idx) != 4 {
                return false;
            }
            let (ring, outside): (Vec<_>, Vec<_>) = mol
                .neighbors(*idx)
                .partition(|&(_, bond)| ring_bonds[bond.0 as usize]);
            ring.len() == 3
                && outside.len() == 1
                && mol.atom(outside[0].0).element == Element::C
                && ring
                    .iter()
                    .chain(outside.iter())
                    .all(|&(_, bond)| mol.bond(bond).order == BondOrder::Single)
        })
        .map(|(idx, _)| idx)
        .collect();
    if is_monosubstituted_adamantane(
        mol,
        ring_bonds,
        hydrogen_bridgeheads,
        &substituted_bridgeheads,
    ) {
        hydrogen_bridgeheads
            .iter()
            .chain(&substituted_bridgeheads)
            .copied()
            .collect()
    } else {
        Vec::new()
    }
}

pub(super) fn additional_3d_ring_centers(mol: &Molecule, mol2_extensions: bool) -> Vec<AtomIdx> {
    // RDKit also records the 3D orientation of the two bridgeheads in a
    // simple bicyclic cage such as norbornane, even when CIP regards its
    // paired ring paths as equivalent. Limit this case to exactly two
    // three-ring-bond carbon bridgeheads with one exocyclic substituent;
    // symmetric multi-bridgehead cages such as unsubstituted adamantane must
    // not acquire arbitrary tetrahedral tags.
    let ring_bonds = chematic_perception::ring_bond_flags(mol);
    let mut additional = Vec::new();
    let bridgeheads: Vec<_> = mol
        .atoms()
        .filter(|(idx, atom)| {
            if atom.element != Element::C || atom.aromatic {
                return false;
            }
            let mut ring_neighbors = 0;
            let mut outside = 0;
            for (neighbor, bond) in mol.neighbors(*idx) {
                if mol.bond(bond).order != BondOrder::Single {
                    return false;
                }
                if ring_bonds[bond.0 as usize] {
                    ring_neighbors += 1;
                } else {
                    // A pendant carbon can replace the bridgehead H without
                    // removing the geometry-defined paired configuration.
                    if mol.atom(neighbor).element == Element::H
                        || (mol2_extensions && mol.atom(neighbor).element == Element::C)
                    {
                        outside += 1;
                    } else {
                        return false;
                    }
                }
            }
            ring_neighbors == 3 && outside + chematic_core::implicit_hcount(mol, *idx) == 1
        })
        .map(|(idx, _)| idx)
        .collect();
    if bridgeheads.len() == 2
        && same_ring_component(mol, &ring_bonds, bridgeheads[0], bridgeheads[1])
    {
        additional.extend(bridgeheads.iter().copied());
    }
    additional.extend(substituted_adamantane_centers(
        mol,
        &ring_bonds,
        &bridgeheads,
    ));
    // Two differently substituted quaternary carbons in one ring can
    // distinguish each other's otherwise identical ring arms. Neither
    // center is found by an atom-at-a-time graph comparison before the
    // partner receives its 3D tag. A lone substituted ring carbon remains
    // achiral, so require a connected pair in the cyclic subgraph.
    let paired_ring_centers: Vec<_> = mol
        .atoms()
        .filter(|(idx, atom)| {
            if atom.element != Element::C || atom.aromatic || mol.degree(*idx) != 4 {
                return false;
            }
            let (ring, outside): (Vec<_>, Vec<_>) = mol
                .neighbors(*idx)
                .partition(|&(_, bond)| ring_bonds[bond.0 as usize]);
            ring.len() == 2
                && outside.len() == 2
                && ring
                    .iter()
                    .chain(outside.iter())
                    .all(|&(_, bond)| mol.bond(bond).order == BondOrder::Single)
                && mol.atom(outside[0].0).element != mol.atom(outside[1].0).element
        })
        .map(|(idx, _)| idx)
        .collect();
    if paired_ring_centers.len() == 2
        && same_ring_component(
            mol,
            &ring_bonds,
            paired_ring_centers[0],
            paired_ring_centers[1],
        )
    {
        additional.extend(paired_ring_centers.iter().copied());
    }
    // A substituted ring CH and a quaternary ring carbon with two distinct
    // exocyclic groups likewise define a paired cis/trans configuration.
    // Neither site need be independently stereogenic when its two ring arms
    // are identical. Require a partner in the same ring component so a
    // singly substituted ring cannot acquire a geometry-dependent tag.
    if paired_ring_centers.len() == 1 {
        let labels = ring_component_labels(mol, &ring_bonds);
        let partner_component = labels[paired_ring_centers[0].0 as usize];
        let ring_ch_sites: Vec<_> = ring_ch_candidates(mol, &ring_bonds)
            .into_iter()
            .filter(|&center| {
                center != paired_ring_centers[0] && labels[center.0 as usize] == partner_component
            })
            .collect();
        if !ring_ch_sites.is_empty() {
            additional.push(paired_ring_centers[0]);
            additional.extend(ring_ch_sites);
        }
    }
    if mol2_extensions {
        additional.extend(mol2_ring_ch_sites(mol, &ring_bonds));
    }
    additional
}

/// Two or more substituted ring CH sites can encode relative cis/trans
/// configuration even when each site's ring paths look equivalent before
/// its partners receive 3D tags (e.g. trimethylcyclohexane).
fn mol2_ring_ch_sites(mol: &Molecule, ring_bonds: &[bool]) -> Vec<AtomIdx> {
    let ring_ch_sites = ring_ch_candidates(mol, ring_bonds);
    let labels = ring_component_labels(mol, ring_bonds);
    let mut counts = std::collections::HashMap::new();
    for &center in &ring_ch_sites {
        *counts.entry(labels[center.0 as usize]).or_insert(0usize) += 1;
    }
    ring_ch_sites
        .into_iter()
        .filter(|center| counts[&labels[center.0 as usize]] >= 2)
        .collect()
}

fn ring_ch_candidates(mol: &Molecule, ring_bonds: &[bool]) -> Vec<AtomIdx> {
    mol.atoms()
        .filter(|(idx, atom)| {
            if atom.element != Element::C || atom.aromatic {
                return false;
            }
            let (ring, outside): (Vec<_>, Vec<_>) = mol
                .neighbors(*idx)
                .partition(|&(_, bond)| ring_bonds[bond.0 as usize]);
            ring.len() == 2
                && ring
                    .iter()
                    .chain(outside.iter())
                    .all(|&(_, bond)| mol.bond(bond).order == BondOrder::Single)
                && outside
                    .iter()
                    .filter(|&&(neighbor, _)| mol.atom(neighbor).element == Element::H)
                    .count()
                    + usize::from(chematic_core::implicit_hcount(mol, *idx))
                    == 1
                && outside
                    .iter()
                    .any(|&(neighbor, _)| mol.atom(neighbor).element != Element::H)
        })
        .map(|(idx, _)| idx)
        .collect()
}

#[cfg(test)]
#[path = "stereo3d_ring_partner_tests.rs"]
mod partner_tests;
