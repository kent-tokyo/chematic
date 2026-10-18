use std::collections::{HashMap, HashSet, VecDeque};

use chematic_core::{Atom, AtomIdx, BondEntry, BondIdx, BondOrder, Molecule, MoleculeBuilder};

#[derive(Default)]
struct NormalizedGroups {
    nitro: HashSet<AtomIdx>,
    azide: HashSet<AtomIdx>,
    sulfoxide: HashSet<AtomIdx>,
    pyridone_atoms: HashSet<AtomIdx>,
    pyridone_bonds: HashSet<BondIdx>,
    pyridone_ring_atoms: HashSet<AtomIdx>,
    pyridone_ring_bonds: HashSet<BondIdx>,
}

struct PyridoneMatch {
    oxygen: AtomIdx,
    bond: BondIdx,
    ring_atoms: Vec<AtomIdx>,
    ring_bonds: Vec<BondIdx>,
}

/// Normalize chemical groups (nitro groups, etc.).
///
/// Transforms:
/// - `[N+](=O)[O-]` → `N(=O)=O` (nitro normalization: N charge 0, O- → double bond)
/// - `R[S](=O)R` → `R[S+]([O-])R` (RDKit-compatible sulfoxide normalization)
/// - ortho/para pyridinium phenoxide → neutral pyridone
///
/// Returns a new molecule with normalized groups.
pub fn normalize_groups(mol: &Molecule) -> Molecule {
    let groups = detect_groups(mol);
    let (mut builder, remap) = copy_normalized_atoms(mol, &groups);
    copy_normalized_bonds(mol, &groups, &remap, &mut builder);

    // Every atom is copied 1:1 in original index order. The identity remap
    // makes verbatim stereo and direction copying correct here.
    builder.copy_stereo_from(mol);
    builder.copy_bond_directions_from(mol);
    builder.copy_stereo_groups_from(mol);
    builder.build()
}

fn detect_groups(mol: &Molecule) -> NormalizedGroups {
    let mut groups = NormalizedGroups::default();
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() == 7 && atom.charge == 1 {
            detect_nitro(mol, idx, &mut groups.nitro);
            detect_azide(mol, idx, &mut groups.azide);
            if let Some(hit) = detect_pyridone(mol, idx) {
                groups.pyridone_atoms.extend([idx, hit.oxygen]);
                groups.pyridone_bonds.insert(hit.bond);
                groups.pyridone_ring_atoms.extend(hit.ring_atoms);
                groups.pyridone_ring_bonds.extend(hit.ring_bonds);
            }
        }
        if atom.element.atomic_number() == 16 {
            detect_sulfoxide(mol, idx, &mut groups.sulfoxide);
        }
    }
    groups
}

fn copy_normalized_atoms(
    mol: &Molecule,
    groups: &NormalizedGroups,
) -> (MoleculeBuilder, HashMap<AtomIdx, AtomIdx>) {
    let mut builder = MoleculeBuilder::new();
    let mut remap = HashMap::new();
    for (idx, atom) in mol.atoms() {
        let mut normalized = atom.clone();
        normalized.charge = normalized_charge(idx, atom, groups);
        if groups.pyridone_ring_atoms.contains(&idx) {
            normalized.aromatic = true;
        }
        remap.insert(idx, builder.add_atom(normalized));
    }
    (builder, remap)
}

fn normalized_charge(idx: AtomIdx, atom: &Atom, groups: &NormalizedGroups) -> i8 {
    let atomic_number = atom.element.atomic_number();
    if groups.pyridone_atoms.contains(&idx)
        || (groups.nitro.contains(&idx) && matches!(atomic_number, 7 | 8))
        || (groups.azide.contains(&idx) && atomic_number == 7)
    {
        0
    } else if groups.sulfoxide.contains(&idx) {
        match atomic_number {
            16 => 1,
            8 => -1,
            _ => atom.charge,
        }
    } else {
        atom.charge
    }
}

fn copy_normalized_bonds(
    mol: &Molecule,
    groups: &NormalizedGroups,
    remap: &HashMap<AtomIdx, AtomIdx>,
    builder: &mut MoleculeBuilder,
) {
    for (bond_idx, bond) in mol.bonds() {
        let order = normalized_bond_order(mol, bond_idx, bond, groups);
        if let (Some(&atom1), Some(&atom2)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(atom1, atom2, order);
        }
    }
}

fn normalized_bond_order(
    mol: &Molecule,
    bond_idx: BondIdx,
    bond: &BondEntry,
    groups: &NormalizedGroups,
) -> BondOrder {
    if groups.pyridone_bonds.contains(&bond_idx) {
        return BondOrder::Double;
    }
    if groups.pyridone_ring_bonds.contains(&bond_idx) {
        return BondOrder::Aromatic;
    }
    if groups.sulfoxide.contains(&bond.atom1)
        && groups.sulfoxide.contains(&bond.atom2)
        && bond.order == BondOrder::Double
    {
        return BondOrder::Single;
    }
    if groups.azide.contains(&bond.atom1)
        && groups.azide.contains(&bond.atom2)
        && bond.order == BondOrder::Single
        && both_elements_are(mol, bond, 7)
    {
        return BondOrder::Double;
    }
    if groups.nitro.contains(&bond.atom1)
        && groups.nitro.contains(&bond.atom2)
        && bond.order == BondOrder::Single
        && is_negative_oxygen_bond(mol, bond)
    {
        return BondOrder::Double;
    }
    bond.order
}

/// RDKit's Cleanup neutralizes the ortho/para 1,2- and 1,4-pyridinium
/// phenoxide resonance forms. Restrict this to an isolated six-membered
/// aromatic C5N ring with one exocyclic O-; fused or ambiguous rings are left
/// unchanged rather than guessing at a bond/charge assignment.
fn detect_pyridone(mol: &Molecule, nitrogen: AtomIdx) -> Option<PyridoneMatch> {
    if mol.neighbors(nitrogen).count() > 3
        || mol
            .neighbors(nitrogen)
            .filter(|(neighbor, _)| mol.atom(*neighbor).element.atomic_number() == 6)
            .count()
            < 2
        || !mol
            .atoms()
            .any(|(_, atom)| atom.element.atomic_number() == 8 && atom.charge == -1)
    {
        return None;
    }
    let assigned =
        (!mol.atom(nitrogen).aromatic).then(|| chematic_perception::assign_aromaticity(mol));
    let aromatic_atom =
        |idx| mol.atom(idx).aromatic || assigned.as_ref().is_some_and(|m| m.is_atom_aromatic(idx));
    let aromatic_bond = |idx| {
        mol.bond(idx).order == BondOrder::Aromatic
            || assigned.as_ref().is_some_and(|m| m.is_bond_aromatic(idx))
    };
    if !aromatic_atom(nitrogen) {
        return None;
    }
    let mut distances = HashMap::from([(nitrogen, 0usize)]);
    let mut queue = VecDeque::from([nitrogen]);
    while let Some(current) = queue.pop_front() {
        for (neighbor, bond_idx) in mol.neighbors(current) {
            if !aromatic_bond(bond_idx) {
                continue;
            }
            if !distances.contains_key(&neighbor) {
                distances.insert(neighbor, distances[&current] + 1);
                queue.push_back(neighbor);
            }
        }
    }
    if distances.len() != 6
        || distances.iter().any(|(&idx, _)| {
            let atom = mol.atom(idx);
            !aromatic_atom(idx)
                || (idx != nitrogen && atom.element.atomic_number() != 6)
                || mol
                    .neighbors(idx)
                    .filter(|(_, bond)| aromatic_bond(*bond))
                    .count()
                    != 2
        })
    {
        return None;
    }
    let (oxygen, bond) = pyridone_oxygen(mol, &distances)?;
    let ring_bonds = mol
        .bonds()
        .filter_map(|(idx, edge)| {
            (distances.contains_key(&edge.atom1)
                && distances.contains_key(&edge.atom2)
                && aromatic_bond(idx))
            .then_some(idx)
        })
        .collect();
    Some(PyridoneMatch {
        oxygen,
        bond,
        ring_atoms: distances.into_keys().collect(),
        ring_bonds,
    })
}

fn pyridone_oxygen(
    mol: &Molecule,
    distances: &HashMap<AtomIdx, usize>,
) -> Option<(AtomIdx, BondIdx)> {
    let mut oxygen = None;
    for (&carbon, &distance) in distances {
        if distance != 1 && distance != 3 {
            continue;
        }
        for (neighbor, bond_idx) in mol.neighbors(carbon) {
            let atom = mol.atom(neighbor);
            if atom.element.atomic_number() == 8
                && atom.charge == -1
                && mol.bond(bond_idx).order == BondOrder::Single
                && mol.neighbors(neighbor).count() == 1
                && oxygen.replace((neighbor, bond_idx)).is_some()
            {
                return None;
            }
        }
    }
    oxygen
}

fn both_elements_are(mol: &Molecule, bond: &BondEntry, atomic_number: u8) -> bool {
    mol.atom(bond.atom1).element.atomic_number() == atomic_number
        && mol.atom(bond.atom2).element.atomic_number() == atomic_number
}

fn is_negative_oxygen_bond(mol: &Molecule, bond: &BondEntry) -> bool {
    [bond.atom1, bond.atom2]
        .into_iter()
        .any(|idx| mol.atom(idx).element.atomic_number() == 8 && mol.atom(idx).charge == -1)
}

/// Mark nitro [N+](=O)[O-] atoms. Aromatic N-oxides retain their
/// charge-separated N+-O- bond; turning it into N=O breaks ring valence.
fn detect_nitro(mol: &Molecule, idx: AtomIdx, nitro_atoms: &mut HashSet<AtomIdx>) {
    let oxygen_neighbors: Vec<_> = mol
        .neighbors(idx)
        .filter(|(neighbor, _)| mol.atom(*neighbor).element.atomic_number() == 8)
        .collect();

    if oxygen_neighbors.len() == 2 {
        let mut has_double_oxygen = false;
        let mut has_single_negative_oxygen = false;
        for (oxygen_idx, bond_idx) in &oxygen_neighbors {
            let oxygen = mol.atom(*oxygen_idx);
            let bond = mol.bond(*bond_idx);
            has_double_oxygen |= bond.order == BondOrder::Double && oxygen.charge == 0;
            if bond.order == BondOrder::Single && oxygen.charge == -1 {
                has_single_negative_oxygen = true;
                nitro_atoms.insert(*oxygen_idx);
            }
        }
        if has_double_oxygen && has_single_negative_oxygen {
            nitro_atoms.insert(idx);
        }
    }
}

/// Mark azide [N-][N+]#N atoms.
fn detect_azide(mol: &Molecule, idx: AtomIdx, azide_atoms: &mut HashSet<AtomIdx>) {
    let nitrogen_neighbors: Vec<_> = mol
        .neighbors(idx)
        .filter(|(neighbor, _)| mol.atom(*neighbor).element.atomic_number() == 7)
        .collect();

    for (nitrogen_idx, bond_idx) in &nitrogen_neighbors {
        let nitrogen = mol.atom(*nitrogen_idx);
        if mol.bond(*bond_idx).order != BondOrder::Triple || nitrogen.charge != 0 {
            continue;
        }
        for (other_idx, other_bond_idx) in &nitrogen_neighbors {
            if other_idx != nitrogen_idx
                && mol.bond(*other_bond_idx).order == BondOrder::Single
                && mol.atom(*other_idx).charge == -1
            {
                azide_atoms.extend([idx, *nitrogen_idx, *other_idx]);
            }
        }
    }
}

/// Mark neutral three-coordinate sulfoxide S=O atoms, but not sulfones.
fn detect_sulfoxide(mol: &Molecule, idx: AtomIdx, sulfoxide_atoms: &mut HashSet<AtomIdx>) {
    let sulfur = mol.atom(idx);
    if sulfur.charge != 0 || mol.neighbors(idx).count() != 3 {
        return;
    }
    let double_oxygens: Vec<_> = mol
        .neighbors(idx)
        .filter(|(oxygen_idx, bond_idx)| {
            let oxygen = mol.atom(*oxygen_idx);
            oxygen.element.atomic_number() == 8
                && oxygen.charge == 0
                && mol.bond(*bond_idx).order == BondOrder::Double
        })
        .map(|(oxygen_idx, _)| oxygen_idx)
        .collect();
    if let [oxygen] = double_oxygens.as_slice() {
        sulfoxide_atoms.extend([idx, *oxygen]);
    }
}
