use std::collections::{HashMap, HashSet};

use chematic_core::{Atom, AtomIdx, BondEntry, BondOrder, Molecule, MoleculeBuilder};

#[derive(Default)]
struct NormalizedGroups {
    nitro: HashSet<AtomIdx>,
    oxide: HashSet<AtomIdx>,
    azide: HashSet<AtomIdx>,
    sulfoxide: HashSet<AtomIdx>,
}

/// Normalize chemical groups (nitro groups, etc.).
///
/// Transforms:
/// - `[N+](=O)[O-]` → `N(=O)=O` (nitro normalization: N charge 0, O- → double bond)
/// - `R[S](=O)R` → `R[S+]([O-])R` (RDKit-compatible sulfoxide normalization)
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
            detect_nitro(mol, idx, atom, &mut groups.nitro, &mut groups.oxide);
            detect_azide(mol, idx, &mut groups.azide);
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
        remap.insert(idx, builder.add_atom(normalized));
    }
    (builder, remap)
}

fn normalized_charge(idx: AtomIdx, atom: &Atom, groups: &NormalizedGroups) -> i8 {
    let atomic_number = atom.element.atomic_number();
    if (groups.nitro.contains(&idx) && matches!(atomic_number, 7 | 8))
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
    for (_, bond) in mol.bonds() {
        let order = normalized_bond_order(mol, bond, groups);
        if let (Some(&atom1), Some(&atom2)) = (remap.get(&bond.atom1), remap.get(&bond.atom2)) {
            let _ = builder.add_bond(atom1, atom2, order);
        }
    }
}

fn normalized_bond_order(mol: &Molecule, bond: &BondEntry, groups: &NormalizedGroups) -> BondOrder {
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

fn both_elements_are(mol: &Molecule, bond: &BondEntry, atomic_number: u8) -> bool {
    mol.atom(bond.atom1).element.atomic_number() == atomic_number
        && mol.atom(bond.atom2).element.atomic_number() == atomic_number
}

fn is_negative_oxygen_bond(mol: &Molecule, bond: &BondEntry) -> bool {
    [bond.atom1, bond.atom2]
        .into_iter()
        .any(|idx| mol.atom(idx).element.atomic_number() == 8 && mol.atom(idx).charge == -1)
}

/// Mark nitro [N+](=O)[O-] and aromatic N-oxide atoms.
fn detect_nitro(
    mol: &Molecule,
    idx: AtomIdx,
    atom: &Atom,
    nitro_atoms: &mut HashSet<AtomIdx>,
    oxide_atoms: &mut HashSet<AtomIdx>,
) {
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
    } else if let Some((oxygen_idx, bond_idx)) = oxygen_neighbors.first() {
        let oxygen = mol.atom(*oxygen_idx);
        let bond = mol.bond(*bond_idx);
        if atom.aromatic && bond.order == BondOrder::Single && oxygen.charge == -1 {
            nitro_atoms.insert(idx);
            oxide_atoms.insert(*oxygen_idx);
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
