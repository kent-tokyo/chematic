//! Bounded acid-priority reionization without independent ionization.

use chematic_core::{AtomIdx, BondOrder, Molecule, implicit_hcount};

use super::clone_molecule;

/// Larger values indicate stronger acids within the supported structural
/// classes. This is deliberately a partial ordering: unfamiliar acid groups
/// cause an abstention rather than a possibly incorrect transfer.
fn oxygen_acid_priority(mol: &Molecule, oxygen: AtomIdx) -> Option<u8> {
    let (neighbor, bond) = mol.neighbors(oxygen).next()?;
    if mol.degree(oxygen) != 1 || mol.bond(bond).order != BondOrder::Single {
        return None;
    }
    let double_oxygens = mol
        .neighbors(neighbor)
        .filter(|&(other, edge)| {
            other != oxygen
                && mol.atom(other).element.atomic_number() == 8
                && mol.bond(edge).order == BondOrder::Double
        })
        .count();
    match mol.atom(neighbor).element.atomic_number() {
        16 if double_oxygens >= 2 => Some(4), // sulfonic acid
        6 if double_oxygens >= 1 => Some(3),  // carboxylic acid
        6 if mol.atom(neighbor).aromatic => Some(2),
        6 => Some(1), // alcohol
        _ => None,
    }
}

fn ionized_priority(mol: &Molecule, idx: AtomIdx) -> Option<u8> {
    let atom = mol.atom(idx);
    if atom.charge != -1 {
        return None;
    }
    match atom.element.atomic_number() {
        8 if atom.hydrogen_count.unwrap_or(0) == 0 => oxygen_acid_priority(mol, idx),
        6 if !atom.aromatic
            && matches!(
                (mol.degree(idx), atom.hydrogen_count),
                (0, Some(3)) | (1, Some(2))
            )
            && mol
                .neighbors(idx)
                .all(|(_, bond)| mol.bond(bond).order == BondOrder::Single) =>
        {
            Some(0)
        }
        _ => None,
    }
}

fn strongest_donor(mol: &Molecule, minimum_priority: u8) -> Option<AtomIdx> {
    let mut donors = Vec::new();
    for (idx, atom) in mol.atoms() {
        if atom.element.atomic_number() != 8 || atom.charge != 0 || implicit_hcount(mol, idx) != 1 {
            continue;
        }
        let Some(priority) = oxygen_acid_priority(mol, idx) else {
            // An unrecognized O-H group may be stronger than the known ones.
            return None;
        };
        if priority > minimum_priority {
            donors.push((idx, priority));
        }
    }
    let strongest = donors.iter().map(|&(_, priority)| priority).max()?;
    let strongest_donors = donors
        .into_iter()
        .filter_map(|(idx, priority)| (priority == strongest).then_some(idx))
        .collect::<Vec<_>>();
    let [donor] = strongest_donors.as_slice() else {
        return None;
    };
    Some(*donor)
}

pub(super) fn reionize(mol: &Molecule) -> Molecule {
    let mut corrected = clone_molecule(mol);
    let start_charge = mol
        .atoms()
        .map(|(_, atom)| i32::from(atom.charge))
        .sum::<i32>();
    let mut correction = 0i32;
    for (idx, atom) in mol.atoms() {
        // RDKit's X0 charge-correction patterns exclude both explicit and
        // implicit hydrogens. In particular `Cl` is HCl while `[Cl]` is bare
        // chlorine, and `[NaH]` is not an isolated sodium atom.
        if atom.charge != 0 || mol.degree(idx) != 0 || implicit_hcount(mol, idx) != 0 {
            continue;
        }
        let charge = match atom.element.atomic_number() {
            3 | 11 | 19 => 1, // isolated Li, Na, K
            12 | 20 => 2,     // isolated Mg, Ca
            17 => -1,         // isolated Cl
            _ => continue,
        };
        corrected.set_charge(idx, charge);
        correction += i32::from(charge);
    }
    if start_charge + correction != 0 {
        for _ in 0..correction.max(0) {
            let Some(donor) = strongest_donor(&corrected, 0) else {
                break;
            };
            corrected.set_charge(donor, -1);
            corrected.set_hydrogen_count(donor, Some(0));
        }
    }

    let mol = &corrected;
    let anions = mol
        .atoms()
        .filter_map(|(idx, atom)| (atom.charge < 0).then_some(idx))
        .collect::<Vec<_>>();
    let [anion] = anions.as_slice() else {
        return clone_molecule(mol);
    };
    let Some(anion_priority) = ionized_priority(mol, *anion) else {
        return clone_molecule(mol);
    };

    let Some(donor) = strongest_donor(mol, anion_priority) else {
        return clone_molecule(mol);
    };

    let mut output = clone_molecule(mol);
    output.set_charge(donor, -1);
    output.set_hydrogen_count(donor, Some(0));
    output.set_charge(*anion, 0);
    output.set_hydrogen_count(
        *anion,
        Some(
            mol.atom(*anion)
                .hydrogen_count
                .unwrap_or(0)
                .saturating_add(1),
        ),
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use chematic_smiles::{canonical_smiles, parse};

    fn output(input: &str) -> String {
        canonical_smiles(&reionize(&parse(input).unwrap()))
    }

    fn expected(input: &str) -> String {
        canonical_smiles(&parse(input).unwrap())
    }

    #[test]
    fn neutral_acids_and_bases_do_not_independently_ionize() {
        for input in [
            "CC(=O)O",
            "Oc1ccccc1",
            "CC(N)C",
            "OC(=O)CCN",
            "[NH3+]CC(=O)O",
            "[Na+].CC(=O)O",
        ] {
            assert_eq!(output(input), expected(input), "{input}");
        }
    }

    #[test]
    fn charged_acid_base_pairs_follow_acid_priority() {
        for (input, product) in [
            ("CC(=O)O.[O-]c1ccccc1", "CC(=O)[O-].Oc1ccccc1"),
            ("CC(=O)O.CC(C)[O-]", "CC(=O)[O-].CC(C)O"),
            ("CS(=O)(=O)O.CC(=O)[O-]", "CS(=O)(=O)[O-].CC(=O)O"),
            ("CC(=O)O.[CH3-]", "CC(=O)[O-].C"),
        ] {
            let actual = output(input);
            assert_eq!(actual, expected(product), "{input}");
            assert_eq!(output(&actual), actual, "must be idempotent: {input}");
        }
    }

    #[test]
    fn weaker_or_ambiguous_donors_do_not_move_protons() {
        for input in [
            "Oc1ccccc1.CC(=O)[O-]",
            "CC(=O)O.CC(=O)O.[CH3-]",
            "O=P(O)(O)C.CC(=O)[O-]",
            "[O-]C(=O)C.[NH3+]C",
        ] {
            assert_eq!(output(input), expected(input), "{input}");
        }
    }

    #[test]
    fn isolated_inorganic_charge_corrections_balance_known_acids() {
        for (input, product) in [
            ("[Li]", "[Li+]"),
            ("[Na]", "[Na+]"),
            ("[K]", "[K+]"),
            ("[Mg]", "[Mg+2]"),
            ("[Ca]", "[Ca+2]"),
            ("[Cl]", "[Cl-]"),
            ("[Na].CC(=O)O", "[Na+].CC(=O)[O-]"),
            ("[K].Oc1ccccc1", "[K+].[O-]c1ccccc1"),
            ("[Na].CC(=O)[O-]", "[Na+].CC(=O)[O-]"),
        ] {
            assert_eq!(output(input), expected(product), "{input}");
        }
        for unchanged in ["Cl", "Cl.Cl.N", "[NaH]", "[KH]"] {
            assert_eq!(output(unchanged), expected(unchanged), "{unchanged}");
        }
    }

    #[test]
    fn proton_transfer_preserves_graph_charge_maps_and_stereo_frame() {
        let source = parse("C[C@H](F)C(=O)[OH:7].[CH3-]").unwrap();
        let result = reionize(&source);
        assert_eq!(source.atom_count(), result.atom_count());
        assert_eq!(source.bond_count(), result.bond_count());
        assert_eq!(
            source
                .atoms()
                .map(|(_, atom)| i32::from(atom.charge))
                .sum::<i32>(),
            result
                .atoms()
                .map(|(_, atom)| i32::from(atom.charge))
                .sum::<i32>()
        );
        for (idx, atom) in source.atoms() {
            assert_eq!(atom.atom_map, result.atom(idx).atom_map);
            assert_eq!(atom.chirality, result.atom(idx).chirality);
            assert_eq!(
                source.stereo_neighbor_order(idx),
                result.stereo_neighbor_order(idx)
            );
        }
        assert_eq!(
            canonical_smiles(&result),
            expected("C[C@H](F)C(=O)[O-:7].C")
        );
    }
}
