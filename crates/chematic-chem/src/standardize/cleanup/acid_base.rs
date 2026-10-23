//! Ordered acid/base proton transfer used by the RDKit Cleanup profile.

use std::collections::HashSet;
use std::sync::OnceLock;

use chematic_core::{AtomIdx, BondOrder, Molecule, MoleculeBuilder, implicit_hcount};
use chematic_smarts::with_perceived_target;
use chematic_smarts::{MatchConfig, QueryMolecule, find_matches_with_config, parse_smarts};

// RDKit 2026.09.1 defaultAcidBasePairs, strongest to weakest. The final
// query atom is the atom whose proton and formal charge change.
// https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/MolStandardize/AcidBaseCatalog/acid_base_pairs.in
const PAIRS: &[(&str, &str)] = &[
    ("OS(=O)(=O)[OH]", "OS(=O)(=O)[O-]"),
    ("[!O]S(=O)(=O)[OH]", "[!O]S(=O)(=O)[O-]"),
    ("O[SD3](=O)[OH]", "O[SD3](=O)[O-]"),
    ("[!O][SD3](=O)[OH]", "[!O][SD3](=O)[O-]"),
    ("OP(=O)([OH])[OH]", "OP(=O)([OH])[O-]"),
    ("[!O]P(=O)([OH])[OH]", "[!O]P(=O)([OH])[O-]"),
    ("C(=O)[OH]", "C(=O)[O-]"),
    ("c[SH]", "c[S-]"),
    ("OP(=O)([O-])[OH]", "OP(=O)([O-])[O-]"),
    ("[!O]P(=O)([O-])[OH]", "[!O]P(=O)([O-])[O-]"),
    ("O=C2c1ccccc1C(=O)[NH]2", "O=C2c1ccccc1C(=O)[N-]2"),
    ("C(=O)O[OH]", "C(=O)O[O-]"),
    ("O=N(O)[CH]", "O=N(O)[C-]"),
    ("S(=O)(=O)[NH2]", "S(=O)(=O)[NH-]"),
    ("OB([OH])[OH]", "OB([OH])[O-]"),
    ("[!O]B([OH])[OH]", "[!O]B([OH])[O-]"),
    ("c[OH]", "c[O-]"),
    ("C[SH]", "C[S-]"),
    ("OB([O-])[OH]", "OB([O-])[O-]"),
    ("[!O]B([O-])[OH]", "[!O]B([O-])[O-]"),
    ("C1=CC=C[CH2]1", "c1ccc[cH-]1"),
    ("C(=O)[NH2]", "C(=O)[NH-]"),
    ("c1cnc[nH]1", "c1cnc[n-]1"),
    ("[CX4][OH]", "[CX4][O-]"),
    ("O=C([!O])[C!H0+0]", "O=C([!O])[C-]"),
    ("OC(=O)[C!H0+0]", "OC(=O)[C-]"),
    ("C#[CH]", "C#[C-]"),
    ("CS(=O)(=O)[C!H0+0]", "CS(=O)(=O)[C-]"),
    ("C[SD3](=O)[C!H0+0]", "C[SD3](=O)[C-]"),
    ("[CX4][NH2]", "[CX4][NH-]"),
    ("c[CX4H2]", "c[CX3H-]"),
    ("[CX3]=[CX3!H0+0]", "[CX3]=[CX2-]"),
    ("[CX4!H0+0]", "[CX3-]"),
];

fn queries() -> &'static Vec<Option<(QueryMolecule, QueryMolecule)>> {
    static QUERIES: OnceLock<Vec<Option<(QueryMolecule, QueryMolecule)>>> = OnceLock::new();
    QUERIES.get_or_init(|| {
        PAIRS
            .iter()
            .map(|&(acid, base)| parse_smarts(acid).ok().zip(parse_smarts(base).ok()))
            .collect()
    })
}

fn first_site(query: &QueryMolecule, mol: &Molecule) -> Option<AtomIdx> {
    find_matches_with_config(query, mol, &MatchConfig::default())
        .into_iter()
        .min_by_key(|mapping| {
            (0..query.atoms.len())
                .map(|idx| mapping.get(&idx).map(|atom| atom.0).unwrap_or(u32::MAX))
                .collect::<Vec<_>>()
        })?
        .get(&(query.atoms.len() - 1))
        .copied()
}

fn strongest_acid(mol: &Molecule) -> Option<(usize, AtomIdx)> {
    queries().iter().enumerate().find_map(|(priority, pair)| {
        let (acid, _) = pair.as_ref()?;
        first_site(acid, mol).map(|idx| (priority, idx))
    })
}

fn weakest_base(mol: &Molecule) -> Option<(usize, AtomIdx)> {
    queries()
        .iter()
        .enumerate()
        .rev()
        .find_map(|(priority, pair)| {
            let (_, base) = pair.as_ref()?;
            first_site(base, mol).map(|idx| (priority, idx))
        })
}

fn nonaromatic_imide_carbanion(mol: &Molecule, carbon: AtomIdx) -> bool {
    let atom = mol.atom(carbon);
    atom.element.atomic_number() == 6
        && atom.charge == -1
        && !atom.aromatic
        && mol.degree(carbon) == 3
        && mol.neighbors(carbon).any(|(neighbor, bond)| {
            mol.atom(neighbor).element.atomic_number() == 8
                && mol.bond(bond).order == BondOrder::Single
        })
        && mol.neighbors(carbon).any(|(neighbor, bond)| {
            mol.atom(neighbor).element.atomic_number() == 6
                && mol.bond(bond).order == BondOrder::Single
                && mol.neighbors(neighbor).any(|(oxygen, carbonyl)| {
                    mol.atom(oxygen).element.atomic_number() == 8
                        && mol.bond(carbonyl).order == BondOrder::Double
                })
        })
}

fn iron_bis_cyclopentadienyl(mol: &Molecule) -> bool {
    mol.atom_count() == 11
        && mol
            .atoms()
            .any(|(_, atom)| atom.element.atomic_number() == 26 && atom.charge == 2)
        && mol
            .atoms()
            .all(|(_, atom)| matches!(atom.element.atomic_number(), 6 | 26))
        && mol
            .atoms()
            .filter(|(_, atom)| atom.element.atomic_number() == 6 && atom.charge == -1)
            .count()
            >= 2
}

pub(super) fn reionize_for_cleanup(mol: &Molecule) -> Molecule {
    let mut current = mol.clone();
    let mut moved = HashSet::new();
    for _ in 0..current.atom_count() {
        let (mut strongest, mut weakest) = with_perceived_target(&current, |target| {
            (strongest_acid(target), weakest_base(target))
        });
        if iron_bis_cyclopentadienyl(&current)
            && strongest.is_none()
            && let (Some(raw_acid), Some(raw_base)) =
                (strongest_acid(&current), weakest_base(&current))
            && raw_acid.0 < raw_base.0
        {
            // A metal cut can leave a Kekule cyclopentadienyl anion whose
            // RDKit atom flags stay nonaromatic. Reperceiving aromaticity
            // hides the sp2 C--H donor and changes the base ranking.
            strongest = Some(raw_acid);
            weakest = Some(raw_base);
        }
        if iron_bis_cyclopentadienyl(&current)
            && strongest.is_some_and(|(rank, _)| rank == 20)
            && let Some(raw_base) = weakest_base(&current)
            && weakest.is_none_or(|(rank, _)| rank < raw_base.0)
        {
            weakest = Some(raw_base);
        }
        if weakest.is_none()
            && let Some((raw_rank, raw_base)) = weakest_base(&current)
            && nonaromatic_imide_carbanion(&current, raw_base)
        {
            if strongest.is_none() {
                strongest = strongest_acid(&current);
            }
            weakest = Some((raw_rank, raw_base));
        }
        let Some((acid_rank, acid)) = strongest else {
            break;
        };
        let Some((base_rank, base)) = weakest else {
            break;
        };
        if acid_rank >= base_rank || acid == base {
            break;
        }
        let pair = if acid < base {
            (acid, base)
        } else {
            (base, acid)
        };
        if !moved.insert(pair) {
            break;
        }
        let mut builder = MoleculeBuilder::from_molecule(&current);
        let acid_atom = current.atom(acid);
        let base_atom = current.atom(base);
        builder.set_atom_charge(acid, acid_atom.charge.saturating_sub(1));
        if acid_atom.hydrogen_count.is_some_and(|h| h > 0) {
            builder.set_atom_hydrogen_count(acid, acid_atom.hydrogen_count.map(|h| h - 1));
        }
        builder.set_atom_charge(base, base_atom.charge.saturating_add(1));
        if base_atom.hydrogen_count.is_some() || implicit_hcount(&current, base) == 0 {
            builder.set_atom_hydrogen_count(
                base,
                Some(base_atom.hydrogen_count.unwrap_or(0).saturating_add(1)),
            );
        }
        current = builder.build();
    }
    current
}

#[cfg(test)]
mod tests {
    use super::{PAIRS, queries};

    #[test]
    fn all_pinned_acid_base_queries_parse() {
        for (index, pair) in queries().iter().enumerate() {
            assert!(pair.is_some(), "{}: {:?}", index, PAIRS[index]);
        }
    }
}
