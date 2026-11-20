//! CIP E/Z of an entire planar cumulene rather than of a constituent bond.

use chematic_core::cumulene::{PlanarCumulene, planar_cumulenes};
use smallvec::SmallVec;

use super::double_bond_cip::highest_stereo_sub_accurate;
use super::{AtomIdx, BondOrder, CipCode, LazyEzRanker, Molecule};

/// E/Z assignment on a chain of at least three consecutive double bonds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CumuleneEzAssignment {
    /// From one terminal atom to the other, inclusive.
    pub atoms: Vec<AtomIdx>,
    pub descriptor: CipCode,
}

/// Assign E/Z to configured odd cumulenes with distinguishable terminal ligands.
pub fn assign_cumulene_ez(mol: &Molecule) -> Vec<CumuleneEzAssignment> {
    let ranker = LazyEzRanker::new(mol, chematic_cip::CipBudget::default_budget());
    planar_cumulenes(mol)
        .into_iter()
        .filter_map(|axis| assign_axis(mol, axis, &ranker))
        .collect()
}

fn assign_axis(
    mol: &Molecule,
    axis: PlanarCumulene,
    ranker: &LazyEzRanker<'_>,
) -> Option<CumuleneEzAssignment> {
    let side = |end: chematic_core::ez_markers::EzFactEnd| {
        let reference = mol.bond(end.reference).other(end.atom)?;
        let substituents: SmallVec<[AtomIdx; 2]> = mol
            .neighbors(end.atom)
            .filter(|&(_, bond)| mol.bond(bond).order != BondOrder::Double)
            .map(|(atom, _)| atom)
            .collect();
        highest_stereo_sub_accurate(mol, end.atom, &substituents, ranker, false, |atom| {
            Some(end.up ^ (atom != reference))
        })
        .map(|(up, _)| up)
    };
    let descriptor = if side(axis.ends[0])? == side(axis.ends[1])? {
        CipCode::Z
    } else {
        CipCode::E
    };
    Some(CumuleneEzAssignment {
        atoms: axis.atoms,
        descriptor,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labeled_axes_survive_writing_and_ranking_at_both_ends() {
        for (text, code, atoms) in [
            (
                r"O[C@H]/1CCCCCCCCCCCCCCCCCCCCCCC[C@H]2CC1=C=C=C/CCCCCCCCCCCCCCCCC2",
                CipCode::E,
                4,
            ), // VS063
            (r"Br\C=C=C=C\Cl", CipCode::E, 4),   // VS118
            (r"C/C=C=C=C\C", CipCode::Z, 4),     // VS135
            (r"C/C=C=C=C/C=C/C", CipCode::E, 4), // VS154
            (
                r"C/1CCCCCCCCCCCCCCCCCCC[C@@H]2CC1=C=C=C/CCCCCCCCCCCCCCCCCCCCC2",
                CipCode::Z,
                4,
            ), // VS164
            (r"F/C=C=C=C=C=C\Br", CipCode::Z, 6),
            (r"F/C(Cl)=C=C=C/Br", CipCode::Z, 4),
            (r"[2H]/C=C=C=C/[1H]", CipCode::Z, 4),
        ] {
            let mol = chematic_smiles::parse(text).unwrap();
            let spellings = [
                text.to_owned(),
                chematic_smiles::write(&mol),
                chematic_smiles::canonical_smiles(&mol),
            ]
            .into_iter()
            .chain((1..=16).map(|seed| chematic_smiles::random_smiles(&mol, seed)));
            for spelling in spellings {
                let reparsed = chematic_smiles::parse(&spelling).unwrap();
                let labels = assign_cumulene_ez(&reparsed);
                assert_eq!(labels.len(), 1, "{text} -> {spelling}");
                assert_eq!(labels[0].descriptor, code, "{text} -> {spelling}");
                assert_eq!(labels[0].atoms.len(), atoms, "{text} -> {spelling}");
            }
        }
    }

    #[test]
    fn unconfigured_and_nondistinguishable_ligands_stay_unassigned() {
        for text in [
            r"F/C=C=C/F",
            r"F/C=C=C=CF",
            r"F/C(F)=C=C=C/Br",
            r"F/C=C=C=C(Br)I",
        ] {
            let mol = chematic_smiles::parse(text).unwrap();
            assert!(assign_cumulene_ez(&mol).is_empty(), "{text}");
        }
    }

    #[test]
    fn canonical_keeps_isomers_distinct_and_converges_across_carriers() {
        let canonical =
            |text| chematic_smiles::canonical_smiles(&chematic_smiles::parse(text).unwrap());
        assert_ne!(canonical(r"F/C=C=C=C/F"), canonical(r"F/C=C=C=C\F"));
        assert_eq!(
            canonical(r"F/C(Cl)=C=C=C(/Br)I"),
            canonical(r"Cl\C(F)=C=C=C(\I)Br")
        );
    }
}
