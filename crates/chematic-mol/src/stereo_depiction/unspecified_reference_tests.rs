use super::super::{AtomIdx, BondIdx, BondOrder, HashMap, Molecule, in_ring_smaller_than};
use super::unspecified_double_bonds;

/// Double bonds outside small rings whose ends each have one or two
/// distinct substituents (an implicit H or a lone pair making up the
/// rest), with no E/Z in `declared`.
fn reference_unspecified_double_bonds(
    mol: &Molecule,
    declared: &HashMap<BondIdx, chematic_core::atom::CipCode>,
) -> Vec<BondIdx> {
    // The substituents of `end` when it can carry E/Z, before the check
    // that two substituents differ (classes are computed only if needed).
    let end_subs = |end: AtomIdx, other: AtomIdx, db: BondIdx| {
        let subs: Vec<AtomIdx> = mol
            .neighbors(end)
            .filter(|&(nb, bi)| nb != other && mol.bond(bi).order != BondOrder::Dative)
            .map(|(nb, _)| nb)
            .collect();
        let cumulated = mol
            .neighbors(end)
            .any(|(_, bi)| bi != db && mol.bond(bi).order == BondOrder::Double);
        let h = chematic_core::implicit_hcount(mol, end) as usize;
        (!cumulated && !subs.is_empty() && subs.len() + h <= 2).then_some(subs)
    };
    let classes = std::cell::OnceCell::new();
    let distinct = |subs: &[AtomIdx]| {
        subs.len() < 2 || {
            let classes =
                classes.get_or_init(|| chematic_smiles::topological_equivalence_classes(mol));
            classes[subs[0].0 as usize] != classes[subs[1].0 as usize]
        }
    };
    let mut out: Vec<BondIdx> = mol
        .bonds()
        .filter(|&(bi, b)| {
            if b.order != BondOrder::Double
                || declared.contains_key(&bi)
                || (mol.atom(b.atom1).aromatic && mol.atom(b.atom2).aromatic)
            {
                return false;
            }
            let (Some(s1), Some(s2)) = (
                end_subs(b.atom1, b.atom2, bi),
                end_subs(b.atom2, b.atom1, bi),
            ) else {
                return false;
            };
            !in_ring_smaller_than(mol, bi, 8) && distinct(&s1) && distinct(&s2)
        })
        .map(|(bi, _)| bi)
        .collect();
    out.sort_by_key(|b| b.0);
    out
}

fn compare(mol: &Molecule, label: &str) {
    use chematic_core::atom::CipCode;
    for selection in 0..3 {
        let declared = mol
            .bonds()
            .filter(|(bi, _)| match selection {
                0 => false,
                1 => true,
                _ => bi.0 % 2 == 0,
            })
            .map(|(bi, _)| {
                (
                    bi,
                    if bi.0 % 2 == 0 {
                        CipCode::E
                    } else {
                        CipCode::Z
                    },
                )
            })
            .collect();
        assert_eq!(
            unspecified_double_bonds(mol, &declared),
            reference_unspecified_double_bonds(mol, &declared),
            "{label}, declared selection={selection}"
        );
    }
}

#[test]
fn whole_classifier_matches_frozen_code_on_all_corpora() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut count = 0;
    for path in [
        "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi",
        "validation/cosmolkit_comparison/corpora/stress-v1.smi",
    ] {
        for (index, line) in std::fs::read_to_string(root.join(path))
            .unwrap()
            .lines()
            .filter(|line| !line.trim().is_empty())
            .enumerate()
        {
            let mol = chematic_smiles::parse(line.split_whitespace().next().unwrap()).unwrap();
            compare(&mol, &format!("{path}:{index}"));
            count += 1;
        }
    }
    assert_eq!(count, 14072);
    eprintln!("frozen classifier: {count} corpus rows, three declaration maps each");
}

#[test]
fn bounded_storage_preserves_high_degree_dative_and_cumulated_eligibility() {
    use chematic_core::{Atom, Element, MoleculeBuilder};
    compare(&MoleculeBuilder::new().build(), "empty");
    let mut count = 0;
    for left in 0..7 {
        for right in 0..7 {
            for dative in 0..4 {
                for cumulated in 0..4 {
                    let mut builder = MoleculeBuilder::new();
                    let a = builder.add_atom(Atom::new(Element::C));
                    let b = builder.add_atom(Atom::new(Element::C));
                    builder.add_bond(a, b, BondOrder::Double).unwrap();
                    for (side, end, degree) in [(0, a, left), (1, b, right)] {
                        for index in 0..degree {
                            let sub = builder.add_atom(Atom::new(if index % 2 == 0 {
                                Element::F
                            } else {
                                Element::CL
                            }));
                            let order = if index == 0 && dative & (1 << side) != 0 {
                                BondOrder::Dative
                            } else if index + 1 == degree && cumulated & (1 << side) != 0 {
                                BondOrder::Double
                            } else {
                                BondOrder::Single
                            };
                            builder.add_bond(end, sub, order).unwrap();
                        }
                    }
                    compare(
                        &builder.build(),
                        &format!("degrees={left},{right};dative={dative};cumulated={cumulated}"),
                    );
                    count += 1;
                }
            }
        }
    }
    assert_eq!(count, 784);
    for smiles in [
        "CC=CC",
        "CC=C(C)C",
        "C=C=C",
        "N=N",
        "C1CC=CCC1",
        "C1CCCCCC=CC1",
        "F/C=C/F",
        "CC(=O)O",
        "C=C(Cl)Br",
        "[H]C(F)=C(Cl)Br",
    ] {
        compare(&chematic_smiles::parse(smiles).unwrap(), smiles);
    }
}
