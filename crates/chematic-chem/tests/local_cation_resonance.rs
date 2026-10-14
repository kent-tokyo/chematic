use std::collections::BTreeSet;

use chematic_chem::{
    LocalResonanceError as Error, LocalResonanceLimits as Limits, LocalResonanceReport as Report,
    LocalResonanceStatus as Status, enumerate_local_cation_resonance as enumerate,
};
use chematic_core::{
    AtomIdx, BondIdx, BondOrder, Molecule, MoleculeBuilder, StereoGroup, StereoGroupKind,
    implicit_hcount,
};
use chematic_smiles::{canonical_smiles, parse};

fn run(smiles: &str) -> Report {
    enumerate(&parse(smiles).unwrap(), &Limits::default()).unwrap()
}

fn strings(report: &Report) -> BTreeSet<String> {
    report.forms.iter().map(canonical_smiles).collect()
}

fn assert_conserved(input: &Molecule, form: &Molecule) {
    assert_eq!(input.atom_count(), form.atom_count());
    assert_eq!(input.bond_count(), form.bond_count());
    assert_eq!(input.total_formula(), form.total_formula());
    let charge = |m: &Molecule| m.atoms().map(|(_, a)| i64::from(a.charge)).sum::<i64>();
    assert_eq!(charge(input), charge(form));
    for (idx, atom) in input.atoms() {
        let next = form.atom(idx);
        assert_eq!(atom.element, next.element);
        assert_eq!(atom.isotope, next.isotope);
        assert_eq!(atom.atom_map, next.atom_map);
        assert_eq!(implicit_hcount(input, idx), implicit_hcount(form, idx));
    }
    for (idx, bond) in input.bonds() {
        let next = form.bond(idx);
        assert_eq!((bond.atom1, bond.atom2), (next.atom1, next.atom2));
    }
}

#[test]
fn rdkit_9458_oxygen_and_sulfur_pairs_are_starting_form_independent() {
    let fixtures: Vec<_> =
        include_str!("../../../validation/rdkit_issues/resonance/local_cation_pairs.smi")
            .lines()
            .filter(|s| !s.starts_with('#') && !s.trim().is_empty())
            .map(|s| s.split_whitespace().next().unwrap())
            .collect();
    assert_eq!(fixtures.len(), 4);
    for pair in fixtures.as_chunks::<2>().0 {
        let expected: BTreeSet<_> = pair
            .iter()
            .map(|s| canonical_smiles(&parse(s).unwrap()))
            .collect();
        for smiles in pair {
            let input = parse(smiles).unwrap();
            let report = enumerate(&input, &Limits::default()).unwrap();
            assert_eq!(report.status, Status::CompleteWithinScope);
            assert_eq!(report.forms.len(), 2, "{smiles}");
            assert_eq!(strings(&report), expected, "{smiles}");
            for form in &report.forms {
                assert_conserved(&input, form);
            }
        }
    }
}

#[test]
fn isotope_atom_map_and_explicit_hydrogen_nodes_survive() {
    for (smiles, count) in [
        ("[13CH3][S:9]C([2H])=[N+:7](C)C", 2),
        // Either neutral N or S can donate into the initial iminium bond.
        ("[H]N([H])C(SC)=[NH2+]", 3),
    ] {
        let input = parse(smiles).unwrap();
        let report = enumerate(&input, &Limits::default()).unwrap();
        assert_eq!(report.forms.len(), count, "{smiles}");
        for form in &report.forms {
            assert_conserved(&input, form);
        }
    }
}

#[test]
fn independent_sites_form_the_cartesian_product_and_reenumerate_stably() {
    let input = parse("CSC(C)=[N+](C)C.COC(C)=[N+](C)C").unwrap();
    let report = enumerate(&input, &Limits::default()).unwrap();
    assert_eq!(report.status, Status::CompleteWithinScope);
    assert_eq!(report.forms.len(), 4);
    let expected = strings(&report);
    assert_eq!(expected.len(), 4);
    for form in &report.forms {
        assert_conserved(&input, form);
        assert_eq!(
            strings(&enumerate(form, &Limits::default()).unwrap()),
            expected
        );
    }
}

#[test]
fn symmetric_indexed_charge_placements_are_retained_without_cycles() {
    let report = run("NC(N)=[NH2+]");
    assert_eq!(report.status, Status::CompleteWithinScope);
    assert_eq!(report.forms.len(), 3);
    assert_eq!(strings(&report).len(), 1);
    assert_eq!(report.transforms_attempted, 6);
}

#[test]
fn exact_form_cap_completes_and_smaller_caps_report_truncation() {
    let input = parse("CSC(C)=[N+](C)C").unwrap();
    let exact = enumerate(
        &input,
        &Limits {
            max_forms: 2,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(exact.status, Status::CompleteWithinScope);
    assert_eq!(exact.transforms_attempted, 2);
    let small = enumerate(
        &input,
        &Limits {
            max_forms: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(small.status, Status::FormLimitReached);
    assert_eq!(small.forms.len(), 1);
    let bounded = enumerate(
        &input,
        &Limits {
            max_transforms: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(bounded.status, Status::TransformLimitReached);
    assert_eq!(bounded.transforms_attempted, 1);
    assert_eq!(bounded.forms.len(), 2);
}

#[test]
fn invalid_limits_and_input_caps_fail_before_enumeration() {
    let input = parse("CCO").unwrap();
    for limits in [
        Limits {
            max_forms: 0,
            ..Default::default()
        },
        Limits {
            max_forms: 1025,
            ..Default::default()
        },
        Limits {
            max_transforms: 0,
            ..Default::default()
        },
        Limits {
            max_transforms: 100_001,
            ..Default::default()
        },
        Limits {
            max_atoms: 0,
            ..Default::default()
        },
        Limits {
            max_atoms: 1025,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            enumerate(&input, &limits),
            Err(Error::InvalidLimit { .. })
        ));
    }
    assert!(matches!(
        enumerate(
            &input,
            &Limits {
                max_atoms: 2,
                ..Default::default()
            }
        ),
        Err(Error::AtomLimitExceeded {
            actual: 3,
            limit: 2
        })
    ));
    assert!(matches!(
        enumerate(&MoleculeBuilder::new().build(), &Limits::default()),
        Err(Error::EmptyMolecule)
    ));
}

#[test]
fn rejection_reasons_have_stable_actionable_messages() {
    let cases = [
        (
            Error::InvalidLimit {
                field: "max_forms",
                value: 0,
                maximum: 1024,
            },
            "max_forms must be in 1..=1024, got 0",
        ),
        (
            Error::EmptyMolecule,
            "local resonance requires a non-empty molecule",
        ),
        (
            Error::AtomLimitExceeded {
                actual: 3,
                limit: 2,
            },
            "local resonance atom count 3 exceeds limit 2",
        ),
        (
            Error::UnsupportedStereo,
            "local cation resonance does not support stereochemical annotations",
        ),
        (
            Error::QueryAtom { atom: AtomIdx(7) },
            "local resonance does not support query atom 7",
        ),
        (
            Error::QueryBond { bond: BondIdx(9) },
            "local resonance does not support query bond 9",
        ),
    ];
    for (error, expected) in cases {
        assert_eq!(error.to_string(), expected);
    }
}

#[test]
fn triple_bonds_are_counted_without_becoming_resonance_moves() {
    let report = run("N#CC=[NH2+]");
    assert_eq!(report.status, Status::CompleteWithinScope);
    assert_eq!(report.forms.len(), 1);
    assert_eq!(report.transforms_attempted, 0);
}

#[test]
fn stereochemical_and_query_inputs_are_rejected_explicitly() {
    for smiles in ["C/C(SC)=[N+](C)/C", "C[C@H](F)CSC(C)=[N+](C)C"] {
        assert!(matches!(
            enumerate(&parse(smiles).unwrap(), &Limits::default()),
            Err(Error::UnsupportedStereo)
        ));
    }
    let mut grouped = parse("CSC(C)=[N+](C)C").unwrap();
    grouped.add_stereo_group(StereoGroup::new(StereoGroupKind::And(1), vec![AtomIdx(0)]));
    assert!(matches!(
        enumerate(&grouped, &Limits::default()),
        Err(Error::UnsupportedStereo)
    ));
    assert!(matches!(
        enumerate(&parse("*C").unwrap(), &Limits::default()),
        Err(Error::QueryAtom { .. })
    ));
    let mut query = parse("CC").unwrap();
    query.set_bond_order(BondIdx(0), BondOrder::QueryAny);
    assert!(matches!(
        enumerate(&query, &Limits::default()),
        Err(Error::QueryBond { .. })
    ));
}

#[test]
fn unrelated_chemistry_does_not_gain_new_charge_or_hydrogens() {
    for smiles in [
        "CCO",
        "COC(C)=NC",
        "CS(=O)C",
        "[O-]C=[NH2+]",
        "C[N+](=O)[O-]",
        "Cn1cc[n+](C)c1",
    ] {
        let input = parse(smiles).unwrap();
        let report = enumerate(&input, &Limits::default()).unwrap();
        assert_eq!(report.status, Status::CompleteWithinScope);
        assert_eq!(report.forms.len(), 1, "{smiles}");
        assert_eq!(report.transforms_attempted, 0, "{smiles}");
        assert_conserved(&input, &report.forms[0]);
    }
}

#[test]
fn reversing_atom_insertion_order_preserves_the_resonance_set() {
    let input = parse("CSC(C)=[N+](C)C").unwrap();
    let mut builder = MoleculeBuilder::new();
    let mut mapping = vec![AtomIdx(0); input.atom_count()];
    for old in (0..input.atom_count()).rev() {
        mapping[old] = builder.add_atom(input.atom(AtomIdx(old as u32)).clone());
    }
    for (_, bond) in input.bonds() {
        builder
            .add_bond(
                mapping[bond.atom1.0 as usize],
                mapping[bond.atom2.0 as usize],
                bond.order,
            )
            .unwrap();
    }
    let reversed = builder.build();
    assert_eq!(
        strings(&enumerate(&input, &Limits::default()).unwrap()),
        strings(&enumerate(&reversed, &Limits::default()).unwrap())
    );
}
