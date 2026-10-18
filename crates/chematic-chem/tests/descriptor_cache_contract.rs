use chematic_chem::{labute_asa, qed, rdkit_tpsa, tpsa};
use chematic_core::derived_cache::DerivedScalar;
use chematic_core::{Atom, AtomIdx, BondIdx, BondOrder, Element, Molecule};

fn values(mol: &Molecule) -> [u64; 4] {
    [
        tpsa(mol).to_bits(),
        rdkit_tpsa(mol).to_bits(),
        labute_asa(mol).to_bits(),
        qed(mol).to_bits(),
    ]
}

fn check_after_edit(mol: &Molecule) {
    let fresh = mol.clone();
    assert_eq!(values(mol), values(&fresh));
    assert_eq!(values(mol), values(&fresh));
}

#[test]
fn descriptor_variants_remain_distinct_and_clones_start_uncached() {
    let mol = chematic_smiles::parse("OP(=O)(O)O").unwrap();
    let first = values(&mol);
    assert_ne!(first[0], first[1]);
    assert_eq!(values(&mol), first);
    let clone = mol.clone();
    for (index, slot) in [
        DerivedScalar::Tpsa,
        DerivedScalar::RdkitTpsa,
        DerivedScalar::LabuteAsa,
        DerivedScalar::Qed,
    ]
    .into_iter()
    .enumerate()
    {
        let computed = std::cell::Cell::new(false);
        clone.derived_scalar(slot, || {
            computed.set(true);
            f64::from_bits(first[index])
        });
        assert!(computed.get(), "clone inherited a scalar cache");
    }
    assert_eq!(values(&clone), first);
}

#[test]
fn atom_bond_and_topology_edits_invalidate_scalar_descriptors() {
    let mut mol = chematic_smiles::parse("CCO").unwrap();
    let ethanol = values(&mol);
    mol.set_element(AtomIdx(2), Element::N);
    check_after_edit(&mol);
    assert_ne!(values(&mol), ethanol);
    mol.set_charge(AtomIdx(2), 1);
    check_after_edit(&mol);
    mol.set_hydrogen_count(AtomIdx(2), Some(3));
    check_after_edit(&mol);
    mol.set_isotope(AtomIdx(2), Some(15));
    check_after_edit(&mol);
    mol.set_bond_order(BondIdx(0), BondOrder::Double);
    check_after_edit(&mol);
    let added = mol.add_atom(Atom::new(Element::O));
    check_after_edit(&mol);
    mol.add_bond(AtomIdx(0), added, BondOrder::Single).unwrap();
    check_after_edit(&mol);
    mol.remove_bond(BondIdx(2));
    check_after_edit(&mol);
    mol.remove_atom(added);
    check_after_edit(&mol);
}

#[test]
fn explicit_hydrogen_and_cleanup_views_keep_repeat_results() {
    for smiles in [
        "[H]O[H]",
        "OCl(=O)(=O)=O",
        "[NH4+]",
        "c1ncccc1",
        "C1CO1",
        "N1CC1",
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        check_after_edit(&mol);
    }
}

#[test]
fn concurrent_readers_observe_the_same_descriptor_values() {
    let mol = std::sync::Arc::new(chematic_smiles::parse("CC(=O)Oc1ccccc1C(=O)O").unwrap());
    let fresh = mol.as_ref().clone();
    let expected = values(&fresh);
    let threads: Vec<_> = (0..8)
        .map(|_| {
            let mol = mol.clone();
            std::thread::spawn(move || (0..16).map(|_| values(&mol)).collect::<Vec<_>>())
        })
        .collect();
    for thread in threads {
        assert!(
            thread
                .join()
                .unwrap()
                .into_iter()
                .all(|value| value == expected)
        );
    }
}

#[test]
fn nested_model_calculations_bypass_an_outer_model_scalar() {
    let mol = chematic_smiles::parse("[H]O[H]").unwrap();
    let count = |view: &Molecule| view.atom_count() as f64;
    let slot = DerivedScalar::RdkitTpsa;
    assert_eq!(
        chematic_smiles::with_rdkit_model_scalar(&mol, slot, count),
        1.0
    );
    let nested = chematic_smiles::with_rdkit_model_molecule(&mol, |_| {
        chematic_smiles::with_rdkit_model_scalar(&mol, slot, count)
    });
    assert_eq!(nested, 3.0);
    assert_eq!(
        chematic_smiles::with_rdkit_model_scalar(&mol, slot, count),
        1.0
    );
}
