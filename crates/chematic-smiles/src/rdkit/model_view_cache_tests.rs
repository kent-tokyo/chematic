//! Frozen pre-cache callback contract, including error fallback and TLS scope.
use super::*;
use chematic_core::{AtomIdx, DerivedSlot, Element};
use std::sync::Arc;

fn frozen_model_callback<R>(
    mol: &Molecule,
    f: impl FnOnce(&Molecule) -> R,
    reuse_aromatic_view: bool,
) -> R {
    if IN_RDKIT_MODEL.with(|c| c.get()) {
        return f(mol);
    }
    let needed = chematic_perception::rdkit_model_may_disagree(mol)
        || mol
            .atoms()
            .any(|(_, a)| !a.wildcard && a.element == chematic_core::Element::H);
    let model = if needed {
        rdkit_model_molecule(mol, reuse_aromatic_view).ok()
    } else {
        None
    };
    match model {
        Some(m) => {
            struct Reset;
            impl Drop for Reset {
                fn drop(&mut self) {
                    IN_RDKIT_MODEL.with(|c| c.set(false));
                }
            }
            IN_RDKIT_MODEL.with(|c| c.set(true));
            let _reset = Reset;
            f(&m)
        }
        None => f(mol),
    }
}

fn snapshot(mol: &Molecule) -> String {
    let atoms: Vec<_> = mol
        .atoms()
        .map(|(i, a)| {
            (
                i,
                a,
                mol.neighbor_slice(i),
                mol.stereo_neighbor_order(i),
                mol.atom_tag(i),
                mol.r_group_label(i),
            )
        })
        .collect();
    let bonds: Vec<_> = mol
        .bonds()
        .map(|(i, b)| (i, b, mol.bond_direction(i)))
        .collect();
    let groups: Vec<_> = mol
        .stereo_groups()
        .iter()
        .enumerate()
        .map(|(i, g)| (g, mol.stereo_group_bond_indices(i)))
        .collect();
    format!(
        "{atoms:?}/{bonds:?}/{groups:?}/{}",
        IN_RDKIT_MODEL.with(|c| c.get())
    )
}

fn compare(mol: &Molecule) {
    for reuse in [false, true] {
        let reference = mol.clone();
        let candidate = mol.clone();
        let expected = frozen_model_callback(&reference, snapshot, reuse);
        for _ in 0..2 {
            assert_eq!(
                with_rdkit_model_molecule_impl(&candidate, snapshot, reuse),
                expected
            );
            assert!(!IN_RDKIT_MODEL.with(|c| c.get()));
        }
    }
}

#[test]
fn cached_model_matches_frozen_callback_on_pinned_corpora() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut count = 0;
    for path in [
        "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi",
        "validation/cosmolkit_comparison/corpora/stress-v1.smi",
    ] {
        let corpus = std::fs::read_to_string(root.join(path)).unwrap();
        for line in corpus.lines().filter(|s| !s.trim().is_empty()) {
            let smiles = line.split_whitespace().next().unwrap();
            compare(&crate::parse_template(smiles).unwrap());
            count += 1;
        }
    }
    assert_eq!(count, 14072);
}

#[test]
fn model_slot_shares_only_scalar_views_and_clears_after_edits() {
    let mut mol = crate::parse("[H][C@](F)(Cl)Br").unwrap();
    let slot = DerivedSlot::RdkitScalarModel;
    with_rdkit_model_molecule(&mol, |_| ());
    assert!(mol.derived_if_computed::<Option<Molecule>>(slot).is_none());
    with_rdkit_model_scalar(
        &mol,
        chematic_core::derived_cache::DerivedScalar::Tpsa,
        |m| m.atom_count() as f64,
    );
    let first = mol.derived_if_computed::<Option<Molecule>>(slot).unwrap();
    assert!(first.is_some());
    with_rdkit_model_scalar(
        &mol,
        chematic_core::derived_cache::DerivedScalar::LabuteAsa,
        |m| m.bond_count() as f64,
    );
    let second = mol.derived_if_computed::<Option<Molecule>>(slot).unwrap();
    assert!(Arc::ptr_eq(&first, &second));
    assert!(
        mol.clone()
            .derived_if_computed::<Option<Molecule>>(slot)
            .is_none()
    );
    mol.set_element(AtomIdx(2), Element::CL);
    assert!(mol.derived_if_computed::<Option<Molecule>>(slot).is_none());
    compare(&mol);
    mol.set_hydrogen_count(AtomIdx(1), Some(1));
    compare(&mol);
    let plain = crate::parse("CCO").unwrap();
    with_rdkit_model_molecule_impl(&plain, |_| (), true);
    assert!(
        plain
            .derived_if_computed::<Option<Molecule>>(slot)
            .is_none()
    );
}

#[test]
fn binding_preparation_only_materializes_required_scalar_models() {
    let slot = DerivedSlot::RdkitScalarModel;
    let mut unusual = crate::parse("[H][C@](F)(Cl)Br").unwrap();
    assert!(prepare_rdkit_scalar_model(&unusual));
    assert!(
        unusual
            .derived_if_computed::<Option<Molecule>>(slot)
            .is_some_and(|model| model.is_some())
    );
    unusual.set_element(AtomIdx(2), Element::CL);
    assert!(
        unusual
            .derived_if_computed::<Option<Molecule>>(slot)
            .is_none()
    );

    let ordinary = crate::parse("CCO").unwrap();
    assert!(!prepare_rdkit_scalar_model(&ordinary));
    assert!(
        ordinary
            .derived_if_computed::<Option<Molecule>>(slot)
            .is_none()
    );
}

#[test]
fn cached_failure_falls_back_and_panic_restores_nested_scope() {
    for s in [
        "[H]O[H]",
        "[2H][C@](F)(Cl)Br",
        "OCl(=O)(=O)=O",
        "[NH4+]",
        "[O]",
        "[CH2]C",
        "[Na+].[Cl-]",
        "c1cccc1",
        "CO(C)C",
        "C1C[N+]2=CC3=CC=CC=C3O[Cu]24OC5=CC=CC=C5C=[N+]14",
    ] {
        compare(&crate::parse_template(s).unwrap());
    }
    let bad = crate::parse_template("[H]O(C)C").unwrap();
    assert!(rdkit_model_molecule(&bad, true).is_err());
    assert_eq!(
        with_rdkit_model_molecule_impl(&bad, snapshot, true),
        snapshot(&bad)
    );
    assert!(
        bad.derived_if_computed::<Option<Molecule>>(DerivedSlot::RdkitScalarModel)
            .unwrap()
            .is_none()
    );
    let mol = crate::parse("[H]O[H]").unwrap();
    let panic = std::panic::catch_unwind(|| {
        with_rdkit_model_molecule_impl(
            &mol,
            |m| {
                assert!(IN_RDKIT_MODEL.with(|c| c.get()));
                assert_eq!(
                    with_rdkit_model_molecule_impl(m, snapshot, true),
                    snapshot(m)
                );
                panic!("callback panic");
            },
            true,
        )
    });
    assert!(panic.is_err());
    assert!(!IN_RDKIT_MODEL.with(|c| c.get()));
    compare(&mol);
}

#[test]
fn concurrent_cold_model_reads_keep_exact_results() {
    let mol = Arc::new(crate::parse("[H][C@](F)(Cl)Br").unwrap());
    let expected = frozen_model_callback(&mol, snapshot, true);
    let threads: Vec<_> = (0..8)
        .map(|_| {
            let mol = Arc::clone(&mol);
            std::thread::spawn(move || with_rdkit_model_molecule_impl(&mol, snapshot, true))
        })
        .collect();
    for t in threads {
        assert_eq!(t.join().unwrap(), expected);
    }
}
