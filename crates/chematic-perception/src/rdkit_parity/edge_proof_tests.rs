//! Exact old shortcut oracle; no chemical rules or fallback are changed.
use super::*;

fn frozen_shortcut(mol: &Molecule) -> ParityShortcut {
    *mol.derived(chematic_core::DerivedSlot::RdkitParityShortcut, || {
        let model_may_disagree = crate::rdkit_model_hook::rdkit_model_may_disagree(mol);
        // Molecules whose view the RDKit-model hook may correct take the
        // full path, where the hook runs.
        if model_may_disagree && crate::rdkit_model_hook::hook_for(mol).is_some() {
            return ParityShortcut::Full;
        }
        // RDKit's sanitization rewrites some graphs (perchlorate,
        // organometallic dative bonds); the view is then never `mol`. Every
        // cleanup candidate is already included in the broader model gate,
        // so ordinary molecules do not need a second full atom scan.
        if model_may_disagree && crate::rdkit_cleanup::rdkit_sanitize_cleanup_needed(mol) {
            return ParityShortcut::Full;
        }
        // The existing complete-preperceived test returns Identity for a
        // forest with no aromatic flags, and for fully saturated graphs.
        // Reuse those proofs before building cyclic-component scratch data.
        // Unusual models retain the full checks and cleanup behavior above.
        if !model_may_disagree
            && !mol.atoms().any(|(_, atom)| atom.aromatic)
            && (mol.bonds().all(|(_, bond)| {
                matches!(
                    bond.order,
                    BondOrder::Single | BondOrder::Up | BondOrder::Down
                )
            }) || mol
                .derived_if_computed::<Vec<bool>>(chematic_core::DerivedSlot::RingBondFlags)
                .is_some_and(|flags| !flags.iter().any(|&cyclic| cyclic)))
        {
            return ParityShortcut::Identity;
        }
        let info = ComponentInfo::new(mol);
        if has_aromatic_macrocycle_bridge_donor(mol) {
            return ParityShortcut::Full;
        }
        match preperceived_kind_with(mol, &info) {
            Some(kind) => kind,
            None if unchanged_without_kekulization_with(mol, &info) => ParityShortcut::Identity,
            None if explicit_kept_by_kekule_verdict(mol, &info) => ParityShortcut::Identity,
            None => ParityShortcut::Full,
        }
    })
}

fn check(mol: &Molecule, label: &str) -> bool {
    let expected = frozen_shortcut(&mol.clone());
    assert_eq!(parity_shortcut(mol), expected, "cold: {label}");
    assert_eq!(parity_shortcut(mol), expected, "warm: {label}");
    assert_eq!(parity_shortcut(&mol.clone()), expected, "clone: {label}");
    let eligible = !crate::rdkit_model_hook::rdkit_model_may_disagree(mol)
        && !has_aromatic_macrocycle_bridge_donor(mol);
    let proven = eligible && preperceived_edges_are_complete(mol);
    if proven {
        assert_eq!(expected, ParityShortcut::Identity, "proof: {label}");
        assert_eq!(
            preperceived_kind_with(mol, &ComponentInfo::new(mol)),
            Some(ParityShortcut::Identity),
            "component proof: {label}"
        );
    }
    proven
}

#[test]
fn edge_proof_matches_frozen_shortcut_on_all_pinned_inputs() {
    let corpora = [
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
        )),
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi"
        )),
    ];
    let mut rows = 0;
    let mut proven = 0;
    let mut slice_proven = [0; 2];
    for (corpus_index, corpus) in corpora.iter().enumerate() {
        for (row, line) in corpus
            .lines()
            .filter(|line| !line.trim().is_empty())
            .enumerate()
        {
            let smiles = line.split_whitespace().next().unwrap();
            let mol = chematic_smiles::parse_template(smiles).unwrap();
            let yes = check(&mol, smiles);
            proven += usize::from(yes);
            if corpus_index == 0 {
                if row < 1000 {
                    slice_proven[0] += usize::from(yes);
                }
                if (5000..6000).contains(&row) {
                    slice_proven[1] += usize::from(yes);
                }
            }
            rows += 1;
        }
    }
    assert_eq!(rows, 14072);
    assert!(proven > 100);
    eprintln!(
        "Edge identity proof: {rows} inputs, {proven} positive proofs, slices {slice_proven:?}; cold/warm/clone decisions match"
    );
}

#[test]
fn edge_proof_matches_all_five_atom_graphs_with_single_or_multiple_bonds() {
    use chematic_core::{Atom, Element, MoleculeBuilder};
    chematic_smiles::register_rdkit_model_hook();
    let edges: Vec<_> = (0..5)
        .flat_map(|a| (a + 1..5).map(move |b| (a, b)))
        .collect();
    let mut rows = 0;
    for mask in 0..1usize << edges.len() {
        for multiple in [false, true] {
            let mut builder = MoleculeBuilder::new();
            for _ in 0..5 {
                builder.add_atom(Atom::new(Element::C));
            }
            for (index, &(a, b)) in edges.iter().enumerate() {
                if mask & (1 << index) != 0 {
                    let order = if multiple && index % 3 == 0 {
                        BondOrder::Double
                    } else {
                        BondOrder::Single
                    };
                    builder.add_bond(AtomIdx(a), AtomIdx(b), order).unwrap();
                }
            }
            check(&builder.build(), "five-atom graph");
            rows += 1;
        }
    }
    assert_eq!(rows, 2048);
    eprintln!("Edge identity graph oracle: {rows} graphs");
}

#[test]
fn mixed_exocyclic_macrocycle_and_mutated_graphs_keep_the_fallback() {
    for smiles in [
        "c1ccccc1",
        "c1ccc2ccccc2c1",
        "c1ccoc1",
        "c1ccsc1",
        "o1cccccccc1",
        "s1cccccccc1",
        "c1ccccc1.C1CCCCC1",
        "c1ccc2c(c1)CCC2",
        "O=C1NNC(=O)N1c1ccccc1",
        "O=C1NC(=O)/C(=C/c2ccccc2)N1",
        "C1=CC=CC=C1",
        "c1ccccc1.[*]",
        "[2H]Oc1ccccc1",
        "[Na+].[O-]c1ccccc1",
    ] {
        let mut mol = chematic_smiles::parse_template(smiles).unwrap();
        check(&mol, smiles);
        if mol.bond_count() != 0 {
            mol.set_bond_order(BondIdx(0), BondOrder::Single);
            check(&mol, "after bond mutation");
        }
        let flags = crate::sssr::ring_bond_flags_shared(&mol);
        let seeded = mol.clone();
        seeded.seed_derived(chematic_core::DerivedSlot::RingBondFlags, flags);
        check(&seeded, "valid independently seeded ring flags");
    }
}
