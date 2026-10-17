use super::*;

fn frozen_parity_shortcut(mol: &Molecule) -> ParityShortcut {
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


fn compare(mol: &Molecule) {
    let expected = mol.clone();
    let actual = mol.clone();
    assert_eq!(parity_shortcut(&actual), frozen_parity_shortcut(&expected));
    assert_eq!(parity_shortcut(&actual), frozen_parity_shortcut(&expected));
}

#[test]
fn all_14072_graphs_and_kekule_forms_keep_the_frozen_shortcut() {
    let corpora = [
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")),
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi")),
    ];
    let mut count = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let mol = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            compare(&mol);
            let kekule = chematic_core::kekulize(&mol).map(|result| chematic_core::apply_kekule(&mol, &result)).unwrap_or_else(|_| mol.clone());
            compare(&kekule);
            count += 1;
        }
    }
    assert_eq!(count, 14072);
}

#[test]
fn fused_mixed_exocyclic_macrocycle_and_unusual_models_keep_fallbacks() {
    for smiles in ["c1ccccc1", "c1ccccc1-c2ccccc2", "CC(=O)Nc1ccccc1", "c1ccccc1C1CCCCC1",
        "c1ccc2ccccc2c1", "C1=CC=CC=C1", "O=C1NNC(=O)N1", "O=C1CCCCC1",
        "c1ccc2c(c1)-c1ccccc1-2", "c1ccccccco1", "c1cccccccs1", "[o+]1ccccc1",
        "[sH]1cccc1", "[se]1cccc1", "[nH]1cccc1", "[H]c1ccccc1", "[Fe]", "*c1ccccc1"] {
        compare(&chematic_smiles::parse_template(smiles).unwrap());
    }
}

#[test]
fn small_graphs_flags_and_hydrogen_overrides_keep_frozen_decisions() {
    use chematic_core::{Atom, Element, MoleculeBuilder};
    let edges: Vec<_>=(0..5).flat_map(|a| ((a+1)..5).map(move |b| (a,b))).collect();
    for mask in 0..(1u32<<edges.len()) {
        for pattern in 0..4 {
            let mut builder=MoleculeBuilder::new();
            for a in 0..5 {
                let mut atom=Atom::new(if a==0 && pattern==3 {Element::S} else {Element::C});
                atom.aromatic=pattern!=0 && (pattern!=2 || a!=0);
                if pattern==3 {atom.hydrogen_count=Some(1);}
                builder.add_atom(atom);
            }
            for (i,&(a,b)) in edges.iter().enumerate() {
                if mask>>i&1==1 {
                    let order=if pattern==0 || (pattern==2 && a==0) {BondOrder::Single} else {BondOrder::Aromatic};
                    builder.add_bond(AtomIdx(a),AtomIdx(b),order).unwrap();
                }
            }
            compare(&builder.build());
        }
    }
}
