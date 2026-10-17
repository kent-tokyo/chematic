from pathlib import Path
import json,hashlib,subprocess
R=Path('/private/tmp/chematic-python-asa-preparation-20261017');O=Path(__file__).resolve().parent
p=R/'crates/chematic-perception/src/rdkit_parity/view.rs';s=p.read_text();start=s.index('fn parity_shortcut(');end=s.index('\n/// Whether the RDKit-parity aromatic view',start);frozen=s[start:end].replace('fn parity_shortcut(', 'fn frozen_parity_shortcut(',1)
marker='        let info = ComponentInfo::new(mol);\n        if has_aromatic_macrocycle_bridge_donor(mol) {';assert s.count(marker)==1
s=s.replace(marker,'''        // Uniform cyclic components already satisfy the complete-preperceived
        // test. Keep the ring flags/components cache, but avoid temporary
        // component facts when no mixed or multiple-bond ring can occur.
        if !model_may_disagree
            && uniform_components::preperceived(mol)
            && !has_aromatic_macrocycle_bridge_donor(mol)
        {
            return ParityShortcut::Identity;
        }
'''+marker)
s+='\n#[path = "uniform_components.rs"]\nmod uniform_components;\n\n#[cfg(test)]\n#[path = "uniform_component_tests.rs"]\nmod uniform_component_tests;\n';p.write_text(s)
(R/'crates/chematic-perception/src/rdkit_parity/uniform_components.rs').write_text('''//! Exact complete-preperceived case without temporary component facts.
use super::{AtomIdx, BondOrder, Molecule};

/// The caller has rejected RDKit-model disagreement. Every cyclic component
/// here is either fully explicit aromatic or saturated with no incident
/// multiple bond. Thus the original complete-preperceived test returns
/// Identity; mixed components, cyclic multiples and exocyclic donor cases
/// keep that original test. Ring cache preparation remains unchanged.
pub(super) fn preperceived(mol: &Molecule) -> bool {
    let cyclic = crate::sssr::ring_bond_flags_shared(mol);
    let non_aromatic_ring_atom = |atom: AtomIdx| {
        !mol.atom(atom).aromatic
            && mol.neighbor_slice(atom).iter().any(|&(_, b)| cyclic[b.0 as usize])
    };
    for (idx, bond) in mol.bonds() {
        let (a, b) = (bond.atom1, bond.atom2);
        let aromatic_a = mol.atom(a).aromatic;
        let aromatic_b = mol.atom(b).aromatic;
        match bond.order {
            BondOrder::Aromatic => {
                if !cyclic[idx.0 as usize] || !aromatic_a || !aromatic_b {
                    return false;
                }
            }
            BondOrder::Single | BondOrder::Up | BondOrder::Down => {
                if cyclic[idx.0 as usize] && (aromatic_a || aromatic_b) {
                    return false;
                }
            }
            BondOrder::Double | BondOrder::Triple => {
                if cyclic[idx.0 as usize]
                    || non_aromatic_ring_atom(a)
                    || non_aromatic_ring_atom(b)
                {
                    return false;
                }
            }
            _ => return false,
        }
    }
    // Preserve ComponentInfo's explicit-representation consistency check.
    !mol.atoms().any(|(idx, atom)| {
        atom.aromatic
            && atom.charge == 0
            && matches!(atom.element.atomic_number(), 8 | 16 | 34 | 52)
            && mol.degree(idx) + chematic_core::implicit_hcount(mol, idx) as usize > 2
    })
}
''')
tests='''use super::*;

FROZEN

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
'''.replace('FROZEN',frozen)
(R/'crates/chematic-perception/src/rdkit_parity/uniform_component_tests.rs').write_text(tests)
(O/'frozen-reference.json').write_text(json.dumps({'base_commit':(O/'baseline-revision.txt').read_text().strip(),'function':'parity_shortcut','original_sha256':hashlib.sha256(frozen.replace('fn frozen_parity_shortcut(', 'fn parity_shortcut(',1).encode()).hexdigest(),'change':'Conservative structural proof of the existing complete-preperceived Identity; full fallback unchanged'},indent=2)+'\n')
