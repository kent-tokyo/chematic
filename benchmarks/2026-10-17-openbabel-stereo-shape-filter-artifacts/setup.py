from pathlib import Path
import json, shutil, subprocess, hashlib
R=Path('/private/tmp/chematic-stereo-shape-filter-20261017'); O=Path(__file__).resolve().parent
p=R/'crates/chematic-perception/src/stereo_validation.rs'; original=p.read_text()
start=original.index('pub fn stereo_centers('); end=original.index('\n/// Summarise',start)
frozen=original[start:end].replace('pub fn stereo_centers(', 'fn frozen_stereo_centers(',1)
base=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip()
B=O/'baseline-source'; B.mkdir()
for d in ['crates','tools']: shutil.copytree(R/d,B/d)
for f in ['Cargo.toml','Cargo.lock','README.md','LICENSE','LICENSE-MIT','LICENSE-APACHE']:
 if (R/f).exists(): shutil.copy2(R/f,B/f)
block='''        if implicit_h > 0 && heavy_nbs.len() == 3 {
            // An implicit H is a fourth group and is distinct from every
            // heavy-atom branch. A ring nitrogen's fourth group is its lone
            // pair, handled by the same branch-distinctness condition.
        } else if heavy_nbs.len() != 4 && !is_ring_nitrogen {
            continue;
        }
'''
assert original.count(block)==1
changed=original.replace(block,'')
marker='        // Compare heavy-atom branches from this center.'
assert changed.count(marker)==1
changed=changed.replace(marker,block+'\n'+marker)
p.write_text(changed+'\n#[cfg(test)]\n#[path = "stereo_shape_tests.rs"]\nmod shape_tests;\n')
tests=r'''
use super::*;
use chematic_core::{Atom, Element, MoleculeBuilder};

FROZEN

fn compare(mol: &Molecule) {
    let before: Vec<_> = mol.atoms().map(|(i,a)| format!("{i:?}:{a:?}:{:?}",mol.stereo_neighbor_order(i))).collect();
    assert_eq!(stereo_centers(mol), frozen_stereo_centers(mol));
    assert_eq!(before, mol.atoms().map(|(i,a)| format!("{i:?}:{a:?}:{:?}",mol.stereo_neighbor_order(i))).collect::<Vec<_>>());
}

#[test]
fn all_14072_pinned_graphs_match_frozen_classifier() {
    let corpora = [
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")),
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi")),
    ];
    let mut count=0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            compare(&chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap());
            count+=1;
        }
    }
    assert_eq!(count,14072);
}

#[test]
fn every_five_atom_graph_and_hydrogen_override_matches() {
    let edges: Vec<_>=(0..5).flat_map(|a| ((a+1)..5).map(move |b| (a,b))).collect();
    for mask in 0..(1u32<<edges.len()) {
        for hydrogen_count in [None, Some(0), Some(1), Some(2), Some(3), Some(4)] {
            let mut b=MoleculeBuilder::new();
            for _ in 0..5 {
                let mut atom=Atom::new(Element::C);
                atom.hydrogen_count=hydrogen_count;
                b.add_atom(atom);
            }
            for (i,&(a,c)) in edges.iter().enumerate() {
                if mask>>i&1==1 { b.add_bond(AtomIdx(a),AtomIdx(c),BondOrder::Single).unwrap(); }
            }
            compare(&b.build());
        }
    }
}

#[test]
fn stereocenter_boundary_groups_and_ring_nitrogen_match() {
    for smiles in ["C", "CCCCCC", "CC(F)CCCCCCCCCCCCCCCC", "c1ccccc1", "N[C@@H](C)C(=O)O",
        "[H][C@](F)(Cl)Br", "FC(Cl)(Br)I", "F[C@](Cl)(Br)I", "F[C@@](Cl)(Br)I",
        "C[N+](F)(Cl)Br", "N1CC2CCC1C2", "C1CC2CCC1C2", "C1C2CC3CC1CC(C2)C3",
        "C[S+](CC)CCC", "CS(=O)CC", "C[Se+](CC)CCC", "C[Se](=O)CC", "[13CH3]C(F)C",
        "[CH2]", "[C]", "[NH4+]", "[H]N1CC2CCC1C2", "[Fe]", "ClC1CC(F)(Br)CC1I"] {
        compare(&chematic_smiles::parse_template(smiles).unwrap());
    }
}
'''.replace('FROZEN',frozen)
(R/'crates/chematic-perception/src/stereo_shape_tests.rs').write_text(tests)
(O/'source-proof.json').write_text(json.dumps({'base_commit':base,'frozen_body_sha256':hashlib.sha256(frozen.replace('fn frozen_stereo_centers(', 'pub fn stereo_centers(',1).encode()).hexdigest(),'transformation':'Move unchanged local group-shape rejection before substituent_signatures_distinct; unchanged branch signature algorithm'},indent=2)+'\n')
print(base)
