from pathlib import Path
import hashlib, json, re, subprocess
ROOT=Path('/private/tmp/chematic-sdf-buffer-scan-20261016')
OUT=Path(__file__).resolve().parent
path='crates/chematic-perception/src/rdkit_model_hook.rs'
base=json.loads((OUT/'measurement-plan.json').read_text())['base']
source=subprocess.check_output(['git','show',base+':'+path],cwd=ROOT,text=True)
names=['rdkit_model_may_disagree_uncached','aromatic_systems_are_huckel','fused_rings_are_huckel']
proof={}
parts=['//! Frozen complete model-gate functions from main 3e7f04b2.\nuse super::*;\n']
for name in names:
    begin=source.index('fn '+name+'(')
    pos=source.index('{',begin); depth=1; end=pos+1
    while depth:
        depth += (source[end]=='{')-(source[end]=='}'); end+=1
    original=source[begin:end]
    frozen=original
    for rename in names: frozen=re.sub(r'\b'+rename+r'\b','frozen_'+rename,frozen)
    restored=frozen
    for rename in names: restored=restored.replace('frozen_'+rename,rename)
    assert restored==original
    proof[name]={'source_sha256':hashlib.sha256(original.encode()).hexdigest(),'only_identifiers_renamed':True}
    parts.append(frozen+'\n')
tests=r'''
fn compare(mol: &Molecule) {
    let baseline = mol.clone();
    let candidate = mol.clone();
    for _ in 0..2 {
        assert_eq!(rdkit_model_may_disagree_uncached(&candidate),
                   frozen_rdkit_model_may_disagree_uncached(&baseline));
        assert_eq!(aromatic_systems_are_huckel(&candidate),
                   frozen_aromatic_systems_are_huckel(&baseline));
    }
}

#[test]
fn all_pinned_model_gate_outcomes_are_unchanged() {
    let corpora = [
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")),
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi")),
    ];
    let mut rows=0;
    let mut inline=0;
    let mut spilled=0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let mol=chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            inline += usize::from(mol.atom_count() <= 64);
            spilled += usize::from(mol.atom_count() > 64);
            compare(&mol);
            rows+=1;
        }
    }
    assert_eq!(rows,14072);
    assert!(inline>0 && spilled>0);
}

#[test]
fn unusual_and_disconnected_graphs_keep_gate_decisions() {
    compare(&chematic_core::MoleculeBuilder::new().build());
    for smiles in ["*", "[H]", "[C]", "[nH]1cccc1", "c1cccc1", "c1ccccccc1",
                   "c1ccccc#1", "C1=CSC#C1", "c1ccccc1.c1ccccc1", "[Na+].[Cl-]",
                   "[O-][Cl](=O)(=O)=O", "N->[Cu]<-N", "[13CH3][NH3+]",
                   "o1cccccccc1", "s1cccccccc1", "O=C1NNC(=O)N1"] {
        compare(&chematic_smiles::parse_template(smiles).unwrap());
    }
}

fn fused_hexagons(rings: usize, spectators: usize) -> Molecule {
    use chematic_core::{Atom, MoleculeBuilder};
    let mut b=MoleculeBuilder::new();
    let add=|b: &mut MoleculeBuilder| {
        let mut atom=Atom::new(Element::C);
        atom.aromatic=true;
        b.add_atom(atom)
    };
    let first: Vec<_>=(0..6).map(|_| add(&mut b)).collect();
    for i in 0..6 { b.add_bond(first[i], first[(i+1)%6], BondOrder::Aromatic).unwrap(); }
    let mut edge=(first[0],first[1]);
    for _ in 1..rings {
        let fresh: Vec<_>=(0..4).map(|_| add(&mut b)).collect();
        let path=[edge.0,fresh[0],fresh[1],fresh[2],fresh[3],edge.1];
        for pair in path.windows(2) { b.add_bond(pair[0],pair[1],BondOrder::Aromatic).unwrap(); }
        edge=(fresh[1],fresh[2]);
    }
    for _ in 0..spectators { b.add_atom(Atom::new(Element::C)); }
    b.build()
}

#[test]
fn fused_bfs_keeps_paths_across_inline_and_spill_boundaries() {
    for rings in [2,8,15,16,17,31,32] {
        for spectators in [0,64] {
            let mol=fused_hexagons(rings,spectators);
            let n=6+4*(rings-1);
            let comp: Vec<_>=(0..n).collect();
            for style in 0..3 {
                let electrons: Vec<_>=(0..mol.atom_count()).map(|i| match style {
                    0=>1, 1=>(i%3) as u32, _=>2,
                }).collect();
                assert_eq!(fused_rings_are_huckel(&mol,&comp,&electrons),
                           frozen_fused_rings_are_huckel(&mol,&comp,&electrons),
                           "rings={rings}, spectators={spectators}, style={style}");
            }
            compare(&mol);
        }
    }
}
'''
parts.append(tests)
dest=ROOT/'crates/chematic-perception/src/rdkit_model_hook_inline_tests.rs'
dest.write_text('\n'.join(parts))
(OUT/'frozen-oracle-proof.json').write_text(json.dumps({'base':base,'file':path,'functions':proof},indent=2)+'\n')
print('Frozen',len(names),'complete source functions')
