from pathlib import Path
s=Path('/private/tmp/chematic-basis-roots-source-20261015');r=Path('/private/tmp/chematic-basis-roots-20261015');p=s/'crates/chematic-perception/src/sssr/basis.rs';old=p.read_text();a=old.index('fn sssr_rings_in_components(');b=old.index('/// Appends the minimum cycle basis',a);reference=old[a:b].replace('fn sssr_rings_in_components(', 'fn reference_rings_in_components(',1)
x=old.index('    let mut edges = vec![0u32; n];',a);y=old.index('    let mut rings: Vec<Vec<AtomIdx>>',x)
new='''    let mut edges = vec![0u32; n];
    let mut atoms = vec![0u32; n];
    let mut root_of = vec![u32::MAX; n];
    // Every non-root participated in a selected cyclic edge. Its final root
    // is unchanged by path compression; count each non-root and its root once.
    for a in 0..n {
        if parent[a] == a as u32 {
            continue;
        }
        let root = find(&mut parent, a as u32) as usize;
        root_of[a] = root as u32;
        atoms[root] += 1;
        if root_of[root] == u32::MAX {
            root_of[root] = root as u32;
            atoms[root] += 1;
        }
    }
    for (bidx, b) in mol.bonds() {
        if cyclic[bidx.0 as usize] && keep[b.atom1.0 as usize] {
            let a = b.atom1.0 as usize;
            // A selected singleton self-edge has no non-root to mark its
            // endpoint. Keep that existing builder/cache boundary as well.
            if root_of[a] == u32::MAX {
                root_of[a] = a as u32;
                atoms[a] = 1;
            }
            edges[root_of[a] as usize] += 1;
        }
    }
'''
t=old[:x]+new+old[y:];needle='.map(|a| is_ring_atom[a] && general[root_of[a] as usize])';assert t.count(needle)==1;t=t.replace(needle,'.map(|a| root_of[a] != u32::MAX && general[root_of[a] as usize])',1);t+='\n#[cfg(test)]\n#[path = "basis_roots_tests.rs"]\nmod basis_roots_tests;\n';p.write_text(t)
test='//! Frozen component-basis oracle; existing pinned 10k and stress corpora only.\nuse super::*;\n\n'+reference+'''
#[test]
fn component_roots_keep_exact_ring_selection_and_order() {
    let corpora = [include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")), include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi"))];
    let mut rows = 0;
    let mut comparisons = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let input = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            let components = ring_components_shared(&input);
            let masks = [vec![true;input.atom_count()],components.label.iter().map(|&root| root % 2 == 0).collect(),components.label.iter().map(|&root| root % 2 != 0).collect()];
            for keep in masks {
                for allow_ranks in [false,true] {
                    // Clones start with independent empty derived caches.
                    let original = input.clone();
                    let actual = input.clone();
                    for _ in 0..2 {
                        let expected = reference_rings_in_components(&original,&keep,false,allow_ranks);
                        let observed = sssr_rings_in_components(&actual,&keep,false,allow_ranks);
                        assert_eq!(observed,expected,"row {rows}, ranked {allow_ranks}");
                        comparisons += 1;
                    }
                }
            }
            rows += 1;
        }
    }
    assert_eq!(rows,14072);
    assert_eq!(comparisons,168864);
}
''';(s/'crates/chematic-perception/src/sssr/basis_roots_tests.rs').write_text(test)
p=s/'CHANGELOG.md';p.write_text(p.read_text().replace('### Performance\n','### Performance\n\n- Reuse final union-find roots while counting selected cyclic components,\n  avoiding repeated edge-root searches and a temporary ring-atom flag buffer\n  without changing the selected ring sequences or canonical tie-breaking.\n',1));print('Component basis roots candidate written; original function frozen; no public API or traversal change')
