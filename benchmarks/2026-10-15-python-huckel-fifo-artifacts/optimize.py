from pathlib import Path
s=Path('/private/tmp/chematic-huckel-fifo-source-20261015');r=Path('/private/tmp/chematic-huckel-fifo-20261015');p=s/'crates/chematic-perception/src/rdkit_model_hook.rs';old=p.read_text();assert old==(r/'baseline-rdkit_model_hook.rs').read_text();start=old.index('fn aromatic_systems_are_huckel(');tail=old[start:];newtail=tail.replace('    let mut stack = Vec::new();','    let mut stack = Vec::new();\n    let mut comp = Vec::new();',1).replace('        let mut comp = Vec::new();','        comp.clear();',1)
newtail=newtail.replace('let mut queue = std::collections::VecDeque::with_capacity(m);','let mut queue = Vec::with_capacity(m);').replace('queue.push_back(u);','queue.push(u);\n            let mut head = 0;').replace("'bfs: while let Some(x) = queue.pop_front() {","'bfs: while head < queue.len() {\n                let x = queue[head];\n                head += 1;").replace('queue.push_back(y);','queue.push(y);')
p.write_text(old[:start]+newtail+'\n#[cfg(test)]\n#[path = "rdkit_model_hook_storage_tests.rs"]\nmod storage_tests;\n')
a=old.index('fn rdkit_model_may_disagree_uncached(');b=old.index('fn has_cyclic_triple_bond(',a);reference=old[a:b]+tail;reference=reference.replace('rdkit_model_may_disagree_uncached','reference_gate').replace('aromatic_systems_are_huckel','reference_aromatic_systems').replace('fused_rings_are_huckel','reference_fused_rings')
test='''//! Differential storage oracle frozen from the immediate ten-candidate baseline.
//! All molecular inputs reuse the existing pinned 10k and stress corpora.
use super::*;

'''+reference+'''
#[test]
fn fixed_inputs_match_original_huckel_gate() {
    let corpora = [
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")),
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi")),
    ];
    let mut rows = 0;
    let mut fused_inputs = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let mol = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            let expected = reference_gate(&mol);
            assert_eq!(aromatic_systems_are_huckel(&mol), reference_aromatic_systems(&mol));
            assert_eq!(rdkit_model_may_disagree_uncached(&mol), expected);
            for _ in 0..2 { assert_eq!(rdkit_model_may_disagree(&mol), expected); }
            fused_inputs += usize::from(mol.bonds().filter(|(_, b)| b.order == BondOrder::Aromatic).count() > mol.atoms().filter(|(_, a)| a.aromatic).count());
            rows += 1;
        }
    }
    assert_eq!(rows,14072);
    assert!(fused_inputs > 0, "exercise fused aromatic breadth-first searches");
}
'''
(s/'crates/chematic-perception/src/rdkit_model_hook_storage_tests.rs').write_text(test)
p=s/'CHANGELOG.md';p.write_text(p.read_text().replace('### Performance\n','### Performance\n\n- Reuse aromatic-component scratch storage and a contiguous FIFO during the\n  RDKit model gate, preserving breadth-first traversal order, electron counts\n  and all aromaticity decisions.\n',1));print('Aromatic component reuse and contiguous FIFO candidate written')
