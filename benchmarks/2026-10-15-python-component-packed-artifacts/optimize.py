from pathlib import Path
import re
s=Path('/private/tmp/chematic-component-packed-source-20261015');r=Path('/private/tmp/chematic-component-packed-20261015')
p=s/'crates/chematic-perception/src/rdkit_parity/view.rs';old=p.read_text();a=old.index('struct ComponentInfo {');b=old.index('/// Some cyclic bond of the component is aromatic.',a);refstruct=old[a:b].replace('struct ComponentInfo','struct ReferenceInfo',1);x=old.index('    fn new(mol: &Molecule) -> Self {',old.index('impl ComponentInfo'));y=old.index('    /// Labels of the cyclic components',x);reffn=old[x:y].replace('fn new(mol: &Molecule) -> Self','fn reference_info(mol: &Molecule) -> ReferenceInfo',1).replace('        Self {','        ReferenceInfo {',1)
t=old.replace('    /// Component label of every atom with a cyclic bond (`u32::MAX` otherwise).\n    root: Vec<u32>,','    /// Adjacent root labels (n), cyclic-bond counts (c), atom counts (c).\n    /// Roots use `u32::MAX` for atoms without a cyclic bond.\n    facts: Vec<u32>,\n    atom_count: usize,',1)
t=t.replace('    /// Per label: cyclic bonds and atoms.\n    edges: Vec<u32>,\n    atoms: Vec<u32>,\n','',1)
t=t.replace('        let mut root = vec![u32::MAX; n];\n        let mut kind = vec![0u8; c];\n        let mut edges = vec![0u32; c];\n        let mut atoms = vec![0u32; c];','        let mut facts = vec![0u32; n + 2 * c];\n        let (root, counts) = facts.split_at_mut(n);\n        root.fill(u32::MAX);\n        let (edges, atoms) = counts.split_at_mut(c);\n        let mut kind = vec![0u8; c];',1)
t=t.replace('            root,\n            kind,\n            edges,\n            atoms,','            facts,\n            atom_count: n,\n            kind,',1)
t=t.replace('(0..self.edges.len()).filter(|&c| self.edges[c] > 0)','(0..self.kind.len()).filter(|&c| self.facts[self.atom_count + c] > 0)',1)
t=re.sub(r'info\.root\[([^\]]+)\]',r'info.facts[\1]',t)
t=t.replace('info.edges[r] != info.atoms[r]','info.facts[n + r] != info.facts[n + info.kind.len() + r]',1)
t=t.replace('candidate.iter_mut().zip(&info.root)','candidate.iter_mut().zip(&info.facts[..n])',1)
assert 'info.root[' not in t and 'info.edges' not in t and 'info.atoms' not in t
assert t!=old
t+='\n#[cfg(test)]\n#[path = "component_packed_tests.rs"]\nmod component_packed_tests;\n';p.write_text(t)
test='//! Frozen component facts; existing pinned 10k and stress inputs only.\nuse super::*;\n\n'+refstruct+reffn+'''
#[test]
fn contiguous_component_facts_match_original() {
    let corpora = [include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")), include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi"))];
    let mut rows = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let input = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            let expected = reference_info(&input);
            for _ in 0..2 {
                let actual = ComponentInfo::new(&input);
                let n = input.atom_count();
                let c = expected.kind.len();
                assert_eq!(actual.atom_count, n);
                assert_eq!(actual.facts.len(), n + 2 * c);
                assert_eq!(actual.cyclic, expected.cyclic);
                assert_eq!(&actual.facts[..n], expected.root.as_slice());
                assert_eq!(actual.kind, expected.kind);
                assert_eq!(&actual.facts[n..n+c], expected.edges.as_slice());
                assert_eq!(&actual.facts[n+c..], expected.atoms.as_slice());
                assert_eq!(actual.roots().collect::<Vec<_>>(), (0..c).filter(|&i|expected.edges[i]>0).collect::<Vec<_>>());
                assert_eq!(actual.has_aromatic_bond, expected.has_aromatic_bond);
                assert_eq!(actual.explicit, expected.explicit);
                assert_eq!(actual.non_ring_aromatic_bond, expected.non_ring_aromatic_bond);
                assert_eq!(actual.any_flagged, expected.any_flagged);
            }
            rows += 1;
        }
    }
    assert_eq!(rows,14072);
}
''';(s/'crates/chematic-perception/src/rdkit_parity/component_packed_tests.rs').write_text(test)
p=s/'CHANGELOG.md';p.write_text(p.read_text().replace('### Performance\n','### Performance\n\n- Combine private aromatic-component root and count buffers into one allocation\n  with the same total number of integer elements, preserving component labels,\n  aromaticity decisions and traversal order.\n',1));print('Contiguous component facts written; n+2c integer elements unchanged')
