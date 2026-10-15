from pathlib import Path
s=Path('/private/tmp/chematic-ring-scratch-source-20261015');r=Path('/private/tmp/chematic-ring-scratch-20261015')
p=s/'crates/chematic-perception/src/sssr/candidates.rs';old=p.read_text();a=old.index('fn ring_bond_flags_uncached(');b=old.index('/// Return an index-aligned flag for every atom',a);reference=old[a:b].replace('fn ring_bond_flags_uncached(', 'fn reference_ring_flags(',1);t=old
replacements={
'let mut dl = vec![(UNSEEN, UNSEEN); atom_count];':'let mut dl: SmallVec<[(u32, u32); 64]> = smallvec::smallvec![(UNSEEN, UNSEEN); atom_count];',
'let mut stack: Vec<(u32, u32, u32)> = Vec::with_capacity(atom_count.min(64));':'let mut stack: SmallVec<[(u32, u32, u32); 64]> = SmallVec::new();',
'let mut visit: Vec<u32> = Vec::with_capacity(atom_count);':'let mut visit: SmallVec<[u32; 64]> = SmallVec::with_capacity(atom_count);'}
for x,y in replacements.items():assert x in t;t=t.replace(x,y,1)
t+='\n#[cfg(test)]\n#[path = "ring_scratch_tests.rs"]\nmod ring_scratch_tests;\n';p.write_text(t)
test='//! Frozen bridge/component oracle; existing pinned 10k and stress inputs only.\nuse super::*;\n\n'+reference+'''
#[test]
fn inline_and_spilled_bridge_scratch_match_original() {
    let corpora = [include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")), include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi"))];
    let mut rows = 0;
    let mut spilled = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let input = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            let (expected_flags, expected_components) = reference_ring_flags(&input);
            let (actual_flags, actual_components) = ring_bond_flags_uncached(&input);
            assert_eq!(actual_flags, expected_flags);
            assert_eq!(actual_components.label, expected_components.label);
            assert_eq!(actual_components.count, expected_components.count);
            spilled += usize::from(input.atom_count() > 64);
            rows += 1;
        }
    }
    assert_eq!(rows,14072);
    assert!(spilled > 0 && spilled < rows);
}
''';(s/'crates/chematic-perception/src/sssr/ring_scratch_tests.rs').write_text(test)
p=s/'crates/chematic-perception/src/rdkit_parity/view.rs';old=p.read_text();a=old.index('struct ComponentInfo {');b=old.index('/// Some cyclic bond of the component is aromatic.',a);refstruct=old[a:b].replace('struct ComponentInfo','struct ReferenceInfo',1);x=old.index('    fn new(mol: &Molecule) -> Self {',old.index('impl ComponentInfo'));y=old.index('    /// Labels of the cyclic components',x);reffn=old[x:y].replace('fn new(mol: &Molecule) -> Self','fn reference_info(mol: &Molecule) -> ReferenceInfo',1).replace('        Self {','        ReferenceInfo {',1)
t=old
for field in ['root','edges','atoms']:t=t.replace(field+': Vec<u32>,',field+': smallvec::SmallVec<[u32; 64]>,',1)
t=t.replace('kind: Vec<u8>,','kind: smallvec::SmallVec<[u8; 64]>,',1)
for field,ty,value,size in [('root','u32','u32::MAX','n'),('kind','u8','0u8','c'),('edges','u32','0u32','c'),('atoms','u32','0u32','c')]:
 oldline=f'let mut {field} = vec![{value}; {size}];';newline=f'let mut {field}: smallvec::SmallVec<[{ty}; 64]> = smallvec::smallvec![{value}; {size}];';assert oldline in t;t=t.replace(oldline,newline,1)
t+='\n#[cfg(test)]\n#[path = "component_scratch_tests.rs"]\nmod component_scratch_tests;\n';p.write_text(t)
test='//! Frozen component facts; existing pinned 10k and stress inputs only.\nuse super::*;\n\n'+refstruct+reffn+'''
#[test]
fn inline_and_spilled_component_facts_match_original() {
    let corpora = [include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")), include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi"))];
    let mut rows = 0;
    let mut spilled = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let input = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            let expected = reference_info(&input);
            for _ in 0..2 {
                let actual = ComponentInfo::new(&input);
                assert_eq!(actual.cyclic, expected.cyclic);
                assert_eq!(actual.root.as_slice(), expected.root.as_slice());
                assert_eq!(actual.kind.as_slice(), expected.kind.as_slice());
                assert_eq!(actual.edges.as_slice(), expected.edges.as_slice());
                assert_eq!(actual.atoms.as_slice(), expected.atoms.as_slice());
                assert_eq!(actual.has_aromatic_bond, expected.has_aromatic_bond);
                assert_eq!(actual.explicit, expected.explicit);
                assert_eq!(actual.non_ring_aromatic_bond, expected.non_ring_aromatic_bond);
                assert_eq!(actual.any_flagged, expected.any_flagged);
            }
            spilled += usize::from(input.atom_count() > 64);
            rows += 1;
        }
    }
    assert_eq!(rows,14072);
    assert!(spilled > 0 && spilled < rows);
}
''';(s/'crates/chematic-perception/src/rdkit_parity/component_scratch_tests.rs').write_text(test)
p=s/'CHANGELOG.md';p.write_text(p.read_text().replace('### Performance\n','### Performance\n\n- Keep small bridge-search and aromatic-component scratch buffers inline,\n  spilling to dynamic storage for larger graphs without changing cyclic\n  components, aromaticity decisions or traversal order.\n',1));print('Inline bridge/component scratch candidate written; public cached vectors unchanged')
