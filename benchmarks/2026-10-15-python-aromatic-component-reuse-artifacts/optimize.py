from pathlib import Path
s=Path('/private/tmp/chematic-aromatic-components-source-20261015');r=Path('/private/tmp/chematic-aromatic-components-20261015')
p=s/'crates/chematic-perception/src/rdkit_parity/typing.rs';old=p.read_text();t=old
start=old.index('    fn new(',old.index('impl ParityTyping'))
end=old.index('    /// RDKit treats',start)
refnew=old[start:end].replace('    fn new(', 'fn reference_typing(',1).replace('    ) -> Self {','    ) -> ParityTyping {',1).replace('        Self {','        ParityTyping {',1).replace('candidate_component_mask(', 'reference_component_mask(')
start=old.index('fn candidate_component_mask(');end=old.index('fn parity_ring_set(',start)
refmask=old[start:end].replace('fn candidate_component_mask(', 'fn reference_component_mask(',1)
t=t.replace('        let mut components = UnionFind::new(atom_count);','        let components = crate::sssr::ring_components_shared(mol);',1)
t=t.replace('                components.union(bond.atom1.0, bond.atom2.0);\n','',1)
t=t.replace('candidate_component_mask(mol, cyclic.as_slice(), &ring_atom, &candidate, components)', 'candidate_component_mask(mol, cyclic.as_slice(), &ring_atom, &candidate, &components.label)',1)
t=t.replace('    mut components: UnionFind,','    labels: &[u32],',1)
t=t.replace('keep_root[components.find(left as u32) as usize] = true;', 'keep_root[labels[left] as usize] = true;',1)
t=t.replace('ring_atom[atom] && keep_root[components.find(atom as u32) as usize]', 'ring_atom[atom] && keep_root[labels[atom] as usize]',1)
t=t.replace('        let components = crate::sssr::ring_components_shared(mol);','''        // The bridge pass also cached exactly the partition obtained by
        // unioning all cyclic bonds. Reuse it for component membership;
        // candidate-only cycle detection below still needs its own forest.
        let components = crate::sssr::ring_components_shared(mol);''',1)
t+='\n#[cfg(test)]\n#[path = "component_typing_tests.rs"]\nmod component_typing_tests;\n'
p.write_text(t)
test='''//! Differential oracle for reuse of the cyclic-component partition.
//! Uses only the existing licensed 10k corpus, with its pinned manifest.
use super::*;

'''+refnew+refmask+'''
fn compare(input: &Molecule) {
    let type_atom = |idx, cyclic: &[bool]| {
        let donor = get_atom_electron_donor_type(input, idx, cyclic);
        (donor, is_atom_candidate_for_aromaticity(input, idx, donor))
    };
    let expected = reference_typing(input, type_atom);
    for _ in 0..2 {
        let actual = ParityTyping::new(input, type_atom);
        assert_eq!(actual.donor, expected.donor);
        assert_eq!(actual.candidate, expected.candidate);
        assert_eq!(actual.keep, expected.keep);
    }
}

#[test]
fn fixed_ten_thousand_typings_match_original_partition() {
    let corpus = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"));
    let mut rows = 0;
    for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
        let input = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
        compare(&input);
        rows += 1;
    }
    assert_eq!(rows, 10000);
}

#[test]
fn individually_seeded_ring_data_preserves_original_partition() {
    let corpus = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"));
    for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
        let input = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
        let flags = crate::sssr::ring_bond_flags_shared(&input);
        let components = crate::sssr::ring_components_shared(&input);
        let flags_only = input.clone();
        flags_only.seed_derived(chematic_core::DerivedSlot::RingBondFlags, flags);
        compare(&flags_only);
        let components_only = input.clone();
        components_only.seed_derived(chematic_core::DerivedSlot::RingComponents, components);
        compare(&components_only);
    }
}
'''
(s/'crates/chematic-perception/src/rdkit_parity/component_typing_tests.rs').write_text(test)
p=s/'CHANGELOG.md';t=p.read_text();t=t.replace('### Performance\n','### Performance\n\n- Reuse cached cyclic-component labels during RDKit-compatible aromaticity\n  typing, avoiding a redundant union-find while preserving candidate-cycle\n  detection, donor types and ring selection.\n',1);p.write_text(t)
print('Candidate prepared')
