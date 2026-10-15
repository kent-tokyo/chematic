from pathlib import Path
s=Path('/private/tmp/chematic-ring-component-source-20261015');r=Path('/private/tmp/chematic-ring-component-20261015')
f='crates/chematic-perception/src/sssr/candidates.rs';p=s/f;(r/'baseline-candidates.rs').write_bytes(p.read_bytes());t=p.read_text();t=t.replace('pub(crate) struct RingComponents {','pub struct RingComponents {',1).replace('pub(crate) fn ring_components_shared(', 'pub fn ring_components_shared(',1)
t=t.replace('#[derive(Debug)]\npub struct RingComponents', '#[derive(Debug)]\n#[non_exhaustive]\npub struct RingComponents', 1)
t=t.replace('/// Memoized [`RingComponents`] of `mol`, computed together with\n/// [`ring_bond_flags_shared`].','/// Shared, memoized [`RingComponents`] of `mol`, computed together with\n/// [`ring_bond_flags_shared`]. Labels identify connected components after\n/// removing bridges and ring-ineligible bonds; they are not stable atom IDs.\n/// Graph mutations invalidate both cached values. A graph with externally\n/// seeded bond flags may compute its component labels on this first call.')
p.write_text(t)
p=s/'crates/chematic-chem/src/descriptors/rings_and_filters.rs';old=p.read_text();t=old
start=t.index('    let n = mol.atom_count();',t.index('fn simple_ring_system_count'))
end=t.index('    let mut rings = 0usize;',start)
t=t[:start]+'''    // The bridge pass that produced `flags` also labels the cyclic
    // components. Reuse those labels instead of rebuilding a union-find.
    let components = chematic_perception::sssr::ring_components_shared(mol);
    let labels = &components.label;
    let mut balance = vec![0i32; components.count as usize];
    let mut on_ring = vec![false; mol.atom_count()];
    let mut roots = Vec::new();
    for (bi, bond) in mol.bonds() {
        if !flags[bi.0 as usize] {
            continue;
        }
        let r = labels[bond.atom1.0 as usize] as usize;
        // Before counting atoms, a zero balance means this is the first
        // cyclic bond in the component. Non-ring singleton labels are never
        // added, so they cannot be mistaken for a simple cycle below.
        if balance[r] == 0 {
            roots.push(r);
        }
        balance[r] -= 1;
        on_ring[bond.atom1.0 as usize] = true;
        on_ring[bond.atom2.0 as usize] = true;
    }
    for (a, &is_on_ring) in on_ring.iter().enumerate() {
        if is_on_ring {
            balance[labels[a] as usize] += 1;
        }
    }
''' +t[end:]
t=t.replace('two_cycle_system_rings(mol, &flags, &mut parent, r)', 'two_cycle_system_rings(mol, &flags, labels, r)').replace('system_rings_rdkit(mol, &flags, &mut parent, r)', 'system_rings_rdkit(mol, &flags, labels, r)')
# Both subordinate routines only use union-find to test component membership.
for fn in ['system_rings_rdkit','two_cycle_system_rings']:
 start=t.index('fn '+fn+'(');end=t.index('\n}\n',start)+3;part=t[start:end]
 part=part.replace('parent: &mut [u32],','labels: &[u32],')
 a=part.index('    fn find(');b=part.index('\n    }\n',a)+7;part=part[:a]+part[b:]
 part=part.replace('find(parent, bond.atom1.0) as usize != root','labels[bond.atom1.0 as usize] as usize != root').replace('find(parent, idx.0) as usize != root','labels[idx.0 as usize] as usize != root')
 t=t[:start]+part+t[end:]
t+='\n#[cfg(test)]\n#[path = "ring_component_count_tests.rs"]\nmod ring_component_count_tests;\n';p.write_text(t)
# Preserve the entire pre-change uncached count and its graph helpers in tests.
a=old.index('fn rdkit_num_rings_uncached(');b=old.index('/// RDKit\'s ring information (',a)
reference=old[a:b]
reference=reference.replace('rdkit_num_rings_uncached','reference_count')
(s/'crates/chematic-chem/src/descriptors/ring_component_count_tests.rs').write_text('//! Differential oracle for cyclic-component reuse.\nuse super::*;\n\n'+reference+'''
fn compare(input: &Molecule) {
    // Fresh clones drop derived caches, then exercise both cold and warm
    // calls. Compare again after graph edits below to cover invalidation.
    let expected = reference_count(&input.clone());
    assert_eq!(rdkit_num_rings(input), expected);
    assert_eq!(rdkit_num_rings(input), expected);
}

#[test]
fn fixed_ten_thousand_ring_counts_match_reference() {
    let corpus = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"));
    let mut rows = 0;
    for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
        let mol = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
        compare(&mol);
        rows += 1;
    }
    assert_eq!(rows, 10000);
}

#[test]
fn edited_graphs_keep_ring_counts_and_invalidate_components() {
    for smiles in ["CCC", "C1CCCCC1", "C1CCCCC1.C1CC1", "C1CCC2(CC1)CCCC2", "C1CC2CCC1C2", "C12C3C4C1C5C2C3C45", "C1=CC2=CC=CC=C2C=C1", "N1CCN->[Cu]1", "[Fe]1(Cl)(Cl)Cl.Cl1"] {
        let mut mol = chematic_smiles::parse_template(smiles).unwrap();
        compare(&mol);
        let before = chematic_perception::sssr::ring_components_shared(&mol);
        let branch = mol.add_atom(chematic_core::Atom::new(Element::C));
        mol.add_bond(AtomIdx(0), branch, BondOrder::Single).unwrap();
        compare(&mol);
        let after = chematic_perception::sssr::ring_components_shared(&mol);
        assert!(!std::sync::Arc::ptr_eq(&before, &after));
        if mol.atom_count() > 3 {
            mol.add_bond(branch, AtomIdx(2), BondOrder::Single).unwrap();
            compare(&mol);
        }
        let bond = mol.bond(BondIdx(0)).clone();
        mol.remove_bond(BondIdx(0));
        compare(&mol);
        mol.add_bond(bond.atom1, bond.atom2, BondOrder::QueryAny).unwrap();
        compare(&mol);
    }
}

#[test]
fn seeded_bond_flags_without_components_keep_the_partition() {
    let mol = chematic_smiles::parse_template("C1CCC2(CC1)CCC2.CC").unwrap();
    let flags = chematic_perception::ring_bond_flags_shared(&mol);
    let copy = mol.clone();
    copy.seed_derived(chematic_core::DerivedSlot::RingBondFlags, flags);
    compare(&copy);
}
''')
p=s/'CHANGELOG.md';t=p.read_text();t=t.replace('### Performance\n','### Performance\n\n- Reuse the bridge pass\'s cached cyclic components in RDKit-compatible ring\n  counts, avoiding a second union-find without changing ring selection or\n  macrocycle and coordination-bond fallbacks.\n',1);p.write_text(t)
print('Candidate edits prepared')
p=s/'crates/chematic-chem/src/descriptors/ring_component_count_tests.rs'
with p.open('a') as output:
 output.write('''
#[test]
fn every_five_atom_simple_graph_matches_reference() {
    let edges: Vec<_> = (0..5).flat_map(|a| (a + 1..5).map(move |b| (a, b))).collect();
    for mask in 0u16..(1 << edges.len()) {
        let mut builder = chematic_core::MoleculeBuilder::new();
        for _ in 0..5 {
            builder.add_atom(chematic_core::Atom::new(Element::C));
        }
        for (i, &(a, b)) in edges.iter().enumerate() {
            if mask & (1 << i) != 0 {
                builder.add_bond(AtomIdx(a), AtomIdx(b), BondOrder::Single).unwrap();
            }
        }
        let mol = builder.build();
        compare(&mol);
    }
}
''')
