from pathlib import Path
r=Path('/private/tmp/chematic-parity-lazy-20261015');s=Path('/private/tmp/chematic-parity-lazy-source-20261015')
assert (r/'baseline-install.log').exists(), 'Build and install unchanged baseline first'
p=s/'crates/chematic-perception/src/aromaticity/model.rs';old=(r/'baseline-model.rs').read_text();assert p.read_text()==old
x=old
x=x.replace('let atom_aromatic: FxHashSet<AtomIdx> = mol\n        .atoms()\n        .filter_map(|(idx, atom)| (model.is_atom_aromatic(idx) || atom.aromatic).then_some(idx))','let atom_aromatic: Vec<bool> = mol\n        .atoms()\n        .map(|(idx, atom)| model.is_atom_aromatic(idx) || atom.aromatic)')
x=x.replace('let mut atom_aromatic: FxHashSet<AtomIdx> = mol\n        .atoms()\n        .filter_map(|(idx, _)| model.is_atom_aromatic(idx).then_some(idx))','let mut atom_aromatic: Vec<bool> = mol\n        .atoms()\n        .map(|(idx, _)| model.is_atom_aromatic(idx))')
x=x.replace('bond_orders[&bidx]', 'bond_orders[bidx.0 as usize]')
x=x.replace('atom_aromatic.insert(bond.atom1);\n            atom_aromatic.insert(bond.atom2);','atom_aromatic[bond.atom1.0 as usize] = true;\n            atom_aromatic[bond.atom2.0 as usize] = true;')
x=x.replace('fn compute_bond_orders(mol: &Molecule, model: &AromaticityModel) -> FxHashMap<BondIdx, BondOrder>', 'fn compute_bond_orders(mol: &Molecule, model: &AromaticityModel) -> Vec<BondOrder>')
x=x.replace('            let order = if model.is_bond_aromatic(bidx) {','            if model.is_bond_aromatic(bidx) {').replace('                bond.order\n            };\n            (bidx, order)','                bond.order\n            }')
x=x.replace('atom_aromatic: &FxHashSet<AtomIdx>,\n    bond_orders: &FxHashMap<BondIdx, BondOrder>,','atom_aromatic: &[bool],\n    bond_orders: &[BondOrder],')
x=x.replace('a.aromatic = atom_aromatic.contains(&idx);','a.aromatic = atom_aromatic[idx.0 as usize];')
x=x.replace('fn compute_bond_orders(mol:', '// Molecule atom/bond iterators enumerate dense indices in order. Keep the\n// normalization facts aligned with those indices instead of hashing them twice.\nfn compute_bond_orders(mol:')
x+='\n#[cfg(test)]\n#[path = "indexed_normalization_tests.rs"]\nmod indexed_normalization_tests;\n';p.write_text(x)
# Preserve the earlier independent storage oracle; only adapt its call to the new private representation.
p=s/'crates/chematic-perception/src/aromaticity/normalization_storage_tests.rs';x=(r/'baseline-normalization-storage-tests.rs').read_text();assert p.read_text()==x
x=x.replace('use super::*;', 'use super::*;\nuse rustc_hash::FxHashMap;')
x=x.replace('let orders = compute_bond_orders(input, model);', 'let orders: FxHashMap<BondIdx, BondOrder> = input\n        .bonds()\n        .map(|(idx, bond)| (idx, if model.is_bond_aromatic(idx) { BondOrder::Aromatic } else { bond.order }))\n        .collect();')
x=x.replace('let actual = finish_molecule_with_flags(input, &flags, &orders);', 'let actual = build_molecule_from_model(input, model);');p.write_text(x)
# Freeze complete baseline builders and finish, independently of the new indexed implementation.
body=old[old.index('pub(crate) fn build_molecule_from_model('):old.index('\n#[cfg(test)]',old.index('pub(crate) fn build_molecule_from_model('))]
prev=(r/'baseline-normalization-storage-tests.rs').read_text();compare=prev[prev.index('fn compare_graph('):prev.index('\nfn compare_normalization')]
compare=compare.replace('        assert_eq!(actual.atom_tag(idx), expected.atom_tag(idx));','        assert_eq!(actual.atom_tag(idx), expected.atom_tag(idx));\n        assert_eq!(actual.smiles_ring_closure_count(idx), expected.smiles_ring_closure_count(idx));')
t='//! Complete immediate-baseline normalization oracle; existing pinned corpora only.\nuse super::*;\nuse rustc_hash::FxHashMap;\n\nmod frozen {\nuse super::*;\n'+body+'\n}\n\n'+compare+'''
fn compare_both(input: &Molecule, model: &AromaticityModel) {
    compare_graph(&frozen::build_molecule_from_model(input, model), &build_molecule_from_model(input, model));
    compare_graph(&frozen::build_molecule_from_model_authoritative(input, model), &build_molecule_from_model_authoritative(input, model));
}

#[test]
fn indexed_normalization_matches_both_original_builders_on_pinned_corpora() {
    let corpora = [
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")),
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi")),
    ];
    let mut rows = 0;
    let mut comparisons = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let input = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            let kekulized = chematic_core::kekulize(&input)
                .map(|k| chematic_core::apply_kekule(&input, &k))
                .unwrap_or_else(|_| input.clone());
            for view in [&input, &kekulized] {
                for algorithm in [AromaticityAlgorithm::Huckel, AromaticityAlgorithm::RdkitLike] {
                    compare_both(view, &assign_aromaticity_ex(view, algorithm));
                    comparisons += 2;
                }
            }
            rows += 1;
        }
    }
    assert_eq!(rows, 14072);
    assert_eq!(comparisons, 112576);
}

#[test]
fn indexed_normalization_preserves_metadata_and_authoritative_fallbacks() {
    let mut inputs = vec![chematic_core::MoleculeBuilder::new().build()];
    for smiles in ["C1=CC=CC=C1", "c1ccccc1", "C1=Cc2ccccc2C1", "O=C1NNC(=O)N1", "[C@H](F)(Cl)C1=CC=CC=C1", "C/C=C/C1=CC=CC=C1", "[*:1]C1=CC=CC=C1", "C1=C=CC=CC1", "[O-][Cl](=O)(=O)=O", "CO(C)[Cu]"] {
        let mut input = chematic_smiles::parse_template(smiles).unwrap();
        input.set_tag(AtomIdx(0), Some(23));
        input.add_stereo_group_with_bonds(
            chematic_core::StereoGroup::new(chematic_core::StereoGroupKind::Or(9), vec![AtomIdx(0)]),
            vec![BondIdx(0)],
        );
        if input.atom(AtomIdx(0)).wildcard { input.set_r_group(AtomIdx(0), chematic_core::RGroupLabel::unnumbered()); }
        inputs.push(input.clone());
        input.set_atom_aromatic(AtomIdx(0), true);
        inputs.push(input.clone());
        let first = input.bond(BondIdx(0)).clone();
        input.remove_bond(BondIdx(0));
        input.add_bond(first.atom1, first.atom2, first.order).unwrap();
        inputs.push(input);
    }
    for input in inputs {
        for (atoms, bonds) in [
            (FxHashSet::default(), FxHashSet::default()),
            (FxHashSet::default(), input.bonds().map(|(idx, _)| idx).collect()),
            (input.atoms().map(|(idx, _)| idx).collect(), FxHashSet::default()),
            (input.atoms().filter_map(|(idx, _)| (idx.0 % 2 == 0).then_some(idx)).collect(), input.bonds().filter_map(|(idx, _)| (idx.0 % 2 == 0).then_some(idx)).collect()),
            ([AtomIdx(u32::MAX)].into_iter().collect(), [BondIdx(u32::MAX)].into_iter().collect()),
        ] {
            compare_both(&input, &AromaticityModel::from_atom_bond_sets(atoms, bonds));
        }
    }
}
'''
(s/'crates/chematic-perception/src/aromaticity/indexed_normalization_tests.rs').write_text(t)
p=s/'CHANGELOG.md';x=p.read_text();x=x.replace('### Performance\n','### Performance\n\n- Use dense atom flags and bond orders during aromatic normalization to avoid\n  temporary hash tables, preserving hydrogen counts and stereo metadata.\n',1);p.write_text(x)
print('Indexed aromatic normalization prototype written; full frozen baseline oracle retained')
