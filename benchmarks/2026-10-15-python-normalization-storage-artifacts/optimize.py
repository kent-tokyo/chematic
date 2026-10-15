from pathlib import Path
s=Path('/private/tmp/chematic-descriptor-profile-source-20261015');r=Path('/private/tmp/chematic-descriptor-profile-20261015')
p=s/'crates/chematic-perception/src/aromaticity/model.rs';old=p.read_text();(r/'baseline-model.rs').write_text(old);(r/'baseline-CHANGELOG.md').write_bytes((s/'CHANGELOG.md').read_bytes())
start=old.index('fn finish_molecule_with_flags(');reference=old[start:].strip();assert reference.endswith('}')
new=old.replace('let pre_h: Vec<Option<u8>> = mol','let mut pre_h: Vec<Option<u8>> = mol')
new=new.replace('let mut builder = MoleculeBuilder::new();','let mut builder = MoleculeBuilder::with_capacity(mol.atom_count(), mol.bond_count());')
a=new.index('    let needs_patch: Vec<');b=new.index('    for (idx, atom) in normalized.atoms()',a)
new=new[:a]+'''    // Reuse the pre-normalization storage as an indexed patch map. Clearing
    // unchanged entries avoids a second allocation and per-atom linear scans.
    let mut needs_patch = false;
    for (index, h) in pre_h.iter_mut().enumerate() {
        *h = h.filter(|&pre| pre != implicit_hcount(&normalized, AtomIdx(index as u32)));
        needs_patch |= h.is_some();
    }
    if !needs_patch {
        return normalized;
    }

    let mut patched =
        MoleculeBuilder::with_capacity(normalized.atom_count(), normalized.bond_count());
'''+new[b:]
new=new.replace('if let Some(&(_, h)) = needs_patch.iter().find(|(pidx, _)| *pidx == idx) {','if let Some(h) = pre_h[idx.0 as usize] {')
new+='\n#[cfg(test)]\n#[path = "normalization_storage_tests.rs"]\nmod normalization_storage_tests;\n'
p.write_text(new)
test='''//! Storage-only normalization oracle frozen from the immediate baseline.
//! All chemical inputs come from the existing licensed and pinned 10k corpus.
use super::*;

'''+reference.replace('fn finish_molecule_with_flags(', 'fn reference_finish(')+'''

fn compare_graph(expected: &Molecule, actual: &Molecule) {
    assert_eq!(actual.atom_count(), expected.atom_count());
    assert_eq!(actual.bond_count(), expected.bond_count());
    assert_eq!(actual.stereo_groups(), expected.stereo_groups());
    assert_eq!(actual.rdkit_bond_order(), expected.rdkit_bond_order());
    for group in 0..expected.stereo_groups().len() {
        assert_eq!(actual.stereo_group_bond_indices(group), expected.stereo_group_bond_indices(group));
    }
    for (idx, atom) in expected.atoms() {
        assert_eq!(actual.atom(idx), atom);
        assert_eq!(actual.atom_tag(idx), expected.atom_tag(idx));
        assert_eq!(actual.r_group_label(idx), expected.r_group_label(idx));
        assert_eq!(actual.neighbors(idx).collect::<Vec<_>>(), expected.neighbors(idx).collect::<Vec<_>>());
        assert_eq!(actual.stereo_neighbor_order(idx), expected.stereo_neighbor_order(idx));
        assert_eq!(chematic_core::implicit_hcount(actual, idx), chematic_core::implicit_hcount(expected, idx));
    }
    for (idx, bond) in expected.bonds() {
        assert_eq!(actual.bond(idx), bond);
        assert_eq!(actual.bond_direction(idx), expected.bond_direction(idx));
        assert_eq!(actual.bond_direction_anchor(idx), expected.bond_direction_anchor(idx));
        assert_eq!(actual.smiles_ring_closure_metadata(idx), expected.smiles_ring_closure_metadata(idx));
    }
}

fn compare_normalization(input: &Molecule, model: &AromaticityModel) -> usize {
    let orders = compute_bond_orders(input, model);
    let flags = input.atoms().filter_map(|(idx, atom)| (model.is_atom_aromatic(idx) || atom.aromatic).then_some(idx)).collect();
    let expected = reference_finish(input, &flags, &orders);
    let actual = finish_molecule_with_flags(input, &flags, &orders);
    compare_graph(&expected, &actual);
    input.atoms().filter(|(idx, atom)| atom.hydrogen_count != actual.atom(*idx).hydrogen_count).count()
}

#[test]
fn fixed_corpus_normalization_matches_original_storage() {
    let corpus = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"));
    let mut rows = 0;
    let mut patches = 0;
    for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
        let input = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
        for algorithm in [AromaticityAlgorithm::Huckel, AromaticityAlgorithm::RdkitLike] {
            patches += compare_normalization(&input, &assign_aromaticity_ex(&input, algorithm));
            if let Ok(kekule) = chematic_core::kekulize(&input) {
                let kekulized = chematic_core::apply_kekule(&input, &kekule);
                patches += compare_normalization(&kekulized, &assign_aromaticity_ex(&kekulized, algorithm));
            }
        }
        rows += 1;
    }
    assert_eq!(rows, 10000);
    assert!(patches > 0, "exercise the indexed hydrogen patch map");
}
'''
(s/'crates/chematic-perception/src/aromaticity/normalization_storage_tests.rs').write_text(test)
ch=s/'CHANGELOG.md';ch.write_text(ch.read_text().replace('### Performance\n','### Performance\n\n- Reserve known atom/bond capacities during aromatic normalization and reuse\n  indexed hydrogen-count storage for corrections, preserving aromaticity,\n  hydrogen counts and the existing metadata-copy behavior.\n',1))
print('Storage-only normalization candidate written')
