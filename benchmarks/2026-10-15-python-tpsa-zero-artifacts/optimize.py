from pathlib import Path
r=Path('/private/tmp/chematic-tpsa-zero-20261015');s=Path('/private/tmp/chematic-tpsa-zero-source-20261015');p=s/'crates/chematic-chem/src/descriptors/tpsa.rs';old=(r/'baseline-tpsa.rs').read_text();assert p.read_text()==old
x=old.replace('fn rdkit_tpsa_impl(mol: &Molecule) -> f64 {\n    let mol_arom', '''fn rdkit_tpsa_impl(mol: &Molecule) -> f64 {
    // This default model sums only N/O contributions. With neither element,
    // the original fold adds only positive zeros, without needing perception.
    if !mol.atoms().any(|(_, atom)| matches!(atom.element.atomic_number(), 7 | 8)) {
        return 0.0;
    }
    let mol_arom''',1)
x+='\n#[cfg(test)]\n#[path = "tpsa_zero_tests.rs"]\nmod tpsa_zero_tests;\n';p.write_text(x)
start=old.index('fn rdkit_tpsa_impl(');end=old.index('\n#[cfg(test)]',start);frozen=old[start:end].replace('fn rdkit_tpsa_impl(', 'fn frozen_default_tpsa(')
t='//! Immediate-baseline TPSA oracle on existing pinned input corpora.\nuse super::*;\n\n'+frozen+'''
fn compare(input: &Molecule) {
    let original = input.clone();
    let candidate = input.clone();
    for _ in 0..2 {
        assert_eq!(rdkit_tpsa_impl(&candidate).to_bits(), frozen_default_tpsa(&original).to_bits());
    }
}

#[test]
fn zero_guard_keeps_exact_cold_and_repeated_tpsa_on_all_pinned_inputs() {
    let corpora = [
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi")),
        include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi")),
    ];
    let mut rows = 0;
    let mut zero_inputs = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let input = chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            zero_inputs += usize::from(!input.atoms().any(|(_, atom)| matches!(atom.element.atomic_number(), 7 | 8)));
            compare(&input);
            let kekulized = chematic_core::kekulize(&input).map(|k| chematic_core::apply_kekule(&input, &k)).unwrap_or_else(|_| input.clone());
            compare(&kekulized);
            rows += 1;
        }
    }
    assert_eq!(rows, 14072);
    assert!(zero_inputs > 0);
}

#[test]
fn default_zero_and_late_nitrogen_keep_original_float_bits() {
    compare(&chematic_core::MoleculeBuilder::new().build());
    for smiles in ["C1=CC=CC=C1", "c1ccccc1", "C1=C=CC=CC1", "[*:1]C1=CC=CC=C1", "[H]", "SP", "[Cu]", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCN", "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCO"] {
        let input = chematic_smiles::parse_template(smiles).unwrap();
        compare(&input);
    }
}
''';(s/'crates/chematic-chem/src/descriptors/tpsa_zero_tests.rs').write_text(t)
p=s/'CHANGELOG.md';x=p.read_text().replace('### Performance\n','### Performance\n\n- Skip aromatic perception for RDKit-default TPSA when its selected model view\n  contains no nitrogen or oxygen, retaining exact zero and all contributing-atom\n  calculations without changing parsing or S/P-inclusive TPSA.\n',1);p.write_text(x)
print('Default TPSA zero-contribution guard written; original N/O fold frozen for exact comparison')
