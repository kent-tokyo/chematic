//! Immediate-baseline TPSA oracle on existing pinned input corpora.
use super::*;

fn frozen_default_tpsa(mol: &Molecule) -> f64 {
    let mol_arom = descriptor_aromaticity(mol);
    // RDKit's default excludes S/P. Compute only N/O instead of allocating a
    // full contribution vector and typing atoms whose values are discarded.
    // Iteration and accumulation order remain atom-index order, exactly as in
    // `tpsa_contributions`, so floating-point results are unchanged.
    mol.atoms().fold(0.0, |acc, (idx, atom)| {
        let contribution = match atom.element.atomic_number() {
            7 => {
                let hydrogens = descriptor_attached_hcount(mol, idx);
                let environment_mol = if use_perceived_nitrogen_environment(mol, &mol_arom, idx) {
                    &mol_arom
                } else {
                    mol
                };
                tpsa_nitrogen(tpsa_environment(environment_mol, idx, hydrogens))
            }
            8 => {
                let hydrogens = descriptor_attached_hcount(&mol_arom, idx);
                tpsa_oxygen(tpsa_environment(&mol_arom, idx, hydrogens))
            }
            _ => 0.0,
        };
        acc + contribution
    })
}

fn compare(input: &Molecule) {
    let original = input.clone();
    let candidate = input.clone();
    for _ in 0..2 {
        assert_eq!(
            rdkit_tpsa_impl(&candidate).to_bits(),
            frozen_default_tpsa(&original).to_bits()
        );
    }
}

#[test]
fn zero_guard_keeps_exact_cold_and_repeated_tpsa_on_all_pinned_inputs() {
    let corpora = [
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
        )),
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi"
        )),
    ];
    let mut rows = 0;
    let mut zero_inputs = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let input =
                chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            zero_inputs += usize::from(
                !input
                    .atoms()
                    .any(|(_, atom)| matches!(atom.element.atomic_number(), 7 | 8)),
            );
            compare(&input);
            let kekulized = chematic_core::kekulize(&input)
                .map(|k| chematic_core::apply_kekule(&input, &k))
                .unwrap_or_else(|_| input.clone());
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
    for smiles in [
        "C1=CC=CC=C1",
        "c1ccccc1",
        "C1=C=CC=CC1",
        "[*:1]C1=CC=CC=C1",
        "[H]",
        "SP",
        "[Cu]",
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCN",
        "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCO",
    ] {
        let input = chematic_smiles::parse_template(smiles).unwrap();
        compare(&input);
    }
}
