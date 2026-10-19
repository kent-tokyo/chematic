//! Differential oracle: the unconditionally prepared assignment pipeline.
use super::*;
const CORPORA: [&str; 2] = [
    include_str!("../../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"),
    include_str!("../../../validation/cosmolkit_comparison/corpora/stress-v1.smi"),
];

fn eager(mol: &Molecule, budget: CipBudget) -> Result<AccurateCipAssignment, CipCompareError> {
    let kekule = prepare_kekule_form(mol).ok();
    let pass1 = assign_all(mol, budget, kekule.as_ref())?;
    let pass2 = crate::resolver::apply_rule4b_pass(mol, budget, kekule.as_ref(), pass1);
    Ok(apply_rule5_pass(mol, budget, kekule.as_ref(), pass2))
}
fn compare(smiles: &str, budget: CipBudget) {
    let mol = if smiles.is_empty() {
        chematic_core::MoleculeBuilder::new().build()
    } else {
        chematic_smiles::parse(smiles).unwrap()
    };
    assert_eq!(
        format!("{:?}", assign_cip_accurate_experimental(&mol, budget)),
        format!("{:?}", eager(&mol, budget)),
        "{smiles} {budget:?}"
    );
}
#[test]
fn prepared_and_skipped_preparation_match_all_corpus_rows() {
    let mut count = 0;
    for text in CORPORA {
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            compare(
                line.split_whitespace().next().unwrap(),
                CipBudget::default_budget(),
            );
            count += 1;
        }
    }
    assert_eq!(count, 14_072);
}
#[test]
fn preparation_guard_preserves_tiny_budgets_and_noncarbon_centres() {
    for smiles in [
        "",
        "c1ccccc1",
        "F/C=C/F",
        "F/C(Cl)=C(/Br)I",
        "C[C@H](N)O",
        "N[P@](F)(Cl)Br",
        "C[S@](=O)CC",
        "[H][C@](F)(Cl)Br",
        "N[Pt@SP1](Cl)(Br)I",
        "[C@](F)(F)(Cl)Br",
    ] {
        for budget in [
            CipBudget::default_budget(),
            CipBudget {
                max_nodes: 0,
                max_depth: 0,
                max_expansions: 0,
            },
            CipBudget {
                max_nodes: 8,
                max_depth: 1,
                max_expansions: 2,
            },
        ] {
            compare(smiles, budget);
        }
    }
}
