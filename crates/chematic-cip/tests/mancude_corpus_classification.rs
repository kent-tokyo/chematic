//! Milestone 3B-0: classify all 98 aromatic/MANCUDE-adjacent corpus cases (96
//! `aromatic_mancude` bucket + 2 `uncharacterized` cases M2.5 found mis-tagged, per
//! `tests/uncharacterized_diagnosis.rs`'s `BucketMisclassified` diagnosis) into
//! multi-label structural tags -- evidence for Milestone 3B-1's design, not a behavior
//! change. Every case gets at least one tag; the test fails loudly if any case falls
//! through unclassified.

mod common;

use std::collections::BTreeMap;

use chematic_cip::CipBudget;
use chematic_cip::digraph_diff::first_divergence;
use chematic_cip::mancude::{MancudeBudget, enumerate_kekule_matchings};
use chematic_core::kekulization::{apply_kekule, kekulize};
use chematic_core::{AtomIdx, BondOrder, Element, implicit_hcount};
use chematic_perception::{RingSystemKind, find_ring_families, find_sssr};

const CORPUS: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../validation/cip_label_corpus.jsonl"
));

#[derive(Default)]
struct CaseReport {
    tags: Vec<&'static str>,
}

fn classify(smiles: &str, atom_idx: u32) -> CaseReport {
    let mut report = CaseReport::default();
    let mol = chematic_smiles::parse(smiles)
        .unwrap_or_else(|e| panic!("corpus SMILES failed to parse: {smiles}: {e:?}"));

    let has_hetero_aromatic = mol
        .atoms()
        .any(|(_, a)| a.aromatic && a.element != Element::C);
    let has_aromatic = mol.atoms().any(|(_, a)| a.aromatic);
    if has_hetero_aromatic {
        report.tags.push("hetero_mancude");
    } else if has_aromatic {
        report.tags.push("hydrocarbon_mancude");
    }

    let sssr = find_sssr(&mol);
    let families = find_ring_families(&mol, &sssr);
    if families.iter().any(|f| f.kind == RingSystemKind::Fused) {
        report.tags.push("fused_ring");
    }
    if families.iter().any(|f| f.kind == RingSystemKind::Bridged) {
        report.tags.push("bridged_ring");
    }

    let has_fully_substituted_ipso = mol.atoms().any(|(idx, a)| {
        a.aromatic && mol.neighbors(idx).count() == 3 && implicit_hcount(&mol, idx) == 0
    });
    if has_fully_substituted_ipso {
        report.tags.push("fully_substituted_ipso");
    }

    let has_exocyclic_multiple_bond = mol.atoms().any(|(idx, a)| {
        a.aromatic
            && mol.neighbors(idx).any(|(_, bidx)| {
                matches!(mol.bond(bidx).order, BondOrder::Double | BondOrder::Triple)
            })
    });
    if has_exocyclic_multiple_bond {
        report.tags.push("exocyclic_multiple_bond_adjacent");
    }

    if let Ok(kekule) = kekulize(&mol) {
        report.tags.push("kekulization_succeeds");
        let kekule_mol = apply_kekule(&mol, &kekule);
        if let Ok(Some(_)) = first_divergence(
            &mol,
            AtomIdx(atom_idx),
            &kekule_mol,
            AtomIdx(atom_idx),
            CipBudget::default_budget(),
        ) {
            report.tags.push("aromatic_vs_kekule_digraph_diverges");
        }
    }

    match enumerate_kekule_matchings(&mol, MancudeBudget::default()) {
        Ok(all) if all.len() > 1 => report.tags.push("multiple_kekule_forms"),
        Ok(_) => {}
        Err(_) => report.tags.push("enumeration_over_budget"),
    }

    report
}

#[test]
fn all_98_mancude_cases_are_classified() {
    let mut tag_tally: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mut unclassified: Vec<String> = Vec::new();
    let mut checked = 0usize;

    for line in CORPUS.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let value: serde_json::Value = serde_json::from_str(line).unwrap();
        let Some(smiles) = value.get("smiles").and_then(|v| v.as_str()) else {
            continue; // manifest line
        };
        let atom_idx = value.get("atom_idx").and_then(|v| v.as_u64()).unwrap() as u32;
        let bucket = value.get("bucket").and_then(|v| v.as_str());
        let modern = value.get("modern").and_then(|v| v.as_str()).unwrap();

        let mol = chematic_smiles::parse(smiles)
            .unwrap_or_else(|e| panic!("corpus SMILES failed to parse: {smiles}: {e:?}"));

        if !common::in_mancude_scope(&mol, atom_idx, bucket, modern) {
            continue;
        }

        checked += 1;
        let report = classify(smiles, atom_idx);
        if report.tags.is_empty() {
            unclassified.push(format!("{smiles} atom {atom_idx}: no tags matched"));
        }
        for &tag in &report.tags {
            *tag_tally.entry(tag).or_default() += 1;
        }
    }

    println!("\n=== mancude corpus classification (Milestone 3B-0) ===");
    for (tag, count) in &tag_tally {
        println!("  {tag:36} {count:3}");
    }
    println!("  {:36} {:3}", "checked", checked);

    assert_eq!(
        checked, 98,
        "expected 96 aromatic_mancude + 2 BucketMisclassified uncharacterized cases"
    );
    assert!(
        unclassified.is_empty(),
        "Milestone 3B-0 gate: every mancude-adjacent case must carry at least one \
         structural tag. {} case(s) did not:\n{}",
        unclassified.len(),
        unclassified.join("\n")
    );
}
