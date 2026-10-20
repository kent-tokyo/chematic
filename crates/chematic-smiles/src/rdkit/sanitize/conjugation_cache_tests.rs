//! Frozen-loop differential and opt-in timing checks for RDKit #9475.

use super::*;
use crate::rdkit::parse::from_chematic;
use std::hint::black_box;
use std::time::Instant;

// The pre-optimization loop, deliberately independent of the new cache.
fn reference(mol: &mut Mol) {
    for b in &mut mol.bonds {
        b.conjugated = b.aromatic;
    }
    for a in 0..mol.atoms.len() {
        if !is_atom_conjug_cand(mol, a) {
            continue;
        }
        let sbo = mol.degree(a) + mol.total_num_hs(a) as usize;
        if !(2..=3).contains(&sbo) {
            continue;
        }
        let degree = mol.atom_bonds[a].len();
        for i1 in 0..degree {
            let b1 = mol.atom_bonds[a][i1];
            if mol.bonds[b1].valence_contrib(a) < 1.5
                || !is_atom_conjug_cand(mol, mol.bonds[b1].other(a))
            {
                continue;
            }
            for i2 in 0..degree {
                let b2 = mol.atom_bonds[a][i2];
                if b1 == b2 {
                    continue;
                }
                let at2 = mol.bonds[b2].other(a);
                let sbo2 = mol.degree(at2) + mol.total_num_hs(at2) as usize;
                if sbo2 > 3 {
                    continue;
                }
                if is_atom_conjug_cand(mol, at2) {
                    mol.bonds[b1].conjugated = true;
                    mol.bonds[b2].conjugated = true;
                }
            }
        }
    }
}

fn prepare(smiles: &str) -> Result<Mol, String> {
    let source = crate::parse(smiles).map_err(|e| e.to_string())?;
    let mut mol = from_chematic(&source).map_err(|e| e.to_string())?;
    remove_hs_and_sanitize(&mut mol).map_err(|e| e.to_string())?;
    Ok(mol)
}

type PreparedModels = Vec<(String, Mol)>;
type RefusedModels = Vec<(String, String)>;

fn corpus() -> (PreparedModels, RefusedModels) {
    let mut accepted = Vec::new();
    let mut refused = Vec::new();
    for line in include_str!("../../../tests/data/rdkit_canonical_smiles.tsv").lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let smiles = line.split('\t').next().unwrap();
        match prepare(smiles) {
            Ok(mol) => accepted.push((smiles.into(), mol)),
            Err(error) => refused.push((smiles.into(), error)),
        }
    }
    (accepted, refused)
}

fn compare(smiles: &str, mol: &Mol) {
    for stale in [false, true] {
        let mut old = mol.clone();
        for bond in &mut old.bonds {
            bond.conjugated = stale;
        }
        let mut new = old.clone();
        reference(&mut old);
        set_conjugation(&mut new);
        // The whole internal model must agree, not only canonical output.
        assert_eq!(format!("{new:?}"), format!("{old:?}"), "{smiles}");
        let once = format!("{new:?}");
        set_conjugation(&mut new);
        assert_eq!(format!("{new:?}"), once, "not idempotent: {smiles}");
    }
}

#[test]
fn frozen_loop_matches_entire_fixture_and_records_refusals() {
    let (models, refused) = corpus();
    assert_eq!(models.len() + refused.len(), 1131);
    assert!(models.len() > 1100);
    for (smiles, mol) in &models {
        compare(smiles, mol);
    }
    eprintln!(
        "conjugation differential: {} prepared, {} refused",
        models.len(),
        refused.len()
    );
    for (smiles, error) in refused {
        eprintln!("refused {smiles}: {error}");
    }
}

#[test]
fn predicate_is_independent_of_assigned_conjugation_flags() {
    for smiles in [
        "CC(=O)N",
        "c1cc[nH]c1",
        "C=C=C",
        "[O-][N+](=O)c1ccccc1",
        "N->[Cu+2]<-N",
    ] {
        let mut mol = prepare(smiles).unwrap();
        for bond in &mut mol.bonds {
            bond.conjugated = false;
        }
        let before: Vec<_> = (0..mol.atoms.len())
            .map(|a| is_atom_conjug_cand(&mol, a))
            .collect();
        for bond in &mut mol.bonds {
            bond.conjugated = true;
        }
        let after: Vec<_> = (0..mol.atoms.len())
            .map(|a| is_atom_conjug_cand(&mol, a))
            .collect();
        assert_eq!(before, after, "{smiles}");
    }
}

#[test]
fn cache_spill_and_high_degree_radical_charge_controls() {
    for smiles in [
        "C".repeat(64),
        "C".repeat(65),
        "C".repeat(256),
        "C=C".repeat(65),
        "[CH2]C=C".into(),
        "[O-]C(=O)C=C".into(),
        "F[S](F)(F)(F)(F)F".into(),
        "[Fe](N)(N)(N)(N)(N)N".into(),
        "C[Si](C)(C)C".into(),
        "[H]C([H])([H])C=C".into(),
    ] {
        let mol = prepare(&smiles).unwrap();
        compare(&smiles, &mol);
    }
    compare("empty", &Mol::default());
}

#[test]
fn single_bond_fast_path_keeps_existing_aromatic_flags() {
    let mut mol = prepare("CCO").unwrap();
    // Aromatic flags and the stored bond type are independent model fields.
    // Even on such a model the fast path must retain the original reset rule.
    mol.bonds[0].aromatic = true;
    compare("single types with an aromatic flag", &mol);
    set_conjugation(&mut mol);
    assert!(mol.bonds[0].conjugated);
    assert!(!mol.bonds[1].conjugated);
}

#[test]
#[ignore = "opt-in release timing; measures only prepared-model set_conjugation"]
fn conjugation_cache_timing() {
    let (models, refused) = corpus();
    eprintln!(
        "timing corpus: {} prepared, {} refused",
        models.len(),
        refused.len()
    );
    for group in ["fixture", "single_only", "multiple"] {
        let mut old: Vec<_> = models
            .iter()
            .filter(|(_, mol)| {
                let single = mol
                    .bonds
                    .iter()
                    .all(|b| matches!(b.bt, BondType::Single | BondType::Dative));
                group == "fixture" || (group == "single_only") == single
            })
            .map(|(_, mol)| mol.clone())
            .collect();
        eprintln!("conjugation_group\t{group}\t{}", old.len());
        let mut new = old.clone();
        for mol in &mut old {
            reference(black_box(mol));
        }
        for mol in &mut new {
            set_conjugation(black_box(mol));
        }
        for round in 0..9 {
            let mut elapsed = [0.0; 2];
            for lane in if round % 2 == 0 { [0, 1] } else { [1, 0] } {
                let start = Instant::now();
                for _ in 0..300 {
                    if lane == 0 {
                        for mol in &mut old {
                            reference(black_box(mol));
                        }
                    } else {
                        for mol in &mut new {
                            set_conjugation(black_box(mol));
                        }
                    }
                }
                elapsed[lane] = start.elapsed().as_secs_f64();
            }
            eprintln!(
                "conjugation_timing\t{group}\t{round}\t{:.9}\t{:.9}",
                elapsed[0], elapsed[1]
            );
        }
        for (old, new) in old.iter().zip(&new) {
            assert_eq!(format!("{new:?}"), format!("{old:?}"));
        }
    }
}
