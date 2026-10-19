//! Frozen pre-change cleanup and original scheduling controls for RDKit #9475.

use super::*;
use crate::rdkit::parse::from_chematic;
use std::{hint::black_box, time::Instant};

fn reference(mol: &mut Mol) -> Result<(), RdkitSmilesError> {
    // nitrogensCleanup
    let mut to_consider = vec![false; mol.atoms.len()];
    for a in 0..mol.atoms.len() {
        if mol.atoms[a].anum != 7 || mol.atoms[a].charge != 0 {
            continue;
        }
        if mol.calc_explicit_valence(a, false)? != 5 {
            continue;
        }
        to_consider[a] = true;
        let arom = mol.atoms[a].aromatic;
        mol.atoms[a].aromatic = false;
        let mut update = false;
        for b in mol.atom_bonds[a].clone() {
            let nbr = mol.bonds[b].other(a);
            if mol.atoms[nbr].anum == 8
                && mol.atoms[nbr].charge == 0
                && mol.bonds[b].bt == BondType::Double
            {
                mol.bonds[b].bt = BondType::Single;
                mol.atoms[a].charge = 1;
                mol.atoms[nbr].charge = -1;
                update = true;
                break;
            }
        }
        mol.atoms[a].aromatic = arom;
        if update {
            mol.calc_explicit_valence(a, false)?;
        }
    }
    for a in 0..mol.atoms.len() {
        if !to_consider[a] {
            continue;
        }
        let arom = mol.atoms[a].aromatic;
        mol.atoms[a].aromatic = false;
        let mut update = false;
        for b in mol.atom_bonds[a].clone() {
            let nbr = mol.bonds[b].other(a);
            if mol.atoms[nbr].anum == 7
                && mol.atoms[nbr].charge == 0
                && mol.bonds[b].bt == BondType::Triple
            {
                mol.bonds[b].bt = BondType::Double;
                mol.atoms[a].charge = 1;
                mol.atoms[nbr].charge = -1;
                update = true;
                break;
            }
        }
        mol.atoms[a].aromatic = arom;
        if update {
            mol.calc_explicit_valence(a, false)?;
        }
    }
    for a in 0..mol.atoms.len() {
        match mol.atoms[a].anum {
            15 => phosphorus_cleanup(mol, a)?,
            17 | 35 | 53 => halogen_cleanup(mol, a)?,
            _ => {}
        }
    }
    Ok(())
}

type Models = Vec<(String, Mol)>;

fn prepare(s: &str) -> Result<Mol, String> {
    let source = crate::parse(s).map_err(|e| e.to_string())?;
    from_chematic(&source).map_err(|e| e.to_string())
}

fn corpus() -> (Models, Vec<(String, String)>) {
    let mut models = Vec::new();
    let mut refused = Vec::new();
    for line in include_str!("../../../tests/data/rdkit_canonical_smiles.tsv").lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let s = line.split('\t').next().unwrap();
        match prepare(s) {
            Ok(m) => models.push((s.into(), m)),
            Err(e) => refused.push((s.into(), e)),
        }
    }
    (models, refused)
}

fn compare(s: &str, input: &Mol) {
    let mut old = input.clone();
    let mut new = input.clone();
    for pass in 0..2 {
        assert_eq!(reference(&mut old), clean_up(&mut new), "{s}, pass {pass}");
        assert_eq!(format!("{old:?}"), format!("{new:?}"), "{s}, pass {pass}");
    }
}

#[test]
fn frozen_cleanup_matches_all_fixture_models_and_retains_refusals() {
    let (models, refused) = corpus();
    assert_eq!(models.len() + refused.len(), 1131);
    assert_eq!(models.len(), 1130);
    assert_eq!(refused.len(), 1);
    for (s, m) in &models {
        compare(s, m);
    }
    eprintln!(
        "cleanup differential: {} prepared, {} refused",
        models.len(),
        refused.len()
    );
    for (s, e) in refused {
        eprintln!("refused {s}: {e}");
    }
}

const CONTROLS: &[&str] = &[
    "CN(=O)=O",
    "N(=O)#N",
    "N#N=N",
    "O=N(=O)N(=O)=O",
    "C[N+](=O)[O-]",
    "C[N-]N",
    "N",
    "[NH4+]",
    "n1ccccc1",
    "[nH]1cccc1",
    "CP(=O)=C",
    "NP(=O)=N",
    "O=Cl(=O)O",
    "O=Br(=O)O",
    "O=I(=O)O",
    "[13CH3]N(=O)=O",
    "F[C@H](Cl)CN(=O)=O",
    "[H]N(=O)=O",
    "[2H]N(=O)=O",
    "[Cu]N(=O)=O",
    "N->[Cu+2]",
    "*N(=O)=O",
];

#[test]
fn charged_aromatic_stereo_isotope_and_neighbour_order_controls() {
    for &s in CONTROLS {
        let m = prepare(s).unwrap_or_else(|e| panic!("{s}: {e}"));
        compare(s, &m);
        // Repeated calls must not reuse discovery state from another molecule.
        let mut changed = m.clone();
        for atom in &mut changed.atoms {
            atom.aromatic = !atom.aromatic;
            atom.explicit_valence = -1;
        }
        compare(s, &changed);
    }
    let mut m = prepare("N(=O)#N").unwrap();
    clean_up(&mut m).unwrap();
    // The central N stays a second-pass candidate after its first-pass charge
    // change; filtering on its *new* charge would skip this N#N rewrite.
    assert_eq!(m.bonds[1].bt, BondType::Double);
    assert_eq!(m.atoms[1].charge, -1);
    assert_eq!(m.atoms[2].charge, -1);
}

#[test]
fn inline_spill_and_large_sparse_molecules_preserve_order() {
    compare("empty", &Mol::default());
    for n in [0, 1, 7, 8, 9, 64] {
        let s = if n == 0 {
            "CC".into()
        } else {
            vec!["N(=O)#N"; n].join(".")
        };
        compare(&s, &prepare(&s).unwrap());
    }
    for n in [64, 256, 2048] {
        let s = format!("{}.N(=O)#N", "C".repeat(n));
        compare(&s, &prepare(&s).unwrap());
    }
}

#[inline(never)]
fn measure(models: &[Mol], f: fn(&mut Mol) -> Result<(), RdkitSmilesError>) -> f64 {
    let mut elapsed = std::time::Duration::ZERO;
    for _ in 0..300 {
        // Parsing and mutable-state restoration are deliberately outside the
        // timer: both functions see identical fresh raw internal models.
        let mut work = models.to_vec();
        let start = Instant::now();
        for m in &mut work {
            black_box(f(black_box(m))).unwrap();
        }
        elapsed += start.elapsed();
        black_box(&work);
    }
    elapsed.as_secs_f64()
}

#[test]
#[ignore = "opt-in prepared-cleanup timing; not an end-to-end benchmark"]
fn cleanup_candidates_timing() {
    let (models, _) = corpus();
    let fixture: Vec<_> = models.into_iter().map(|(_, m)| m).collect();
    let controls: Vec<_> = CONTROLS.iter().map(|s| prepare(s).unwrap()).collect();
    for (group, models) in [("fixture", fixture), ("controls", controls)] {
        eprintln!("cleanup_group\t{group}\t{}", models.len());
        for m in &models {
            compare(group, m);
        }
        // Warm both code paths before alternating their measurement order.
        black_box(measure(&models, reference));
        black_box(measure(&models, clean_up));
        for round in 0..9 {
            let (old, new) = if round % 2 == 0 {
                (measure(&models, reference), measure(&models, clean_up))
            } else {
                let new = measure(&models, clean_up);
                (measure(&models, reference), new)
            };
            eprintln!("cleanup_timing\t{group}\t{round}\t{old:.9}\t{new:.9}");
        }
    }
}
