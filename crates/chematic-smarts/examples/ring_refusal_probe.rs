//! Diagnose the six pinned `[RN]` refusal structures without changing matching.
//! Usage: cargo run -p chematic-smarts --example ring_refusal_probe -- CORPUS

use chematic_core::{AtomIdx, Molecule};
use chematic_perception::{RingSet, find_sssr, find_symmetrized_sssr_with_diagnostics_bounded};
use chematic_smarts::perceived::with_perceived_target;
use chematic_smarts::rdkit_ring_model::{
    RdkitRingModelBudget, build_rdkit_parity_ring_model, build_shared_symmetrized_ring_model,
};
use serde_json::{Value, json};

const ROWS: [usize; 6] = [9, 23, 28, 29, 30, 34];

fn atom_sets(rings: &RingSet) -> Vec<Vec<u32>> {
    let mut sets: Vec<Vec<u32>> = rings
        .rings()
        .iter()
        .map(|ring| {
            let mut atoms: Vec<u32> = ring.iter().map(|atom| atom.0).collect();
            atoms.sort_unstable();
            atoms
        })
        .collect();
    sets.sort();
    sets
}

fn bond_sets(mol: &Molecule, rings: &RingSet) -> Vec<Vec<[u32; 2]>> {
    let mut sets: Vec<Vec<[u32; 2]>> = rings
        .rings()
        .iter()
        .map(|ring| {
            let mut bonds: Vec<[u32; 2]> = ring
                .iter()
                .zip(ring.iter().cycle().skip(1))
                .take(ring.len())
                .map(|(a, b)| {
                    assert!(mol.bond_between(*a, *b).is_some(), "ring edge must exist");
                    [a.0.min(b.0), a.0.max(b.0)]
                })
                .collect();
            bonds.sort_unstable();
            bonds
        })
        .collect();
    sets.sort();
    sets
}

fn probe(mol: &Molecule) -> Result<Value, Box<dyn std::error::Error>> {
    let base = find_sssr(mol);
    let budget = RdkitRingModelBudget::default();
    let sym = find_symmetrized_sssr_with_diagnostics_bounded(mol, Some(budget.max_candidates));
    let approximate = build_rdkit_parity_ring_model(mol, &base, &budget)?;
    let shared = build_shared_symmetrized_ring_model(mol, &base, &budget)?;
    let counts = |model: &chematic_smarts::rdkit_ring_model::RdkitParityRingModel| {
        (0..mol.atom_count())
            .map(|i| model.ring_count(AtomIdx(i as u32)))
            .collect::<Vec<_>>()
    };
    Ok(json!({
        "atom_count": mol.atom_count(),
        "elements": mol.atoms().map(|(_, atom)| atom.element.symbol()).collect::<Vec<_>>(),
        "charges": mol.atoms().map(|(_, atom)| atom.charge).collect::<Vec<_>>(),
        "aromatic_cations": mol.atoms().filter(|(_, atom)| atom.aromatic && atom.charge > 0).count(),
        "base_ring_atom_sets": atom_sets(&base),
        "base_ring_bond_sets": bond_sets(mol, &base),
        "symmetrized_ring_atom_sets": atom_sets(sym.rings()),
        "symmetrized_ring_bond_sets": bond_sets(mol, sym.rings()),
        "symmetrized_status": format!("{:?}", sym.status()),
        "approximate_extra_ring_count": approximate.extra_ring_count(),
        "approximate_ring_counts": counts(&approximate),
        "shared_extra_ring_count": shared.extra_ring_count(),
        "shared_ring_counts": counts(&shared),
    }))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("usage: ring_refusal_probe CORPUS")?;
    let corpus = std::fs::read_to_string(path)?;
    let lines: Vec<_> = corpus.lines().collect();
    for row in ROWS {
        let smiles = *lines
            .get(row)
            .ok_or("corpus is shorter than a selected row")?;
        let mol = chematic_smiles::parse(smiles)?;
        let diagnostics = with_perceived_target(&mol, probe)?;
        println!(
            "{}",
            json!({"input_index":row,"smiles":smiles,"diagnostics":diagnostics})
        );
    }
    Ok(())
}
