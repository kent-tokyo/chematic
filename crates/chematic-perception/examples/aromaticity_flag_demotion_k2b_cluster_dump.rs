//! K2b (fix/aromaticity-flag-demotion-k2b) cluster diagnosis: for each
//! given SMILES, dumps everything needed to classify WHY chematic's default
//! Huckel model doesn't confirm a ring RDKit itself calls aromatic --
//! raw SSSR, the augmented ring set, each ring's Pass 1 electron trace
//! (per-atom contribution + `ContributionReason`), the converged
//! `ring_classifications()` (Pass 1 + Pass 2), and the final
//! `apply_aromaticity()` atom/bond flags. Diagnostic only, calls only
//! existing public API.
//!
//! Run:
//! ```text
//! cargo run -p chematic-perception --release \
//!     --example aromaticity_flag_demotion_k2b_cluster_dump \
//!     -- path/to/smiles_list.txt > dump.jsonl
//! ```

use std::fs;

use chematic_core::{BondIdx, BondOrder, Molecule};
use chematic_perception::{
    apply_aromaticity, assign_aromaticity, augmented_ring_set, find_sssr, trace_ring_pi_electrons,
};
use rustc_hash::FxHashSet;
use serde_json::json;

/// Whole-molecule set of bonds lying on ANY ring in `rings` -- mirrors
/// `aromaticity.rs`'s own private `ring_bond_set` (not exported), needed
/// here so `trace_ring_pi_electrons` can distinguish a genuine exocyclic
/// double bond from a ring-fusion bond into a different ring (K2b fix).
fn all_ring_bonds_of(mol: &Molecule, rings: &[Vec<chematic_core::AtomIdx>]) -> FxHashSet<BondIdx> {
    let mut out = FxHashSet::default();
    for ring in rings {
        for i in 0..ring.len() {
            if let Some((bidx, _)) = mol.bond_between(ring[i], ring[(i + 1) % ring.len()]) {
                out.insert(bidx);
            }
        }
    }
    out
}

fn main() {
    let path = std::env::args().nth(1).expect("usage: <smiles_list.txt>");
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));

    for line in text.lines() {
        let smiles = line.trim();
        if smiles.is_empty() || smiles.starts_with('#') {
            continue;
        }
        let raw = match chematic_smiles::parse(smiles) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("PARSE_FAIL {smiles:?}: {e}");
                continue;
            }
        };
        let kek_mol = match chematic_core::kekulize(&raw) {
            Ok(k) => chematic_core::apply_kekule(&raw, &k),
            Err(e) => {
                eprintln!("KEKULIZE_FAIL {smiles:?}: {e}");
                continue;
            }
        };

        let sssr = find_sssr(&kek_mol);
        let raw_sssr: Vec<Vec<u32>> = sssr
            .rings()
            .iter()
            .map(|r| r.iter().map(|a| a.0).collect())
            .collect();
        let augmented = augmented_ring_set(&kek_mol, sssr.rings());
        let augmented_dump: Vec<Vec<u32>> = augmented
            .iter()
            .map(|r| r.iter().map(|a| a.0).collect())
            .collect();

        // Pass 1 trace (empty aromatic_context) for every augmented ring.
        let empty_ctx: FxHashSet<chematic_core::AtomIdx> = FxHashSet::default();
        let all_ring_bonds = all_ring_bonds_of(&kek_mol, &augmented);
        let pass1_traces: Vec<_> = augmented
            .iter()
            .map(|ring| {
                let trace = trace_ring_pi_electrons(
                    &kek_mol,
                    ring,
                    &empty_ctx,
                    Default::default(),
                    &all_ring_bonds,
                );
                json!({
                    "ring": ring.iter().map(|a| a.0).collect::<Vec<u32>>(),
                    "total_pi": trace.total,
                    "atoms": trace.atoms.iter().map(|a| json!({
                        "idx": a.atom_idx.0,
                        "contribution": a.contribution,
                        "reason": format!("{:?}", a.reason),
                    })).collect::<Vec<_>>(),
                })
            })
            .collect();

        let model = assign_aromaticity(&kek_mol);
        let ring_classifications: Vec<_> = model
            .ring_classifications()
            .iter()
            .map(|(ring, cls, count)| {
                json!({
                    "ring": ring.iter().map(|a| a.0).collect::<Vec<u32>>(),
                    "classification": format!("{cls:?}"),
                    "pi_count": count,
                })
            })
            .collect();

        let applied = apply_aromaticity(&kek_mol);
        let atoms: Vec<_> = applied
            .atoms()
            .map(|(idx, atom)| {
                json!({
                    "idx": idx.0,
                    "element": atom.element.symbol(),
                    "charge": atom.charge,
                    "final_aromatic": atom.aromatic,
                    "model_aromatic": model.is_atom_aromatic(idx),
                })
            })
            .collect();
        let bonds: Vec<_> = applied
            .bonds()
            .map(|(idx, bond)| {
                json!({
                    "idx": idx.0,
                    "a1": bond.atom1.0,
                    "a2": bond.atom2.0,
                    "final_aromatic": bond.order == BondOrder::Aromatic,
                    "kekulized_order": format!("{:?}", kek_mol.bond(idx).order),
                })
            })
            .collect();

        let row = json!({
            "smiles": smiles,
            "raw_sssr": raw_sssr,
            "augmented_ring_set": augmented_dump,
            "pass1_traces": pass1_traces,
            "ring_classifications_converged": ring_classifications,
            "huckel_model_aromatic_atom_count": model.aromatic_atom_count(),
            "atoms": atoms,
            "bonds": bonds,
        });
        println!("{row}");
    }
}
