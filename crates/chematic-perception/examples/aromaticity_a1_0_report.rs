//! Aromaticity-A1-0: fused-component characterization report.
//!
//! Diagnostic only -- does not change `assign_aromaticity_ex`'s behavior.
//! Reads `validation/aromaticity_a1_0_corpus.jsonl` (three buckets:
//! `false_positive`, `false_negative`, `negative_control` -- see
//! `scripts/gen_aromaticity_a1_0_corpus.py` for provenance), and for every
//! molecule:
//!
//! 1. builds fused candidate components from the *augmented* ring list
//!    (the same ring list `assign_aromaticity_ex`'s Pass 1/Pass 2 actually
//!    iterate -- not raw SSSR, see `find_ring_families_over`'s doc comment);
//! 2. traces each ring atom's pi-electron contribution reason
//!    (`trace_ring_pi_electrons`), in two contexts: intrinsic (empty,
//!    Pass-1-equivalent) and full (the model's converged aromatic_context,
//!    an observational Pass-2 upper bound -- not a literal iteration replay);
//! 3. reports each ring's electron total and Hückel verdict under both
//!    contexts, alongside the real engine's final per-atom verdict
//!    (`assign_aromaticity_ex(mol).is_atom_aromatic`);
//! 4. reports the component's cycle rank (independent ring count) and
//!    topology kind (Simple/Fused/Spiro/Bridged).
//!
//! Emits one JSONL row per (molecule, ring, atom) to stdout. RDKit's own
//! per-atom aromatic verdict is *not* added here (RDKit bindings are
//! Python-only in this project) -- `scripts/aromaticity_a1_0_diagnosis.py`
//! joins this output against RDKit and runs the false-positive/
//! false-negative polarization check.
//!
//! Run:
//! ```text
//! cargo run -p chematic-perception --release --example aromaticity_a1_0_report \
//!     -- validation/aromaticity_a1_0_corpus.jsonl \
//!     > validation/results/aromaticity_a1_0_trace.jsonl
//! ```

use std::fs;

use chematic_core::{AtomIdx, BondIdx, Molecule};
use chematic_perception::{
    AromaticityAlgorithm, assign_aromaticity_ex, augmented_ring_set, exhaustive_aromaticity_oracle,
    find_ring_families_over, find_sssr, trace_ring_pi_electrons,
};
use rustc_hash::FxHashSet;
use serde_json::{Value, json};

/// Whole-molecule set of bonds lying on ANY ring in `rings` -- mirrors
/// `aromaticity.rs`'s own private `ring_bond_set` (not exported), needed
/// here so `trace_ring_pi_electrons` can distinguish a genuine exocyclic
/// double bond from a ring-fusion bond into a different ring (K2b fix).
fn all_ring_bonds_of(mol: &Molecule, rings: &[Vec<AtomIdx>]) -> FxHashSet<BondIdx> {
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

/// Same 4n+2 rule as the crate's private `classify_ring_aromaticity` --
/// duplicated here (not exposed as a second public API) since it's a
/// one-line, stable rule and this is a diagnostic-only example.
fn is_huckel_aromatic(pi: Option<u32>) -> bool {
    matches!(pi, Some(p) if p >= 2 && (p - 2) % 4 == 0)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let corpus_path = args
        .get(1)
        .map(String::as_str)
        .unwrap_or("validation/aromaticity_a1_0_corpus.jsonl");

    let corpus = fs::read_to_string(corpus_path)
        .unwrap_or_else(|e| panic!("failed to read {corpus_path}: {e}"));

    let algo = AromaticityAlgorithm::RdkitLike;
    let mut rows_written = 0usize;
    let mut molecules_seen = 0usize;
    let mut parse_failures = 0usize;

    for line in corpus.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let case: Value = serde_json::from_str(line).expect("valid corpus JSON line");
        let bucket = case["bucket"].as_str().unwrap_or("unknown").to_string();
        let case_id = case["case_id"].as_str().unwrap_or("").to_string();
        let smiles = case["smiles"].as_str().unwrap_or("").to_string();

        // Kekulize every corpus molecule uniformly, matching the exact
        // input representation the pinned corpus this data comes from uses
        // (`mol_kekulized` in aromaticity.rs's test module) -- this isolates
        // Pass 1/Pass 2's own Hückel-counting logic from the SEPARATE,
        // already-diagnosed parser aromatic-flag behavior (see
        // docs/rdkit_compat.md's SMARTS-A0 section). Note: purine's known
        // false negative only reproduces on Kekulized input -- feeding its
        // corpus SMILES (aromatic lowercase) through the real production
        // path (parse -> apply_aromaticity, no explicit kekulize) does NOT
        // reproduce it, since that path never routes through the
        // exocyclic-to-heteroatom 0π carbon rule at all (it requires an
        // explicit `BondOrder::Double`, which aromatic-form parsing never
        // produces). This representation-dependence is itself a real A1-0
        // finding, not a corpus bug -- see docs/rfcs/aromaticity_a1_rfc.md.
        let mol = match chematic_smiles::parse(&smiles) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("SKIP {case_id} ({smiles}): parse error: {e}");
                parse_failures += 1;
                continue;
            }
        };
        let mol = match chematic_core::kekulize(&mol) {
            Ok(k) => chematic_core::apply_kekule(&mol, &k),
            Err(e) => {
                eprintln!("SKIP {case_id} ({smiles}): kekulize error: {e}");
                parse_failures += 1;
                continue;
            }
        };
        molecules_seen += 1;

        let model = assign_aromaticity_ex(&mol, algo);
        let final_context: FxHashSet<AtomIdx> = mol
            .atoms()
            .map(|(idx, _)| idx)
            .filter(|&idx| model.is_atom_aromatic(idx))
            .collect();
        let empty_context: FxHashSet<AtomIdx> = FxHashSet::default();

        // Aromaticity-A1-1a: the exhaustive-candidate reference oracle,
        // computed once per molecule (whole-molecule, not per-ring).
        let (oracle_atoms, _oracle_bonds) = exhaustive_aromaticity_oracle(&mol, algo);

        let sssr = find_sssr(&mol);
        let rings = augmented_ring_set(&mol, sssr.rings());
        let families = find_ring_families_over(&mol, &rings);
        let all_ring_bonds = all_ring_bonds_of(&mol, &rings);

        for (component_id, family) in families.iter().enumerate() {
            let cycle_rank = family.ring_indices.len();
            let kind = format!("{:?}", family.kind);

            for &ring_idx in &family.ring_indices {
                let ring = &rings[ring_idx];
                let intrinsic =
                    trace_ring_pi_electrons(&mol, ring, &empty_context, algo, &all_ring_bonds);
                let context =
                    trace_ring_pi_electrons(&mol, ring, &final_context, algo, &all_ring_bonds);
                let ring_aromatic_intrinsic = is_huckel_aromatic(intrinsic.total);
                let ring_aromatic_context = is_huckel_aromatic(context.total);

                for (a_intrinsic, a_context) in intrinsic.atoms.iter().zip(context.atoms.iter()) {
                    debug_assert_eq!(a_intrinsic.atom_idx, a_context.atom_idx);
                    let atom_idx = a_intrinsic.atom_idx;
                    let element = mol.atom(atom_idx).element.symbol();

                    let row = json!({
                        "bucket": bucket,
                        "case_id": case_id,
                        "smiles": smiles,
                        "component_id": component_id,
                        "cycle_rank": cycle_rank,
                        "ring_system_kind": kind,
                        "ring_idx": ring_idx,
                        "ring_size": ring.len(),
                        "atom_idx": atom_idx.0,
                        "element": element,
                        "candidate_intrinsic": a_intrinsic.reason.is_eligible(),
                        "contribution_intrinsic": a_intrinsic.contribution,
                        "reason_intrinsic": format!("{:?}", a_intrinsic.reason),
                        "candidate_context": a_context.reason.is_eligible(),
                        "contribution_context": a_context.contribution,
                        "reason_context": format!("{:?}", a_context.reason),
                        "ring_electron_total_intrinsic": intrinsic.total,
                        "ring_electron_total_context": context.total,
                        "ring_aromatic_intrinsic": ring_aromatic_intrinsic,
                        "ring_aromatic_context": ring_aromatic_context,
                        "current_engine_atom_aromatic": model.is_atom_aromatic(atom_idx),
                        "oracle_atom_aromatic": oracle_atoms.contains(&atom_idx),
                    });
                    println!("{row}");
                    rows_written += 1;
                }
            }
        }
    }

    eprintln!(
        "molecules: {molecules_seen} (parse failures: {parse_failures}), rows: {rows_written}"
    );
}
