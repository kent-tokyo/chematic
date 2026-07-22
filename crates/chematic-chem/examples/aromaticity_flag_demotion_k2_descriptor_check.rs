//! fix/aromaticity-flag-demotion-k2: downstream cross-check (read-only) --
//! computes the 14 descriptor values `docs/descriptor_census_rfc.md` §8
//! tags `depends-on-aromaticity` for the same 5000-molecule corpus PR #137's
//! census used, so a before/after diff can be taken without re-running the
//! full Python `descriptor_census.py` pipeline. Does not modify
//! `chematic-chem` production code -- calls only its existing public API.
//!
//! Run once per side of the fix (see the K2 PR description for how the
//! before/after dumps were produced) and diff the two JSONL outputs by
//! `smiles`.
//!
//! ```text
//! cargo run -p chematic-chem --release \
//!     --example aromaticity_flag_demotion_k2_descriptor_check \
//!     -- scripts/descriptor_census_corpus.smi > /tmp/desc_after.jsonl
//! ```

use std::fs;

use chematic_chem::descriptors::{
    aromatic_ring_count, egan_passes, fsp3, hba_count, hba_count_lipinski, hbd_count,
    num_aliphatic_heterocycles, num_aliphatic_rings, num_aromatic_heterocycles,
    num_saturated_heterocycles, num_saturated_rings, pfizer_3_75_passes, tpsa, tpsa_per_atom,
};

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "scripts/descriptor_census_corpus.smi".to_string());
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("failed to read {path}: {e}"));

    let mut n = 0usize;
    let mut n_fail = 0usize;
    for line in text.lines() {
        let smiles = line.trim();
        if smiles.is_empty() {
            continue;
        }
        n += 1;
        let mol = match chematic_smiles::parse(smiles) {
            Ok(m) => m,
            Err(e) => {
                n_fail += 1;
                eprintln!("PARSE_FAIL: {smiles:?}: {e}");
                continue;
            }
        };

        let tpa: Vec<String> = tpsa_per_atom(&mol)
            .iter()
            .map(|v| format!("{v:.6}"))
            .collect();
        println!(
            "{{\"smiles\":{smiles:?},\"hbd_count\":{},\"hba_count\":{},\"tpsa\":{:.6},\"fsp3\":{:.6},\"aromatic_ring_count\":{},\"hba_count_lipinski\":{},\"num_aliphatic_rings\":{},\"num_saturated_rings\":{},\"num_aromatic_heterocycles\":{},\"num_aliphatic_heterocycles\":{},\"num_saturated_heterocycles\":{},\"egan_passes\":{},\"pfizer_3_75_passes\":{},\"tpsa_per_atom\":[{}]}}",
            hbd_count(&mol),
            hba_count(&mol),
            tpsa(&mol),
            fsp3(&mol),
            aromatic_ring_count(&mol),
            hba_count_lipinski(&mol),
            num_aliphatic_rings(&mol),
            num_saturated_rings(&mol),
            num_aromatic_heterocycles(&mol),
            num_aliphatic_heterocycles(&mol),
            num_saturated_heterocycles(&mol),
            egan_passes(&mol),
            pfizer_3_75_passes(&mol),
            tpa.join(",")
        );
    }
    eprintln!("dumped {} molecules, {n_fail} parse failures", n - n_fail);
}
