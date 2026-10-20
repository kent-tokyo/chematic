//! Emit source-built RDKit-compatible ring counts for a SMILES corpus.

use std::io::{self, BufRead};

use chematic_chem::rdkit_num_rings;
use serde_json::json;

fn main() {
    for (index, line) in io::stdin().lock().lines().enumerate() {
        let line = line.expect("read SMILES row");
        let smiles = line.split_whitespace().next().expect("nonempty SMILES row");
        let mol = chematic_smiles::parse(smiles).expect("parse corpus SMILES");
        let relevant_rings = chematic_perception::relevant_cycles(&mol, 1_000_000).ok();
        let legacy_rings = chematic_perception::rdkit_sssr_ring_order(&mol);
        println!(
            "{}",
            json!({
                "index": index,
                "smiles": smiles,
                "num_rings": rdkit_num_rings(&mol),
                "relevant_cycles": relevant_rings.as_ref().map(Vec::len),
                "legacy_sssr": legacy_rings.as_ref().map(Vec::len),
                "legacy_max_ring_size": legacy_rings.as_ref().and_then(|rings| rings.iter().map(Vec::len).max()),
                "relevant_max_ring_size": relevant_rings.as_ref().and_then(|rings| rings.iter().map(Vec::len).max()),
            })
        );
    }
}
