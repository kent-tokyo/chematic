//! Time `find_match_atom_sets_rdkit_parity` over a SMILES corpus with the
//! cross-engine benchmark's 31 SMARTS queries (profiling aid).
//!
//! ```text
//! cargo run -p chematic-smarts --release --example rdkit_parity_profile -- \
//!     scripts/chembl_accuracy_corpus_4999.smi validation/rdkit_rebaseline_smarts_queries.json
//! ```

use chematic_smarts::{RdkitParityConfig, find_match_atom_sets_rdkit_parity, parse_smarts};
use std::time::{Duration, Instant};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let corpus = std::fs::read_to_string(&args[1]).expect("corpus");
    let queries_json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&args[2]).expect("queries")).expect("json");
    let reps: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1);
    let queries: Vec<String> = queries_json["queries"]
        .as_array()
        .expect("queries array")
        .iter()
        .map(|q| q.as_str().expect("query").to_string())
        .collect();
    let compiled: Vec<_> = queries
        .iter()
        .map(|q| parse_smarts(q).expect("smarts"))
        .collect();
    let smiles: Vec<&str> = corpus
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .collect();
    let config = RdkitParityConfig {
        use_rdkit_parity_aromaticity: true,
        ..RdkitParityConfig::default()
    };
    let mut per_query = vec![Duration::ZERO; compiled.len()];
    let mut parse = Duration::ZERO;
    let mut total_matches = 0usize;
    let start = Instant::now();
    for _ in 0..reps {
        for smi in &smiles {
            let t = Instant::now();
            let Ok(mol) = chematic_smiles::parse(smi) else {
                continue;
            };
            parse += t.elapsed();
            for (i, q) in compiled.iter().enumerate() {
                let t = Instant::now();
                if let Ok((m, _)) = find_match_atom_sets_rdkit_parity(q, &mol, &config) {
                    total_matches += m.len();
                }
                per_query[i] += t.elapsed();
            }
        }
    }
    let total = start.elapsed();
    println!(
        "total {:.1} ms, parse {:.1} ms, matches {total_matches}",
        ms(total),
        ms(parse)
    );
    let mut order: Vec<usize> = (0..compiled.len()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(per_query[i]));
    for i in order {
        println!("{:8.1} ms  {}", ms(per_query[i]), queries[i]);
    }
}

fn ms(d: Duration) -> f64 {
    d.as_secs_f64() * 1000.0
}
