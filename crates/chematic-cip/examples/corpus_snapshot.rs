//! Snapshot `assign_cip_accurate_experimental` (or its stable
//! `..._without_mancude` reference point) across every stereocenter in a full SMILES
//! corpus, for full-corpus regression/accuracy verification -- the methodology Milestone
//! 3B-1b used ad hoc (see `docs/rfcs/cip_accurate_rfc.md`), formalized here as a checked-in,
//! reusable tool rather than a one-off scratchpad script, since every future full-corpus
//! CIP gate (Milestone 4 included) needs the same rigor again.
//!
//! `SMILES.csv` itself is not checked into this repo (same convention as
//! `scripts/bench5k.py`/`scripts/cip_ground_truth.py`) -- this is a `cargo run --example`
//! tool for local, on-demand use against a corpus file the user supplies, not a
//! `scripts/check.sh`/CI step.
//!
//! Usage:
//!   cargo run -p chematic-cip --release --example corpus_snapshot -- \
//!     --baseline|--candidate <SMILES.csv> <out.tsv>
//!
//! Output: one TSV row per candidate stereocenter (`smiles\tatom_idx\tvalue`), where
//! `value` is `R`/`S`/`E`/`Z`, `skip:tied`/`skip:budget`/`skip:not4`, or `ERR\t<message>`.
//! Row selection matches `scripts/cip_ground_truth.py`'s corpus convention (skip the CSV
//! header line) but gates on `chematic`'s own chirality flag, not RDKit's -- see
//! `scripts/cip_accurate_full_corpus_report.py`'s module docs for why that set can differ
//! slightly from the oracle's own stereocenter selection, and why that's a distinct cause
//! from an actual correctness disagreement.

use std::env;
use std::fs;
use std::io::Write;
use std::time::Instant;

use chematic_cip::{
    CipBudget, SkipReason, assign_cip_accurate_experimental,
    assign_cip_accurate_experimental_without_mancude,
};
use chematic_core::{AtomIdx, Chirality, CipCode};

fn code_str(c: CipCode) -> &'static str {
    match c {
        CipCode::R => "R",
        CipCode::S => "S",
        CipCode::E => "E",
        CipCode::Z => "Z",
        CipCode::LowerR => "r",
        CipCode::LowerS => "s",
    }
}

fn skip_str(r: SkipReason) -> &'static str {
    match r {
        SkipReason::NotFourSubstituents => "skip:not4",
        SkipReason::Tied => "skip:tied",
        SkipReason::BudgetExceeded => "skip:budget",
        SkipReason::OracleUnstable => "skip:oracle-unstable",
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mode = args.get(1).map(String::as_str).unwrap_or_else(|| {
        panic!("usage: corpus_snapshot --baseline|--candidate <SMILES.csv> <out.tsv>")
    });
    let candidate = match mode {
        "--baseline" => false,
        "--candidate" => true,
        other => panic!("first arg must be --baseline or --candidate, got {other}"),
    };
    let csv_path = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| format!("{}/Downloads/SMILES.csv", env::var("HOME").unwrap()));
    let out_path = args
        .get(3)
        .cloned()
        .unwrap_or_else(|| "snapshot.tsv".to_string());

    let content = fs::read_to_string(&csv_path).expect("read SMILES.csv");
    let smis: Vec<&str> = content
        .lines()
        .skip(1)
        .filter(|l| !l.trim().is_empty())
        .collect();
    let input_rows = smis.len();

    let mut lines: Vec<String> = Vec::new();
    let mut parsed_rows = 0usize;
    let mut molecules_with_chirality = 0usize;
    let mut total_nanos: u128 = 0;
    let mut per_molecule_nanos: Vec<u128> = Vec::new();

    for smi in &smis {
        let mol = match chematic_smiles::parse(smi) {
            Ok(m) => m,
            Err(_) => continue,
        };
        parsed_rows += 1;
        let has_chirality =
            (0..mol.atom_count()).any(|i| mol.atom(AtomIdx(i as u32)).chirality != Chirality::None);
        if !has_chirality {
            continue;
        }
        molecules_with_chirality += 1;

        let budget = CipBudget::default_budget();
        let start = Instant::now();
        let result = if candidate {
            assign_cip_accurate_experimental(&mol, budget)
        } else {
            assign_cip_accurate_experimental_without_mancude(&mol, budget)
        };
        let elapsed = start.elapsed().as_nanos();
        total_nanos += elapsed;
        per_molecule_nanos.push(elapsed);

        match result {
            Ok(result) => {
                for (idx, code) in &result.assignments {
                    lines.push(format!("{smi}\t{}\t{}", idx.0, code_str(*code)));
                }
                for (idx, reason) in &result.skipped {
                    lines.push(format!("{smi}\t{}\t{}", idx.0, skip_str(*reason)));
                }
            }
            Err(e) => {
                lines.push(format!("{smi}\tERR\t{e}"));
            }
        }
    }

    lines.sort();
    let mut f = fs::File::create(&out_path).expect("create output");
    for l in &lines {
        writeln!(f, "{l}").unwrap();
    }

    per_molecule_nanos.sort_unstable();
    let median_ns = per_molecule_nanos
        .get(per_molecule_nanos.len() / 2)
        .copied()
        .unwrap_or(0);
    let p95_idx = (per_molecule_nanos.len() * 95) / 100;
    let p95_ns = per_molecule_nanos
        .get(p95_idx.min(per_molecule_nanos.len().saturating_sub(1)))
        .copied()
        .unwrap_or(0);

    eprintln!(
        "mode={} input_rows={input_rows} parsed_rows={parsed_rows} \
         molecules_with_chirality={molecules_with_chirality} rows_written={} \
         total_ms={:.1} median_us={:.2} p95_us={:.2}",
        if candidate { "candidate" } else { "baseline" },
        lines.len(),
        total_nanos as f64 / 1_000_000.0,
        median_ns as f64 / 1_000.0,
        p95_ns as f64 / 1_000.0,
    );
}
