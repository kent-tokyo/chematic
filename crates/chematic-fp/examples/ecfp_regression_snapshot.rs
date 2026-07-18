//! ECFP RDKit-invariant-mode PR: non-regression snapshot, run once before
//! this PR's `ecfp.rs` `EcfpInvariantMode` refactor and once after. Uses only
//! stable, pre-existing public API (no `EcfpInvariantMode`/`*_rdkit_invariants`
//! reference) so the exact same file compiles against both revisions.
//!
//! Every existing entry point (`ecfp4`, `ecfp6`, `ecfp` with chirality,
//! `ecfp_with_bitinfo` fp + origins, `morgan_fp_counts`) must be byte-for-byte
//! identical before/after, since none of these callers changed which
//! invariant mode they use.
//!
//! Usage:
//! ```text
//! cargo run -p chematic-fp --release --example ecfp_regression_snapshot \
//!     -- <SMILES.csv> <out.tsv>
//! ```

use chematic_fp::{EcfpConfig, ecfp, ecfp_with_bitinfo, ecfp4, ecfp6, morgan_fp_counts};
use chematic_smiles::parse;
use std::fs;
use std::io::Write;

fn bits_string(fp: &chematic_fp::BitVec2048) -> String {
    (0..2048)
        .filter(|&i| fp.get(i))
        .map(|i| i.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn row(smi: &str) -> Option<String> {
    let mol = parse(smi).ok()?;

    let fp4 = ecfp4(&mol);
    let fp6 = ecfp6(&mol);
    let fp4_chiral = ecfp(
        &mol,
        &EcfpConfig {
            use_chirality: true,
            ..EcfpConfig::default()
        },
    );
    let (bi_fp, bi_info) = ecfp_with_bitinfo(&mol, &EcfpConfig::default());
    let mut bi_origins: Vec<(usize, Vec<(u32, u32)>)> = bi_info.into_iter().collect();
    bi_origins.sort_by_key(|(bit, _)| *bit);
    for (_, envs) in &mut bi_origins {
        envs.sort_unstable();
    }
    let bi_origins_str = bi_origins
        .iter()
        .map(|(bit, envs)| format!("{bit}:{envs:?}"))
        .collect::<Vec<_>>()
        .join(";");

    let counts = morgan_fp_counts(&mol, 2);
    let mut counts_vec: Vec<(u64, u32)> = counts.into_iter().collect();
    counts_vec.sort_unstable();
    let counts_str = counts_vec
        .iter()
        .map(|(h, c)| format!("{h}:{c}"))
        .collect::<Vec<_>>()
        .join(",");

    Some(format!(
        "{smi}\t{}\t{}\t{}\t{}\t{}\t{}",
        bits_string(&fp4),
        bits_string(&fp6),
        bits_string(&fp4_chiral),
        bits_string(&bi_fp),
        bi_origins_str,
        counts_str,
    ))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let csv_path = args
        .get(1)
        .cloned()
        .unwrap_or_else(|| panic!("usage: ecfp_regression_snapshot <SMILES.csv> <out.tsv>"));
    let out_path = args
        .get(2)
        .cloned()
        .unwrap_or_else(|| "snapshot.tsv".to_string());

    let content = fs::read_to_string(&csv_path).unwrap_or_else(|e| panic!("read {csv_path}: {e}"));

    let mut lines: Vec<String> = Vec::new();
    let mut parse_fail = 0usize;
    for line in content.lines() {
        let smi = line.trim();
        if smi.is_empty() {
            continue;
        }
        match row(smi) {
            Some(r) => lines.push(r),
            None => parse_fail += 1,
        }
    }

    let mut f = fs::File::create(&out_path).unwrap_or_else(|e| panic!("create {out_path}: {e}"));
    for l in &lines {
        writeln!(f, "{l}").unwrap();
    }

    eprintln!(
        "input_lines={} parse_fail={parse_fail} rows={} out={out_path}",
        content.lines().filter(|l| !l.trim().is_empty()).count(),
        lines.len(),
    );
}
