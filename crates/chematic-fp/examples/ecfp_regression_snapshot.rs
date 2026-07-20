//! Non-regression snapshot for every existing `chematic-fp` Morgan/ECFP
//! entry point -- run once before a PR touching this crate's internals and
//! once after; the two runs must be byte-for-byte identical whenever the PR
//! claims not to change any existing function's output. Originally written
//! for PR #110's `EcfpInvariantMode` refactor (hence the narrow original
//! set); extended for Phase B (`rdkit_morgan_ecfp4_experimental`, PR #124
//! follow-up) to cover every existing public fingerprint function in
//! `ecfp.rs`, since that PR reuses internals
//! ([`chematic_fp`]'s own `rdkit_morgan_hash` module) that are adjacent to,
//! but must not affect, this set.
//!
//! Covers: `ecfp4`, `ecfp6`, `ecfp` with chirality, `ecfp_with_bitinfo` fp +
//! origins, `morgan_fp_counts`, `ecfp4_rdkit_invariants`,
//! `ecfp6_rdkit_invariants`, `ecfp4_rdkit_environment_experimental`,
//! `ecfp6_rdkit_environment_experimental`,
//! `ecfp_with_bitinfo_rdkit_environment_experimental`.
//!
//! Usage:
//! ```text
//! cargo run -p chematic-fp --release --example ecfp_regression_snapshot \
//!     -- <SMILES.csv> <out.tsv>
//! ```

use chematic_fp::{
    EcfpConfig, ecfp, ecfp_with_bitinfo, ecfp_with_bitinfo_rdkit_environment_experimental, ecfp4,
    ecfp4_rdkit_environment_experimental, ecfp4_rdkit_invariants, ecfp6,
    ecfp6_rdkit_environment_experimental, ecfp6_rdkit_invariants, morgan_fp_counts,
};
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

    let rdkit_inv4 = ecfp4_rdkit_invariants(&mol);
    let rdkit_inv6 = ecfp6_rdkit_invariants(&mol);
    let rdkit_env4 = ecfp4_rdkit_environment_experimental(&mol);
    let rdkit_env6 = ecfp6_rdkit_environment_experimental(&mol);
    let (rdkit_env_bi_fp, rdkit_env_bi_info) =
        ecfp_with_bitinfo_rdkit_environment_experimental(&mol, &EcfpConfig::default());
    let mut rdkit_env_bi_origins: Vec<(usize, Vec<(u32, u32)>)> =
        rdkit_env_bi_info.into_iter().collect();
    rdkit_env_bi_origins.sort_by_key(|(bit, _)| *bit);
    for (_, envs) in &mut rdkit_env_bi_origins {
        envs.sort_unstable();
    }
    let rdkit_env_bi_origins_str = rdkit_env_bi_origins
        .iter()
        .map(|(bit, envs)| format!("{bit}:{envs:?}"))
        .collect::<Vec<_>>()
        .join(";");

    Some(format!(
        "{smi}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        bits_string(&fp4),
        bits_string(&fp6),
        bits_string(&fp4_chiral),
        bits_string(&bi_fp),
        bi_origins_str,
        counts_str,
        bits_string(&rdkit_inv4),
        bits_string(&rdkit_inv6),
        bits_string(&rdkit_env4),
        bits_string(&rdkit_env6),
        bits_string(&rdkit_env_bi_fp),
        rdkit_env_bi_origins_str,
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
