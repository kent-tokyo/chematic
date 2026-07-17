//! Canonical-Stereo-D0: for each corpus molecule, compare the legacy CIP
//! element multiset of the direct parse against the multiset after one
//! `canonical_smiles` round trip (parse -> canonical_smiles -> parse).
//! Reports counts of: stable (identical), "gained" (round trip produced an
//! element that wasn't there before -- the bug this PR fixes: a stashed
//! aromatic direction becoming a literal Up/Down order and suddenly
//! visible to `assign_ez`), "lost" (an element disappeared), and any other
//! change.
//!
//! Verification recipe used for the D0 PR: on `main` (before the fix),
//! 4,994/5,000 stable, 6 "gained-only" -- all 6 are the bug this PR fixes,
//! confirmed by hand. On the fix branch: 5,000/5,000 stable, 0 gained, 0
//! lost, 0 changed.
//!
//! Run:
//! ```text
//! cargo run -p chematic-chem --release --example canonical_stereo_d0_roundtrip_audit \
//!     -- ~/Downloads/SMILES.csv
//! ```

use chematic_chem::assign_cip;
use chematic_smiles::{canonical_smiles, parse};
use std::fs;

fn sorted_codes(mol: &chematic_core::Molecule) -> Vec<String> {
    let mut codes: Vec<String> = assign_cip(mol)
        .assignments
        .iter()
        .map(|(_, c)| format!("{c:?}"))
        .collect();
    codes.sort();
    codes
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| panic!("usage: canonical_stereo_d0_roundtrip_audit <smiles_file>"));
    let content = fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {path}: {e}"));

    let mut stable = 0;
    let mut gained = 0;
    let mut lost = 0;
    let mut changed_other = 0;
    let mut parse_fail = 0;

    for (i, line) in content.lines().enumerate() {
        let smi = line.trim();
        if smi.is_empty() {
            continue;
        }
        let mol = match parse(smi) {
            Ok(m) => m,
            Err(_) => {
                parse_fail += 1;
                continue;
            }
        };
        let before = sorted_codes(&mol);
        let c1 = canonical_smiles(&mol);
        let mol2 = match parse(&c1) {
            Ok(m) => m,
            Err(_) => {
                parse_fail += 1;
                continue;
            }
        };
        let after = sorted_codes(&mol2);

        if before == after {
            stable += 1;
            continue;
        }

        let before_set: std::collections::HashSet<&String> = before.iter().collect();
        let after_set: std::collections::HashSet<&String> = after.iter().collect();
        let only_gained = after_set.difference(&before_set).count() > 0
            && before_set.difference(&after_set).count() == 0
            && before.len() < after.len();
        let only_lost = before_set.difference(&after_set).count() > 0
            && after_set.difference(&before_set).count() == 0
            && after.len() < before.len();

        if only_gained {
            gained += 1;
            println!("GAINED\t{i}\t{smi}\tbefore={before:?}\tafter={after:?}");
        } else if only_lost {
            lost += 1;
            println!("LOST\t{i}\t{smi}\tbefore={before:?}\tafter={after:?}");
        } else {
            changed_other += 1;
            println!("CHANGED\t{i}\t{smi}\tbefore={before:?}\tafter={after:?}");
        }
    }

    eprintln!("\n=== summary ===");
    eprintln!("stable:        {stable}");
    eprintln!("gained only:   {gained}");
    eprintln!("lost only:     {lost}");
    eprintln!("changed other: {changed_other}");
    eprintln!("parse_fail:    {parse_fail}");
}
