//! SMARTS-R2: dump chematic's default AND opt-in RDKit-parity SMARTS
//! match-sets across a molecule corpus, one JSON row per molecule, for
//! cross-checking against a live RDKit oracle in
//! `scripts/rdkit_ring_parity_diagnosis.py`.
//!
//! Diagnostic only. Calls only this crate's existing public API
//! (`chematic_smarts::{find_matches, find_matches_rdkit_parity, parse_smarts}`)
//! exactly as an external caller would; does not modify any production
//! module. See `crates/chematic-smarts/src/rdkit_ring_model.rs`'s module doc
//! comment for the root-cause/design rationale this diagnosis is measuring,
//! and `docs/rdkit_compat.md`'s "SMARTS-R2" section for the write-up this
//! feeds.
//!
//! Run:
//! ```text
//! cargo run -p chematic-smarts --release --example rdkit_parity_dump -- \
//!     ~/Downloads/SMILES.csv > validation/results/rdkit_ring_parity_dump.jsonl
//! /tmp/chematic-smartsC-venv/bin/python scripts/rdkit_ring_parity_diagnosis.py
//! ```

use chematic_smarts::{RdkitParityConfig, find_matches, find_matches_rdkit_parity, parse_smarts};
use serde_json::json;
use std::io::Write;

/// Pattern families named in the task spec: ring-membership (`R`/`RN`/`r`/
/// `rN`/`k`/`kN`), aromatic atoms, ring-bond topology (`@`/`!@`), recursive/
/// boolean SMARTS, plus the original 16-pattern `rdkit_compat_diff.py` list
/// (kept so this diagnosis's numbers are directly comparable to the
/// existing SMARTS-R0/R1 figures for the patterns they share).
const PATTERNS: &[&str] = &[
    // -- existing 16-pattern list (SMARTS-R0/R1 baseline) --
    "[OH]",
    "c",
    "[#7]",
    "C=O",
    "[NX3;H2,H1;!$(NC=O)]",
    "[r5]",
    "[r6]",
    "c1ccccc1",
    "[CX4]",
    "[!#6;!#1]",
    "C(=O)[OH]",
    "[nH]",
    "[#6]=[#6]",
    "[F,Cl,Br,I]",
    "[OX2H]",
    "[#16]",
    // -- new: ring-membership-count, this track's actual target --
    "[R]",
    "[R0]",
    "[R1]",
    "[R2]",
    "[R3]",
    // -- new: ring-bond-count --
    "[x2]",
    "[x3]",
    // -- new: any-ring-of-size-N (kN), for regression coverage --
    "[k5]",
    "[k6]",
    // -- new: ring-bond topology --
    "*@*",
    "*!@*",
    // -- new: recursive + boolean SMARTS --
    "[$(C=O)]",
    "[N,O]",
    "[C;!R]",
    "[C;R;!$(C=O)]",
];

/// Hand-built structure classes deliberately under-represented in a random
/// 5,000-molecule drug-like corpus: fused/bridged/spiro/cage systems and the
/// specific scaffolds `docs/rdkit_compat.md` already names (SMARTS-A0's
/// bridgehead-N reproducer, SMARTS-R0/R2's adamantane-type cage class).
/// (id, smiles). SMILES for the cages are literature/PubChem canonical
/// forms, independently re-verified against a live RDKit oracle (see
/// `scripts/rdkit_ring_parity_diagnosis.py`'s alignment check, which would
/// fail loudly if any of these didn't round-trip).
const HAND_CORPUS: &[(&str, &str)] = &[
    ("benzene", "c1ccccc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("azulene", "c1ccc2cccc-2cc1"),
    ("indolizine", "c1ccn2ccccc12"),
    ("purine", "c1ncc2[nH]cnc2n1"),
    ("quinazoline", "c1ccc2ncncc2c1"),
    ("quinoxaline", "c1ccc2nccnc2c1"),
    ("decalin", "C1CCC2CCCCC2C1"),
    ("spiro44nonane", "C1CCC2(C1)CCCC2"),
    ("spiro55undecane", "C1CCCC2(C1)CCCCC2"),
    ("norbornane", "C1CC2CCC1C2"),
    ("bicyclo222octane", "C1CC2CCC1CC2"),
    ("adamantane", "C1C2CC3CC1CC(C2)C3"),
    ("cubane", "C12C3C4C1C5C4C3C25"),
    ("prismane", "C12C3C4C1C1C2C3C41"),
    (
        "dodecahedrane",
        "C12C3C4C5C1C6C7C2C8C3C9C4C1C5C6C2C7C8C9C12",
    ),
    ("morphinan_core", "C1CC2CCCC3C2C1CCN3"),
    ("steroid_core", "C1CCC2CCC3C(CCC4CCCCC34)C2C1"),
    // SMARTS-A0's bridgehead-N over-aromatization reproducer (bare core,
    // `docs/rdkit_compat.md`) -- included so this run also directly measures
    // that named residual, not just the ring-count one.
    ("bridgehead_n_smarts_a0", "C1=Cc2ccccc2C2=NCCCN12"),
    ("quinuclidine", "C1CN2CCC1CC2"),
    ("coronene", "c1cc2ccc3ccc4ccc5ccc1c1c2c3c4c51"),
];

fn atom_elements(mol: &chematic_core::Molecule) -> Vec<String> {
    (0..mol.atom_count())
        .map(|i| {
            mol.atom(chematic_core::AtomIdx(i as u32))
                .element
                .symbol()
                .to_string()
        })
        .collect()
}

fn match_set_json(
    matches: &[rustc_hash::FxHashMap<usize, chematic_core::AtomIdx>],
) -> Vec<Vec<u32>> {
    let mut sets: Vec<Vec<u32>> = matches
        .iter()
        .map(|m| {
            let mut v: Vec<u32> = m.values().map(|a| a.0).collect();
            v.sort_unstable();
            v
        })
        .collect();
    sets.sort();
    sets
}

fn dump_one(id: &str, smiles: &str, stdout: &mut impl Write) {
    let Ok(mol) = chematic_smiles::parse(smiles) else {
        return;
    };
    let compiled: Vec<(&str, Option<_>)> = PATTERNS
        .iter()
        .map(|p| (*p, parse_smarts(p).ok()))
        .collect();

    let mut per_pattern = serde_json::Map::new();
    for (pat, q) in &compiled {
        let Some(query) = q else {
            per_pattern.insert((*pat).to_string(), json!({"parse_error": true}));
            continue;
        };
        let default_matches = find_matches(query, &mol);
        let parity_result = find_matches_rdkit_parity(query, &mol, &RdkitParityConfig::default());
        let entry = match parity_result {
            Ok((parity_matches, budget_exhausted)) => json!({
                "default": match_set_json(&default_matches),
                "parity": match_set_json(&parity_matches),
                "parity_budget_exhausted": budget_exhausted,
            }),
            Err(e) => json!({
                "default": match_set_json(&default_matches),
                "parity_error": format!("{e:?}"),
            }),
        };
        per_pattern.insert((*pat).to_string(), entry);
    }

    let row = json!({
        "id": id,
        "smiles": smiles,
        "atom_elements": atom_elements(&mol),
        "patterns": per_pattern,
    });
    writeln!(stdout, "{row}").expect("stdout write failed");
}

fn main() {
    let stdout = std::io::stdout();
    let mut out = stdout.lock();

    for (id, smi) in HAND_CORPUS {
        dump_one(id, smi, &mut out);
    }

    let path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| shellexpand_home("~/Downloads/SMILES.csv"));
    if let Ok(contents) = std::fs::read_to_string(&path) {
        for (i, line) in contents.lines().enumerate() {
            let smi = line.trim();
            if smi.is_empty() {
                continue;
            }
            dump_one(&format!("corpus_{i}"), smi, &mut out);
        }
    } else {
        eprintln!("warning: could not read corpus file {path}; only HAND_CORPUS was dumped");
    }
}

/// Minimal `~` expansion (no shellexpand dependency for one path).
fn shellexpand_home(path: &str) -> String {
    if let Some(rest) = path.strip_prefix("~/")
        && let Ok(home) = std::env::var("HOME")
    {
        return format!("{home}/{rest}");
    }
    path.to_string()
}
