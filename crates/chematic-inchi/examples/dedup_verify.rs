//! Verification driver for `chematic_inchi::dedup`, used to check:
//!
//! 1. Determinism across independent process runs (run this binary 5 times
//!    as 5 separate OS processes and diff the output -- NOT an in-process
//!    loop, which would not catch process-level nondeterminism such as
//!    HashMap iteration order or a mis-seeded RNG).
//! 2. A worst-of-N permutation sweep over the residual-row fixtures.
//! 3. An end-to-end `group_candidates` pass over a slice of the project's
//!    standard 5,000-molecule corpus, reporting group counts (not identity
//!    claims by itself -- see the module docs).
//!
//! Usage: `cargo run -p chematic-inchi --features native-inchi --example
//! dedup_verify -- /path/to/SMILES.csv`
//!
//! Every line is printed in a fixed, sorted order so the full stdout can be
//! diffed byte-for-byte across independent process invocations.

#![cfg(feature = "native-inchi")]

use chematic_inchi::dedup::{
    DedupRelation, IdentityPolicy, compare, compare_molecules, deduplicate_verified,
    group_candidates,
};
use chematic_smiles::{canonical_smiles, parse, random_smiles, random_smiles_vect};

const POLICIES: [IdentityPolicy; 4] = [
    IdentityPolicy::StandardInchiString,
    IdentityPolicy::StandardInchiKey,
    IdentityPolicy::StereoIgnored,
    IdentityPolicy::IsotopeIgnored,
];

fn mol(s: &str) -> chematic_core::Molecule {
    parse(s).unwrap_or_else(|e| panic!("parse {s:?}: {e}"))
}

fn report_pair(label: &str, a: &str, b: &str) {
    let ma = mol(a);
    let mb = mol(b);
    for policy in POLICIES {
        println!(
            "PAIR {label} policy={policy:?} -> {:?}",
            compare_molecules(&ma, &mb, policy)
        );
    }
}

fn main() {
    // --- Fixture pairs (fixed, sorted output order) -------------------------
    report_pair("atom_renumbered", "CC(=O)Oc1ccccc1C(=O)O", &{
        let a = mol("CC(=O)Oc1ccccc1C(=O)O");
        random_smiles(&a, 42)
    });
    report_pair("alternate_spelling", "CCO", "OCC");
    report_pair("ez", "C/C=C/C", "C/C=C\\C");
    report_pair("enantiomer", "N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O");
    report_pair(
        "diastereomer",
        "OC(=O)[C@H](O)[C@@H](O)C(=O)O",
        "OC(=O)[C@H](O)[C@H](O)C(=O)O",
    );
    report_pair("isotopologue_heavy", "CC", "[13CH3]C");
    report_pair("isotopologue_h", "C", "[2H]C([2H])([2H])[2H]");
    report_pair("protonation_state", "CC(=O)O", "CC(=O)[O-]");
    report_pair("tautomer", "O=c1cccc[nH]1", "Oc1ccccn1");
    report_pair("disconnected_salt", "CC(=O)[O-].[Na+]", "[Na+].CC(=O)[O-]");
    report_pair(
        "residual_relabel_only",
        "OC(=O)[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c1/c(c(c1O)O)=N/CCCCC",
        r"OC(=O)[C@H](Cc2ccc(NC(c3c(Cl)cncc3Cl)=O)cc2)/N=c\1c(/c(c1O)O)=N/CCCCC",
    );

    // Synthetic collision control (low-level `compare`, injected key).
    {
        let benzene = mol("c1ccccc1");
        let methanol = mol("CO");
        for policy in POLICIES {
            let rel = compare(
                "SYNTHETIC-SHARED-KEY",
                &benzene,
                "SYNTHETIC-SHARED-KEY",
                &methanol,
                policy,
            );
            println!("PAIR synthetic_collision policy={policy:?} -> {rel:?}");
        }
    }

    // --- Worst-of-30 permutation sweep on the residual-row fixture ----------
    let reference = mol("CC(=O)Oc1ccccc1C(=O)O");
    let respellings = random_smiles_vect(&reference, 30, 1000);
    let mut flips = 0usize;
    for s in &respellings {
        let variant = mol(s);
        let rel = compare_molecules(&reference, &variant, IdentityPolicy::StandardInchiString);
        if !matches!(
            rel,
            DedupRelation::VerifiedDuplicate | DedupRelation::CanonicalSplit
        ) {
            flips += 1;
        }
    }
    println!(
        "WORST_OF_30 respellings={} flips_to_distinct_or_collision={}",
        respellings.len(),
        flips
    );

    // --- Optional: corpus candidate-grouping + batch reconciliation pass ----
    if let Some(path) = std::env::args().nth(1) {
        let content = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let lines: Vec<&str> = content
            .lines()
            .skip(1)
            .filter(|l| !l.trim().is_empty())
            .collect();
        let n_lines = lines.len();
        // Parallel `mol_lines`/`mols`: `mols[i]` always corresponds to
        // `mol_lines[i]`'s SMILES text, kept in lockstep even though
        // `filter_map` skips parse failures (a plain `lines[idx]` lookup
        // below would silently misalign once any line fails to parse).
        let mut mol_lines: Vec<&str> = Vec::with_capacity(lines.len());
        let mols: Vec<chematic_core::Molecule> = lines
            .iter()
            .filter_map(|l| {
                let trimmed = l.trim();
                let m = parse(trimmed).ok();
                if m.is_some() {
                    mol_lines.push(trimmed);
                }
                m
            })
            .collect();
        let n_parse_failures = n_lines - mols.len();

        let groups = group_candidates(&mols);
        let singleton_groups = groups.iter().filter(|g| g.len() == 1).count();
        let multi_groups = groups.len() - singleton_groups;
        println!(
            "CORPUS n_lines={} n_parse_failures={} n_mols={} n_groups={} singleton_groups={} multi_member_groups={}",
            n_lines,
            n_parse_failures,
            mols.len(),
            groups.len(),
            singleton_groups,
            multi_groups
        );
        // Fixed-order digest: sort canonical keys, join, print length only
        // (not the full corpus text) so output stays short and diffable.
        let mut keys: Vec<String> = mols.iter().map(canonical_smiles).collect();
        keys.sort();
        let joined = keys.join("\n");
        println!("CORPUS canonical_keys_total_len={}", joined.len());

        // --- Batch reconciliation (`deduplicate_verified`), one pass per
        // policy -- reports the real pass/fail/unavailable counts over the
        // whole corpus, not just the fixture pairs above.
        for policy in POLICIES {
            let report = deduplicate_verified(&mols, policy);
            let n_verified =
                mols.len() - report.verification_unavailable.len() - report.invalid_molecules.len();
            let group_members: usize = report.groups.iter().map(|g| g.members.len()).sum();
            println!(
                "CORPUS_DEDUP policy={:?} n_verified={} n_verification_unavailable={} n_invalid_molecules={} n_groups={} group_members={} n_canonical_splits={} n_canonical_collisions={}",
                policy,
                n_verified,
                report.verification_unavailable.len(),
                report.invalid_molecules.len(),
                report.groups.len(),
                group_members,
                report.canonical_splits.len(),
                report.canonical_collisions.len(),
            );
            // Diagnostic: identify exactly which corpus molecule(s) failed
            // verification and why (re-running the public `standard_inchi`
            // directly, not reaching into `dedup`'s private `Verify` enum),
            // sorted for stable output.
            let mut unavailable = report.verification_unavailable.clone();
            unavailable.sort_unstable();
            for idx in unavailable {
                let smiles = mol_lines[idx];
                let reason = match chematic_inchi::standard_inchi(&mols[idx]) {
                    Ok(s) => format!("unexpectedly Ok({s:?}) on retry"),
                    Err(e) => format!("{e:?}"),
                };
                println!(
                    "CORPUS_DEDUP_UNAVAILABLE policy={policy:?} idx={idx} smiles={smiles:?} reason={reason}"
                );
            }
            // Diagnostic: print each verified group's member SMILES + their
            // canonical-SMILES keys, so a `CanonicalSplit` group found live
            // in the corpus (not just the fixed fixture pairs) is fully
            // inspectable from stdout.
            for group in &report.groups {
                let members: Vec<String> = group
                    .members
                    .iter()
                    .map(|&i| format!("{i}:{:?}", mol_lines[i]))
                    .collect();
                println!(
                    "CORPUS_DEDUP_GROUP policy={policy:?} members=[{}]",
                    members.join(", ")
                );
            }
        }
    }
}
