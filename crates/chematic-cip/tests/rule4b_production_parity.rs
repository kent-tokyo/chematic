//! Milestone 4B-2: production-path parity for the ported Rule 4b engine
//! (`src/auxiliary.rs`/`src/rule4b.rs`/`src/resolver.rs`, wired into
//! `assign_cip_accurate_experimental` via `assign.rs::assign_cip_accurate_experimental`).
//!
//! Reruns the exact 4 oracle corpora (+ mirrors) the diagnostic reference engine
//! (`examples/rule4b_bottom_up.rs`) validated 72/72 against, but through the real
//! production entry point (`assign_cip_accurate_experimental`, not a hand-built
//! digraph) -- this is the parity gate: if the mechanical port introduced any
//! divergence from the validated reference, it surfaces here, not just in the
//! full-corpus report.
//!
//! `ROWS_ORIGINAL` (quinic acid / galloyl ester molecules) contain aromatic gallic-acid
//! rings, so running them through `assign_cip_accurate_experimental` (which always
//! attempts `prepare_kekule_form`/`CipDigraph::new_with_mancude` when it succeeds) is
//! this suite's MANCUDE/aromatic coverage for the new Rule 4b pass -- a code path the
//! diagnostic reference engine itself never exercised (it used plain `CipDigraph::new`
//! throughout). `aromatic_content_present_in_rows_original` checks this claim directly
//! rather than assuming it.
//!
//! A dedicated "Rule 5 didn't regress" test is deliberately not duplicated here: the
//! full-corpus gate (`cip_accurate_full_corpus_report.py`'s independent per-row oracle
//! recompute) already covers every corpus row, including the two known Rule-5
//! (lowercase `r`/`s`) rows, more thoroughly than a hardcoded 2-row unit test would.

use chematic_cip::{CipBudget, assign_cip_accurate_experimental};
use chematic_core::{AtomIdx, CipCode};

const VS196_SMILES: &str = "Cl[C@H]1[C@H]([C@H]([C@@H]([C@@H]([C@H]1Cl)Cl)Cl)Cl)Cl";
const VS196_ROWS: &[(u32, CipCode)] = &[
    (1, CipCode::R),
    (2, CipCode::S),
    (3, CipCode::S),
    (4, CipCode::R),
    (5, CipCode::S),
    (6, CipCode::R),
];

const VS197_SMILES: &str = "Cl[C@H]1[C@@H]([C@H]([C@@H]([C@@H]([C@H]1Cl)Cl)Cl)Cl)Cl";
const VS197_ROWS: &[(u32, CipCode)] = &[
    (1, CipCode::R),
    (2, CipCode::R),
    (3, CipCode::R),
    (4, CipCode::R),
    (5, CipCode::S),
    (6, CipCode::S),
];

const ROWS_ORIGINAL: &[(&str, u32, &str, CipCode)] = &[
    (
        "O=C(O[C@@H]1C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1OC(=O)c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1",
        5,
        "tri-galloyl-a",
        CipCode::S,
    ),
    (
        "O=C(O[C@@H]1C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1OC(=O)c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1",
        35,
        "tri-galloyl-b",
        CipCode::S,
    ),
    (
        "O=C(O[C@H]1[C@H](O)C[C@](O)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
        3,
        "quinic-a",
        CipCode::S,
    ),
    (
        "O=C(O[C@H]1[C@H](O)C[C@](O)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
        7,
        "quinic-b",
        CipCode::S,
    ),
    (
        "O=C(O[C@H]1[C@H](O)C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
        3,
        "mono-galloyl-a",
        CipCode::S,
    ),
    (
        "O=C(O[C@H]1[C@H](O)C[C@](OC(=O)c2cc(O)c(O)c(O)c2)(C(=O)O)C[C@H]1O)c1cc(O)c(O)c(O)c1",
        7,
        "mono-galloyl-b",
        CipCode::S,
    ),
    (
        "O=C(Oc1c(O)cc(C(=O)O[C@H]2[C@H](OC(=O)c3cc(O)c(O)c(O)c3)C[C@](O)(C(=O)O)C[C@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
        11,
        "tetra-galloyl-a",
        CipCode::S,
    ),
    (
        "O=C(Oc1c(O)cc(C(=O)O[C@H]2[C@H](OC(=O)c3cc(O)c(O)c(O)c3)C[C@](O)(C(=O)O)C[C@H]2OC(=O)c2cc(O)c(O)c(O)c2)cc1O)c1cc(O)c(O)c(O)c1",
        26,
        "tetra-galloyl-b",
        CipCode::S,
    ),
];

const ROWS_DISCRIMINATING: &[(&str, u32, &str, CipCode)] = &[
    (
        "C[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C",
        7,
        "indep-ref-1",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)C",
        7,
        "indep-ref-2",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@@H](O)C",
        7,
        "indep-ref-3",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C",
        7,
        "indep-ref-4",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)C",
        7,
        "indep-ref-5",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)C",
        7,
        "indep-ref-6",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)C",
        7,
        "indep-ref-7",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)C",
        7,
        "indep-ref-8",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)C",
        7,
        "deep-div-1",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)C",
        7,
        "deep-div-2",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)C",
        7,
        "deep-div-3",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)C",
        7,
        "deep-div-4",
        CipCode::S,
    ),
    (
        "C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C",
        7,
        "deep-div-5",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C",
        7,
        "deep-div-6",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)C",
        7,
        "deep-div-7",
        CipCode::R,
    ),
    (
        "C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C",
        7,
        "deep-div-8",
        CipCode::R,
    ),
];

fn mirror_smiles(smiles: &str) -> String {
    let mut out = String::with_capacity(smiles.len());
    let mut chars = smiles.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '@' {
            if chars.peek() == Some(&'@') {
                chars.next();
                out.push('@');
            } else {
                out.push_str("@@");
            }
        } else {
            out.push(c);
        }
    }
    out
}

fn mirror_code(c: CipCode) -> CipCode {
    match c {
        CipCode::R => CipCode::S,
        CipCode::S => CipCode::R,
        other => other,
    }
}

fn assert_case(smiles: &str, atom_idx: u32, case_id: &str, oracle: CipCode) {
    let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
    let budget = CipBudget::default_budget();
    let result = assign_cip_accurate_experimental(&mol, budget).expect("assignment succeeds");
    let got = result
        .assignments
        .iter()
        .find(|(idx, _)| *idx == AtomIdx(atom_idx))
        .map(|(_, code)| *code);
    let skip_reason = result
        .skipped
        .iter()
        .find(|(idx, _)| *idx == AtomIdx(atom_idx))
        .map(|(_, reason)| reason);
    assert_eq!(
        got,
        Some(oracle),
        "{case_id} (atom {atom_idx}): expected {oracle:?}, got {got:?} (skip_reason={skip_reason:?})"
    );
}

#[test]
fn vs196_production_parity() {
    for &(atom_idx, oracle) in VS196_ROWS {
        assert_case(VS196_SMILES, atom_idx, "VS196", oracle);
    }
}

#[test]
fn vs196_mirror_production_parity() {
    let mirrored = mirror_smiles(VS196_SMILES);
    for &(atom_idx, oracle) in VS196_ROWS {
        assert_case(&mirrored, atom_idx, "VS196-mirror", mirror_code(oracle));
    }
}

#[test]
fn vs197_production_parity() {
    for &(atom_idx, oracle) in VS197_ROWS {
        assert_case(VS197_SMILES, atom_idx, "VS197", oracle);
    }
}

#[test]
fn vs197_mirror_production_parity() {
    let mirrored = mirror_smiles(VS197_SMILES);
    for &(atom_idx, oracle) in VS197_ROWS {
        assert_case(&mirrored, atom_idx, "VS197-mirror", mirror_code(oracle));
    }
}

#[test]
fn rows_original_production_parity() {
    for &(smiles, atom_idx, case_id, oracle) in ROWS_ORIGINAL {
        assert_case(smiles, atom_idx, case_id, oracle);
    }
}

#[test]
fn rows_original_mirror_production_parity() {
    for &(smiles, atom_idx, case_id, oracle) in ROWS_ORIGINAL {
        let mirrored = mirror_smiles(smiles);
        assert_case(&mirrored, atom_idx, case_id, mirror_code(oracle));
    }
}

#[test]
fn rows_discriminating_production_parity() {
    for &(smiles, atom_idx, case_id, oracle) in ROWS_DISCRIMINATING {
        assert_case(smiles, atom_idx, case_id, oracle);
    }
}

#[test]
fn rows_discriminating_mirror_production_parity() {
    for &(smiles, atom_idx, case_id, oracle) in ROWS_DISCRIMINATING {
        let mirrored = mirror_smiles(smiles);
        assert_case(&mirrored, atom_idx, case_id, mirror_code(oracle));
    }
}

/// Confirms `ROWS_ORIGINAL`'s molecules genuinely contain aromatic atoms -- the basis
/// for this suite's claim that `rows_original_production_parity` exercises the
/// `CipDigraph::new_with_mancude` path, not just plain `CipDigraph::new`.
#[test]
fn aromatic_content_present_in_rows_original() {
    for &(smiles, _, case_id, _) in ROWS_ORIGINAL {
        let mol = chematic_smiles::parse(smiles).expect("valid SMILES");
        let has_aromatic = (0..mol.atom_count()).any(|i| mol.atom(AtomIdx(i as u32)).aromatic);
        assert!(
            has_aromatic,
            "{case_id}: expected an aromatic ring (gallic acid ester), found none -- \
             MANCUDE-path coverage claim would be false"
        );
    }
}
