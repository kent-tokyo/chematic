//! Coordination-chemistry compatibility benchmark for anticancer platinum
//! complexes -- see `validation/platinum/FEASIBILITY.md` for the full
//! writeup. Reads `validation/platinum/pt_corpus.jsonl` (hand-authored
//! corpus with documented provenance; see the corpus file itself and the
//! FEASIBILITY doc for why this is not simply re-derived from PubChem's own
//! SMILES/InChI) and exercises: SMILES parse, formula/mass/charge, canonical
//! SMILES identity (the cisplatin/transplatin killer-benchmark check), MOL
//! V3000 round-trip, `validate_valence`, and panic/determinism smoke checks
//! for ECFP4 and a simple Pt-ligand SMARTS match.
//!
//! Writes one JSONL row per corpus entry to the path given as the first CLI
//! argument (defaults to stdout only).
//!
//! Run: `cargo run --release -p chematic-mol --example platinum_benchmark -- validation/results/platinum_baseline_chematic.jsonl`

use chematic_core::{AtomIdx, Molecule};
use serde_json::{Value, json};
use std::collections::BTreeMap;
use std::fs;
use std::panic::{self, AssertUnwindSafe};

/// Element-count multiset (atomic_number -> count), implicit H included.
/// Deliberately data (not a Hill-notation string) -- comparison against
/// `formula_expected` (parsed the same way) doesn't need string formatting.
fn formula_counts(mol: &Molecule) -> BTreeMap<u8, u32> {
    let mut counts: BTreeMap<u8, u32> = BTreeMap::new();
    for (idx, atom) in mol.atoms() {
        if atom.wildcard {
            continue;
        }
        let an = atom.element.atomic_number();
        if an != 1 {
            *counts.entry(an).or_insert(0) += 1;
            let h = chematic_core::valence::implicit_hcount(mol, idx) as u32;
            if h > 0 {
                *counts.entry(1).or_insert(0) += h;
            }
        } else {
            *counts.entry(1).or_insert(0) += 1;
        }
    }
    counts
}

/// Parse a Hill-notation-ish "expected" formula string (e.g. "Cl2H6N2Pt")
/// into the same (atomic_number -> count) shape as `formula_counts`, so the
/// two can be compared directly without needing chematic's own formula
/// writer (there isn't a public one outside chematic-wasm as of this
/// benchmark -- see FEASIBILITY.md).
fn parse_expected_formula(formula: &str) -> BTreeMap<u8, u32> {
    let counts = chematic_chem::formula::parse_formula(formula).unwrap_or_default();
    let mut out = BTreeMap::new();
    for (symbol, count) in counts {
        if let Some(el) = chematic_core::Element::from_symbol(&symbol) {
            out.insert(el.atomic_number(), count);
        }
    }
    out
}

fn net_charge(mol: &Molecule) -> i64 {
    mol.atoms().map(|(_, a)| a.charge as i64).sum()
}

fn pt_coordination_number(mol: &Molecule) -> Option<usize> {
    mol.atoms()
        .find(|(_, a)| a.element == chematic_core::Element::PT)
        .map(|(idx, _)| mol.degree(idx))
}

fn main() {
    let out_path = std::env::args().nth(1);
    let corpus_path = "validation/platinum/pt_corpus.jsonl";
    let content = fs::read_to_string(corpus_path)
        .unwrap_or_else(|e| panic!("failed to read {corpus_path}: {e}"));

    let mut rows: Vec<Value> = Vec::new();

    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let entry: Value = serde_json::from_str(line).expect("corpus row must be valid JSON");
        let id = entry["id"].as_str().unwrap().to_string();
        let smiles = entry["smiles_dative"].as_str().unwrap().to_string();
        let formula_expected = entry["formula_expected"].as_str().unwrap().to_string();
        let charge_expected = entry["charge_expected"].as_i64().unwrap();
        let coordination_expected = entry["pt_coordination_number"].as_i64();

        let mut row = json!({
            "id": id,
            "smiles_input": smiles,
            "formula_expected": formula_expected,
            "charge_expected": charge_expected,
        });

        // --- 4.1 Parsing ---
        let parsed = chematic_smiles::parse(&smiles);
        let mol = match parsed {
            Ok(m) => m,
            Err(e) => {
                row["parse_ok"] = json!(false);
                row["parse_error"] = json!(format!("{e:?}"));
                rows.push(row);
                continue;
            }
        };
        row["parse_ok"] = json!(true);
        row["atom_count"] = json!(mol.atom_count());
        row["bond_count"] = json!(mol.bond_count());

        // --- 4.3 Formula / charge / mass ---
        let counts = formula_counts(&mol);
        let expected_counts = parse_expected_formula(&formula_expected);
        row["formula_counts"] = json!(
            counts
                .iter()
                .map(|(k, v)| (k.to_string(), *v))
                .collect::<BTreeMap<_, _>>()
        );
        row["formula_matches_expected"] = json!(counts == expected_counts);
        row["net_charge"] = json!(net_charge(&mol));
        row["charge_matches_expected"] = json!(net_charge(&mol) == charge_expected);
        row["molecular_weight"] = json!(chematic_chem::descriptors::molecular_weight(&mol));
        row["exact_mass"] = json!(chematic_chem::descriptors::exact_mass(&mol));

        // --- 4.4 Coordination topology ---
        let cn = pt_coordination_number(&mol);
        row["pt_coordination_number_observed"] = json!(cn);
        row["pt_coordination_number_matches_expected"] = json!(
            coordination_expected
                .map(|e| Some(e as usize) == cn)
                .unwrap_or(true)
        );
        row["connected_components"] = json!(count_components(&mol));

        // --- 4.6 Canonicalization ---
        let canonical_result =
            panic::catch_unwind(AssertUnwindSafe(|| chematic_smiles::canonical_smiles(&mol)));
        match canonical_result {
            Ok(c) => row["canonical_smiles"] = json!(c),
            Err(_) => row["canonical_smiles_panicked"] = json!(true),
        }

        // --- validate_valence (opt-in sanitizer check; not on the default
        // parse path, but relevant to 4.1's "unsupported vs malformed") ---
        let valence_errors = chematic_core::validate_valence(&mol);
        row["valence_errors"] = json!(
            valence_errors
                .iter()
                .map(|e| format!("{e}"))
                .collect::<Vec<_>>()
        );

        // --- 4.2 Round-trip through MOL V3000 ---
        let coords: Vec<(f64, f64)> = vec![(0.0, 0.0); mol.atom_count()];
        let metadata = chematic_mol::MolMetadata::default();
        let mol_block = chematic_mol::mol3000::write_mol_v3000(&mol, &metadata, &coords);
        match chematic_mol::mol3000::parse_mol_v3000(&mol_block) {
            Ok((rt_mol, _meta)) => {
                let rt_counts = formula_counts(&rt_mol);
                row["mol_v3000_roundtrip_ok"] = json!(true);
                row["mol_v3000_roundtrip_formula_preserved"] = json!(rt_counts == counts);
                row["mol_v3000_roundtrip_charge_preserved"] =
                    json!(net_charge(&rt_mol) == net_charge(&mol));
                row["mol_v3000_roundtrip_coordination_preserved"] =
                    json!(pt_coordination_number(&rt_mol) == cn);
                let dative_before = mol
                    .bonds()
                    .filter(|(_, b)| b.order == chematic_core::BondOrder::Dative)
                    .count();
                let dative_after = rt_mol
                    .bonds()
                    .filter(|(_, b)| b.order == chematic_core::BondOrder::Dative)
                    .count();
                row["mol_v3000_roundtrip_dative_bonds_before"] = json!(dative_before);
                row["mol_v3000_roundtrip_dative_bonds_after"] = json!(dative_after);
                row["mol_v3000_roundtrip_dative_preserved"] = json!(dative_before == dative_after);
            }
            Err(e) => {
                row["mol_v3000_roundtrip_ok"] = json!(false);
                row["mol_v3000_roundtrip_error"] = json!(format!("{e:?}"));
            }
        }

        // --- 4.7 Fingerprint: panic + determinism smoke check only (per
        // FEASIBILITY.md scope decision -- bit-exact RDKit parity is not
        // required for metal-containing structures) ---
        let fp1 = panic::catch_unwind(AssertUnwindSafe(|| chematic_fp::ecfp4(&mol)));
        let fp2 = panic::catch_unwind(AssertUnwindSafe(|| chematic_fp::ecfp4(&mol)));
        match (fp1, fp2) {
            (Ok(a), Ok(b)) => {
                row["ecfp4_ok"] = json!(true);
                row["ecfp4_deterministic"] = json!(a == b);
            }
            _ => {
                row["ecfp4_ok"] = json!(false);
            }
        }

        // --- 4.8 Substructure: does a simple Pt-Cl SMARTS panic/find matches? ---
        let smarts_result = panic::catch_unwind(AssertUnwindSafe(|| {
            let query = chematic_smarts::parse_smarts("[#78]~[#17]").ok()?;
            Some(chematic_smarts::find_matches(&query, &mol).len())
        }));
        match smarts_result {
            Ok(Some(n)) => {
                row["smarts_pt_cl_ok"] = json!(true);
                row["smarts_pt_cl_match_count"] = json!(n);
            }
            Ok(None) => row["smarts_pt_cl_ok"] = json!(false),
            Err(_) => row["smarts_pt_cl_panicked"] = json!(true),
        }

        rows.push(row);
    }

    // --- 6. Killer benchmark: cisplatin/transplatin canonical identity ---
    let cisplatin_canon = rows
        .iter()
        .find(|r| r["id"] == "cisplatin")
        .and_then(|r| r["canonical_smiles"].as_str());
    let transplatin_canon = rows
        .iter()
        .find(|r| r["id"] == "transplatin")
        .and_then(|r| r["canonical_smiles"].as_str());
    let killer_benchmark = json!({
        "cisplatin_canonical": cisplatin_canon,
        "transplatin_canonical": transplatin_canon,
        "cis_trans_distinguished": cisplatin_canon.is_some()
            && transplatin_canon.is_some()
            && cisplatin_canon != transplatin_canon,
    });

    let output = json!({
        "rows": rows,
        "killer_benchmark_cisplatin_transplatin": killer_benchmark,
    });

    if let Some(path) = out_path {
        let mut text = String::new();
        for row in &rows {
            text.push_str(&row.to_string());
            text.push('\n');
        }
        fs::write(&path, &text).unwrap_or_else(|e| panic!("failed to write {path}: {e}"));
        let summary_path = path.replace(".jsonl", "_summary.json");
        fs::write(
            &summary_path,
            serde_json::to_string_pretty(&killer_benchmark).unwrap(),
        )
        .unwrap_or_else(|e| panic!("failed to write {summary_path}: {e}"));
        eprintln!("wrote {} rows to {path}", rows.len());
        eprintln!("killer benchmark summary written to {summary_path}");
    } else {
        println!("{}", serde_json::to_string_pretty(&output).unwrap());
    }
}

fn count_components(mol: &Molecule) -> usize {
    let n = mol.atom_count();
    if n == 0 {
        return 0;
    }
    let mut seen = vec![false; n];
    let mut components = 0;
    for start in 0..n {
        if seen[start] {
            continue;
        }
        components += 1;
        let mut stack = vec![AtomIdx(start as u32)];
        while let Some(cur) = stack.pop() {
            if seen[cur.0 as usize] {
                continue;
            }
            seen[cur.0 as usize] = true;
            for (nb, _) in mol.neighbors(cur) {
                if !seen[nb.0 as usize] {
                    stack.push(nb);
                }
            }
        }
    }
    components
}
