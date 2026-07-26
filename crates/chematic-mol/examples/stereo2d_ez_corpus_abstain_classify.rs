//! Reconciles the "RDKit-resolved but chematic-abstained" bond population
//! (the 346 bonds in `bond_chematic_abstained`, from
//! `scripts/stereo2d_ez_corpus_diagnosis.py`'s InChI `/b`-layer per-bond
//! comparison over the standard 4,999-molecule corpus) against chematic's
//! own internal per-bond classification, so the abstention reasons can be
//! reported with a real denominator instead of assumed.
//!
//! Reads `validation/results/stereo2d_ez_corpus_diagnosis_summary.json`'s
//! `chematic_abstained_detail` list: one entry per abstained bond, each
//! carrying the original RDKit-generated `mol_block` and a native (0-based,
//! chematic `AtomIdx`-order -- chematic preserves the MOL file's atom order
//! exactly, so no separate mapping is needed on this side) atom-index pair,
//! reverse-mapped from InChI's own canonical numbering via the AuxInfo `/N:`
//! layer (see that script's `inchi_b_layer_and_auxinfo_n`/
//! `native_pair_from_inchi_key` and its own `_self_test_auxinfo_mapping`,
//! which verified the mapping recipe on two independent hand-built fixtures
//! before trusting it at corpus scale).
//!
//! For each entry, re-parses the SAME `mol_block` through the real reader
//! (`chematic_mol::read_mol_with_diagnostics`), looks up the bond via
//! `Molecule::bond_between`, and classifies it:
//!   - present in `report.ez_diagnostics` -> tally by its
//!     `EzDirectionRejectionReason` variant.
//!   - not present at all (silently `NotRequested`) -> sub-classify by
//!     re-deriving the same two conditions
//!     `stereo2d_ez_direction::{classify_double_bond, resolve_end}` check
//!     before ever reaching geometry: is the bond itself aromatic (per a
//!     fresh, non-mutating `assign_aromaticity` query -- the production
//!     module's own mechanism, not a new one), or do the two substituents
//!     at one end compare topologically equal (per
//!     `cip_priority::compare_branches`, the same oracle `resolve_end` uses
//!     for its `NonStereogenic` check). Anything not positively identified
//!     goes into `other_unclassified` with its own SMILES/bond recorded --
//!     never silently folded into another bucket to make the sum look
//!     cleaner.
//!
//! Run:
//! ```text
//! cargo run -p chematic-mol --example stereo2d_ez_corpus_abstain_classify --release
//! ```

use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs;

use chematic_core::{AtomIdx, BondOrder, Molecule};
use chematic_perception::assign_aromaticity;
use chematic_perception::cip_priority::compare_branches;
use serde_json::{Value, json};

/// Non-double-bond neighbors of `end`, excluding `other_end` -- mirrors
/// `stereo2d_ez_direction::substituents` exactly (that helper is private to
/// its module; reimplemented here rather than exposed as `pub`, since this
/// is validation tooling, not a change to the production algorithm).
fn substituents(mol: &Molecule, end: AtomIdx, other_end: AtomIdx) -> Vec<AtomIdx> {
    mol.neighbors(end)
        .filter(|&(nb, bidx)| nb != other_end && mol.bond(bidx).order != BondOrder::Double)
        .map(|(nb, _)| nb)
        .collect()
}

fn main() {
    let summary_path = "validation/results/stereo2d_ez_corpus_diagnosis_summary.json";
    let text = fs::read_to_string(summary_path).unwrap_or_else(|e| {
        panic!(
            "failed to read {summary_path}: {e} -- run \
             scripts/stereo2d_ez_corpus_diagnosis.py first"
        )
    });
    let root: Value = serde_json::from_str(&text).expect("parse summary json");
    let details = root["chematic_abstained_detail"]
        .as_array()
        .expect("chematic_abstained_detail missing or not an array -- regenerate the summary with the updated script");

    let total = details.len();
    let mut reason_counts: HashMap<String, usize> = HashMap::new();
    let mut reason_samples: HashMap<String, Vec<String>> = HashMap::new();
    let mut not_requested_aromatic = 0usize;
    let mut not_requested_equivalent_substituents = 0usize;
    let mut lost_in_canonicalization = 0usize;
    let mut lost_in_canonicalization_detail: Vec<Value> = Vec::new();
    let mut other_unclassified: Vec<Value> = Vec::new();

    for entry in details {
        let mol_block = entry["mol_block"].as_str().expect("mol_block field");
        let pair = entry["native_atom_pair"]
            .as_array()
            .expect("native_atom_pair field (should never be null -- the diagnosis script gates on 0 unmapped)");
        let a = AtomIdx(pair[0].as_u64().expect("native_atom_pair[0]") as u32);
        let b = AtomIdx(pair[1].as_u64().expect("native_atom_pair[1]") as u32);

        let report = match chematic_mol::read_mol_with_diagnostics(mol_block) {
            Ok(r) => r,
            Err(e) => {
                other_unclassified.push(json!({
                    "reason": "reparse_failed",
                    "error": format!("{e:?}"),
                    "smiles": entry["smiles"],
                }));
                continue;
            }
        };
        let mol = &report.mol;

        let Some((bond_idx, bond)) = mol.bond_between(a, b) else {
            other_unclassified.push(json!({
                "reason": "bond_not_found_after_reparse",
                "smiles": entry["smiles"],
                "native_atom_pair": entry["native_atom_pair"],
            }));
            continue;
        };
        if bond.order != BondOrder::Double {
            other_unclassified.push(json!({
                "reason": "resolved_bond_is_not_double",
                "smiles": entry["smiles"],
                "native_atom_pair": entry["native_atom_pair"],
                "actual_order": format!("{:?}", bond.order),
            }));
            continue;
        }

        if let Some(diag) = report.ez_diagnostics.iter().find(|d| d.bond == bond_idx) {
            let key = format!("{:?}", diag.reason);
            *reason_counts.entry(key.clone()).or_insert(0) += 1;
            let samples = reason_samples.entry(key).or_default();
            if samples.len() < 3 {
                samples.push(entry["smiles"].as_str().unwrap_or("").to_string());
            }
            continue;
        }

        // Silently NotRequested (no diagnostic at all): sub-classify by the
        // same two structural checks the production module itself applies
        // before ever touching geometry.
        let aromaticity = assign_aromaticity(mol);
        if aromaticity.is_bond_aromatic(bond_idx) {
            not_requested_aromatic += 1;
            continue;
        }

        let a1 = bond.atom1;
        let a2 = bond.atom2;
        let subs_a1 = substituents(mol, a1, a2);
        let subs_a2 = substituents(mol, a2, a1);
        let equivalent_at_a1 = subs_a1.len() == 2
            && compare_branches(mol, a1, subs_a1[0], subs_a1[1]) == Ordering::Equal;
        let equivalent_at_a2 = subs_a2.len() == 2
            && compare_branches(mol, a2, subs_a2[0], subs_a2[1]) == Ordering::Equal;
        if equivalent_at_a1 || equivalent_at_a2 {
            not_requested_equivalent_substituents += 1;
            continue;
        }

        // Did `stereo2d_ez_direction` actually assign a carrier direction
        // for THIS double bond (on either end's substituent bond), even
        // though no `ez_diagnostics` entry exists? If so, this bond was
        // never abstained by the perception module at all -- any loss is
        // happening downstream, in `chematic_smiles::canonical_smiles`'s
        // own carrier/grouping logic (verified empirically, not assumed:
        // found via individual inspection of a sample of this exact shape --
        // `crates/chematic-perception/src/stereo2d_ez_direction.rs` itself
        // is untouched by this discovery, and this is NOT fixed here, per
        // this addendum's explicit scope). A carrier bond is any of the
        // double bond's own (<=2) substituent bonds at either end.
        let carrier_candidates: Vec<_> = subs_a1
            .iter()
            .map(|&s| mol.bond_between(a1, s).unwrap().0)
            .chain(subs_a2.iter().map(|&s| mol.bond_between(a2, s).unwrap().0))
            .collect();
        let assigned_carrier = carrier_candidates
            .iter()
            .any(|&bidx| mol.bond_direction(bidx).is_some());
        if assigned_carrier {
            let write_out = chematic_smiles::write(mol);
            let canonical_out = chematic_smiles::canonical_smiles(mol);
            let write_has_token = write_out.contains('/') || write_out.contains('\\');
            lost_in_canonicalization += 1;
            lost_in_canonicalization_detail.push(json!({
                "smiles": entry["smiles"],
                "native_atom_pair": entry["native_atom_pair"],
                "write_has_directional_token_somewhere": write_has_token,
                "canonical_smiles": canonical_out,
            }));
            continue;
        }

        other_unclassified.push(json!({
            "reason": "unclassified_not_requested",
            "smiles": entry["smiles"],
            "native_atom_pair": entry["native_atom_pair"],
            "subs_a1_count": subs_a1.len(),
            "subs_a2_count": subs_a2.len(),
        }));
    }

    let rejected_total: usize = reason_counts.values().sum();
    let not_requested_total = not_requested_aromatic + not_requested_equivalent_substituents;
    let sum =
        rejected_total + not_requested_total + lost_in_canonicalization + other_unclassified.len();

    println!(
        "RDKit-resolved but chematic-abstained (population: chematic_abstained_detail from \
         scripts/stereo2d_ez_corpus_diagnosis.py's bond_chematic_abstained, i.e. rv real sign, \
         cv is None/\"?\"): {total}"
    );
    let mut sorted_reasons: Vec<(&String, &usize)> = reason_counts.iter().collect();
    sorted_reasons.sort_by_key(|(a, _)| *a);
    for (reason, count) in &sorted_reasons {
        let empty = Vec::new();
        let samples = reason_samples.get(*reason).unwrap_or(&empty);
        println!("  Rejected({reason}): {count}  (e.g. {samples:?})");
    }
    for reason in [
        "CarrierConflict",
        "MissingCoordinate",
        "NonFiniteCoordinate",
        "DegenerateGeometry",
        "ExplicitlyUnspecified",
        "UnsupportedTopology",
    ] {
        if !reason_counts.contains_key(reason) {
            println!("  Rejected({reason}): 0");
        }
    }
    println!(
        "  NotRequested (sub-classified where possible): {not_requested_total}\n\
        \x20   aromatic Kekule bond: {not_requested_aromatic}\n\
        \x20   topologically-equivalent substituents: {not_requested_equivalent_substituents}"
    );
    println!(
        "  Assigned by stereo2d_ez_direction but lost in canonical_smiles (NOT a \
         stereo2d_ez_direction rejection -- see doc comment above; not fixed in this \
         addendum, production algorithm code unchanged): {lost_in_canonicalization}"
    );
    println!("  Other/unclassified: {}", other_unclassified.len());
    println!("  sum: {sum} (population: {total})");

    let result = json!({
        "population": total,
        "rejected_by_reason": reason_counts,
        "rejected_by_reason_samples": reason_samples,
        "not_requested_total": not_requested_total,
        "not_requested_aromatic": not_requested_aromatic,
        "not_requested_equivalent_substituents": not_requested_equivalent_substituents,
        "lost_in_canonicalization_count": lost_in_canonicalization,
        "lost_in_canonicalization_detail": lost_in_canonicalization_detail,
        "other_unclassified_count": other_unclassified.len(),
        "other_unclassified": other_unclassified,
        "sum": sum,
        "sum_matches_population": sum == total,
    });
    let out_path = "validation/results/stereo2d_ez_corpus_abstain_classify_summary.json";
    fs::write(out_path, serde_json::to_string_pretty(&result).unwrap())
        .expect("write classify summary");
    println!("wrote {out_path}");

    if sum != total {
        eprintln!(
            "FATAL: sum ({sum}) != population ({total}) -- the reconciliation does not close, \
             investigate before trusting this table"
        );
        std::process::exit(1);
    }
    if !other_unclassified.is_empty() {
        eprintln!(
            "NOTE: {} bond(s) landed in Other/unclassified -- each one is listed individually \
             in {out_path}, not just counted.",
            other_unclassified.len()
        );
    }
}
