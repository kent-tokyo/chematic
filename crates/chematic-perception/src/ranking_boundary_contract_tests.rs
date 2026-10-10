use crate::rdkit_canon::*;
use serde_json::Value;

#[test]
fn raw_cis_trans_fragment_and_symbol_rankings_match_pinned_rdkit() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-2026.03.1-ranking-boundary.json"
    ))
    .unwrap();
    assert_eq!(fixture["rdkit_version"], "2026.03.1");
    let mut failures = Vec::new();
    for row in fixture["rows"].as_array().unwrap() {
        let atoms: Vec<_> = row["atoms"]
            .as_array()
            .unwrap()
            .iter()
            .map(|a| RdkitRankAtom {
                atomic_num: a["atomic_num"].as_u64().unwrap() as u32,
                isotope: a["isotope"].as_u64().unwrap() as u32,
                formal_charge: a["formal_charge"].as_i64().unwrap() as i32,
                atom_map: a["atom_map"].as_i64().unwrap() as i32,
                chiral_tag: a["chiral_tag"].as_u64().unwrap() as u8,
                total_num_hs: a["total_num_hs"].as_u64().unwrap() as u32,
                num_rings: a["num_rings"].as_u64().unwrap() as u32,
                ring_stereo: a["ring_stereo"].as_bool().unwrap(),
            })
            .collect();
        let bonds: Vec<_> = row["bonds"]
            .as_array()
            .unwrap()
            .iter()
            .map(|b| RdkitRankBond {
                begin: b["begin"].as_u64().unwrap() as u32,
                end: b["end"].as_u64().unwrap() as u32,
                bond_type: b["bond_type"].as_u64().unwrap() as u32,
                stereo: b["stereo"].as_u64().unwrap() as u32,
                stereo_atoms: b["stereo_atoms"]
                    .as_array()
                    .map(|s| (s[0].as_u64().unwrap() as u32, s[1].as_u64().unwrap() as u32)),
            })
            .collect();
        let mut actual = match row["kind"].as_str().unwrap() {
            "molecule" => rdkit_rank_mol_atoms(&atoms, &bonds),
            "fragment" => {
                let selected = |name: &str| {
                    row[name]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|x| x.as_bool().unwrap())
                        .collect::<Vec<_>>()
                };
                rdkit_rank_fragment_atoms(
                    &atoms,
                    &bonds,
                    &selected("atoms_in_play"),
                    &selected("bonds_in_play"),
                    row["include_chirality"].as_bool().unwrap(),
                )
            }
            "symbols" => {
                let symbols = row["atom_symbols"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|s| s.as_str().unwrap().to_string())
                    .collect::<Vec<_>>();
                rdkit_rank_fragment_atoms_with_symbols(
                    &atoms,
                    &bonds,
                    &symbols,
                    &vec![String::new(); bonds.len()],
                )
            }
            other => panic!("unknown ranking fixture {other}"),
        };
        if row["kind"] == "fragment" {
            // Python masks unused C++ ranks to -1. Only atoms in play are
            // part of the fragment ranking contract; mirror that output mask.
            for (rank, selected) in actual
                .iter_mut()
                .zip(row["atoms_in_play"].as_array().unwrap())
            {
                if !selected.as_bool().unwrap() {
                    *rank = u32::MAX;
                }
            }
        }
        let expected: Vec<u32> = serde_json::from_value(row["ranks"].clone()).unwrap();
        if actual != expected {
            failures.push(format!("{}: {actual:?} != {expected:?}", row["label"]));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
