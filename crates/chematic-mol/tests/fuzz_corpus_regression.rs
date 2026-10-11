//! Replay the checked-in fuzz corpora as deterministic regression tests.
//!
//! `cargo fuzz` remains the discovery lane. This test makes every minimized
//! input found there part of the ordinary workspace test and coverage lanes.

use std::fs;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::{Path, PathBuf};

fn corpus_files(name: &str) -> Vec<PathBuf> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fuzz/corpus")
        .join(name);
    let mut files: Vec<_> = fs::read_dir(&root)
        .unwrap_or_else(|error| panic!("read {}: {error}", root.display()))
        .map(|entry| entry.expect("read corpus entry").path())
        .filter(|path| path.is_file())
        .collect();
    files.sort();
    files
}

fn replay(name: &str, minimum_files: usize, mut run: impl FnMut(&[u8])) {
    let files = corpus_files(name);
    assert!(
        files.len() >= minimum_files,
        "{name} corpus unexpectedly shrank: {} < {minimum_files}",
        files.len()
    );
    for path in files {
        let data =
            fs::read(&path).unwrap_or_else(|error| panic!("read {}: {error}", path.display()));
        let result = catch_unwind(AssertUnwindSafe(|| run(&data)));
        assert!(result.is_ok(), "parser panicked for {}", path.display());
    }
}

#[test]
fn formats_corpus_is_panic_free() {
    replay("formats", 5_000, |data| {
        let Some((&selector, payload)) = data.split_first() else {
            return;
        };
        let Ok(input) = std::str::from_utf8(payload) else {
            return;
        };
        let input = &input[..input.len().min(1 << 20)];

        match selector % 22 {
            0 => {
                let _ = chematic_chem::parse_formula(input);
            }
            1 => {
                let _ = chematic_chem::parse_condensed(input);
            }
            2 => {
                let _ = chematic_smarts::parse_smarts(input);
            }
            3 => {
                let _ = chematic_smarts::parse_cxsmarts(input);
            }
            4 => {
                let _ = chematic_rxn::parse_reaction(input);
            }
            5 => {
                let _ = chematic_rxn::parse_reaction_query(input);
            }
            6 => {
                let _ = chematic_mol::parse_cml(input);
            }
            7 => {
                let _ = chematic_mol::parse_cif(input);
            }
            8 => {
                let _ = chematic_mol::parse_cjson(input);
            }
            9 => {
                let _ = chematic_mol::parse_cube(input);
            }
            10 => {
                let _ = chematic_mol::parse_ket(input);
            }
            11 => {
                let _ = chematic_mol::parse_pdbqt(input);
            }
            12 => {
                let _ = chematic_mol::parse_pqr(input);
            }
            13 => {
                let _ = chematic_mol::parse_gjf(input);
            }
            14 => {
                let _ = chematic_mol::parse_gaussian_log(input);
            }
            15 => {
                let _ = chematic_mol::parse_mmcif(input);
            }
            16 => {
                let _ = chematic_mol::parse_mol2(input);
            }
            17 => {
                let _ = chematic_mol::parse_mrv(input);
            }
            18 => {
                let _ = chematic_mol::parse_opendx(input);
            }
            19 => {
                let _ = chematic_mol::parse_orca_input(input);
            }
            20 => {
                let _ = chematic_mol::parse_orca_output(input);
            }
            21 => {
                let _ = chematic_3d::parse_pdb_atoms_with_limits(
                    input,
                    &chematic_3d::PdbParseLimits {
                        max_input_bytes: 1 << 20,
                        max_atoms: 10_000,
                        max_line_bytes: 4_096,
                        max_models: 256,
                    },
                );
            }
            _ => unreachable!(),
        }
    });
}

#[test]
fn moljson_corpus_is_panic_free() {
    let limits = chematic_mol::MolJsonParseLimits {
        max_input_bytes: 1 << 20,
        max_json_depth: 64,
        max_array_items: 1_024,
        max_string_bytes: 100_000,
        max_atoms: 10_000,
        max_bonds: 20_000,
    };
    replay("moljson", 1_100, |data| {
        let Ok(input) = std::str::from_utf8(data) else {
            return;
        };
        let _ = chematic_mol::parse_moljson_with_limits(input, &limits);
    });
}

#[test]
fn smiles_corpus_is_panic_free_and_canonical_output_reparses() {
    replay("smiles", 1_500, |data| {
        let Ok(input) = std::str::from_utf8(data) else {
            return;
        };
        if let Ok(molecule) = chematic_smiles::parse(input) {
            let canonical = chematic_smiles::canonical_smiles(&molecule);
            if !canonical.is_empty() {
                let _ = chematic_smiles::parse(&canonical);
            }
        }
    });
}

#[test]
fn xyz_corpus_is_panic_free() {
    let limits = chematic_3d::XyzParseLimits {
        max_input_bytes: 1 << 20,
        max_atoms: 10_000,
        max_line_bytes: 4_096,
    };
    replay("xyz", 500, |data| {
        let Ok(input) = std::str::from_utf8(data) else {
            return;
        };
        let _ = chematic_3d::parse_xyz_with_limits(input, limits);
    });
}
