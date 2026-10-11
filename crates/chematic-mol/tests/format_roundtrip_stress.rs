//! Exercise every general-purpose molecular writer on the checked-in stress
//! corpus.  Exact textual parity belongs to the format-specific fixtures;
//! this lane protects the public conversion boundary from panics and ensures
//! that emitted records remain readable.

use std::fs;
use std::path::Path;

#[test]
fn stress_corpus_writes_and_reparses_every_general_molecular_format() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../validation/cosmolkit_comparison/corpora/stress-v1.smi");
    let corpus = fs::read_to_string(&path).unwrap();
    let mut checked = 0usize;

    for (row, line) in corpus.lines().enumerate() {
        let smiles = line.split_whitespace().next().unwrap_or("");
        if smiles.is_empty() || smiles.starts_with('#') {
            continue;
        }
        let Ok(mol) = chematic_smiles::parse(smiles) else {
            continue;
        };
        let metadata = chematic_mol::MolMetadata::default();

        let v2000 = chematic_mol::write_mol(&mol, &metadata);
        assert!(chematic_mol::parse_mol(&v2000).is_ok(), "V2000 row {row}");

        let v3000 = chematic_mol::write_mol_v3000(&mol, &metadata, &[]);
        assert!(
            chematic_mol::parse_mol_v3000(&v3000).is_ok(),
            "V3000 row {row}"
        );

        let cml = chematic_mol::write_cml(&mol, None);
        assert!(chematic_mol::parse_cml(&cml).is_ok(), "CML row {row}");

        let cjson = chematic_mol::write_cjson(&mol, &[]);
        assert!(chematic_mol::parse_cjson(&cjson).is_ok(), "CJSON row {row}");

        let moljson = chematic_mol::write_moljson(&mol);
        assert!(
            chematic_mol::parse_moljson(&moljson).is_ok(),
            "MolJSON row {row}"
        );

        let cdxml = chematic_mol::write_cdxml(&mol, &[]);
        assert!(chematic_mol::parse_cdxml(&cdxml).is_ok(), "CDXML row {row}");

        let mol2 = chematic_mol::write_mol2(&mol, &[]);
        assert!(chematic_mol::parse_mol2(&mol2).is_ok(), "MOL2 row {row}");

        checked += 1;
    }
    assert!(
        checked >= 4_000,
        "stress corpus unexpectedly shrank: {checked}"
    );
}
