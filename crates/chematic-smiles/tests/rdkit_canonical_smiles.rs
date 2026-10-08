//! `rdkit_canonical_smiles` against RDKit 2026.03.1 on a fixture corpus of
//! hand-picked edge cases (stereo, rings, charges, isotopes, metals,
//! fragments, explicit hydrogens) and RDKit-randomized / Kekulé spellings of
//! them. `<NONE>` marks inputs RDKit's `MolFromSmiles` rejects; those must
//! be errors.

const FIXTURES: &str = include_str!("data/rdkit_canonical_smiles.tsv");

#[test]
fn matches_rdkit_on_fixture_corpus() {
    let mut failures = Vec::new();
    let mut n = 0;
    for line in FIXTURES.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let (input, expected) = line.split_once('\t').expect("tab-separated");
        n += 1;
        let got = chematic_smiles::parse(input)
            .map_err(|e| e.to_string())
            .and_then(|m| chematic_smiles::rdkit_canonical_smiles(&m).map_err(|e| e.to_string()));
        let ok = match (&got, expected) {
            (Err(_), "<NONE>") => true,
            (Ok(s), e) => s == e,
            _ => false,
        };
        if !ok {
            failures.push(format!("{input}: expected {expected}, got {got:?}"));
        }
    }
    assert!(n > 400, "fixture corpus too small: {n}");
    assert!(
        failures.is_empty(),
        "{} failures:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
