//! `rdkit_mol_from_inchi` against RDKit 2026.03.1:
//! `Chem.MolToSmiles(Chem.MolFromInchi(inchi))`.

#![cfg(feature = "native-inchi")]

use chematic_inchi::rdkit_mol_from_inchi;
use chematic_smiles::rdkit_canonical_smiles;

const CASES: &[(&str, &str)] = &[
    ("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", "CCO"),
    (
        "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1",
        "C[C@H](N)C(=O)O",
    ),
    (
        "InChI=1S/C4H6O2/c1-2-3-4(5)6/h2-3H,1H3,(H,5,6)/b3-2+",
        "C/C=C/C(=O)O",
    ),
    ("InChI=1S/C5H5NO/c7-6-4-2-1-3-5-6/h1-5H", "[O-][n+]1ccccc1"),
    (
        "InChI=1S/C2H4O2.Na/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1",
        "CC(=O)[O-].[Na+]",
    ),
    ("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i1+1", "[13CH3]CO"),
    (
        "InChI=1S/C5H11NO2/c1-6(2,3)4-5(7)8/h4H2,1-3H3/p+1",
        "C[N+](C)(C)CC(=O)O",
    ),
    (
        "InChI=1S/C8H7N/c1-2-4-8-7(3-1)5-6-9-8/h1-6,9H",
        "c1ccc2[nH]ccc2c1",
    ),
    // Pentavalent nitro N repaired by the InChI clean-up rules.
    ("InChI=1S/CH3NO2/c1-2(3)4/h1H3", "C[N+](=O)[O-]"),
    ("InChI=1S/CH4O3S/c1-5(2,3)4/h1H3,(H,2,3,4)", "CS(=O)(=O)O"),
    ("InChI=1S/CH4O/c1-2/h2H,1H3/i1D3", "[2H]C([2H])([2H])O"),
];

#[test]
fn inchi_read_matches_rdkit_mol_from_inchi() {
    for &(inchi, expected) in CASES {
        let mol = rdkit_mol_from_inchi(inchi).unwrap_or_else(|e| panic!("{inchi}: {e}"));
        let got = rdkit_canonical_smiles(&mol).unwrap_or_else(|e| panic!("{inchi}: {e:?}"));
        assert_eq!(got, expected, "{inchi}");
    }
}

#[test]
fn invalid_inchi_is_an_error() {
    assert!(rdkit_mol_from_inchi("InChI=1S/garbage").is_err());
    assert!(rdkit_mol_from_inchi("not an inchi").is_err());
}

#[test]
fn official_inchi_reader_boundary_regressions_match_pinned_rdkit() {
    let mut mismatches = Vec::new();
    let mut count = 0;
    for line in
        include_str!("../../../validation/rdkit-2026.03.1-inchi-reader-boundary.tsv").lines()
    {
        let (inchi, expected) = line.split_once('\t').unwrap();
        let result = rdkit_mol_from_inchi(inchi)
            .map_err(|e| e.to_string())
            .and_then(|m| rdkit_canonical_smiles(&m).map_err(|e| e.to_string()));
        let agrees = match &result {
            Ok(actual) => actual == expected,
            Err(_) => expected == "<NONE>",
        };
        if !agrees {
            mismatches.push(format!("{inchi}: expected {expected}, got {result:?}"));
        }
        count += 1;
    }
    assert_eq!(count, 37);
    assert!(mismatches.is_empty(), "{}", mismatches.join("\n"));
}
