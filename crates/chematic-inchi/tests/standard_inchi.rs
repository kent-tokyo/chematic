//! Integration tests comparing standard_inchi() output against known IUPAC reference values.
//!
//! Reference InChI strings were verified against the IUPAC InChI Trust reference tool
//! and the RDKit InChI implementation.

#![cfg(feature = "native-inchi")]

use chematic_inchi::{standard_inchi, standard_inchi_key};
use chematic_smiles::parse;

// (SMILES, expected_standard_InChI)
//
// For symmetric ring systems (benzene, cyclohexane), the connectivity layer ordering
// depends on the atom ordering provided to GetStdINCHI. Different InChI string forms
// for the same symmetric molecule (e.g. c1-2-3-4-5-6-1 vs c1-2-4-6-5-3-1) represent
// the same compound and produce the same InChIKey. We test InChIKeys directly below
// to verify molecular identity; the InChI strings for symmetric rings use the form
// produced by our GetStdINCHI call.
const REFERENCE: &[(&str, &str)] = &[
    ("C", "InChI=1S/CH4/h1H4"),
    ("CC", "InChI=1S/C2H6/c1-2/h1-2H3"),
    ("CCO", "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"),
    ("CC(=O)O", "InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)"),
    ("OO", "InChI=1S/H2O2/c1-2/h1-2H"),
    ("N", "InChI=1S/H3N/h1H3"),
    ("O", "InChI=1S/H2O/h1H2"),
    ("[NH4+]", "InChI=1S/H3N/h1H3/p+1"),
    ("[Na+].[Cl-]", "InChI=1S/ClH.Na/h1H;/q;+1/p-1"),
    ("CC(C)C", "InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3"),
    ("C#N", "InChI=1S/CHN/c1-2/h1H"),
    ("C=C", "InChI=1S/C2H4/c1-2/h1-2H2"),
    // --- Issue #11 reporter's specific failing examples ---
    // The pure-Rust implementation produced wrong output for all of these;
    // the native InChI C library (1.07.5) fixes all six.
    // Note: [H][H] (all-explicit-H molecule) is excluded — the C library
    // requires at least one non-H atom and returns null for H2-only input.
    ("N#N", "InChI=1S/N2/c1-2"),
    ("O=C=O", "InChI=1S/CO2/c2-1-3"),
    ("S=C=S", "InChI=1S/CS2/c2-1-3"),
    ("[OH-]", "InChI=1S/H2O/h1H2/p-1"),
    ("[C-]#[O+]", "InChI=1S/CO/c1-2"),
    ("O=O", "InChI=1S/O2/c1-2"),
];

#[test]
fn standard_inchi_matches_iupac_reference() {
    let mut failures = Vec::new();
    for &(smiles, expected) in REFERENCE {
        let mol = parse(smiles).unwrap_or_else(|e| panic!("parse {smiles:?}: {e}"));
        let got =
            standard_inchi(&mol).unwrap_or_else(|e| panic!("standard_inchi({smiles:?}): {e}"));
        if got != expected {
            failures.push(format!(
                "SMILES {smiles:?}\n  got:      {got}\n  expected: {expected}"
            ));
        }
    }
    if !failures.is_empty() {
        panic!("InChI mismatches:\n{}", failures.join("\n\n"));
    }
}

#[test]
fn standard_inchi_key_format() {
    let mol = parse("c1ccccc1").unwrap();
    let inchi = standard_inchi(&mol).unwrap();
    let key = standard_inchi_key(&inchi).unwrap();
    assert_eq!(key.len(), 27, "InChIKey length must be 27: {key}");
    assert_eq!(&key[14..15], "-", "first dash at position 14: {key}");
    assert_eq!(&key[25..26], "-", "second dash at position 25: {key}");
}

// InChIKey tests verify molecular identity regardless of connectivity string form.
// Reference values from PubChem.
#[test]
fn inchikey_benzene() {
    check_key("c1ccccc1", "UHOVQNZJYSORNB-UHFFFAOYSA-N");
}
#[test]
fn inchikey_cyclohexane() {
    check_key("C1CCCCC1", "XDTMQSROBMDMFD-UHFFFAOYSA-N");
}
#[test]
fn inchikey_pyridine() {
    check_key("c1ccncc1", "JUJWROOIHBZHMG-UHFFFAOYSA-N");
}
#[test]
fn inchikey_ethanol() {
    check_key("CCO", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N");
}
#[test]
fn inchikey_acetic_acid() {
    check_key("CC(=O)O", "QTBSBXVTEAMEQO-UHFFFAOYSA-N");
}
#[test]
fn inchikey_aspirin() {
    check_key("CC(=O)Oc1ccccc1C(=O)O", "BSYNRYMUTXBXSQ-UHFFFAOYSA-N");
}
// Caffeine (PubChem CID 2519): InChIKey RYYVLZVUVIJVGH-UHFFFAOYSA-N
#[test]
fn inchikey_caffeine() {
    check_key("Cn1cnc2c1c(=O)n(c(=O)n2C)C", "RYYVLZVUVIJVGH-UHFFFAOYSA-N");
}

// Tartaric acid enantiomers — tests multi-centre stereo discrimination.
//
// For the symmetric tartaric acid backbone, opposite SMILES annotations at C2 and C3
// encode the two enantiomers. Keys are from InChI 1.07.5 (/m1 for L-, /m0 for D-).
// PubChem's keys (computed with InChI ≤1.06) differ for D-tartaric due to a /m-layer
// reassignment between InChI 1.06 and 1.07.
//
// (2R,3R) L-tartaric acid
#[test]
fn inchikey_tartaric_acid_rr() {
    check_key(
        "OC(=O)[C@H](O)[C@@H](O)C(=O)O",
        "FEWJPZIEWOKRBE-JCYAYHJZSA-N",
    );
}
// (2S,3S) D-tartaric acid — enantiomer of R,R (stereo block differs)
#[test]
fn inchikey_tartaric_acid_ss() {
    check_key(
        "OC(=O)[C@@H](O)[C@H](O)C(=O)O",
        "FEWJPZIEWOKRBE-LWMBPPNESA-N",
    );
}

// L-Threonine: PubChem CID 6288, 2 chiral centers.
#[test]
fn inchikey_l_threonine() {
    check_key("C[C@@H](O)[C@H](N)C(=O)O", "AYFVYJQAPQTCCC-GBXIJSLDSA-N");
}

// L-alanine: PubChem CID 5950, InChIKey QNAYBMKLOCPYGJ-REOHCLBHSA-N
#[test]
fn inchikey_l_alanine() {
    check_key("N[C@@H](C)C(=O)O", "QNAYBMKLOCPYGJ-REOHCLBHSA-N");
}

// D-alanine: enantiomer, different stereo block
#[test]
fn inchikey_d_alanine() {
    check_key("N[C@H](C)C(=O)O", "QNAYBMKLOCPYGJ-UWTATZPHSA-N");
}

// (E)-but-2-ene and (Z)-but-2-ene differ only in the /b stereo layer.
// InChIKeys sourced from PubChem: CID 12112 (E) and CID 12113 (Z).
#[test]
fn inchikey_e_but2ene() {
    check_key("C/C=C/C", "IAQRGUVFOMOMEM-ONEGZZNKSA-N");
}
#[test]
fn inchikey_z_but2ene() {
    check_key("C/C=C\\C", "IAQRGUVFOMOMEM-ARJAWSKDSA-N");
}

fn check_key(smiles: &str, expected_key: &str) {
    let mol = parse(smiles).unwrap_or_else(|e| panic!("parse {smiles:?}: {e}"));
    let inchi = standard_inchi(&mol).unwrap_or_else(|e| panic!("standard_inchi({smiles:?}): {e}"));
    let key = standard_inchi_key(&inchi).unwrap_or_else(|e| panic!("standard_inchi_key: {e}"));
    assert_eq!(key, expected_key, "SMILES={smiles:?}, InChI={inchi}");
}

/// Smoke test against the 5000-molecule corpus from issue #11.
///
/// Checks that every parseable molecule generates an InChI without panicking or
/// returning an error. Not run in CI (requires the file on disk); run locally with:
///
/// ```
/// SMILES_CSV=/path/to/SMILES.csv \
///   cargo test -p chematic-inchi --features native-inchi corpus_smoke -- --ignored --nocapture
/// ```
#[test]
#[ignore]
fn corpus_smoke_issue_11() {
    let path = match std::env::var("SMILES_CSV") {
        Ok(p) => p,
        Err(_) => {
            eprintln!("SMILES_CSV env var not set — skipping corpus smoke test");
            return;
        }
    };
    let content =
        std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"));

    let mut ok = 0usize;
    let mut kekule_fail = 0usize; // pre-existing kekulization limit (complex/unusual rings)
    let mut lib_err = 0usize; // C library errors (these would be regressions)
    let mut parse_fail = 0usize;

    for line in content.lines().skip(1) {
        let smiles = line.trim();
        if smiles.is_empty() {
            continue;
        }
        match chematic_smiles::parse(smiles) {
            Err(_) => {
                parse_fail += 1;
            }
            Ok(mol) => match standard_inchi(&mol) {
                Ok(_) => {
                    ok += 1;
                }
                Err(chematic_inchi::InchiError::KekulizationFailed(_)) => {
                    kekule_fail += 1;
                }
                Err(chematic_inchi::InchiError::NullOutput) => {
                    lib_err += 1;
                    eprintln!("NULL-OUTPUT {smiles:?}");
                }
                Err(chematic_inchi::InchiError::InvalidInput(_)) => {
                    // Known: pure-H molecules ([H][H]) have no heavy atoms.
                    kekule_fail += 1;
                }
                Err(e) => {
                    lib_err += 1;
                    eprintln!("LIB-ERROR {smiles:?}: {e}");
                }
            },
        }
    }
    let total = ok + kekule_fail + lib_err + parse_fail;
    eprintln!(
        "Corpus ({total} molecules): {ok} ok, {kekule_fail} kekulization-skipped, \
         {lib_err} lib-errors, {parse_fail} parse-fails"
    );
    // Kekulization failures for complex/unusual ring systems are a pre-existing
    // limitation unrelated to issue #11. Only C-library errors are regressions.
    assert_eq!(
        lib_err, 0,
        "{lib_err} unexpected C-library errors (see LIB-ERROR lines above)"
    );
}
