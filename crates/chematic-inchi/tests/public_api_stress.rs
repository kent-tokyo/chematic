//! Public InChI and identity APIs over unusual but parseable structures.

use chematic_inchi::dedup::{
    AtomMapPolicy, GraphStrictness, IdentityPolicy, IndexedGraphRelationMode,
    compare_indexed_graph_relation, compare_molecules,
    compare_molecules_with_accurate_cip_preflight, deduplicate_verified,
    deduplicate_verified_with_accurate_cip_preflight, group_candidates,
};
use chematic_inchi::{
    InchiParseLimits, build_inchi_index, inchi, inchi_key, parse_inchi, parse_inchi_with_limits,
    rdkit_mol_from_inchi, standard_inchi, standard_inchi_key,
};
use chematic_smiles::parse;

const STRESS: &str = include_str!("../../../validation/cosmolkit_comparison/corpora/stress-v1.smi");

#[test]
fn inchi_roundtrip_and_identity_surface_handles_special_chemistry() {
    let mut molecules = Vec::new();
    for line in STRESS.lines() {
        if molecules.len() == 384 {
            break;
        }
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        if mol.atom_count() == 0 || mol.atom_count() > 100 {
            continue;
        }

        let approx = inchi(&mol);
        let _ = build_inchi_index(&mol);
        let _ = inchi_key(&approx);
        let _ = parse_inchi(&approx);
        let _ = parse_inchi_with_limits(&approx, &InchiParseLimits::default());

        if let Ok(standard) = standard_inchi(&mol) {
            let _ = standard_inchi_key(&standard);
            let _ = rdkit_mol_from_inchi(&standard);
            let _ = parse_inchi(&standard);
        }
        molecules.push(mol);
    }
    assert_eq!(molecules.len(), 384);

    let candidates = group_candidates(&molecules);
    assert!(!candidates.is_empty());
    for chunk in molecules.chunks(2).take(96) {
        if let [left, right] = chunk {
            for policy in [
                IdentityPolicy::StandardInchiString,
                IdentityPolicy::StandardInchiKey,
                IdentityPolicy::StereoIgnored,
                IdentityPolicy::IsotopeIgnored,
            ] {
                let _ = compare_molecules(left, right, policy);
                let _ = compare_molecules_with_accurate_cip_preflight(left, right, policy);
            }
            let _ = compare_indexed_graph_relation(
                left,
                right,
                IndexedGraphRelationMode {
                    graph_strictness: GraphStrictness::ChemicalGraphExact,
                    atom_maps: AtomMapPolicy::Ignore,
                },
            );
        }
    }
    let subset = &molecules[..64];
    let _ = deduplicate_verified(subset, IdentityPolicy::StandardInchiString);
    let _ = deduplicate_verified_with_accurate_cip_preflight(
        subset,
        IdentityPolicy::StandardInchiString,
    );
}

#[test]
fn inchi_parser_rejects_malformed_and_bounded_inputs_without_panicking() {
    let malformed = [
        "",
        "InChI=",
        "InChI=1S/",
        "InChI=1S/C2H6/c1-2/h1-2H3/q",
        "InChI=1S/C2H6/c1-999/h1-2H3",
        "InChI=1S/C2H6/c1-(2/h1-2H3",
        "InChI=1S/999999999999C",
        "not-an-inchi",
    ];
    for text in malformed {
        let _ = parse_inchi(text);
        let limits = InchiParseLimits {
            max_input_bytes: 16,
            ..InchiParseLimits::default()
        };
        let _ = parse_inchi_with_limits(text, &limits);
        let _ = rdkit_mol_from_inchi(text);
        let _ = standard_inchi_key(text);
    }
}

#[test]
fn approximate_inchi_roundtrips_the_full_special_chemistry_corpus() {
    let mut parsed = 0usize;
    let mut reparsed = 0usize;
    for line in STRESS.lines() {
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        if mol.atom_count() == 0 || mol.atom_count() > 160 {
            continue;
        }
        parsed += 1;
        let value = inchi(&mol);
        let _ = inchi_key(&value);
        let _ = build_inchi_index(&mol);
        if parse_inchi(&value).is_ok() {
            reparsed += 1;
        }
    }
    assert!(parsed > 3_500);
    assert!(reparsed > 3_000);
}

#[test]
fn standard_inchi_reader_handles_diverse_native_outputs() {
    let mut generated = 0usize;
    let mut reconstructed = 0usize;
    for line in STRESS.lines() {
        if generated == 1_024 {
            break;
        }
        let smiles = line.split_whitespace().next().unwrap_or("");
        let Ok(mol) = parse(smiles) else { continue };
        if mol.atom_count() == 0 || mol.atom_count() > 120 {
            continue;
        }
        let Ok(value) = standard_inchi(&mol) else {
            continue;
        };
        generated += 1;
        let _ = standard_inchi_key(&value);
        if rdkit_mol_from_inchi(&value).is_ok() {
            reconstructed += 1;
        }
    }
    assert_eq!(generated, 1_024);
    assert!(reconstructed > 900);
}

#[test]
fn parser_layer_boundary_matrix_is_typed_and_panic_free() {
    let cases = [
        "InChI=1S/C2/c1=2",
        "InChI=1S/C2/c1#2",
        "InChI=1S/C4/c1-4(2)3",
        "InChI=1S/C5/c1-5(2)(3)4",
        "InChI=1S/C2H6/c1-2/h",
        "InChI=1S/C2H6/c1-2/h,1H3,,2H3",
        "InChI=1S/C2H6/c1-2/h1H0,2H255",
        "InChI=1S/C3/c1-2-3/q1+1,2-3-1",
        "InChI=1S/C3/c1-2-3/q1-2+1",
        "InChI=1S/C3/c1-2-3/qinvalid+1",
        "InChI=1S/C2/c1-2/i1/13C,2/2H",
        "InChI=1S/C2/c1-2/i,1/foo,2/999C",
        "InChI=1S/C4/c1-2-3-4/b2-3+,3-2-,bad+",
        "InChI=1S/C4/c1-2-3-4/t1+,2-,bad+",
        "InChI=1S/C4/c1-2-3-4/m",
        "InChI=1S/C4/c1-2-3-4/m1,,2/s",
        "InChI=1S/C999999999999999999999999999999999999",
        "InChI=1S/Xx2",
        "InChI=1S/c1-2",
    ];
    for input in cases {
        let _ = parse_inchi(input);
        let _ = rdkit_mol_from_inchi(input);
    }
}
