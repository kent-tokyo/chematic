//! ChemDraw-saved CDXML files from RDKit's test corpus (BSD-3, see
//! `validation/rdkit_cdxml/README.md`): abbreviation expansion (#1316),
//! generic labels and query bonds (#1317), and writer round trips (#1319).
use std::path::{Path, PathBuf};

use chematic_core::{BondOrder, Molecule};
use chematic_mol::{
    CdxmlError, CdxmlParseLimits, CdxmlParseOptions, CdxmlParseReport, parse_cdxml_all,
    parse_cdxml_all_with_report, write_cdxml,
};

fn corpus_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../validation/rdkit_cdxml")
}

fn read(name: &str) -> String {
    std::fs::read_to_string(corpus_dir().join(name)).unwrap()
}

fn corpus() -> Vec<String> {
    fn walk(dir: &Path, base: &Path, out: &mut Vec<String>) {
        for entry in std::fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                walk(&path, base, out);
            } else if path.extension().is_some_and(|ext| ext == "cdxml") {
                let relative = path.strip_prefix(base).unwrap();
                out.push(relative.to_string_lossy().replace('\\', "/"));
            }
        }
    }
    let mut names = Vec::new();
    walk(&corpus_dir(), &corpus_dir(), &mut names);
    names.sort();
    names
}

fn intentionally_invalid(name: &str) -> bool {
    name.starts_with("bad-") || name == "malformed.cdxml"
}

fn report(input: &str) -> Result<CdxmlParseReport, CdxmlError> {
    parse_cdxml_all_with_report(
        input,
        &CdxmlParseOptions::default(),
        &CdxmlParseLimits::default(),
    )
}

fn canonical(mol: &Molecule) -> String {
    chematic_smiles::canonical_smiles(mol)
}

/// Canonical SMILES with aromaticity perceived, so a Kekulé ChemDraw ring
/// and RDKit's aromatic spelling compare equal.
fn normalized(mol: &Molecule) -> String {
    let mut mol = mol.clone();
    chematic_perception::aromatize(&mut mol);
    canonical(&mol)
}

fn normalized_smiles(smiles: &str) -> String {
    // RDKit writes R groups as `[*:n]`; chematic reads them as plain `*`.
    let mut cleaned = String::new();
    let mut rest = smiles;
    while let Some(start) = rest.find("[*:") {
        cleaned.push_str(&rest[..start]);
        cleaned.push_str("[*]");
        let tail = &rest[start..];
        rest = &tail[tail.find(']').unwrap() + 1..];
    }
    cleaned.push_str(rest);
    normalized(&chematic_smiles::parse(&cleaned).unwrap())
}

fn sorted(mut values: Vec<String>) -> Vec<String> {
    values.sort();
    values
}

#[test]
fn every_valid_corpus_file_parses_with_the_report_api() {
    let names = corpus();
    assert!(names.len() >= 60, "corpus missing: {}", names.len());
    for name in names.iter().filter(|name| !intentionally_invalid(name)) {
        let parsed = report(&read(name)).unwrap_or_else(|error| panic!("{name}: {error}"));
        assert!(!parsed.molecules.is_empty(), "{name}");
    }
}

#[test]
fn strict_api_fails_only_on_unrepresentable_features() {
    for name in corpus().iter().filter(|name| !intentionally_invalid(name)) {
        if let Err(error) = parse_cdxml_all(&read(name)) {
            assert!(
                matches!(error, CdxmlError::UnsupportedPseudoatomLabel(_))
                    || error.unsupported_bond_order() == Some("hydrogen"),
                "{name}: {error}"
            );
        }
    }
}

#[test]
fn matches_rdkit_expected_structures() {
    let expected: serde_json::Value = serde_json::from_str(&read("expected_smiles.json")).unwrap();
    for (name, smiles) in expected.as_object().unwrap() {
        let got = parse_cdxml_all(&read(name))
            .unwrap_or_else(|error| panic!("{name}: {error}"))
            .iter()
            .map(|(mol, _)| normalized(mol))
            .collect();
        let want = smiles
            .as_array()
            .unwrap()
            .iter()
            .map(|smiles| normalized_smiles(smiles.as_str().unwrap()))
            .collect();
        assert_eq!(sorted(got), sorted(want), "{name}");
    }
}

#[test]
fn parse_write_parse_keeps_canonical_smiles() {
    for name in corpus().iter().filter(|name| !intentionally_invalid(name)) {
        let Ok(first) = parse_cdxml_all(&read(name)) else {
            continue;
        };
        for (mol, coords) in &first {
            let written = write_cdxml(mol, coords);
            let again = parse_cdxml_all(&written).unwrap_or_else(|e| panic!("{name}: {e}"));
            assert_eq!(again.len(), 1, "{name}");
            assert_eq!(canonical(&again[0].0), canonical(mol), "{name}");
        }
    }
}

#[test]
fn hashed_wedges_are_written_directionally() {
    let (mol, coords) = parse_cdxml_all(&read("stereo1.cdxml")).unwrap().remove(0);
    let written = write_cdxml(&mol, &coords);
    assert!(written.contains("Display=\"WedgedHashBegin\""));
    assert!(!written.contains("Display=\"Hash\""));
}

#[test]
fn wavy_alkene_stays_unspecified_after_round_trip() {
    for (mol, coords) in parse_cdxml_all(&read("wavy.cdxml")).unwrap() {
        let again = parse_cdxml_all(&write_cdxml(&mol, &coords)).unwrap();
        assert_eq!(canonical(&again[0].0), canonical(&mol));
    }
}

#[test]
fn abbreviations_expand_and_keep_their_labels() {
    let parsed = report(&read("reaction-with-boc.cdxml")).unwrap();
    let boc = parsed
        .abbreviations
        .iter()
        .find(|abbreviation| abbreviation.label.as_deref() == Some("Boc"))
        .expect("Boc abbreviation is reported");
    let (mol, _) = &parsed.molecules[boc.molecule_index];
    // Boc = C(=O)OC(CH3)3: 7 heavy atoms.
    assert_eq!(boc.atoms.len(), 7);
    let mut elements = boc
        .atoms
        .iter()
        .map(|&atom| mol.atom(atom).element.symbol())
        .collect::<Vec<_>>();
    elements.sort_unstable();
    assert_eq!(elements, ["C", "C", "C", "C", "C", "O", "O"]);
    assert!(parsed.diagnostics.is_empty(), "{:?}", parsed.diagnostics);
}

#[test]
fn multi_attachment_abbreviations_follow_bond_ordering() {
    let mols = parse_cdxml_all(&read("stereo5.cdxml")).unwrap();
    assert_eq!(mols.len(), 1);
    let smiles = canonical(&mols[0].0);
    assert!(!smiles.contains('*'), "{smiles}");
    assert!(smiles.contains('@'), "stereo kept: {smiles}");
}

const MISSING_CONNECTION_POINT: &str = r#"<CDXML><page><fragment>
<n id="1" p="0 0"/>
<n id="2" p="10 0" NodeType="Nickname"><fragment>
<n id="3" p="10 0" Element="8"/><n id="4" p="20 0"/>
<b id="5" B="3" E="4"/>
</fragment><t><s>OMe</s></t></n>
<n id="6" p="-10 0"/>
<b id="7" B="1" E="2"/><b id="8" B="6" E="2"/>
</fragment></page></CDXML>"#;

#[test]
fn unexpandable_abbreviation_is_a_typed_error_naming_it() {
    let error = match parse_cdxml_all(MISSING_CONNECTION_POINT) {
        Ok(_) => panic!("ambiguous attachment accepted"),
        Err(error) => error,
    };
    assert_eq!(error.unexpandable_abbreviation(), Some("OMe"));
    assert!(
        error
            .to_string()
            .starts_with("cannot expand CDXML abbreviation \"OMe\"")
    );
    assert!(!matches!(error, CdxmlError::UnknownAtomRef(_)));
}

#[test]
fn single_atom_abbreviation_without_connection_point_attaches_to_that_atom() {
    let input = r#"<CDXML><page><fragment>
<n id="1" p="0 0" Element="8"/>
<n id="2" p="10 0" NodeType="Fragment"><fragment><n id="3" p="10 0"/></fragment></n>
<b id="4" B="1" E="2"/>
</fragment></page></CDXML>"#;
    let (mol, _) = parse_cdxml_all(input).unwrap().remove(0);
    assert_eq!(canonical(&mol), "OC");
}

#[test]
fn generic_query_labels_refuse_by_default_and_become_wildcards_in_the_report() {
    let input = read("queries/qatom_MH.cdxml");
    assert!(matches!(
        parse_cdxml_all(&input),
        Err(CdxmlError::UnsupportedPseudoatomLabel(label)) if label == "MH"
    ));
    let parsed = report(&input).unwrap();
    let (mol, _) = &parsed.molecules[0];
    assert_eq!(mol.atoms().filter(|(_, atom)| atom.wildcard).count(), 1);
    let diagnostic = &parsed.diagnostics[0];
    assert_eq!(diagnostic.code, "generic_query_atom_as_wildcard");
    assert!(
        diagnostic.message.contains("\"MH\""),
        "{}",
        diagnostic.message
    );
    assert!(diagnostic.object_id.is_some());
}

#[test]
fn query_bond_orders_map_to_core_query_bonds() {
    for (name, order) in [
        ("anybond.cdxml", BondOrder::QueryAny),
        ("queries/furan_sa.cdxml", BondOrder::QuerySingleOrAromatic),
        ("queries/furan_da.cdxml", BondOrder::QueryDoubleOrAromatic),
        ("queries/furan_sd.cdxml", BondOrder::QuerySingleOrDouble),
    ] {
        let mols = parse_cdxml_all(&read(name)).unwrap_or_else(|e| panic!("{name}: {e}"));
        let (mol, coords) = &mols[0];
        assert!(
            mol.bonds().any(|(_, bond)| bond.order == order),
            "{name}: {order:?}"
        );
        let again = parse_cdxml_all(&write_cdxml(mol, coords)).unwrap();
        assert!(
            again[0].0.bonds().any(|(_, bond)| bond.order == order),
            "{name}"
        );
    }
}

#[test]
fn unrepresentable_bond_is_dropped_with_a_diagnostic_in_the_report() {
    let input = read("queries/qbond_hydrogen.cdxml");
    assert_eq!(
        parse_cdxml_all(&input)
            .err()
            .and_then(|e| e.unsupported_bond_order().map(str::to_string))
            .as_deref(),
        Some("hydrogen")
    );
    let parsed = report(&input).unwrap();
    assert_eq!(parsed.molecules.len(), 1);
    assert_eq!(parsed.diagnostics.len(), 1);
    assert_eq!(parsed.diagnostics[0].code, "unsupported_bond_order_dropped");
}

fn nickname_doc(inner: &str, outer_extra: &str, bonds: &str) -> String {
    format!(
        r#"<CDXML><page><fragment>
<n id="1" p="0 0"/>
<n id="2" p="10 0" NodeType="Nickname"{outer_extra}><fragment>{inner}</fragment><t><s>Me&amp;Co</s></t></n>
{bonds}
</fragment></page></CDXML>"#
    )
}

fn expansion_error(input: &str) -> String {
    match parse_cdxml_all(input) {
        Ok(_) => panic!("expansion accepted"),
        Err(error) => {
            assert_eq!(error.unexpandable_abbreviation(), Some("Me&Co"), "{error}");
            error.to_string()
        }
    }
}

#[test]
fn connection_point_must_have_exactly_one_bond_to_a_known_atom() {
    let two_bonds = nickname_doc(
        r#"<n id="3" p="10 0"/><n id="4" p="20 0"/><n id="9" NodeType="ExternalConnectionPoint"/>
<b B="3" E="9"/><b B="4" E="9"/>"#,
        "",
        r#"<b B="1" E="2"/>"#,
    );
    assert!(expansion_error(&two_bonds).contains("exactly one bond"));
    let unknown = nickname_doc(
        r#"<n id="3" p="10 0"/><n id="9" NodeType="ExternalConnectionPoint"/><b B="9" E="42"/>"#,
        "",
        r#"<b B="1" E="2"/>"#,
    );
    assert!(expansion_error(&unknown).contains("unknown atom 42"));
}

#[test]
fn bond_ordering_selects_numbered_connection_points() {
    let inner = r#"<n id="3" p="10 0" Element="8"/><n id="4" p="20 0" Element="7"/>
<n id="8" NodeType="ExternalConnectionPoint" ExternalConnectionNum="2"/>
<n id="9" NodeType="ExternalConnectionPoint" ExternalConnectionNum="1"/>
<b B="3" E="4"/><b B="3" E="8"/><b B="4" E="9"/>"#;
    let bonds = r#"<n id="5" p="30 0" Element="9"/>
<b id="20" B="1" E="2"/><b id="21" B="5" E="2"/>"#;
    let input = nickname_doc(inner, r#" BondOrdering="21 20""#, bonds);
    let (mol, _) = parse_cdxml_all(&input).unwrap().remove(0);
    // Bond 21 (F) takes connection 1 (N); bond 20 (C) takes connection 2 (O).
    assert_eq!(
        canonical(&mol),
        canonical(&chematic_smiles::parse("CONF").unwrap())
    );
    let missing = nickname_doc(inner, r#" BondOrdering="21 0 20""#, bonds);
    assert!(expansion_error(&missing).contains("no external connection point for bond 20"));
}

#[test]
fn report_lists_label_free_and_unexpanded_nodes() {
    let input = r#"<CDXML><page><fragment>
<n id="1" p="0 0"/>
<n id="2" p="10 0" NodeType="Fragment"><fragment><n id="3" p="10 0"/></fragment></n>
<n id="4" p="20 0" NodeType="Nickname"/>
<n id="5" p="30 0" NodeType="ExternalConnectionPoint"/>
<n p="40 0"><fragment/><fragment><n id="6" p="40 0"/></fragment></n>
<b B="1" E="2"/><b B="2" E="4"/><b B="4" E="5"/>
</fragment></page></CDXML>"#;
    let parsed = report(input).unwrap();
    assert_eq!(parsed.abbreviations.len(), 1);
    assert_eq!(parsed.abbreviations[0].label, None);
    assert_eq!(parsed.abbreviations[0].node_id, "2");
    let codes = parsed
        .diagnostics
        .iter()
        .map(|diagnostic| diagnostic.code)
        .collect::<Vec<_>>();
    assert_eq!(
        codes,
        [
            "nested_fragment_ignored",
            "connection_point_dropped",
            "unexpanded_abbreviation"
        ]
    );
    assert_eq!(canonical(&parsed.molecules[0].0), "CCC");
}

#[test]
fn unlabelled_abbreviation_error_names_the_node() {
    let input = r#"<CDXML><page><fragment>
<n id="1" p="0 0"/>
<n id="2" p="10 0" NodeType="Nickname"><fragment><n id="3"/><n id="4"/><b B="3" E="4"/></fragment></n>
<b B="1" E="2"/>
</fragment></page></CDXML>"#;
    let error = parse_cdxml_all(input).err().unwrap();
    assert_eq!(error.unexpandable_abbreviation(), Some("node 2"));
    assert_eq!(
        CdxmlError::MissingBondEndpoint.unexpandable_abbreviation(),
        None
    );
    assert_eq!(
        CdxmlError::InvalidDocument("x".into()).unexpandable_abbreviation(),
        None
    );
}

#[test]
fn expansion_respects_the_atom_limit() {
    let input = nickname_doc(
        r#"<n id="3"/><n id="4"/><n id="5"/><n id="9" NodeType="ExternalConnectionPoint"/>
<b B="3" E="4"/><b B="4" E="5"/><b B="3" E="9"/>"#,
        "",
        r#"<b B="1" E="2"/>"#,
    );
    let limits = CdxmlParseLimits {
        max_atoms: 3,
        ..CdxmlParseLimits::default()
    };
    let error = chematic_mol::parse_cdxml_all_with_limits(&input, &limits)
        .err()
        .unwrap();
    assert!(matches!(
        error,
        CdxmlError::ResourceLimit {
            resource: "atoms",
            ..
        }
    ));
}

#[test]
fn unrecognized_generic_label_is_reported_not_guessed() {
    let input = r#"<CDXML><page><fragment>
<n id="1" p="0 0" NodeType="GenericNickname" GenericNickname="Foo"/>
<n id="2" p="10 0"/><b B="1" E="2"/>
</fragment></page></CDXML>"#;
    assert!(parse_cdxml_all(input).is_err());
    let parsed = report(input).unwrap();
    assert_eq!(parsed.diagnostics[0].code, "unrecognized_generic_label");
    assert_eq!(parsed.diagnostics[0].object_id.as_deref(), Some("1"));
}
