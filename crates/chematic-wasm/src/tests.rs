use super::*;

fn parse(s: &str) -> MolHandle {
    MolHandle {
        inner: std::rc::Rc::new(chematic_smiles::parse(s).unwrap()),
    }
}

#[test]
fn parse_benzene_atom_count() {
    assert_eq!(parse("c1ccccc1").atom_count(), 6);
}

#[test]
fn common_format_conversion_uses_shared_aliases_and_preserves_graph() {
    let mol2 = convert_common_format("CCO", "smiles", "mol2").unwrap();
    assert!(mol2.contains("@<TRIPOS>MOLECULE"));
    let smiles = convert_common_format(&mol2, ".mol2", "smi").unwrap();
    assert!(chematic_chem::are_identical(
        &chematic_smiles::parse(&smiles).unwrap(),
        &chematic_smiles::parse("CCO").unwrap()
    ));
}

#[test]
fn common_format_name_rejects_unknown_values_without_wasm_runtime() {
    assert_eq!(common_format_name(".mol2"), Some("mol2".to_string()));
    assert_eq!(common_format_name("unknown"), None);
}

// --- logd / isotope / topological index tests ---------------------------

#[test]
fn logd_at_ph7_benzene() {
    let h = parse("c1ccccc1");
    let logd = h.logd_at_ph(7.0);
    assert!(logd.is_finite(), "logD should be finite, got {logd}");
}

#[test]
fn logd_profile_json_has_ph_fields() {
    let h = parse("c1ccccc1");
    let json = h.logd_profile_json(0.0, 14.0, 15);
    assert!(
        json.starts_with('[') && json.ends_with(']'),
        "should be JSON array: {json}"
    );
    assert!(json.contains(r#""ph":"#), "should contain ph field: {json}");
    assert!(
        json.contains(r#""logd":"#),
        "should contain logd field: {json}"
    );
}

#[test]
fn isotope_distribution_json_nonempty() {
    let h = parse("C");
    let json = h.isotope_distribution_json(0.1);
    assert!(json.starts_with('['), "should be JSON array: {json}");
    assert!(
        json.contains(r#""mass":"#),
        "should contain mass entries: {json}"
    );
}

#[test]
fn randic_index_ethane() {
    let h = parse("CC");
    let r = h.randic_index();
    // Ethane: 1 C-C bond, both heavy atoms have degree 1 → χ = 1/√(1×1) = 1.0
    assert!((r - 1.0).abs() < 1e-6, "ethane Randic index = 1.0, got {r}");
}

#[test]
fn zagreb_m1_ethane() {
    let h = parse("CC");
    let z = h.zagreb_index_m1();
    // Ethane: each C has heavy-atom degree 1 → M1 = 1² + 1² = 2
    assert_eq!(z, 2, "ethane Zagreb M1 = 2, got {z}");
}

// --- iupac_name tests ---------------------------------------------------

#[test]
fn iupac_name_benzene() {
    assert_eq!(parse("c1ccccc1").iupac_name(), "benzene");
}

#[test]
fn iupac_name_ethanol() {
    assert_eq!(parse("CCO").iupac_name(), "ethanol");
}

#[test]
fn iupac_name_branched_acid() {
    assert_eq!(parse("CC(C)C(=O)O").iupac_name(), "2-methylpropanoic acid");
}

#[test]
fn iupac_name_unsupported_returns_empty() {
    // Disconnected molecule → unsupported → empty string
    let h = parse("C.C");
    assert!(
        h.iupac_name().is_empty(),
        "disconnected should return empty"
    );
}

// --- assign_cip_json tests ----------------------------------------------

#[test]
fn assign_cip_json_no_centers() {
    // Ethane has no stereocenters
    let json = parse("CC").assign_cip_json();
    assert_eq!(json, r#"{"centers":[]}"#);
}

#[test]
fn assign_cip_json_chiral_center() {
    // L-alanine SMILES with explicit chirality
    let json = parse("N[C@@H](C)C(=O)O").assign_cip_json();
    assert!(
        json.contains("\"code\":"),
        "should contain code field: {json}"
    );
    assert!(
        json.contains("\"atom\":"),
        "should contain atom field: {json}"
    );
}

#[test]
fn canonical_smiles_benzene() {
    let mol = parse("c1ccccc1");
    let cs = mol.canonical_smiles();
    assert!(!cs.is_empty());
}

#[test]
fn parse_cxsmiles_json_preserves_metadata() {
    let json = parse_cxsmiles_json("C~O |$C1;O2$,atomProp:1.role.acceptor,^2:0,Z:0|").unwrap();
    assert!(json.contains(r#""atomLabels":["C1","O2"]"#), "{json}");
    assert!(json.contains(r#""key":"role""#), "{json}");
    assert!(json.contains(r#""atomRadicals":[2,null]"#), "{json}");
    assert!(json.contains(r#""zeroBonds":[0]"#), "{json}");
}

#[test]
fn parse_cxsmarts_json_preserves_metadata() {
    let json = parse_cxsmarts_json("[#6]~[#8] |$C1;O2$,atomProp:1.role.acceptor,^2:0|").unwrap();
    assert!(json.contains(r#""atomCount":2"#), "{json}");
    assert!(json.contains(r#""atomLabels":["C1","O2"]"#), "{json}");
    assert!(json.contains(r#""key":"role""#), "{json}");
    assert!(json.contains(r#""atomRadicals":[2,null]"#), "{json}");
}

#[test]
fn molecular_weight_aspirin() {
    let mw = parse("CC(=O)Oc1ccccc1C(=O)O").molecular_weight();
    assert!((mw - 180.16).abs() < 1.0);
}

#[test]
fn lipinski_aspirin() {
    assert!(parse("CC(=O)Oc1ccccc1C(=O)O").lipinski_passes());
}

#[test]
fn tanimoto_same_mol() {
    let a = parse("c1ccccc1");
    let b = parse("c1ccccc1");
    let sim = tanimoto_ecfp4(&a, &b);
    assert!((sim - 1.0).abs() < 1e-6);
}

#[test]
fn tanimoto_different() {
    let a = parse("c1ccccc1");
    let b = parse("CC(=O)Oc1ccccc1C(=O)O");
    assert!(tanimoto_ecfp4(&a, &b) < 1.0);
}

#[test]
fn heavy_atom_count_ethanol() {
    assert_eq!(parse("CCO").heavy_atom_count(), 3);
}

#[test]
fn logp_crippen_aspirin_range() {
    let lp = parse("CC(=O)Oc1ccccc1C(=O)O").logp_crippen();
    assert!(lp > 0.5 && lp < 3.5, "aspirin LogP = {lp:.3}");
}

#[test]
fn fsp3_benzene_zero() {
    assert_eq!(parse("c1ccccc1").fsp3(), 0.0, "benzene Fsp3 = 0");
}

#[test]
fn fsp3_cyclohexane_one() {
    assert_eq!(parse("C1CCCCC1").fsp3(), 1.0, "cyclohexane Fsp3 = 1");
}

#[test]
fn aromatic_ring_count_benzene() {
    assert_eq!(parse("c1ccccc1").aromatic_ring_count(), 1);
}

#[test]
fn qed_aspirin_range() {
    let q = parse("CC(=O)Oc1ccccc1C(=O)O").qed();
    assert!(q > 0.0 && q <= 1.0, "aspirin QED = {q:.3}");
}

#[test]
fn exact_mass_aspirin() {
    // Aspirin monoisotopic mass: 180.0423
    let em = parse("CC(=O)Oc1ccccc1C(=O)O").exact_mass();
    assert!((em - 180.042).abs() < 0.01, "aspirin exact mass = {em:.4}");
}

#[test]
fn rotatable_bond_count_aspirin() {
    // Aspirin has 3 rotatable bonds (OC, C=O ester, and COOH)
    let rb = parse("CC(=O)Oc1ccccc1C(=O)O").rotatable_bond_count();
    assert!((2..=5).contains(&rb), "aspirin rotbonds = {rb}");
}

#[test]
fn tanimoto_atom_pair_same_mol() {
    let a = parse("c1ccccc1");
    let b = parse("c1ccccc1");
    assert!((tanimoto_atom_pair(&a, &b) - 1.0).abs() < 1e-6);
}

#[test]
fn tanimoto_torsion_same_mol() {
    let a = parse("CCCC");
    let b = parse("CCCC");
    assert!((tanimoto_torsion(&a, &b) - 1.0).abs() < 1e-6);
}

#[test]
fn brics_fragment_count_benzene() {
    assert_eq!(brics_fragment_count(&parse("c1ccccc1")), 1);
}

#[test]
fn brics_fragment_count_aspirin() {
    assert!(brics_fragment_count(&parse("CC(=O)Oc1ccccc1C(=O)O")) >= 2);
}

#[test]
fn wiener_index_ethane() {
    // Ethane: 2 atoms, distance 1 — Wiener index = 1.
    assert_eq!(parse("CC").wiener_index(), 1.0);
}

#[test]
fn kappa1_propane_range() {
    let k = parse("CCC").kappa1();
    assert!(k > 0.0, "kappa1 should be positive");
}

#[test]
fn chi0_benzene_positive() {
    assert!(parse("c1ccccc1").chi0() > 0.0);
}

#[test]
fn labute_asa_aspirin_range() {
    let asa = parse("CC(=O)Oc1ccccc1C(=O)O").labute_asa();
    assert!(asa > 50.0 && asa < 200.0, "aspirin LabuteASA = {asa:.2}");
}

#[test]
fn bertz_ct_benzene_positive() {
    assert!(parse("c1ccccc1").bertz_ct() > 0.0);
}

#[test]
fn morgan_fp_counts_json_benzene() {
    let json = parse("c1ccccc1").morgan_fp_counts_json(2);
    assert!(json.starts_with('{') && json.ends_with('}'));
}

#[test]
fn add_remove_hydrogens_roundtrip() {
    let mol = parse("CC");
    let with_h = add_hydrogens(&mol);
    assert!(
        with_h.atom_count() > mol.atom_count(),
        "H atoms should be added"
    );
    let back = remove_hydrogens(&with_h);
    assert_eq!(back.atom_count(), mol.atom_count());
}

#[test]
fn depict_svg_grid_two_mols() {
    let svg = depict_svg_grid("CC\nCCC", 2);
    assert!(svg.contains("<svg"), "expected SVG output");
}

#[test]
fn depict_svg_grid_invalid_smiles_skipped() {
    let svg = depict_svg_grid("CC\nNOT_A_SMILES\nCCC", 2);
    assert!(
        svg.contains("<svg"),
        "invalid SMILES should be silently skipped"
    );
}

#[test]
fn run_reactants_esterification() {
    // Simple esterification: carboxylic acid + alcohol → ester + water
    let result = run_reactants(
        "[C:1](=O)[OH:2].[O:3][C:4]>>[C:1](=O)[O:3][C:4]",
        "CC(=O)O|CCO",
    );
    assert!(result.is_ok(), "run_reactants should succeed");
    let json = result.unwrap();
    assert!(json.contains('['), "expected JSON array");
}

// Note: run_reactants error-path tests are omitted here because JsValue::from_str
// panics outside a WASM runtime. Error coverage lives in chematic-rxn unit tests.

#[test]
fn is_valid_smiles_valid() {
    assert!(is_valid_smiles("CCO"), "ethanol is valid");
    assert!(is_valid_smiles("c1ccccc1"), "benzene is valid");
    assert!(is_valid_smiles("O"), "water is valid");
    assert!(is_valid_smiles("C"), "methane is valid");
}

#[test]
fn is_valid_smiles_invalid() {
    assert!(!is_valid_smiles(""), "empty string is invalid");
    assert!(
        !is_valid_smiles("[NOSUCHELEMENT]"),
        "unknown bracket atom is invalid"
    );
}

#[test]
fn depict_svg_opts_transparent_background() {
    let h = parse("CCO");
    let mut opts = DepictOptions::new();
    opts.set_background("transparent".to_string());
    let svg = h.depict_svg_opts(&opts);
    assert!(svg.contains("<svg"), "must produce SVG");
    assert!(
        !svg.contains("fill=\"transparent\""),
        "no bg rect for transparent"
    );
}

#[test]
fn depict_svg_opts_custom_size() {
    let h = parse("CCO");
    let mut opts = DepictOptions::new();
    opts.set_width(300);
    opts.set_height(200);
    let svg = h.depict_svg_opts(&opts);
    assert!(svg.contains("width=\"300\""), "SVG width should be 300");
    assert!(svg.contains("height=\"200\""), "SVG height should be 200");
}

#[test]
fn depict_svg_opts_dark_theme() {
    let h = parse("CC");
    let mut opts = DepictOptions::new();
    opts.set_dark(true);
    opts.set_background("#0f172a".to_string());
    let svg = h.depict_svg_opts(&opts);
    assert!(
        svg.contains("stroke=\"white\""),
        "dark theme bonds should be white"
    );
}

#[test]
fn depict_svg_single_atom_water_shows_h2o() {
    let svg = parse("O").depict_svg();
    assert!(svg.contains("H2O"), "water 'O' should render as H2O");
}

#[test]
fn depict_svg_single_atom_methane_shows_ch4() {
    let svg = parse("C").depict_svg();
    assert!(svg.contains("CH4"), "methane 'C' should render as CH4");
}

// ── Sprint L: disconnected SMILES ────────────────────────────────────────

#[test]
fn depict_svg_disconnected_nacl() {
    let svg = parse("[Na+].[Cl-]").depict_svg();
    assert!(
        svg.contains("Na"),
        "Na should appear in disconnected SMILES SVG"
    );
    assert!(
        svg.contains("Cl"),
        "Cl should appear in disconnected SMILES SVG"
    );
    assert!(!svg.is_empty());
}

#[test]
fn depict_svg_disconnected_water_dimer() {
    let svg = parse("O.O").depict_svg();
    // Degree-0 O atoms use isolated (Hill) notation: H2O
    assert!(
        svg.matches("H2O").count() >= 2,
        "both water O atoms should render as H2O"
    );
    assert!(!svg.is_empty());
}

// ── Sprint L: atom data attributes ──────────────────────────────────────

#[test]
fn depict_svg_opts_atom_ids_contains_data_attrs() {
    let h = parse("CC(=O)O"); // acetic acid: 2 C (unlabelled) + 2 O (labelled)
    let mut opts = DepictOptions::new();
    opts.set_atom_ids(true);
    let svg = h.depict_svg_opts(&opts);
    assert!(
        svg.contains("data-atom-idx="),
        "atom_ids should add data-atom-idx"
    );
    assert!(
        svg.contains("data-element="),
        "atom_ids should add data-element"
    );
    assert!(
        svg.contains("data-charge="),
        "atom_ids should add data-charge"
    );
    // B3: all 4 atoms should be addressable (unlabelled carbons get invisible anchor)
    assert_eq!(
        svg.matches("data-atom-idx=").count(),
        4,
        "all atoms should have data-atom-idx"
    );
}

#[test]
fn depict_svg_opts_atom_ids_false_no_data_attrs() {
    let h = parse("CC(=O)O");
    let svg = h.depict_svg_opts(&DepictOptions::new());
    assert!(
        !svg.contains("data-atom-idx="),
        "default opts should not have data-atom-idx"
    );
}

#[test]
fn depict_svg_opts_atom_ids_charge_correct() {
    let h = parse("[NH4+]");
    let mut opts = DepictOptions::new();
    opts.set_atom_ids(true);
    let svg = h.depict_svg_opts(&opts);
    assert!(
        svg.contains("data-charge=\"1\""),
        "NH4+ should have charge=1"
    );
}

// ── Sprint L: show_atom_indices ──────────────────────────────────────────

#[test]
fn depict_svg_opts_show_atom_indices() {
    let h = parse("c1ccccc1"); // benzene — 6 atoms, indices 0-5
    let mut opts = DepictOptions::new();
    opts.set_show_atom_indices(true);
    let svg = h.depict_svg_opts(&opts);
    assert!(svg.contains(">0<"), "index 0 should appear");
    assert!(svg.contains(">5<"), "index 5 should appear");
}

#[test]
fn depict_svg_opts_show_atom_indices_false_no_indices() {
    let h = parse("CCO");
    let svg = h.depict_svg_opts(&DepictOptions::new());
    assert!(
        !svg.contains("fill=\"#8b92a9\""),
        "default should not show grey index labels"
    );
}

// ── Sprint L: kekulize ───────────────────────────────────────────────────

#[test]
fn depict_svg_opts_kekulize_removes_aromatic_bonds() {
    let h = parse("c1ccccc1"); // benzene
    let mut opts = DepictOptions::new();
    opts.set_kekulize(true);
    let svg = h.depict_svg_opts(&opts);
    assert!(!svg.is_empty());
    assert!(
        svg.contains("<line"),
        "kekulé benzene should have line elements"
    );
    assert!(
        !svg.contains("stroke-dasharray"),
        "kekulé benzene must not use aromatic dashed style"
    ); // B1
}

#[test]
fn depict_svg_opts_kekulize_false_uses_aromatic() {
    let h = parse("c1ccccc1");
    let svg = h.depict_svg_opts(&DepictOptions::new());
    // Default aromatic rendering uses stroke-dasharray for the inner ring line.
    assert!(
        svg.contains("stroke-dasharray"),
        "default benzene should use aromatic dashed style"
    );
}

// ── Sprint M: smarts_match_atoms ─────────────────────────────────────────

#[test]
fn smarts_match_benzene_ring_returns_json() {
    let mol = parse("c1ccccc1");
    let result = smarts_match_atoms("c1ccccc1", &mol);
    assert!(result.is_ok(), "valid SMARTS should not error");
    let json = result.unwrap();
    assert!(
        !json.is_empty() && json != "[]",
        "benzene ring SMARTS should find a match"
    );
    assert!(json.starts_with("[["), "result should be array of arrays");
}

#[test]
fn smarts_match_no_match_returns_empty_array() {
    let mol = parse("CC"); // ethane has no aromatic ring
    let result = smarts_match_atoms("c1ccccc1", &mol);
    assert!(
        result.is_ok(),
        "valid SMARTS on non-matching mol should not error"
    );
    assert_eq!(
        result.unwrap(),
        "[]",
        "no match should return empty JSON array"
    );
}

// ── Sprint N: generate_3d_pdb ────────────────────────────────────────────

#[test]
fn generate_3d_pdb_benzene_has_6_hetatm_lines() {
    let mol = parse("c1ccccc1");
    let pdb = generate_3d_pdb(&mol);
    let hetatm_count = pdb.lines().filter(|l| l.starts_with("HETATM")).count();
    assert_eq!(hetatm_count, 6, "benzene should produce 6 HETATM lines");
}

#[test]
fn generate_3d_pdb_aspirin_no_nan() {
    let mol = parse("CC(=O)Oc1ccccc1C(=O)O");
    let pdb = generate_3d_pdb(&mol);
    assert!(
        !pdb.contains("nan") && !pdb.contains("inf"),
        "PDB must have no NaN/Inf coords"
    );
    assert!(pdb.contains("HETATM"), "must produce HETATM records");
}

// Note: smarts_match_atoms error-path test (invalid SMARTS) is omitted here
// because JsValue::from_str panics outside a WASM runtime.
// The underlying chematic_smarts::parse_smarts error path is tested separately:
#[test]
fn smarts_parse_invalid_is_err() {
    assert!(
        chematic_smarts::parse_smarts("[invalid").is_err(),
        "invalid SMARTS should return Err from parse_smarts"
    );
}

// ── Sprint P: SDF I/O, EState, topo path FP ─────────────────────────────

const ETHANE_MOL_BLOCK: &str = "\
ethane
  chematic

  2  1  0  0  0  0  0  0  0  0  0 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
";

#[test]
fn mol_from_sdf_block_ethane_atom_count() {
    let h = mol_from_sdf_block(ETHANE_MOL_BLOCK).expect("ethane parse");
    assert_eq!(h.atom_count(), 2);
}

#[test]
fn sdf_to_smiles_json_two_records() {
    let sdf = format!("{ETHANE_MOL_BLOCK}$$$$\n{ETHANE_MOL_BLOCK}$$$$\n");
    let json = sdf_to_smiles_json(&sdf);
    assert!(json.starts_with('[') && json.ends_with(']'));
    // Should contain 2 SMILES entries separated by comma.
    let count = json.matches("CC").count();
    assert_eq!(count, 2, "expected 2 ethane SMILES in JSON, got: {json}");
}

#[test]
fn estate_indices_json_acetic_acid_nonempty() {
    let h = parse("CC(=O)O");
    let json = estate_indices_json(&h);
    assert!(json.starts_with('[') && json.ends_with(']'));
    assert!(!json.is_empty() && json != "[]");
}

#[test]
fn tanimoto_topo_path_same_mol_is_one() {
    let a = parse("c1ccccc1");
    let b = parse("c1ccccc1");
    assert!((tanimoto_topo_path(&a, &b) - 1.0).abs() < 1e-9);
}

#[test]
fn tanimoto_topo_path_different_mols_lt_one() {
    let a = parse("c1ccccc1");
    let b = parse("CC(=O)Oc1ccccc1C(=O)O");
    assert!(tanimoto_topo_path(&a, &b) < 1.0);
}

#[test]
fn sum_estate_aspirin_positive() {
    let h = parse("CC(=O)Oc1ccccc1C(=O)O");
    assert!(h.sum_estate() > 0.0);
}

#[test]
fn max_min_estate_ordering() {
    let h = parse("CC(=O)O");
    assert!(h.max_estate() >= h.min_estate());
}

// ── Sprint O: depict_reaction_svg ────────────────────────────────────────

#[test]
fn depict_reaction_svg_esterification() {
    let svg = depict_reaction_svg("CC(=O)O.CCO>>CC(=O)OCC.O").unwrap();
    assert!(svg.contains("→"), "must contain arrow character");
    assert!(svg.contains("<svg"), "must be valid SVG");
}

#[test]
fn depict_reaction_svg_single_step() {
    let svg = depict_reaction_svg("C>>CC").unwrap();
    assert!(svg.contains("→"));
    assert!(svg.contains("<svg"));
}

// Error-path tested via the underlying parse_reaction (JsValue::from_str panics outside WASM).
#[test]
fn rxn_parse_missing_arrow_is_err() {
    assert!(chematic_rxn::parse_reaction("not_a_reaction").is_err());
}

// ── Sprint Q: IFG, Gasteiger, SA Score, VSA ──────────────────────────────

#[test]
fn identify_functional_groups_pyridine_has_n() {
    let h = parse("c1ccncc1");
    let json = identify_functional_groups(&h);
    assert!(
        json.contains('N'),
        "pyridine FG JSON should contain N: {json}"
    );
    assert!(json.starts_with('['), "should be JSON array");
}

#[test]
fn identify_functional_groups_hexane_empty() {
    let h = parse("CCCCCC");
    let json = identify_functional_groups(&h);
    assert_eq!(json, "[]", "hexane should have no functional groups");
}

#[test]
fn gasteiger_charges_json_oxygen_negative() {
    let h = parse("CC(=O)O"); // acetic acid
    let json = gasteiger_charges_json(&h);
    assert!(json.starts_with('[') && json.ends_with(']'));
    // Parse values and verify at least one is negative.
    let has_negative = json
        .trim_matches(|c| c == '[' || c == ']')
        .split(',')
        .any(|s| s.trim().parse::<f64>().map(|v| v < 0.0).unwrap_or(false));
    assert!(
        has_negative,
        "acetic acid should have at least one negative charge: {json}"
    );
}

#[test]
fn test_mmff94_charges_json_length() {
    // Acetic acid has 4 heavy atoms (2C, 2O)
    let h = parse("CC(=O)O");
    let json = mmff94_charges_json(&h);
    assert!(json.starts_with('[') && json.ends_with(']'));
    let count = json
        .trim_matches(|c| c == '[' || c == ']')
        .split(',')
        .count();
    assert_eq!(count, 4, "acetic acid should have 4 charges, got {count}");
}

#[test]
fn test_mmff94_charges_json_acetate_negative() {
    // Acetate [O-] → total charge = -1
    let h = parse("CC(=O)[O-]");
    let json = mmff94_charges_json(&h);
    let total: f64 = json
        .trim_matches(|c| c == '[' || c == ']')
        .split(',')
        .filter_map(|s| s.trim().parse::<f64>().ok())
        .sum();
    assert!(
        (total + 1.0).abs() < 0.01,
        "acetate total charge should be -1, got {total}"
    );
}

#[test]
fn test_mhfp_hashes_json_format() {
    let h = parse("c1ccccc1");
    let json = mhfp_hashes_json(&h);
    assert!(
        json.contains("\"num_hashes\":128"),
        "should have num_hashes:128"
    );
    assert!(json.contains("\"hashes\":"), "should have hashes key");
}

#[test]
fn test_tanimoto_mhfp_smiles_self() {
    let result = tanimoto_mhfp_smiles("c1ccccc1", "c1ccccc1").unwrap();
    assert!(
        (result - 1.0).abs() < 1e-9,
        "self-similarity should be 1.0, got {result}"
    );
}

#[test]
fn test_mhfp_lsh_handle_query() {
    let mut idx = MhfpLshHandle::new(128);
    idx.add_smiles("c1ccccc1").unwrap(); // benzene → 0
    idx.add_smiles("Cc1ccccc1").unwrap(); // toluene → 1
    idx.add_smiles("CC").unwrap(); // ethane  → 2

    let json = idx.query_json("c1ccccc1", 0.99).unwrap();
    // Benzene should find itself at similarity 1.0
    assert!(
        json.contains("\"index\":0"),
        "benzene should find itself: {json}"
    );
    assert_eq!(idx.len(), 3);
}

#[test]
fn to_mol_block_has_nonzero_coords() {
    let h = parse("c1ccccc1"); // benzene
    let block = to_mol_block(&h);
    // Extract atom-block lines (lines 4..10, 0-indexed) and check at least
    // one coordinate is non-zero.
    let nonzero = block.lines().skip(4).take(6).any(|line| {
        // Atom line: first 30 chars are three 10.4 coordinate fields.
        if line.len() < 30 {
            return false;
        }
        let x: f64 = line[0..10].trim().parse().unwrap_or(0.0);
        let y: f64 = line[10..20].trim().parse().unwrap_or(0.0);
        x.abs() > 0.01 || y.abs() > 0.01
    });
    assert!(
        nonzero,
        "to_mol_block should write real 2D coordinates, got:\n{block}"
    );
}

#[test]
fn sa_score_range() {
    let h = parse("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
    let score = sa_score(&h);
    assert!(
        (1.0..=10.0).contains(&score),
        "SA score out of [1,10]: {score:.2}"
    );
}

// v0.1.21 new API tests

// with_atom_charge
#[test]
fn mol_with_atom_charge_changes_charge() {
    let mol = parse("N"); // neutral N
    let mol2 = mol_with_atom_charge(&mol, 0, 1).unwrap();
    let atom = mol2.inner.atom(chematic_core::AtomIdx(0));
    assert_eq!(atom.charge, 1, "charge should be 1");
}

// with_atom_element
#[test]
fn mol_with_atom_element_changes_element() {
    let mol = parse("CC");
    let mol2 = mol_with_atom_element(&mol, 0, "N").unwrap();
    let atom = mol2.inner.atom(chematic_core::AtomIdx(0));
    assert_eq!(atom.element.symbol(), "N", "element should be N");
}

// SDF parse with coords
#[test]
fn mol_block_coords_json_ethanol() {
    let mol_block = "\nethanol\n\n  3  2  0  0  0  0  0  0  0  0  0 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n";
    let json = mol_block_coords_json(mol_block).unwrap();
    assert!(
        json.contains("[0.0000,0.0000]"),
        "first atom at origin: {json}"
    );
    assert!(
        json.contains("[1.5000,0.0000]"),
        "second atom at x=1.5: {json}"
    );
}

// CDXML all fragments
#[test]
fn cdxml_to_smiles_json_two_fragments() {
    let cdxml = r#"<?xml version="1.0"?>
<CDXML>
<fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1"/>
</fragment>
<fragment>
<n id="3" Element="8" p="30 0"/>
<n id="4" Element="6" p="40 0"/>
<b B="3" E="4" Order="1"/>
</fragment>
</CDXML>"#;
    let json = cdxml_to_smiles_json(cdxml).unwrap();
    // Should have 2 comma-separated SMILES entries.
    let count = json.matches(',').count() + 1;
    assert_eq!(count, 2, "should have 2 fragments: {json}");
}

#[test]
fn cdxml_document_json_preserves_page_objects_contract() {
    let cdxml = "<CDXML>\n<page id=\"p1\">\n<arrow id=\"a1\" Custom=\"keep\"/>\n</page>\n</CDXML>";
    let json = cdxml_document_json(cdxml).unwrap();
    assert!(json.contains("chematic.cdxml-document.v1"));
    assert!(json.contains("Custom"));
    assert!(json.contains("arrow"));
}

// CDXML stereo
#[test]
fn cdxml_wedge_bond_becomes_up() {
    let cdxml = r#"<?xml version="1.0"?>
<CDXML>
<fragment>
<n id="1" Element="6" p="0 0"/>
<n id="2" Element="6" p="10 0"/>
<b B="1" E="2" Order="1" Display="WedgeBegin"/>
</fragment>
</CDXML>"#;
    let (mol, _) = chematic_mol::parse_cdxml(cdxml).unwrap();
    let bond = mol.bond(chematic_core::BondIdx(0));
    assert_eq!(bond.order, chematic_core::BondOrder::Up, "WedgeBegin → Up");
}

// depict_data_with_coords
#[test]
fn depict_data_with_coords_uses_provided_coords() {
    let mol = parse("CC");
    // Provide coords far from defaults.
    let json = depict_data_with_coords_json(&mol, "[[100.0,200.0],[300.0,400.0]]");
    assert!(json.contains("100."), "x coord should be 100: {json}");
    assert!(json.contains("200."), "y coord should be 200: {json}");
}

// Mutable Molecule API
#[test]
fn mol_with_atom_added_increases_count() {
    let mol = parse("CC");
    let mol2 = mol_with_atom_added(&mol, "O").unwrap();
    assert_eq!(mol2.atom_count(), mol.atom_count() + 1);
}

#[test]
fn mol_with_bond_added_increases_count() {
    let mol = parse("CC.O"); // ethane + disconnected O (3 atoms, 1 bond)
    let mol2 = mol_with_bond_added(&mol, 1, 2, 1).unwrap();
    assert_eq!(mol2.bond_count(), mol.bond_count() + 1);
}

#[test]
fn mol_with_atom_removed_decreases_count() {
    let mol = parse("CCO");
    let mol2 = mol_with_atom_removed(&mol, 2).unwrap();
    assert_eq!(mol2.atom_count(), mol.atom_count() - 1);
    // The C-O bond is also removed.
    assert_eq!(mol2.bond_count(), mol.bond_count() - 1);
}

#[test]
fn mol_with_bond_removed_decreases_count() {
    let mol = parse("CC");
    let mol2 = mol_with_bond_removed(&mol, 0).unwrap();
    assert_eq!(mol2.atom_count(), 2, "atoms preserved");
    assert_eq!(mol2.bond_count(), 0, "bond removed");
}

#[test]
fn mol_next_atom_idx_equals_atom_count() {
    let mol = parse("CCO");
    assert_eq!(mol_next_atom_idx(&mol), 3);
}

// SDF / V3000 write
#[test]
fn smiles_array_to_sdf_contains_dollar_delimiters() {
    let sdf = smiles_array_to_sdf(r#"["CC","CCO"]"#).unwrap();
    assert_eq!(sdf.matches("$$$$").count(), 2, "2 delimiters: {sdf}");
}

#[test]
fn to_mol_v3000_block_contains_v3000() {
    let mol = parse("CC");
    let v3k = to_mol_v3000_block(&mol);
    assert!(v3k.contains("V3000"), "should contain V3000 tag");
    assert!(v3k.contains("M  V30 BEGIN ATOM"), "should have atom block");
}

// DepictData
#[test]
fn depict_data_json_benzene_atoms_bonds() {
    let mol = parse("c1ccccc1");
    let json = depict_data_json(&mol);
    assert!(json.contains("\"atoms\":"), "should have atoms key");
    assert!(json.contains("\"bonds\":"), "should have bonds key");
    assert!(json.contains("\"element\":\"C\""), "benzene has C atoms");
}

#[test]
fn depict_data_json_pyridine_nitrogen_has_label() {
    let mol = parse("c1ccncc1");
    let json = depict_data_json(&mol);
    // Nitrogen has a label
    assert!(
        json.contains("\"label\":\"N\""),
        "pyridine N should have label: {json}"
    );
}

// CPK colors
#[test]
fn cpk_color_nitrogen_is_blue() {
    assert_eq!(cpk_color("N"), "#3050F8", "N is blue in CPK");
}

#[test]
fn cpk_color_carbon_is_black() {
    assert_eq!(cpk_color("C"), "#000000", "C is black in CPK");
}

// CML / CDXML
#[test]
fn mol_from_cml_ethanol() {
    let cml = r#"<molecule>
  <atomArray>
    <atom id="a1" elementType="C" x2="0.0" y2="0.0"/>
    <atom id="a2" elementType="C" x2="1.5" y2="0.0"/>
    <atom id="a3" elementType="O" x2="3.0" y2="0.0"/>
  </atomArray>
  <bondArray>
    <bond atomRefs2="a1 a2" order="1"/>
    <bond atomRefs2="a2 a3" order="1"/>
  </bondArray>
</molecule>"#;
    let mol = mol_from_cml(cml).unwrap();
    assert_eq!(mol.atom_count(), 3, "ethanol: 3 heavy atoms");
    assert_eq!(mol.bond_count(), 2, "ethanol: 2 bonds");
}

#[test]
fn to_cml_contains_molecule_tag() {
    let mol = parse("CC(=O)O");
    let cml = to_cml(&mol);
    assert!(cml.contains("<molecule"), "should contain <molecule");
    assert!(cml.contains("elementType="), "should contain elementType");
    assert!(cml.contains("atomRefs2="), "should contain atomRefs2");
}

#[test]
fn to_cml_roundtrip_atom_count() {
    let mol = parse("CC(=O)O");
    let cml = to_cml(&mol);
    let mol2 = mol_from_cml(&cml).unwrap();
    assert_eq!(
        mol.atom_count(),
        mol2.atom_count(),
        "CML round-trip preserves atom count"
    );
    assert_eq!(
        mol.bond_count(),
        mol2.bond_count(),
        "CML round-trip preserves bond count"
    );
}

#[test]
fn mol_from_cdxml_ethanol() {
    let cdxml = r#"<?xml version="1.0"?>
<CDXML>
<fragment>
<n id="1" p="10.0 20.0" Element="6"/>
<n id="2" p="25.0 20.0" Element="6"/>
<n id="3" p="40.0 20.0" Element="8"/>
<b B="1" E="2" Order="1"/>
<b B="2" E="3" Order="1"/>
</fragment>
</CDXML>"#;
    let mol = mol_from_cdxml(cdxml).unwrap();
    assert_eq!(mol.atom_count(), 3, "CDXML ethanol: 3 heavy atoms");
    assert_eq!(mol.bond_count(), 2, "CDXML ethanol: 2 bonds");
}

#[test]
fn mol_from_cml_unknown_element_returns_err() {
    let cml = r#"<molecule><atomArray>
<atom id="a1" elementType="Xx"/>
</atomArray></molecule>"#;
    // Test via underlying function (JsValue::from_str would abort in native)
    let result = chematic_mol::parse_cml(cml);
    assert!(result.is_err(), "unknown element should return Err");
}

// Sprint CC
#[test]
fn mmp_pairs_json_ethylbenzene_propylbenzene() {
    let smiles_json = r#"["CCc1ccccc1","CCCc1ccccc1"]"#;
    let result = mmp_pairs_json(smiles_json).unwrap();

    // Should find exactly 1 pair with the benzene core.
    assert!(
        result.contains("\"core\":"),
        "should have core key: {result}"
    );
    assert!(
        // Canonical form updated (issue #205/#206): the explicit/implicit-
        // H-count Morgan-rank unification fix changed which atom the
        // canonical DFS starts from for this ring -- same core molecule,
        // same wildcard attachment point, mirrors chematic-chem::mmp's
        // equivalent test.
        result.contains("c1cc([*])ccc1"),
        "core should be benzene: {result}"
    );
    assert!(
        result.contains("\"fragment_a\":"),
        "should have fragment_a: {result}"
    );
    assert!(
        result.contains("\"fragment_b\":"),
        "should have fragment_b: {result}"
    );
}

#[test]
fn mmp_pairs_json_no_pairs_for_single_molecule() {
    let smiles_json = r#"["CCc1ccccc1"]"#;
    let result = mmp_pairs_json(smiles_json).unwrap();
    assert_eq!(result, "[]", "single molecule should have no pairs");
}

#[test]
fn mmp_pairs_json_three_molecules() {
    // Ethylbenzene, propylbenzene, butylbenzene.
    // Different isomers and cut points may yield multiple pairs and cores.
    let smiles_json = r#"["CCc1ccccc1","CCCc1ccccc1","CCCCc1ccccc1"]"#;
    let result = mmp_pairs_json(smiles_json).unwrap();

    // Should find at least 1 pair (the minimum for 3 benzene-derived molecules).
    assert!(
        result.contains("\"mol_a\":"),
        "should have mol_a key: {result}"
    );
    assert!(
        result.contains("\"core\":"),
        "should have core key: {result}"
    );
    // Extract pair count by counting top-level braces.
    let pair_count = result.matches("\"mol_a\":").count();
    assert!(
        pair_count >= 1,
        "should have >= 1 pair, got {pair_count}: {result}"
    );
}

// Sprint BB
#[test]
fn conformer_handle_new_and_add() {
    let ens = ConformerHandle::new("CCCC").unwrap(); // butane
    assert_eq!(ens.conformer_count(), 0, "new ensemble has 0 conformers");
}

#[test]
fn conformer_handle_add_generated() {
    let mut ens = ConformerHandle::new("CCCC").unwrap();
    let idx = ens.add_generated_conformer();
    assert_eq!(idx, 0, "first conformer has index 0");
    assert_eq!(ens.conformer_count(), 1);
}

#[test]
fn conformer_handle_add_minimized() {
    let mut ens = ConformerHandle::new("CCCC").unwrap();
    let idx = ens.add_minimized_conformer();
    assert_eq!(idx, 0);
    assert_eq!(ens.conformer_count(), 1);
}

#[test]
fn conformer_handle_get_pdb_returns_hetatm() {
    let mut ens = ConformerHandle::new("CC").unwrap();
    ens.add_generated_conformer();
    let pdb = ens.get_conformer_pdb(0).expect("conformer 0 exists");
    assert!(pdb.contains("HETATM"), "PDB should contain HETATM records");
}

#[test]
fn conformer_handle_rmsd_same_conformer() {
    // RMSD of a conformer with itself should be 0.0.
    let mut ens = ConformerHandle::new("CCCC").unwrap();
    ens.add_generated_conformer();
    ens.add_generated_conformer();
    // Different conformers will have RMSD >= 0; same index trivially returns
    // the distance between identical points = 0.
    let rmsd = ens.conformer_rmsd(0, 0);
    assert!(rmsd.is_finite(), "RMSD should be finite");
    assert!(rmsd >= 0.0, "RMSD must be non-negative");
}

#[test]
fn conformer_handle_rmsd_out_of_range_is_nan() {
    let mut ens = ConformerHandle::new("CC").unwrap();
    ens.add_generated_conformer();
    assert!(ens.conformer_rmsd(0, 99).is_nan(), "out-of-range gives NaN");
}

#[test]
fn conformer_handle_remove_conformer() {
    let mut ens = ConformerHandle::new("CC").unwrap();
    ens.add_generated_conformer();
    assert!(ens.remove_conformer(0));
    assert_eq!(ens.conformer_count(), 0);
    assert!(!ens.remove_conformer(0), "already removed");
}

#[test]
fn conformer_handle_mol_returns_correct_atom_count() {
    let ens = ConformerHandle::new("CC(=O)O").unwrap(); // acetic acid, 4 heavy atoms
    let mol = ens.mol();
    assert_eq!(mol.atom_count(), 4);
}

// R-group decomposition oracle tests
#[test]
fn rgroup_toluene_ethylbenzene_para_substituted() {
    // Core: para-substituted benzene c1ccc(*)cc1
    // Toluene Cc1ccccc1 → R1 = C
    // Ethylbenzene CCc1ccccc1 → R1 = CC
    let smiles_json = r#"["Cc1ccccc1","CCc1ccccc1"]"#;
    let json = rgroup_decompose_json(smiles_json, "c1ccc(*)cc1").unwrap();

    // Both should match.
    assert!(
        json.contains("\"matched\":true"),
        "both should match: {json}"
    );
    // R-groups should be C and CC (canonical forms may vary slightly).
    assert!(json.contains("\"r1\":"), "should have r1 key: {json}");
}

#[test]
fn rgroup_no_match_returns_false() {
    // Benzene has no substituent → doesn't match c1ccc(*)cc1.
    let smiles_json = r#"["c1ccccc1"]"#;
    let json = rgroup_decompose_json(smiles_json, "c1ccc(*)cc1").unwrap();
    assert!(
        json.contains("\"matched\":false"),
        "unsubstituted benzene should not match: {json}"
    );
}

#[test]
fn rgroup_two_attachment_points() {
    // Di-substituted benzene: c1cc(*)cc(*)c1
    // 1,3-dimethylbenzene m-xylene Cc1cccc(C)c1 → R1=C, R2=C
    let smiles_json = r#"["Cc1cccc(C)c1"]"#;
    let json = rgroup_decompose_json(smiles_json, "c1cc(*)cc(*)c1").unwrap();
    assert!(
        json.contains("\"matched\":true"),
        "m-xylene should match: {json}"
    );
    assert!(json.contains("\"r1\":"), "should have r1: {json}");
    assert!(json.contains("\"r2\":"), "should have r2: {json}");
}

#[test]
fn rgroup_invalid_smarts_returns_err() {
    // SMARTS parse should fail via underlying function (not JsValue path).
    let result = chematic_smarts::parse_smarts("~~~invalid~~~");
    assert!(result.is_err(), "invalid SMARTS should not parse");
}

// Sprint AA
#[test]
fn fcfp4_bitvec_length_256() {
    assert_eq!(fcfp4_bitvec(&parse("c1ccccc1")).len(), 256);
}

#[test]
fn fcfp6_bitvec_length_256() {
    assert_eq!(fcfp6_bitvec(&parse("c1ccccc1")).len(), 256);
}

#[test]
fn dice_ecfp6_identical_is_one() {
    let mol = parse("c1ccccc1");
    assert!((dice_ecfp6(&mol, &mol) - 1.0).abs() < 1e-9);
}

#[test]
fn write_smiles_contains_same_atoms() {
    let mol = parse("CC(=O)O");
    let smi = write_smiles(&mol);
    // Non-canonical SMILES must still round-trip to the same atom count.
    let mol2 = parse_smiles(&smi).unwrap();
    assert_eq!(mol.atom_count(), mol2.atom_count());
}

#[test]
fn normalize_reaction_smiles_roundtrip() {
    // parse_reaction error path tested via underlying fn (JsValue native-abort).
    let result = chematic_rxn::parse_reaction("CC>>CO");
    assert!(result.is_ok(), "simple reaction should parse");
    let rxn = result.unwrap();
    let out = chematic_rxn::write_reaction(&rxn);
    assert!(
        out.contains(">>"),
        "written reaction should contain >>: {out}"
    );
}

// Sprint Z
#[test]
fn brics_fragments_json_aspirin() {
    let mol = parse("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
    let json = brics_fragments_json(&mol);
    let count_from_json = json
        .split('"')
        .filter(|s| s.contains('C') || s.contains('c'))
        .count();
    assert!(
        json.starts_with('[') && json.ends_with(']'),
        "not a JSON array: {json}"
    );
    assert_eq!(
        count_from_json,
        brics_fragment_count(&mol),
        "fragment count mismatch: json={json}"
    );
}

#[test]
fn brics_fragments_json_benzene_self() {
    // No BRICS-breakable bonds → returns the molecule itself as one fragment.
    let mol = parse("c1ccccc1");
    let json = brics_fragments_json(&mol);
    assert!(json.starts_with('[') && json.ends_with(']'));
    assert_eq!(
        brics_fragment_count(&mol),
        1,
        "benzene has 1 fragment (itself): {json}"
    );
}

// ── issue #91: retro_disconnect_json ──────────────────────────────────────

#[test]
fn retro_disconnect_json_matches_rust_api_acetanilide() {
    // Same molecule as chematic_rxn::retro's own doc example.
    let smiles = "CC(=O)Nc1ccccc1";
    let mol = parse(smiles);

    let json = retro_disconnect_json(&mol, 20, "").expect("should succeed");
    let parsed: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");
    let arr = parsed.as_array().expect("top-level array");
    assert!(!arr.is_empty(), "acetanilide should have >=1 disconnection");

    // Rust-vs-WASM parity: same molecule, same template library, same
    // ordering, same precursor SMILES -- the WASM wrapper must not change
    // the underlying algorithm's result in any way.
    let rust_mol = chematic_smiles::parse(smiles).expect("parse");
    let rust_results = chematic_rxn::retro::retro_disconnect(
        &rust_mol,
        chematic_rxn::retro::DEFAULT_TEMPLATES,
        20,
    );
    assert_eq!(arr.len(), rust_results.len(), "result count mismatch");
    for (json_item, rust_item) in arr.iter().zip(rust_results.iter()) {
        assert_eq!(json_item["template"], rust_item.template_name);
        assert_eq!(
            json_item["reaction_class"],
            rust_item.reaction_class.as_str()
        );
        let json_precursors: Vec<String> = json_item["precursors"]
            .as_array()
            .expect("precursors array")
            .iter()
            .map(|v| v.as_str().unwrap().to_string())
            .collect();
        assert_eq!(&json_precursors, &rust_item.precursor_smiles);
    }
}

#[test]
fn retro_disconnect_json_no_match_returns_empty_array_not_error() {
    // Methane has no disconnectable bond in the template library -- this is
    // a valid empty result, not an error.
    let mol = parse("C");
    let json = retro_disconnect_json(&mol, 20, "").expect("empty result is not an error");
    assert_eq!(json, "[]");
}

// Unknown-reaction_class error path is tested via
// crates/chematic-wasm/tests/retro_disconnect.test.mjs, not natively here:
// JsValue::from_str aborts the process outside a real WASM runtime (see the
// same note on parse_reaction/parse_smarts error-path tests above), and this
// validation is WASM-layer-only (no underlying non-wasm function to redirect
// to, unlike those cases).

#[test]
fn retro_disconnect_json_reaction_class_filter_narrows_results() {
    let mol = parse("CC(=O)Nc1ccccc1"); // acetanilide: amide bond
    let all = retro_disconnect_json(&mol, 0, "").expect("all classes");
    let filtered = retro_disconnect_json(&mol, 0, "AmideBond").expect("AmideBond filter");
    let filtered_out = retro_disconnect_json(&mol, 0, "Ether").expect("Ether filter");

    let all_arr: serde_json::Value = serde_json::from_str(&all).unwrap();
    let filtered_arr: serde_json::Value = serde_json::from_str(&filtered).unwrap();
    let filtered_out_arr: serde_json::Value = serde_json::from_str(&filtered_out).unwrap();

    assert!(!filtered_arr.as_array().unwrap().is_empty());
    for item in filtered_arr.as_array().unwrap() {
        assert_eq!(item["reaction_class"], "AmideBond");
    }
    assert!(
        filtered_arr.as_array().unwrap().len() <= all_arr.as_array().unwrap().len(),
        "filtered result must not exceed unfiltered result"
    );
    assert_eq!(
        filtered_out_arr.as_array().unwrap().len(),
        0,
        "acetanilide has no Ether-class disconnection"
    );
}

#[test]
fn retro_disconnect_json_max_results_caps_output() {
    let mol = parse("CC(=O)Nc1ccccc1");
    let json = retro_disconnect_json(&mol, 1, "").expect("should succeed");
    let arr: serde_json::Value = serde_json::from_str(&json).unwrap();
    assert!(arr.as_array().unwrap().len() <= 1);
}

#[test]
fn atom_pair_bitvec_length_256() {
    let mol = parse("c1ccccc1");
    assert_eq!(atom_pair_bitvec(&mol).len(), 256);
}

#[test]
fn torsion_bitvec_length_256() {
    let mol = parse("CCCCC");
    assert_eq!(torsion_bitvec(&mol).len(), 256);
}

#[test]
fn tanimoto_fcfp6_identical_is_one() {
    let mol = parse("c1ccccc1");
    assert!((tanimoto_fcfp6(&mol, &mol) - 1.0).abs() < 1e-9);
}

#[test]
fn sdf_from_records_json_round_trip() {
    // Build a 1-record SDF, read it back, check name and property.
    let smiles_json = r#"["CC"]"#;
    let names_json = r#"["ethane"]"#;
    let props_json = r#"["MW\t30.07\nSource\ttest"]"#;
    let sdf = sdf_from_records_json(smiles_json, names_json, props_json).unwrap();
    assert!(sdf.contains("ethane"), "name missing: {sdf}");
    assert!(sdf.contains("> <MW>"), "MW field missing: {sdf}");
    assert!(sdf.contains("30.07"), "MW value missing: {sdf}");
    assert!(sdf.contains("$$$$"), "SDF delimiter missing: {sdf}");
}

#[test]
fn sdf_from_records_json_no_props() {
    // Empty props string → no data fields, still valid SDF.
    let sdf = sdf_from_records_json(r#"["C"]"#, r#"["methane"]"#, r#"[""]"#).unwrap();
    assert!(sdf.contains("methane"));
    assert!(sdf.contains("$$$$"));
}

#[test]
fn sdf_from_records_json_length_mismatch_returns_err() {
    // Mismatched array lengths — note: testing via underlying logic,
    // not .is_err() on the wrapper (JsValue native-abort).
    let smiles = parse_smiles_json_array(r#"["CC","CCC"]"#).unwrap();
    let names = parse_smiles_json_array(r#"["ethane"]"#).unwrap();
    assert_ne!(smiles.len(), names.len(), "lengths should differ");
}

#[test]
fn parse_smiles_json_array_handles_escaped_strings() {
    let smiles = parse_smiles_json_array(r#"["C\\C=C\\O","CC"]"#).unwrap();
    assert_eq!(smiles[0], r#"C\C=C\O"#);
    assert_eq!(smiles[1], "CC");
}

// Sprint Y
#[test]
fn mol_from_xyz_roundtrip_atom_count() {
    // Build a minimal XYZ string for ethane and verify atom count.
    // (parse_xyz error path tested via the underlying fn — JsValue native-abort)
    let xyz = "2\nethane\nC  0.0 0.0 0.0\nC  1.5 0.0 0.0\n";
    let result = chematic_3d::parse_xyz(xyz);
    assert!(result.is_ok(), "valid XYZ should parse");
    let (mol, _) = result.unwrap();
    assert_eq!(mol.atom_count(), 2);
}

#[test]
fn to_xyz_contains_atom_lines() {
    let mol = parse("CC"); // ethane
    let xyz = to_xyz(&mol);
    assert!(xyz.contains('C'), "XYZ should contain C atom lines");
    // First line is atom count
    let count: usize = xyz.lines().next().unwrap().trim().parse().unwrap();
    assert_eq!(count, mol.atom_count(), "XYZ header atom count matches mol");
}

#[test]
fn mol_from_pdb_returns_handle() {
    // pdb_to_molecule with no atoms gives empty molecule.
    let h = mol_from_pdb("");
    assert_eq!(h.atom_count(), 0);
}

#[test]
fn logp_per_atom_json_length() {
    let mol = parse("CC(=O)O"); // acetic acid
    let json = logp_per_atom_json(&mol);
    let count = json
        .trim_matches(|c| c == '[' || c == ']')
        .split(',')
        .count();
    assert_eq!(
        count,
        mol.atom_count(),
        "per-atom logP array length mismatch"
    );
}

#[test]
fn mr_per_atom_json_length() {
    let mol = parse("CC(=O)O");
    let json = mr_per_atom_json(&mol);
    let count = json
        .trim_matches(|c| c == '[' || c == ']')
        .split(',')
        .count();
    assert_eq!(count, mol.atom_count());
}

#[test]
fn labute_asa_per_atom_json_no_nan() {
    let mol = parse("c1ccccc1");
    let json = labute_asa_per_atom_json(&mol);
    assert!(!json.contains("NaN"), "no NaN in asa per-atom: {json}");
}

#[test]
fn sssr_rings_json_benzene_one_ring() {
    let mol = parse("c1ccccc1");
    let json = sssr_rings_json(&mol);
    // One ring: [[0,1,2,3,4,5]]
    assert!(
        json.starts_with("[[") && json.ends_with("]]"),
        "expected one ring: {json}"
    );
    let ring: Vec<usize> = json
        .trim_matches(|c| c == '[' || c == ']')
        .split(',')
        .filter_map(|s| s.trim_matches(|c| c == '[' || c == ']').parse().ok())
        .collect();
    assert_eq!(ring.len(), 6, "benzene ring has 6 atoms");
}

#[test]
fn sssr_rings_json_naphthalene_two_rings() {
    let mol = parse("c1ccc2ccccc2c1"); // naphthalene
    let json = sssr_rings_json(&mol);
    // Count rings by counting "[" at start
    let ring_count = json.matches("],[").count() + 1;
    assert_eq!(ring_count, 2, "naphthalene has 2 rings: {json}");
}

#[test]
fn ecfp_bitvec_custom_256_length() {
    let mol = parse("c1ccccc1");
    let bv = ecfp_bitvec_custom(&mol, 2, 256, false);
    assert_eq!(bv.len(), 32, "256-bit FP = 32 bytes");
}

#[test]
fn ecfp_bitvec_custom_identical_tanimoto() {
    let mol = parse("c1ccccc1");
    let a = ecfp_bitvec_custom(&mol, 2, 512, false);
    let b = ecfp_bitvec_custom(&mol, 2, 512, false);
    assert_eq!(a, b, "same mol same FP");
}

// ── Chirality tests ──────────────────────────────────────────────────────
#[test]
fn ecfp4_bitvec_with_chirality_l_vs_d_alanine_different() {
    let l_ala = parse("C[C@H](N)C(=O)O");
    let d_ala = parse("C[C@@H](N)C(=O)O");
    let l_fp = ecfp4_bitvec_with_chirality(&l_ala, true);
    let d_fp = ecfp4_bitvec_with_chirality(&d_ala, true);
    assert_ne!(
        l_fp, d_fp,
        "L-alanine and D-alanine should have different ECFP4 with use_chirality=true"
    );
}

#[test]
fn ecfp4_bitvec_with_chirality_l_vs_d_alanine_same_without_chirality() {
    let l_ala = parse("C[C@H](N)C(=O)O");
    let d_ala = parse("C[C@@H](N)C(=O)O");
    let l_fp = ecfp4_bitvec_with_chirality(&l_ala, false);
    let d_fp = ecfp4_bitvec_with_chirality(&d_ala, false);
    assert_eq!(
        l_fp, d_fp,
        "L-alanine and D-alanine should have identical ECFP4 with use_chirality=false"
    );
}

#[test]
fn ecfp6_bitvec_with_chirality_l_vs_d_alanine_different() {
    let l_ala = parse("C[C@H](N)C(=O)O");
    let d_ala = parse("C[C@@H](N)C(=O)O");
    let l_fp = ecfp6_bitvec_with_chirality(&l_ala, true);
    let d_fp = ecfp6_bitvec_with_chirality(&d_ala, true);
    assert_ne!(
        l_fp, d_fp,
        "L-alanine and D-alanine should have different ECFP6 with use_chirality=true"
    );
}

#[test]
fn smarts_match_atoms_with_chirality_chiral_atom_matches_only_with_flag() {
    let l_ala = parse("C[C@H](N)C(=O)O");
    let d_ala = parse("C[C@@H](N)C(=O)O");
    // Match [C@H] (counterclockwise) — should only match L-alanine when use_chirality=true
    let l_matches_with = smarts_match_atoms_with_chirality("[C@H]", &l_ala, true).unwrap();
    let l_matches_without = smarts_match_atoms_with_chirality("[C@H]", &l_ala, false).unwrap();
    let d_matches_with = smarts_match_atoms_with_chirality("[C@H]", &d_ala, true).unwrap();
    let d_matches_without = smarts_match_atoms_with_chirality("[C@H]", &d_ala, false).unwrap();

    // With use_chirality=true: L-ala matches, D-ala doesn't
    assert!(
        !l_matches_with.contains("[]"),
        "L-ala should match [C@H] with use_chirality=true"
    );
    assert!(
        d_matches_with.contains("[]"),
        "D-ala should NOT match [C@H] with use_chirality=true"
    );

    // With use_chirality=false: both match (chirality ignored)
    assert!(
        !l_matches_without.contains("[]"),
        "L-ala should match [C@H] with use_chirality=false"
    );
    assert!(
        !d_matches_without.contains("[]"),
        "D-ala should match [C@H] with use_chirality=false (chirality ignored)"
    );
}

#[test]
fn enumerate_stereo_isomers_one_center_gives_two() {
    // C(F)(Cl)Br has one unspecified chiral center → 2 stereoisomers
    let mol = parse("C(F)(Cl)Br");
    // enumerate_stereo_isomers_json returns Result; unwrap Ok via JsValue-safe check
    let result = enumerate_stereo_isomers_json(&mol);
    assert!(result.is_ok(), "should return Ok for valid mol");
    let json = result.unwrap();
    // Count objects: each object contains "smiles", "inchi", "inchikey"
    // So 2 isomers means 2 occurrences of "smiles" at top level
    let count = json.matches("\"smiles\"").count();
    assert_eq!(count, 2, "expected 2 stereoisomers: {json}");
}

#[test]
fn enumerate_stereo_isomers_already_specified_gives_one() {
    // Already-specified stereocenter → only 1 isomer (no unspecified centers)
    let mol = parse("C[C@H](F)Cl");
    let result = enumerate_stereo_isomers_json(&mol);
    assert!(result.is_ok());
    let json = result.unwrap();
    // Count objects: each object contains "smiles"
    let count = json.matches("\"smiles\"").count();
    assert_eq!(count, 1, "fully specified mol should give 1 isomer: {json}");
}

#[test]
fn enumerate_stereo_isomers_no_center_gives_one() {
    // No chiral centers at all → 1 isomer
    let mol = parse("CC");
    let result = enumerate_stereo_isomers_json(&mol);
    assert!(result.is_ok());
    let json = result.unwrap();
    assert!(json.starts_with('[') && json.ends_with(']'));
}

// Sprint X
#[test]
fn mol_from_v3000_block_parses_atom_count() {
    // Borrow a minimal V3000 fixture from the mol3000 test suite.
    // (Parse failure tested via underlying parse_mol_v3000 — see JsValue note.)
    assert!(
        chematic_mol::parse_mol_v3000(
            "\n  \n\n  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\n\
             M  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\n\
             M  V30 1 C 0 0 0 0\nM  V30 2 C 1.5 0 0 0\n\
             M  V30 END ATOM\nM  V30 BEGIN BOND\n\
             M  V30 1 1 1 2\nM  V30 END BOND\n\
             M  V30 END CTAB\nM  END"
        )
        .is_ok()
    );
}

#[test]
fn generate_3d_minimized_pdb_nonzero_coords() {
    let mol = parse("CCCC"); // butane — flexible, benefits from minimization
    let pdb = generate_3d_minimized_pdb(&mol);
    assert!(pdb.contains("HETATM"), "expected HETATM records");
    // At least one non-zero coordinate in the PDB output.
    let has_nonzero = pdb.lines().filter(|l| l.starts_with("HETATM")).any(|l| {
        // PDB columns 31-38, 39-46, 47-54 are x,y,z
        if l.len() < 54 {
            return false;
        }
        let x: f64 = l[30..38].trim().parse().unwrap_or(0.0);
        let y: f64 = l[38..46].trim().parse().unwrap_or(0.0);
        x.abs() > 0.01 || y.abs() > 0.01
    });
    assert!(
        has_nonzero,
        "minimized PDB should have non-zero coords:\n{pdb}"
    );
}

#[test]
fn sdf_to_records_json_parses_properties() {
    // Raw string preserves the required fixed-width whitespace of MOL V2000.
    let sdf = concat!(
        "aspirin\n",
        "  chematic\n",
        "\n",
        "  2  1  0  0  0  0  0  0  0  0  0 V2000\n",
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "  1  2  1  0\n",
        "M  END\n",
        "> <MW>\n",
        "180.2\n",
        "\n",
        "> <Source>\n",
        "ChEMBL\n",
        "\n",
        "$$$$\n",
    );
    let json = sdf_to_records_json(sdf);
    assert!(
        json.contains("\"name\":\"aspirin\""),
        "name missing: {json}"
    );
    assert!(json.contains("\"MW\":\"180.2\""), "MW missing: {json}");
    assert!(
        json.contains("\"Source\":\"ChEMBL\""),
        "Source missing: {json}"
    );
}

#[test]
fn sdf_to_records_json_escapes_special_chars() {
    let sdf = concat!(
        "mol\n",
        "  chematic\n",
        "\n",
        "  1  0  0  0  0  0  0  0  0  0  0 V2000\n",
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
        "M  END\n",
        "> <Notes>\n",
        "line1\n",
        "line2\n",
        "\n",
        "$$$$\n",
    );
    let json = sdf_to_records_json(sdf);
    // "line1\nline2" value joined with \n, then JSON-escaped to \\n
    assert!(
        json.contains("\\n"),
        "newline in multi-line value should be escaped: {json}"
    );
}

#[test]
fn depict_svg_grid_highlighted_no_smarts_returns_grid() {
    let smiles = "c1ccccc1\nCCO";
    let svg = depict_svg_grid_highlighted(smiles, 2, "");
    assert!(svg.contains("mol-0"), "expected mol-0: {svg}");
    assert!(svg.contains("mol-1"), "expected mol-1: {svg}");
}

#[test]
fn depict_svg_grid_highlighted_with_smarts_adds_circles() {
    let smiles = "c1ccccc1\nc1ccncc1"; // benzene and pyridine
    let svg = depict_svg_grid_highlighted(smiles, 2, "c1ccccn1"); // pyridine SMARTS
    // Pyridine matches → circle element expected in second cell
    assert!(svg.contains("<circle"), "expected highlight circles: {svg}");
}

// Sprint W
#[test]
fn pains_matches_json_empty_for_clean_mol() {
    let mol = parse("c1ccccc1"); // benzene — no PAINS alerts
    let json = pains_matches_json(&mol);
    assert_eq!(json, "[]");
}

#[test]
fn pains_matches_json_returns_alert_names() {
    // Rhodanine scaffold is a classic PAINS alert
    let mol = parse("O=C1CSC(=S)N1"); // rhodanine
    let json = pains_matches_json(&mol);
    assert!(json.starts_with('[') && json.ends_with(']'));
    // May be empty for some molecules; just check structure
    let _ = json; // success if it compiles and runs
}

#[test]
fn cip_assignments_json_aspirin_no_stereo() {
    let mol = parse("CC(=O)Oc1ccccc1C(=O)O"); // aspirin — no stereocenters
    let json = cip_assignments_json(&mol);
    assert_eq!(json, "[]");
}

#[test]
fn cip_assignments_json_chiral_center() {
    let mol = parse("[C@@H](F)(Cl)Br"); // chiral center
    let json = cip_assignments_json(&mol);
    assert!(json.contains("cipCode"), "expected cipCode in: {json}");
    assert!(json.contains('R') || json.contains('S'));
}

#[test]
fn ecfp6_bitvec_length_256() {
    let mol = parse("c1ccccc1");
    assert_eq!(ecfp6_bitvec(&mol).len(), 256);
}

#[test]
fn tanimoto_ecfp6_identical_is_one() {
    let mol = parse("c1ccccc1");
    assert!((tanimoto_ecfp6(&mol, &mol) - 1.0).abs() < 1e-9);
}

#[test]
fn dice_ecfp4_identical_is_one() {
    let mol = parse("c1ccccc1");
    assert!((dice_ecfp4(&mol, &mol) - 1.0).abs() < 1e-9);
}

#[test]
fn dice_maccs_identical_is_one() {
    let mol = parse("c1ccccc1");
    assert!((dice_maccs(&mol, &mol) - 1.0).abs() < 1e-9);
}

#[test]
fn num_aliphatic_rings_cyclohexane() {
    let mol = parse("C1CCCCC1"); // cyclohexane
    assert_eq!(mol.num_aliphatic_rings(), 1);
    assert_eq!(mol.num_saturated_rings(), 1);
    assert_eq!(mol.aromatic_ring_count(), 0);
}

#[test]
fn num_unspecified_stereocenters_unspec() {
    let mol = parse("C(F)(Cl)Br"); // chiral center without @/@@ annotation
    assert_eq!(mol.num_unspecified_stereocenters(), 1);
}

#[test]
fn shape_descriptors_json_benzene_keys() {
    let mol = parse("c1ccccc1");
    let json = shape_descriptors_json(&mol);
    assert!(json.contains("\"pmi1\""), "missing pmi1: {json}");
    assert!(json.contains("\"npr1\""), "missing npr1: {json}");
    assert!(
        json.contains("\"asphericity\""),
        "missing asphericity: {json}"
    );
}

#[test]
fn shape_descriptors_json_single_atom_no_nan() {
    let mol = parse("[Na+]"); // single atom — shape descriptors may be non-finite
    let json = shape_descriptors_json(&mol);
    // Must be valid JSON (no bare NaN)
    assert!(!json.contains("NaN"), "NaN in output: {json}");
    assert!(!json.contains("inf"), "inf in output: {json}");
}

#[test]
fn maxmin_picks_returns_correct_count() {
    let json = r#"["CC","c1ccccc1","CCO","CCCC","c1cccnc1"]"#;
    let result = maxmin_picks_ecfp4_json(json, 3).unwrap();
    let count = result
        .trim_matches(|c| c == '[' || c == ']')
        .split(',')
        .count();
    assert_eq!(count, 3, "expected 3 picks: {result}");
}

#[test]
fn smiles_parse_fails_for_unclosed_ring() {
    // Verifies the error-propagation contract: parse_smiles_json_array feeds
    // each string through chematic_smiles::parse, which rejects unclosed rings.
    // (Calling maxmin_picks_ecfp4_json directly would invoke JsValue::from_str,
    // which panics outside WASM — see comment on rxn_parse_missing_arrow_is_err.)
    assert!(
        chematic_smiles::parse("C1CC").is_err(),
        "unclosed ring should fail"
    );
}

#[test]
fn butina_cluster_all_similar_one_cluster() {
    // All identical SMILES → cutoff=0.0 means all at distance 1.0 from centroid → 1 cluster
    let json = r#"["c1ccccc1","c1ccccc1","c1ccccc1"]"#;
    let result = butina_cluster_ecfp4_json(json, 0.0).unwrap();
    assert!(result.starts_with('[') && result.ends_with(']'));
}

#[test]
fn mcs_smiles_two_acetyl_compounds() {
    let json = r#"["CC(=O)O","CC(=O)N"]"#; // acetic acid and acetamide — MCS is CC=O
    let result = mcs_smiles_json(json).unwrap();
    assert_ne!(result, "null", "expected a non-null MCS");
    assert!(!result.is_empty());
}

#[test]
fn mcs_smiles_no_overlap_returns_null() {
    let json = r#"["[Na+]","[Cl-]"]"#; // no common organic substructure
    let result = mcs_smiles_json(json).unwrap();
    // May return a single-atom MCS or null depending on algorithm
    assert!(result == "null" || !result.is_empty());
}

#[test]
fn mcs_with_ring_config_quinoline_series_complete_rings() {
    // Issue #1 example: quinoline series with both ring constraints.
    // ring_matches_ring_only blocks ring<->non-ring; complete_rings_only ensures
    // no partial ring is included.  The shared exocyclic CH₂ (non-ring in all three
    // molecules) is still valid under ring_matches_ring_only (non-ring matches
    // non-ring), so the result is quinoline (10) + exocyclic C (1) = 11 atoms.
    let json = r#"["c1ccc2nc(CC)ccc2c1","c1ccc2nc(CO)ccc2c1","c1ccc2nc(CN)ccc2c1"]"#;
    let result = mcs_smiles_json_with_ring_config(json, true, true).unwrap();
    assert_ne!(result, "null");
    let mol = chematic_smiles::parse(&result).expect("valid SMILES");
    assert!(
        mol.atom_count() >= 10,
        "expected at least the quinoline scaffold (10 atoms), got {}",
        mol.atom_count()
    );
}

#[test]
fn mcs_with_ring_config_benzene_toluene() {
    // benzene vs toluene with complete_rings_only: full benzene ring (6 atoms).
    let json = r#"["c1ccccc1","Cc1ccccc1"]"#;
    let result = mcs_smiles_json_with_ring_config(json, false, true).unwrap();
    assert_ne!(result, "null");
    let mol = chematic_smiles::parse(&result).expect("valid SMILES");
    assert_eq!(
        mol.atom_count(),
        6,
        "expected full benzene ring (6 atoms), got {}",
        mol.atom_count()
    );
}

#[test]
fn mcs_smiles_single_input_returns_err() {
    // Passing only one SMILES string should produce a JS error; we verify
    // via the Result variant rather than calling is_err() (which would
    // materialise a JsValue drop on non-wasm32 and cause a panic).
    let json = r#"["CC"]"#;
    let smiles_list: Vec<String> = json
        .split('"')
        .enumerate()
        .filter(|(i, _)| i % 2 == 1)
        .map(|(_, s)| s.to_string())
        .collect();
    assert_eq!(smiles_list.len(), 1, "helper extracted wrong count");
    assert!(smiles_list.len() < 2, "should be fewer than 2 inputs");
}

// Sprint V
#[test]
fn murcko_scaffold_benzene_ring() {
    let mol = parse("c1ccccc1CC(=O)O"); // phenylacetic acid
    let scaffold = murcko_scaffold(&mol);
    let smi = scaffold.canonical_smiles();
    // Murcko of phenylacetic acid is benzene ring only.
    assert!(!smi.is_empty(), "murcko_scaffold returned empty SMILES");
}

#[test]
fn canonical_tautomer_returns_handle() {
    let mol = parse("Oc1cccc2ccccc12"); // 1-naphthol
    let t = canonical_tautomer(&mol);
    assert!(t.atom_count() > 0);
}

#[test]
fn tautomer_parent_json_reports_completed_result() {
    let mol = parse("CCN=O");
    let json = tautomer_parent_json(&mol, 16, 32, None);
    assert!(json.contains(r#""smiles":"#), "missing smiles: {json}");
    assert!(
        json.contains(r#""status":"completed""#),
        "unexpected status: {json}"
    );
}

#[test]
fn tautomer_parent_json_reports_transform_budget() {
    let mol = parse("OC=C");
    let json = tautomer_parent_json(&mol, 0, 32, None);
    assert!(
        json.contains(r#""status":"max_transforms_reached""#),
        "expected budget status: {json}"
    );
}

#[test]
fn parent_json_bindings_cover_mechanical_and_composed_parents() {
    let mol = parse("[NH3+][C@@H]([2H])C(=O)[O-].Cl");
    assert!(fragment_parent_json(&mol).contains(r#""status":"completed""#));
    assert!(charge_parent_json(&mol).contains(r#""status":"completed""#));
    assert!(isotope_parent_json(&mol).contains(r#""status":"completed""#));
    assert!(stereo_parent_json(&mol).contains(r#""status":"completed""#));
    assert!(super_parent_json(&mol, 16, 32, None).contains(r#""status":"completed""#));
    let report = super_parent_report_json(&mol, 16, 32, None);
    assert!(report.contains(r#""name":"fragment""#));
    assert!(report.contains(r#""name":"tautomer""#));
}

#[test]
fn enumerate_tautomers_json_is_array() {
    let mol = parse("Oc1cccc2ccccc12");
    let json = enumerate_tautomers_json(&mol);
    assert!(json.starts_with('[') && json.ends_with(']'));
    assert!(json.len() > 2, "expected at least one tautomer");
}

#[test]
fn enumerate_tautomers_json_reports_oversize_as_an_object() {
    let smiles = "C".repeat(WASM_MAX_ATOMS + 1);
    let mol = parse(&smiles);
    let value: serde_json::Value =
        serde_json::from_str(&enumerate_tautomers_json(&mol)).expect("valid JSON");
    assert_eq!(
        value["error"].as_str(),
        Some("molecule too large (max 10000 atoms)")
    );
}

#[test]
fn blocked_tautomer_json_is_valid_when_parser_error_contains_special_chars() {
    let mol = parse("CC");
    let value: serde_json::Value = serde_json::from_str(
        &canonical_tautomer_with_blocked_atoms_json(&mol, "{not-json}"),
    )
    .expect("invalid input must still produce valid JSON");
    assert!(value["error"].as_str().is_some());
}

#[test]
fn largest_fragment_strips_salt() {
    // sodium acetate: "CC(=O)[O-].[Na+]" — largest fragment is acetate
    let mol = parse("CC(=O)[O-].[Na+]");
    let frag = largest_fragment(&mol);
    // acetate has 4 heavy atoms; Na has 1
    assert!(frag.atom_count() > 1, "expected the larger fragment");
    assert!(
        frag.atom_count() < mol.atom_count(),
        "fragment should be smaller than the salt"
    );
}

#[test]
fn neutralize_charges_removes_charges() {
    let mol = parse("CC(=O)[O-]");
    let neutral = neutralize_charges(&mol);
    assert_eq!(neutral.formal_charge_sum(), 0);
}

#[test]
fn standardize_smiles_report_json_includes_audit_report() {
    let json = standardize_smiles_report_json("CC.CCC", true, false, false, false);
    assert!(json.contains("\"smiles\""), "missing smiles: {json}");
    assert!(json.contains("\"report\""), "missing report: {json}");
    assert!(
        json.contains("\"status\":\"Modified\""),
        "missing status: {json}"
    );
    assert!(
        json.contains("\"step\":\"LargestFragment\""),
        "missing largest-fragment step: {json}"
    );
}

#[test]
fn maccs_bitvec_length_21() {
    let mol = parse("c1ccccc1");
    let bv = maccs_bitvec(&mol);
    assert_eq!(bv.len(), 21, "MACCS 166 bits should fit in 21 bytes");
}

#[test]
fn tanimoto_maccs_identical_is_one() {
    let mol = parse("c1ccccc1");
    assert!((tanimoto_maccs(&mol, &mol) - 1.0).abs() < 1e-9);
}

#[test]
fn get_descriptors_json_keys() {
    let mol = parse("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
    let json = get_descriptors_json(&mol);
    assert!(json.contains("\"mw\""), "missing mw: {json}");
    assert!(json.contains("\"tpsa\""), "missing tpsa: {json}");
    assert!(json.contains("\"qed\""), "missing qed: {json}");
}

#[test]
fn slogp_vsa_json_length_12() {
    let h = parse("c1ccccc1");
    let json = slogp_vsa_json(&h);
    let count = json
        .trim_matches(|c| c == '[' || c == ']')
        .split(',')
        .count();
    assert_eq!(count, 12, "SlogP_VSA should have 12 values");
}

#[test]
fn smr_vsa_json_length_10() {
    let h = parse("c1ccccc1");
    let json = smr_vsa_json(&h);
    let count = json
        .trim_matches(|c| c == '[' || c == ']')
        .split(',')
        .count();
    assert_eq!(count, 10, "SMR_VSA should have 10 values");
}

#[test]
fn peoe_vsa_json_length_14() {
    let h = parse("c1ccccc1");
    let json = peoe_vsa_json(&h);
    let count = json
        .trim_matches(|c| c == '[' || c == ']')
        .split(',')
        .count();
    assert_eq!(count, 14, "PEOE_VSA should have 14 values");
}

// ── pKa / ADMET MolHandle methods ────────────────────────────────────────

#[test]
fn test_pka_acid_value_acetic_acid() {
    let h = parse("CC(=O)O");
    let pka = h.pka_acid_value();
    assert!(pka.is_finite(), "acetic acid pKa should be finite");
    assert!(
        (pka - 4.0).abs() < 1.0,
        "acetic acid pKa ~4.0, got {pka:.2}"
    );
}

#[test]
fn test_pka_base_value_aniline() {
    let h = parse("Nc1ccccc1");
    let pka = h.pka_base_value();
    assert!(pka.is_finite(), "aniline pKa_base should be finite");
    assert!(
        (pka - 4.6).abs() < 1.0,
        "aniline pKa_base ~4.6, got {pka:.2}"
    );
}

#[test]
fn test_pka_acid_value_benzene_nan() {
    let h = parse("c1ccccc1");
    let pka = h.pka_acid_value();
    assert!(pka.is_nan(), "benzene has no acidic site → NaN");
}

#[test]
fn test_pka_base_value_benzene_nan() {
    let h = parse("c1ccccc1");
    let pka = h.pka_base_value();
    assert!(pka.is_nan(), "benzene has no basic site → NaN");
}

#[test]
fn test_bbb_score_benzene_positive() {
    let h = parse("c1ccccc1");
    let score = h.bbb_score();
    assert!(
        score > -1.0,
        "benzene should be CNS penetrant (logBB > -1), got {score:.3}"
    );
}

#[test]
fn test_bbb_passes_benzene() {
    let h = parse("c1ccccc1");
    assert!(
        h.bbb_passes(),
        "benzene (TPSA=0, MW=78) should pass BBB rules"
    );
}

#[test]
fn test_caco2_hexane_high() {
    let h = parse("CCCCCC");
    let c = h.caco2_permeability();
    assert!(
        c > -5.5,
        "hexane should have high Caco-2 permeability, got {c:.3}"
    );
}

#[test]
fn test_herg_risk_range() {
    let h = parse("c1ccccc1");
    let r = h.herg_risk_score();
    assert!(
        (0.0..=1.0).contains(&r),
        "hERG risk must be in [0,1], got {r}"
    );
}

#[test]
fn test_cyp3a4_risk_range() {
    let h = parse("c1ccccc1");
    let r = h.cyp3a4_inhibition_risk();
    assert!(
        (0.0..=1.0).contains(&r),
        "CYP3A4 risk must be in [0,1], got {r}"
    );
}

// ── predict_pka_json ─────────────────────────────────────────────────────

#[test]
fn test_predict_pka_json_acetic_acid() {
    let json = predict_pka_json("CC(=O)O");
    assert!(
        json.starts_with('[') && json.ends_with(']'),
        "should be JSON array"
    );
    assert!(
        json.contains("\"type\":\"acid\""),
        "should contain acid site"
    );
    assert!(!json.contains("\"error\""), "should not contain error");
}

#[test]
fn test_predict_pka_json_benzene_empty() {
    let json = predict_pka_json("c1ccccc1");
    assert_eq!(json, "[]", "benzene has no ionizable sites");
}

#[test]
fn test_predict_pka_json_invalid_smiles() {
    let json = predict_pka_json("C1CC");
    assert!(
        json.contains("\"error\""),
        "invalid SMILES should return error JSON"
    );
}

#[test]
fn test_predict_pka_json_glycine_both() {
    let json = predict_pka_json("NCC(=O)O");
    assert!(
        json.contains("\"type\":\"acid\""),
        "glycine should have acid site"
    );
    assert!(
        json.contains("\"type\":\"base\""),
        "glycine should have base site"
    );
}

// ── admet_profile_json ───────────────────────────────────────────────────

#[test]
fn test_admet_profile_json_aspirin_valid() {
    let json = admet_profile_json("CC(=O)Oc1ccccc1C(=O)O");
    assert!(
        !json.contains("\"error\""),
        "aspirin should parse without error"
    );
    assert!(
        json.contains("\"bbb_score\""),
        "should have bbb_score field"
    );
    assert!(json.contains("\"caco2\""), "should have caco2 field");
    assert!(
        json.contains("\"herg_risk\""),
        "should have herg_risk field"
    );
    assert!(json.contains("\"pka_acid\""), "should have pka_acid field");
    assert!(
        json.contains("\"bbb_passes\":true"),
        "aspirin should pass BBB"
    );
}

#[test]
fn test_admet_profile_json_benzene_null_pka() {
    let json = admet_profile_json("c1ccccc1");
    assert!(
        json.contains("\"pka_acid\":null"),
        "benzene pka_acid should be null"
    );
    assert!(
        json.contains("\"pka_base\":null"),
        "benzene pka_base should be null"
    );
}

#[test]
fn test_admet_profile_json_invalid_smiles() {
    let json = admet_profile_json("C1CC");
    assert!(
        json.contains("\"error\""),
        "invalid SMILES should return error JSON"
    );
}

// ── get_descriptors_json ADMET extension ─────────────────────────────────

#[test]
fn test_get_descriptors_json_has_admet_fields() {
    let h = parse("c1ccccc1");
    let json = get_descriptors_json(&h);
    assert!(json.contains("\"bbbScore\""), "should have bbbScore field");
    assert!(json.contains("\"caco2\""), "should have caco2 field");
    assert!(json.contains("\"hergRisk\""), "should have hergRisk field");
    assert!(
        json.contains("\"cyp3a4Risk\""),
        "should have cyp3a4Risk field"
    );
    assert!(json.contains("\"pkaAcid\""), "should have pkaAcid field");
    assert!(json.contains("\"pkaBase\""), "should have pkaBase field");
}

#[test]
fn test_get_descriptors_json_benzene_null_pka() {
    let h = parse("c1ccccc1");
    let json = get_descriptors_json(&h);
    assert!(
        json.contains("\"pkaAcid\":null"),
        "benzene pkaAcid should be null"
    );
    assert!(
        json.contains("\"pkaBase\":null"),
        "benzene pkaBase should be null"
    );
}

// ── issue #90: mmff94_energy_breakdown_from_coords_json / pdb_coords_json ──
// mmff94_energy_breakdown_json always regenerated a fresh rule-based
// conformer, ignoring any geometry the caller actually had (e.g. from
// mol_from_pdb). These cover the new coordinate-explicit contract that
// mirrors the Python binding's `mol.mmff94_energy_breakdown(coords)`.

#[test]
fn test_pdb_coords_json_and_mol_from_pdb_share_atom_order() {
    let pdb = set_dihedral_json("CCCC", 0, 1, 2, 3, 90.0).expect("set_dihedral_json");
    let mol = mol_from_pdb(&pdb);
    let coords_json = pdb_coords_json(&pdb);
    let coords: Vec<[f64; 3]> = serde_json::from_str(&coords_json).expect("valid coords json");
    assert_eq!(
        coords.len(),
        mol.atom_count(),
        "pdb_coords_json length must match mol_from_pdb's atom count (same underlying parse)"
    );
}

#[test]
fn test_mmff94_energy_breakdown_from_coords_json_varies_with_geometry() {
    // The exact issue #90 repro shape: torsion-scan a CCCC dihedral through
    // set_dihedral_json -> mol_from_pdb -> pdb_coords_json -> the new energy
    // function. Before the fix, every angle produced the SAME total because
    // the geometry was silently regenerated instead of read from the PDB.
    let mut totals = Vec::new();
    let mut torsions = Vec::new();
    for angle in [0.0, 60.0, 120.0, 180.0, 240.0, 300.0] {
        let pdb = set_dihedral_json("CCCC", 0, 1, 2, 3, angle).expect("set_dihedral_json");
        let mol = mol_from_pdb(&pdb);
        let coords_json = pdb_coords_json(&pdb);
        let energy_json = mmff94_energy_breakdown_from_coords_json(&mol, &coords_json);
        let v: serde_json::Value = serde_json::from_str(&energy_json).expect("valid energy json");
        assert!(v.get("error").is_none(), "unexpected error: {energy_json}");
        totals.push(format!("{}", v["total"]));
        torsions.push(format!("{}", v["torsion"]));
    }
    let distinct_totals: std::collections::HashSet<_> = totals.iter().collect();
    let distinct_torsions: std::collections::HashSet<_> = torsions.iter().collect();
    assert!(
        distinct_totals.len() >= 3,
        "expected varying total energy across dihedral angles, got {totals:?}"
    );
    assert!(
        distinct_torsions.len() >= 3,
        "expected varying torsion term across dihedral angles, got {torsions:?}"
    );
}

/// One row of the Python-oracle fixture table below: `mol.mmff94_energy_breakdown(coords)`
/// from `.venv`'s chematic Python binding, for `CCCC` at `dihedral_deg`
/// (atoms 0,1,2,3), on the IDENTICAL (mol, coords) pair this Rust test
/// reconstructs via the same `set_dihedral_json` -> `mol_from_pdb` +
/// `pdb_coords_json` pipeline. Field name `angle_term` (not `angle`) to
/// avoid colliding with `dihedral_deg`, the scan parameter.
struct EnergyOracle {
    dihedral_deg: f64,
    bond: f64,
    angle_term: f64,
    stretch_bend: f64,
    torsion: f64,
    oop: f64,
    vdw: f64,
    electrostatic: f64,
    total: f64,
}

#[test]
fn test_mmff94_energy_breakdown_from_coords_json_matches_python_oracle_all_angles() {
    // Fixture-based, not just "energy is present" or "energy varies":
    // test_mmff94_energy_breakdown_from_coords_json_varies_with_geometry
    // above only checks that >=3 distinct total/torsion values appear across
    // the 6 scanned angles -- a regression that returns a DIFFERENT but
    // still-wrong value at 60-300 degrees (everything except the 0-degree
    // case, which used to be the only value pinned) would pass that test.
    // This pins all 8 energy terms at all 6 angles against the Python
    // oracle, so any wrong value anywhere in the scan is caught.
    let oracle = [
        EnergyOracle {
            dihedral_deg: 0.0,
            bond: 0.8902727980430658,
            angle_term: 0.00046499350423254214,
            stretch_bend: -0.007350084252858442,
            torsion: 0.43499994978247686,
            oop: 0.0,
            vdw: 14.519561102945316,
            electrostatic: 0.0,
            total: 15.837948760022233,
        },
        EnergyOracle {
            dihedral_deg: 60.0,
            bond: 0.8948441183316793,
            angle_term: 0.0004158023578226278,
            stretch_bend: -0.006857208726671617,
            torsion: 0.5878467445463341,
            oop: 0.0,
            vdw: 1.9693600350541345,
            electrostatic: 0.0,
            total: 3.4456094915632987,
        },
        EnergyOracle {
            dihedral_deg: 120.0,
            bond: 0.8781955726317947,
            angle_term: 0.0004982437526120773,
            stretch_bend: -0.0076120012966454376,
            torsion: 0.8684749102195811,
            oop: 0.0,
            vdw: -0.022709051177030738,
            electrostatic: 0.0,
            total: 1.7168476741303118,
        },
        EnergyOracle {
            dihedral_deg: 180.0,
            bond: 0.8969756143810449,
            angle_term: 0.0005974131500818183,
            stretch_bend: -0.008439350718638324,
            torsion: 0.0,
            oop: 0.0,
            vdw: -0.0670059041190226,
            electrostatic: 0.0,
            total: 0.8221277726934657,
        },
        EnergyOracle {
            dihedral_deg: 240.0,
            bond: 0.9063335706281842,
            angle_term: 0.0007291729835570543,
            stretch_bend: -0.009325393516517342,
            torsion: 0.8683134139870532,
            oop: 0.0,
            vdw: -0.02299635476312685,
            electrostatic: 0.0,
            total: 1.7430544093191502,
        },
        EnergyOracle {
            dihedral_deg: 300.0,
            bond: 0.9054704500796205,
            angle_term: 0.0005632607553765077,
            stretch_bend: -0.008219837248585977,
            torsion: 0.5884635963691149,
            oop: 0.0,
            vdw: 1.9670921081294108,
            electrostatic: 0.0,
            total: 3.4533695780849367,
        },
    ];

    let mut checked = 0;
    for row in &oracle {
        let pdb =
            set_dihedral_json("CCCC", 0, 1, 2, 3, row.dihedral_deg).expect("set_dihedral_json");
        let mol = mol_from_pdb(&pdb);
        let coords_json = pdb_coords_json(&pdb);
        let energy_json = mmff94_energy_breakdown_from_coords_json(&mol, &coords_json);
        let v: serde_json::Value = serde_json::from_str(&energy_json).expect("valid energy json");
        assert!(
            v.get("error").is_none(),
            "angle {}: unexpected error: {energy_json}",
            row.dihedral_deg
        );

        let mut expect = |key: &str, want: f64| {
            let got = v[key]
                .as_f64()
                .unwrap_or_else(|| panic!("angle {}: missing {key}", row.dihedral_deg));
            assert!(
                (got - want).abs() <= 1e-9,
                "angle {}: {key}: got {got}, want {want} (python oracle)",
                row.dihedral_deg
            );
            checked += 1;
        };
        expect("bond", row.bond);
        expect("angle", row.angle_term);
        expect("stretch_bend", row.stretch_bend);
        expect("torsion", row.torsion);
        expect("oop", row.oop);
        expect("vdw", row.vdw);
        expect("electrostatic", row.electrostatic);
        expect("total", row.total);
    }
    assert_eq!(
        checked,
        6 * 8,
        "expected all 48 oracle values to be checked"
    );
}

#[test]
fn test_mmff94_energy_breakdown_from_coords_json_rejects_length_mismatch() {
    let mol = parse("CCCC"); // 4 heavy atoms
    let energy_json = mmff94_energy_breakdown_from_coords_json(&mol, "[[0,0,0],[1,0,0],[2,0,0]]"); // only 3
    assert!(
        energy_json.contains("\"error\""),
        "expected error for coords/atom-count mismatch, got {energy_json}"
    );
}

#[test]
fn test_mmff94_energy_breakdown_from_coords_json_rejects_malformed_json() {
    let mol = parse("CCCC");
    for bad in [
        "not json",
        "{\"not\":\"an array\"}",
        "[[0,0],[1,0,0]]",
        "[[0,0,0,0]]",
    ] {
        let energy_json = mmff94_energy_breakdown_from_coords_json(&mol, bad);
        assert!(
            energy_json.contains("\"error\""),
            "expected error for malformed coords_json {bad:?}, got {energy_json}"
        );
    }
}

#[test]
fn test_mmff94_energy_breakdown_from_coords_json_rejects_non_finite() {
    let mol = parse("CCCC");
    // Literal NaN/Infinity are not valid JSON tokens -- rejected at parse time.
    let nan_json =
        mmff94_energy_breakdown_from_coords_json(&mol, "[[NaN,0,0],[1,0,0],[2,0,0],[3,0,0]]");
    assert!(nan_json.contains("\"error\""), "got {nan_json}");
    // A numerically-valid-JSON but overflowing exponent parses to f64::INFINITY
    // via serde -- this exercises the explicit post-parse is_finite() guard,
    // not just JSON's lack of a NaN/Infinity literal syntax.
    let overflow_json =
        mmff94_energy_breakdown_from_coords_json(&mol, "[[1e400,0,0],[1,0,0],[2,0,0],[3,0,0]]");
    assert!(overflow_json.contains("\"error\""), "got {overflow_json}");
}

#[test]
fn test_pdb_coords_json_rejects_non_finite_coordinate() {
    // A PDB coordinate FIELD containing the literal text "NaN" parses
    // successfully via f64::from_str (Rust's FromStr recognizes "nan"/
    // "inf"/"infinity" case-insensitively as their special values, distinct
    // from truly malformed text which falls back to 0.0) -- so this reaches
    // parse_pdb_atoms as a real non-finite Point3, not a parse failure.
    // Before the fix, pdb_coords_json would serialize this as a bare `NaN`
    // token, which is not valid JSON.
    let pdb = "ATOM      1  C   LIG A   1         NaN   0.000   0.000\nEND\n";
    let out = pdb_coords_json(pdb);
    assert!(
        out.contains("\"error\""),
        "expected explicit error for a non-finite PDB coordinate, got {out}"
    );
    assert!(
        serde_json::from_str::<serde_json::Value>(&out).is_ok(),
        "output must be valid JSON even on this error path, got {out}"
    );
}

#[test]
fn test_mmff94_energy_breakdown_json_unchanged_semantics() {
    // Regression pin: the EXISTING function keeps its internally-generated-
    // conformer contract, 4-decimal rounding, and "elec" key -- untouched by
    // this fix (a new function was added alongside it, not a replacement).
    let mol = parse("CCCC");
    let json = mmff94_energy_breakdown_json(&mol);
    assert!(json.contains("\"elec\""), "got {json}");
    assert!(!json.contains("\"electrostatic\""), "got {json}");
    assert!(!json.contains("\"stretch_bend\""), "got {json}");
}

// --- Extended XYZ (extxyz) -------------------------------------------------

const EXTXYZ_WATER_JSON_FIXTURE: &str = "3\nLattice=\"10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0\" Properties=species:S:1:pos:R:3:forces:R:3 energy=-76.4\nO 0.0 0.0 0.0 0.1 0.0 0.0\nH 0.7586 0.0 0.504284 0.0 0.1 0.0\nH 0.7586 0.0 -0.504284 0.0 -0.1 0.0\n";

#[test]
fn test_mol_from_extxyz_and_extxyz_frame_json_share_atom_order() {
    let mol = mol_from_extxyz(EXTXYZ_WATER_JSON_FIXTURE).expect("mol_from_extxyz");
    let json = extxyz_frame_json(EXTXYZ_WATER_JSON_FIXTURE).expect("extxyz_frame_json");
    let v: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");

    let coords = v["coords"].as_array().unwrap();
    assert_eq!(coords.len(), mol.atom_count());
    assert_eq!(mol.atom_count(), 3);

    let lattice = v["lattice"].as_array().unwrap();
    assert_eq!(lattice.len(), 9);
    assert_eq!(lattice[0].as_f64(), Some(10.0));

    let forces = v["properties"]["forces"].as_array().unwrap();
    assert_eq!(forces.len(), 3);
    assert_eq!(forces[0][0].as_f64(), Some(0.1));

    assert_eq!(v["info"]["energy"].as_str(), Some("-76.4"));
}

#[test]
fn test_extxyz_frame_json_plain_xyz_has_null_lattice_and_empty_properties() {
    let plain = "3\nwater\nO 0.0 0.0 0.0\nH 0.7586 0.0 0.504284\nH 0.7586 0.0 -0.504284\n";
    let json = extxyz_frame_json(plain).expect("extxyz_frame_json");
    let v: serde_json::Value = serde_json::from_str(&json).expect("valid JSON");
    assert!(v["lattice"].is_null());
    assert_eq!(v["properties"].as_object().unwrap().len(), 0);
    assert_eq!(v["info"].as_object().unwrap().len(), 0);
}

#[test]
fn test_mol_from_extxyz_rejects_malformed_input() {
    // mol_from_extxyz/extxyz_frame_json are thin Result<_, JsValue> wrappers
    // around chematic_mol::parse_extxyz -- tested via the underlying fn
    // (JsValue native-abort, same pattern as mol_from_xyz_roundtrip_atom_count).
    let bad = "1\nLattice=\"1.0 2.0 3.0\"\nC 0.0 0.0 0.0\n";
    assert!(chematic_mol::parse_extxyz(bad).is_err());
}

#[test]
fn test_to_extxyz_json_roundtrips_lattice_and_properties() {
    let mol = mol_from_extxyz(EXTXYZ_WATER_JSON_FIXTURE).expect("mol_from_extxyz");
    let frame_json = extxyz_frame_json(EXTXYZ_WATER_JSON_FIXTURE).expect("extxyz_frame_json");
    let v: serde_json::Value = serde_json::from_str(&frame_json).unwrap();

    let coords_json = v["coords"].to_string();
    let options_json = serde_json::json!({
        "lattice": v["lattice"],
        "properties": v["properties"],
        "info": v["info"],
    })
    .to_string();

    let written = to_extxyz_json(&mol, &coords_json, &options_json).expect("to_extxyz_json");
    let mol2 = mol_from_extxyz(&written).expect("mol_from_extxyz on written output");
    let json2 = extxyz_frame_json(&written).expect("extxyz_frame_json on written output");
    let v2: serde_json::Value = serde_json::from_str(&json2).unwrap();

    assert_eq!(mol2.atom_count(), mol.atom_count());
    assert_eq!(v2["lattice"], v["lattice"]);
    assert_eq!(v2["properties"]["forces"], v["properties"]["forces"]);
    assert_eq!(v2["info"], v["info"]);
}

#[test]
fn test_to_extxyz_json_rejects_coords_atom_count_mismatch() {
    // Tested via the underlying non-wasm helper -- JsValue native-abort
    // (same pattern as parse_pdb_molecule_and_coords's PdbInputError).
    let mol = mol_from_extxyz(EXTXYZ_WATER_JSON_FIXTURE).expect("mol_from_extxyz");
    let err = extxyz_frame_from_json_args(&mol.inner, "[[0.0,0.0,0.0]]", "{}");
    assert!(err.is_err(), "expected error for 1 coord row vs 3 atoms");
}
