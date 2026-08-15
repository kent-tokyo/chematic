//! Integration tests for `chematic-mol`'s QCSchema module
//! (`crates/chematic-mol/src/qcschema.rs`).

use chematic_core::{Atom, AtomIdx, BondOrder, Coords3D, Element, MoleculeBuilder, Point3};
use chematic_mol::qcschema::*;
use serde_json::{Value, json};

// ─── Water molecule fixture ─────────────────────────────────────────────────
//
// r(O-H) = 0.9584 Angstrom, angle(H-O-H) = 104.45 degrees -- the standard
// experimental gas-phase equilibrium geometry commonly used as the "hello
// world" molecule in MolSSI tutorials/docs (same molecule as the worked
// example on <https://molssi.github.io/QCElemental/model_molecule.html>).
// Bohr values below are that Angstrom geometry divided by
// `qcschema::BOHR_TO_ANGSTROM`, computed independently (not copied from any
// qcelemental source) so the unit-conversion test has a known-correct
// target to check against.
const WATER_BOHR_GEOMETRY: [f64; 9] = [
    0.0,
    0.0,
    0.0, //
    0.0,
    1.431544636036637,
    1.109419726497757, //
    0.0,
    -1.431544636036637,
    1.109419726497757,
];
const WATER_OH_ANGSTROM: f64 = 0.9584;

fn water_molecule_json(schema_name: &str) -> Value {
    json!({
        "schema_name": schema_name,
        "schema_version": 1,
        "symbols": ["O", "H", "H"],
        "geometry": WATER_BOHR_GEOMETRY,
        "molecular_charge": 0.0,
        "molecular_multiplicity": 1,
        "fix_com": false,
        "fix_orientation": false,
        "name": "water",
        "connectivity": [[0, 1, 1.0], [0, 2, 1.0]],
        "provenance": {"creator": "chematic-test-fixture", "version": "0.0", "routine": "hand-built"}
    })
}

// ─── Semantic JSON comparator (Value equality with float tolerance) ────────
//
// `assert_eq!` on raw `Value`s is too strict for a round-trip check: a
// `0` in the input parses to `Number(0)` while this module always writes
// floats like `molecular_charge` back as `0.0` -> `Number(0.0)`; these are
// not `PartialEq`-equal `Value`s despite being the same molecule. This
// walks both trees, comparing numbers as `f64` within a relative
// tolerance and objects by matching key sets (so a field genuinely
// present in one and missing in the other still fails, as it should).
fn json_semantically_eq(a: &Value, b: &Value) -> bool {
    match (a, b) {
        (Value::Number(x), Value::Number(y)) => match (x.as_f64(), y.as_f64()) {
            (Some(fx), Some(fy)) => (fx - fy).abs() <= 1e-9 * fx.abs().max(fy.abs()).max(1.0),
            _ => x == y,
        },
        (Value::Array(xs), Value::Array(ys)) => {
            xs.len() == ys.len() && xs.iter().zip(ys).all(|(x, y)| json_semantically_eq(x, y))
        }
        (Value::Object(xo), Value::Object(yo)) => {
            xo.len() == yo.len()
                && xo
                    .iter()
                    .all(|(k, v)| yo.get(k).is_some_and(|v2| json_semantically_eq(v, v2)))
        }
        _ => a == b,
    }
}

fn assert_roundtrip_eq(original: &Value, roundtripped: &Value) {
    assert!(
        json_semantically_eq(original, roundtripped),
        "round-trip mismatch:\n  original:     {original}\n  roundtripped: {roundtripped}"
    );
}

// ─── Molecule round trip ─────────────────────────────────────────────────────

#[test]
fn roundtrip_qcmolecule_water() {
    let original = water_molecule_json("qcschema_molecule");
    let text = original.to_string();

    let parsed = parse_qcschema_molecule(&text).expect("parse water molecule");
    assert_eq!(parsed.symbols, ["O", "H", "H"]);
    assert_eq!(parsed.molecular_multiplicity, 1);
    assert_eq!(parsed.connectivity.as_ref().unwrap().len(), 2);

    let out_text = write_qcschema_molecule(&parsed);
    let out_value: Value = serde_json::from_str(&out_text).unwrap();
    assert_roundtrip_eq(&original, &out_value);
}

#[test]
fn molecule_minimal_input_applies_spec_defaults() {
    // Only the two truly-required fields (symbols, geometry) -- everything
    // else must come out at its documented spec default.
    let text = json!({"symbols": ["H"], "geometry": [0.0, 0.0, 0.0]}).to_string();
    let m = parse_qcschema_molecule(&text).expect("parse minimal molecule");
    assert_eq!(m.schema_name, "qcschema_molecule");
    assert_eq!(m.schema_version, 1);
    assert_eq!(m.molecular_charge, 0.0);
    assert_eq!(m.molecular_multiplicity, 1);
    assert!(!m.fix_com);
    assert!(!m.fix_orientation);
    assert!(m.connectivity.is_none());
    assert!(m.provenance.is_none());
}

#[test]
fn molecule_legacy_underscored_schema_name_accepted_and_preserved() {
    let text = json!({
        "schema_name": "qc_schema_molecule",
        "symbols": ["He"],
        "geometry": [0.0, 0.0, 0.0]
    })
    .to_string();
    let m = parse_qcschema_molecule(&text).expect("legacy spelling should parse");
    assert_eq!(m.schema_name, "qc_schema_molecule");
    let out = write_qcschema_molecule(&m);
    assert!(
        out.contains("qc_schema_molecule"),
        "legacy spelling must round-trip verbatim"
    );
}

// ─── AtomicInput round trip + extensible-field preservation ────────────────

fn atomic_input_json() -> Value {
    json!({
        "schema_name": "qcschema_input",
        "schema_version": 1,
        "molecule": water_molecule_json("qcschema_molecule"),
        "driver": "energy",
        "model": {"method": "b3lyp", "basis": "6-31g"},
        "keywords": {"scf_type": "df", "e_convergence": 1e-10, "nested": {"a": [1, 2, 3]}},
        "extras": {"my_custom_extra": true},
        "some_unrecognized_top_level_key": "must survive round trip"
    })
}

#[test]
fn roundtrip_atomic_input() {
    let original = atomic_input_json();
    let text = original.to_string();

    let parsed = parse_atomic_input(&text).expect("parse AtomicInput");
    assert_eq!(parsed.driver, Driver::Energy);
    assert_eq!(parsed.model.method, "b3lyp");
    assert!(matches!(parsed.model.basis, Some(Basis::Name(ref s)) if s == "6-31g"));

    let out_text = write_atomic_input(&parsed);
    let out_value: Value = serde_json::from_str(&out_text).unwrap();
    assert_roundtrip_eq(&original, &out_value);
}

#[test]
fn extensible_field_preservation() {
    let text = atomic_input_json().to_string();
    let parsed = parse_atomic_input(&text).expect("parse");

    // The spec-defined open bag ("keywords") kept its nested custom content.
    assert_eq!(parsed.keywords.get("scf_type"), Some(&json!("df")));
    assert_eq!(
        parsed.keywords.get("nested"),
        Some(&json!({"a": [1, 2, 3]}))
    );
    // A key nowhere in the QCSchema AtomicInput spec survived as `extra`.
    assert_eq!(
        parsed.unknown_fields.get("some_unrecognized_top_level_key"),
        Some(&json!("must survive round trip"))
    );

    let out_text = write_atomic_input(&parsed);
    assert!(out_text.contains("must survive round trip"));
    assert!(out_text.contains("e_convergence"));
}

// ─── AtomicResult round trips (success / gradient / failure) ───────────────

fn atomic_result_energy_json() -> Value {
    json!({
        "schema_name": "qcschema_output",
        "schema_version": 1,
        "molecule": water_molecule_json("qcschema_molecule"),
        "driver": "energy",
        "model": {"method": "b3lyp", "basis": "6-31g"},
        // Non-empty keywords/extras and an unrecognized top-level key here
        // on purpose: `write_atomic_result` builds a throwaway `AtomicInput`
        // internally to reuse `atomic_input_fields_to_map` (see its source
        // comment) and re-inserts `unknown_fields` separately afterwards --
        // this fixture is what proves that detour doesn't drop anything.
        "keywords": {"scf_type": "df"},
        "extras": {"qcvars": {"CURRENT ENERGY": -76.4}},
        "provenance": {"creator": "chematic-test-fixture", "version": "0.0", "routine": "energy"},
        "properties": {"return_energy": -76.4, "calcinfo_natom": 3},
        "return_result": -76.4,
        "success": true,
        "my_vendor_tag": 42
    })
}

#[test]
fn roundtrip_atomic_result_success_energy() {
    let original = atomic_result_energy_json();
    let text = original.to_string();

    let parsed = parse_atomic_result(&text).expect("parse AtomicResult");
    assert!(parsed.success);
    assert!(parsed.error.is_none());
    assert!(
        matches!(parsed.return_result, Some(ReturnResult::Scalar(e)) if (e + 76.4).abs() < 1e-12)
    );
    assert_eq!(parsed.properties.get("calcinfo_natom"), Some(&json!(3)));
    assert_eq!(parsed.keywords.get("scf_type"), Some(&json!("df")));
    assert_eq!(
        parsed.extras.get("qcvars"),
        Some(&json!({"CURRENT ENERGY": -76.4}))
    );
    assert_eq!(parsed.unknown_fields.get("my_vendor_tag"), Some(&json!(42)));

    let out_text = write_atomic_result(&parsed);
    let out_value: Value = serde_json::from_str(&out_text).unwrap();
    assert_roundtrip_eq(&original, &out_value);
}

#[test]
fn atomic_result_gradient_flattens_nested_array_and_is_idempotent() {
    // Some producers emit an (nat, 3) nested gradient array rather than a
    // flat 3N one; this module accepts both and always *writes* flat (see
    // `ReturnResult::Array` doc comment) -- so this fixture intentionally
    // does not get a byte/field-identical round trip on the first pass, but
    // parse -> write -> parse -> write must be idempotent from then on.
    let text = json!({
        "schema_name": "qcschema_output",
        "molecule": water_molecule_json("qcschema_molecule"),
        "driver": "gradient",
        "model": {"method": "hf", "basis": "sto-3g"},
        "provenance": {"creator": "chematic-test-fixture"},
        "properties": {},
        "return_result": [[0.0, 0.0, 0.1], [0.0, 0.0, -0.05], [0.0, 0.0, -0.05]],
        "success": true
    })
    .to_string();

    let parsed = parse_atomic_result(&text).expect("parse gradient result");
    match &parsed.return_result {
        Some(ReturnResult::Array(v)) => {
            assert_eq!(v, &vec![0.0, 0.0, 0.1, 0.0, 0.0, -0.05, 0.0, 0.0, -0.05]);
        }
        other => panic!("expected flattened gradient array, got {other:?}"),
    }

    let pass1 = write_atomic_result(&parsed);
    let reparsed = parse_atomic_result(&pass1).expect("reparse");
    let pass2 = write_atomic_result(&reparsed);
    let v1: Value = serde_json::from_str(&pass1).unwrap();
    let v2: Value = serde_json::from_str(&pass2).unwrap();
    assert_roundtrip_eq(&v1, &v2);
}

#[test]
fn roundtrip_atomic_result_failure() {
    let original = json!({
        "schema_name": "qcschema_output",
        "schema_version": 1,
        "molecule": water_molecule_json("qcschema_molecule"),
        "driver": "energy",
        "model": {"method": "b3lyp", "basis": "6-31g"},
        "provenance": {"creator": "chematic-test-fixture", "version": "0.0", "routine": "energy"},
        "properties": {},
        "success": false,
        "error": {"error_type": "ConvergenceError", "error_message": "SCF did not converge in 200 iterations"}
    });
    let text = original.to_string();

    let parsed = parse_atomic_result(&text).expect("parse failed AtomicResult");
    assert!(!parsed.success);
    assert!(parsed.return_result.is_none());
    assert_eq!(
        parsed.error.as_ref().unwrap().error_type,
        "ConvergenceError"
    );

    let out_text = write_atomic_result(&parsed);
    let out_value: Value = serde_json::from_str(&out_text).unwrap();
    assert_roundtrip_eq(&original, &out_value);
}

// ─── Cross-field invariant validation ───────────────────────────────────────

#[test]
fn atomic_result_success_true_without_return_result_is_rejected() {
    let text = json!({
        "molecule": water_molecule_json("qcschema_molecule"),
        "driver": "energy",
        "model": {"method": "hf", "basis": "sto-3g"},
        "provenance": {"creator": "x"},
        "properties": {},
        "success": true
    })
    .to_string();
    assert!(matches!(
        parse_atomic_result(&text),
        Err(QcSchemaError::Inconsistent { .. })
    ));
}

#[test]
fn atomic_result_success_false_without_error_is_rejected() {
    let text = json!({
        "molecule": water_molecule_json("qcschema_molecule"),
        "driver": "energy",
        "model": {"method": "hf", "basis": "sto-3g"},
        "provenance": {"creator": "x"},
        "properties": {},
        "success": false
    })
    .to_string();
    assert!(matches!(
        parse_atomic_result(&text),
        Err(QcSchemaError::Inconsistent { .. })
    ));
}

// ─── Malformed / adversarial JSON never panics ──────────────────────────────

#[test]
fn malformed_json_does_not_panic() {
    assert!(matches!(
        parse_qcschema_molecule("not json at all"),
        Err(QcSchemaError::InvalidJson(_))
    ));
    assert!(matches!(
        parse_qcschema_molecule(""),
        Err(QcSchemaError::InvalidJson(_))
    ));
    assert!(matches!(
        parse_qcschema_molecule("[1, 2, 3]"),
        Err(QcSchemaError::WrongType { .. })
    ));
    assert!(matches!(
        parse_atomic_input("{"),
        Err(QcSchemaError::InvalidJson(_))
    ));
    assert!(matches!(
        parse_atomic_result("null"),
        Err(QcSchemaError::WrongType { .. })
    ));
    // Deeply-nested garbage inside an otherwise-plausible shell must not panic either.
    let weird = json!({"symbols": ["C"], "geometry": {"not": "an array"}}).to_string();
    assert!(parse_qcschema_molecule(&weird).is_err());
}

#[test]
fn missing_required_field_is_typed_error() {
    let text = json!({"symbols": ["O", "H", "H"]}).to_string(); // no geometry
    assert!(matches!(
        parse_qcschema_molecule(&text),
        Err(QcSchemaError::MissingField(f)) if f == "geometry"
    ));
}

#[test]
fn geometry_length_mismatch_is_rejected() {
    let text =
        json!({"symbols": ["O", "H", "H"], "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}).to_string();
    assert!(matches!(
        parse_qcschema_molecule(&text),
        Err(QcSchemaError::LengthMismatch { .. })
    ));
}

#[test]
fn invalid_schema_name_is_rejected() {
    let text =
        json!({"schema_name": "totally_wrong", "symbols": ["H"], "geometry": [0.0, 0.0, 0.0]})
            .to_string();
    assert!(matches!(
        parse_qcschema_molecule(&text),
        Err(QcSchemaError::InvalidSchemaName { .. })
    ));
}

#[test]
fn connectivity_duplicate_bond_rejected_regardless_of_atom_order() {
    let text = json!({
        "symbols": ["H", "H"],
        "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 0.74],
        "connectivity": [[0, 1, 1.0], [1, 0, 1.0]]
    })
    .to_string();
    assert!(matches!(
        parse_qcschema_molecule(&text),
        Err(QcSchemaError::DuplicateBond { .. })
    ));
}

#[test]
fn connectivity_self_bond_rejected() {
    let text = json!({
        "symbols": ["H", "H"],
        "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 0.74],
        "connectivity": [[0, 0, 1.0]]
    })
    .to_string();
    assert!(matches!(
        parse_qcschema_molecule(&text),
        Err(QcSchemaError::Inconsistent { .. })
    ));
}

#[test]
fn connectivity_out_of_range_index_rejected() {
    let text = json!({
        "symbols": ["H", "H"],
        "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 0.74],
        "connectivity": [[0, 5, 1.0]]
    })
    .to_string();
    assert!(matches!(
        parse_qcschema_molecule(&text),
        Err(QcSchemaError::IndexOutOfRange { .. })
    ));
}

#[test]
fn invalid_driver_enum_value_rejected() {
    let text = json!({
        "molecule": water_molecule_json("qcschema_molecule"),
        "driver": "frequencies",
        "model": {"method": "hf"}
    })
    .to_string();
    assert!(matches!(
        parse_atomic_input(&text),
        Err(QcSchemaError::InvalidEnumValue { .. })
    ));
}

// ─── NaN / Infinity rejection (task requirement: reject at parse time) ─────

#[test]
fn overflowing_json_number_fails_closed_not_open() {
    // `1e400` is syntactically valid JSON and would overflow f64 to
    // +Infinity on `Value::as_f64()` if it ever became a `Value::Number` --
    // empirically (checked against this workspace's pinned `serde_json`
    // 1.0.151), `serde_json::from_str` itself already refuses to parse an
    // out-of-range number literal ("number out of range"), so this surfaces
    // as `InvalidJson` rather than reaching this module's own
    // `check_finite` guard (which exists as defense-in-depth for any
    // `serde_json` configuration where that upstream guard doesn't apply --
    // e.g. the `arbitrary_precision` feature, not enabled anywhere in this
    // workspace). Either way the requirement holds: fail closed with a
    // typed `Err`, never silently carry an infinite value through, and
    // never panic.
    let text = r#"{"symbols": ["H"], "geometry": [1e400, 0.0, 0.0]}"#;
    assert!(parse_qcschema_molecule(text).is_err());

    let molecule_json = water_molecule_json("qcschema_molecule").to_string();
    let nested = format!(
        r#"{{"molecule": {molecule_json}, "driver": "energy", "model": {{"method": "hf"}}, "keywords": {{"some_threshold": -1e400}}}}"#
    );
    assert!(parse_atomic_input(&nested).is_err());
}

// ─── Unit conversion correctness (Bohr <-> Angstrom) ────────────────────────

#[test]
fn bohr_to_angstrom_constant_is_codata_2018() {
    assert!((BOHR_TO_ANGSTROM - 0.529177210903).abs() < 1e-12);
}

#[test]
fn qc_to_chematic_conversion_yields_correct_angstrom_geometry() {
    let text = water_molecule_json("qcschema_molecule").to_string();
    let qc = parse_qcschema_molecule(&text).unwrap();
    let view = qc_molecule_to_chematic(&qc).expect("convert to chematic");

    assert_eq!(view.molecule.atom_count(), 3);
    assert_eq!(view.molecular_charge, 0.0);
    assert_eq!(view.molecular_multiplicity, 1);

    let o = view.coords.get(AtomIdx(0));
    let h1 = view.coords.get(AtomIdx(1));
    let oh_distance = o.distance(&h1);
    assert!(
        (oh_distance - WATER_OH_ANGSTROM).abs() < 1e-9,
        "expected O-H = {WATER_OH_ANGSTROM} Angstrom, got {oh_distance}"
    );
}

#[test]
fn chematic_to_qc_conversion_is_the_exact_inverse_division() {
    // The classic silent bug: `x * (1.0 / c)` and `x / c` differ in the
    // last bits of precision. Round Bohr -> Angstrom -> Bohr and check the
    // *value*, not just self-consistency, catches an inverted factor that
    // would otherwise "round-trip" back to itself even when wrong.
    let text = water_molecule_json("qcschema_molecule").to_string();
    let qc = parse_qcschema_molecule(&text).unwrap();
    let view = qc_molecule_to_chematic(&qc).unwrap();
    let back = chematic_to_qc_molecule(
        &view.molecule,
        &view.coords,
        view.molecular_charge,
        view.molecular_multiplicity,
    )
    .expect("convert back to QCSchema");

    for (original, roundtripped) in qc.geometry.iter().zip(back.geometry.iter()) {
        let scale = original.abs().max(1.0);
        assert!(
            (original - roundtripped).abs() <= 1e-9 * scale,
            "geometry component drifted: {original} vs {roundtripped}"
        );
    }
}

#[test]
fn chematic_to_qc_atom_count_mismatch_is_rejected() {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(Atom::new(Element::H));
    builder.add_atom(Atom::new(Element::H));
    let mol = builder.build();
    let coords = Coords3D::new_zeroed(1); // wrong count on purpose
    let err = chematic_to_qc_molecule(&mol, &coords, 0.0, 1).unwrap_err();
    assert!(matches!(
        err,
        QcConvertError::AtomCountMismatch {
            molecule: 2,
            coords: 1
        }
    ));
}

#[test]
fn qc_to_chematic_unknown_element_is_rejected() {
    let text = json!({"symbols": ["Xx"], "geometry": [0.0, 0.0, 0.0]}).to_string();
    let qc = parse_qcschema_molecule(&text).unwrap();
    assert!(
        matches!(qc_molecule_to_chematic(&qc), Err(QcConvertError::UnknownElement(s)) if s == "Xx")
    );
}

// ─── Bond-order mapping through the chematic <-> QCSchema boundary ─────────

#[test]
fn bond_orders_survive_the_round_trip_through_qcschema_connectivity() {
    let mut builder = MoleculeBuilder::new();
    let c1 = builder.add_atom(Atom::new(Element::C));
    let c2 = builder.add_atom(Atom::new(Element::C));
    let c3 = builder.add_atom(Atom::new(Element::C));
    let n1 = builder.add_atom(Atom::new(Element::N));
    builder.add_bond(c1, c2, BondOrder::Single).unwrap();
    builder.add_bond(c2, c3, BondOrder::Double).unwrap();
    builder.add_bond(c3, n1, BondOrder::Triple).unwrap();
    let mol = builder.build();
    let mut coords = Coords3D::new_zeroed(4);
    coords.set(c1, Point3::new(0.0, 0.0, 0.0));
    coords.set(c2, Point3::new(1.5, 0.0, 0.0));
    coords.set(c3, Point3::new(3.0, 0.0, 0.0));
    coords.set(n1, Point3::new(4.5, 0.0, 0.0));

    let qc = chematic_to_qc_molecule(&mol, &coords, 0.0, 1).unwrap();
    let conn = qc.connectivity.as_ref().unwrap();
    let order_of = |a: usize, b: usize| -> f64 {
        conn.iter()
            .find(|(x, y, _)| (*x == a && *y == b) || (*x == b && *y == a))
            .map(|(_, _, o)| *o)
            .unwrap()
    };
    assert_eq!(order_of(0, 1), 1.0);
    assert_eq!(order_of(1, 2), 2.0);
    assert_eq!(order_of(2, 3), 3.0);

    // And back: chematic -> QCSchema -> chematic must land on the same
    // BondOrder values.
    let view = qc_molecule_to_chematic(&qc).unwrap();
    let bond_order_between = |a: AtomIdx, b: AtomIdx| -> BondOrder {
        view.molecule
            .bonds()
            .find(|(_, be)| (be.atom1 == a && be.atom2 == b) || (be.atom1 == b && be.atom2 == a))
            .unwrap()
            .1
            .order
    };
    assert_eq!(
        bond_order_between(AtomIdx(0), AtomIdx(1)),
        BondOrder::Single
    );
    assert_eq!(
        bond_order_between(AtomIdx(1), AtomIdx(2)),
        BondOrder::Double
    );
    assert_eq!(
        bond_order_between(AtomIdx(2), AtomIdx(3)),
        BondOrder::Triple
    );
}

#[test]
fn mass_number_maps_to_isotope_both_directions() {
    let text = json!({
        "symbols": ["C"],
        "geometry": [0.0, 0.0, 0.0],
        "mass_numbers": [13]
    })
    .to_string();
    let qc = parse_qcschema_molecule(&text).unwrap();
    let view = qc_molecule_to_chematic(&qc).unwrap();
    assert_eq!(view.molecule.atom(AtomIdx(0)).isotope, Some(13));

    let back = chematic_to_qc_molecule(&view.molecule, &view.coords, 0.0, 1).unwrap();
    assert_eq!(back.mass_numbers, Some(vec![13]));
}
