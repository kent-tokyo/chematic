//! Invalid records must produce typed, field-specific diagnostics, and resource
//! limits must be enforced before constructing oversized graphs.
use crate::*;
use serde_json::json;

fn v2000(atoms: &[&str], bonds: &[&str], properties: &str) -> String {
    format!(
        "contract\n  chematic\n\n{:>3}{:>3}  0  0  0  0  0  0  0  0  0 V2000\n{}{}{}M  END\n",
        atoms.len(),
        bonds.len(),
        atoms.iter().map(|a| format!("{a}\n")).collect::<String>(),
        bonds.iter().map(|b| format!("{b}\n")).collect::<String>(),
        properties
    )
}
const C: &str = "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0";

#[test]
fn mol2000_rejects_malformed_charges_and_rgroup_properties() {
    for property in [
        "M  CHG",
        "M  CHG x",
        "M  CHG 1",
        "M  CHG 1 0 1",
        "M  CHG 1 2 1",
        "M  CHG 1 x 1",
        "M  CHG 1 1 128",
        "M  CHG 1 1 nope",
        "M  RGP",
        "M  RGP x",
        "M  RGP 1",
        "M  RGP 1 0 1",
        "M  RGP 1 2 1",
        "M  RGP 1 x 1",
        "M  RGP 1 1 0",
        "M  RGP 1 1 10000",
        "M  RGP 1 1 nope",
    ] {
        let input = v2000(&[C], &[], &format!("{property}\n"));
        let error = parse_mol(&input).err().unwrap();
        assert!(
            matches!(error, MolParseError::InvalidPropertyLine { .. }),
            "{property}: {error}"
        );
        assert!(error.to_string().contains("line"), "{error}");
    }
}

#[test]
fn mol2000_rejects_invalid_atom_and_bond_records() {
    for bond in [
        "  0  1  1",
        "  1  0  1",
        "  1  3  1",
        "  x  2  1",
        "  1  x  1",
        "  1  2  x",
        "  1",
    ] {
        let error = parse_mol(&v2000(&[C, C], &[bond], ""))
            .err()
            .unwrap_or_else(|| panic!("accepted bond {bond}"));
        assert!(
            matches!(error, MolParseError::InvalidBondLine { .. }),
            "{bond}: {error}"
        );
    }
    for (range, replacement) in [
        (20..30, "       NaN"),
        (20..30, "  Infinity"),
        (20..30, "     wrong"),
        (31..34, "Qq "),
        (31..34, "R0 "),
    ] {
        let mut atom = C.to_owned();
        atom.replace_range(range, replacement);
        let error = parse_mol(&v2000(&[&atom], &[], "")).err().unwrap();
        assert!(
            matches!(
                error,
                MolParseError::InvalidAtomLine { .. } | MolParseError::UnknownElement { .. }
            ),
            "{error}"
        );
    }
    assert!(matches!(
        parse_mol(&v2000(&["short"], &[], "")),
        Err(MolParseError::InvalidAtomLine { .. })
    ));
    let error = parse_mol(&v2000(&[C, C], &["  1  2  1", "  1  2  1"], ""))
        .err()
        .unwrap();
    assert!(error.to_string().contains("duplicate"), "{error}");
}

fn v3000(atom: &str, bond: &str) -> String {
    format!(
        "contract\n  chematic\n\n  0  0  0  0  0  0  0  0  0  0  0 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 {atom}\nM  V30 2 C 1.5 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 {bond}\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n"
    )
}

#[test]
fn mol3000_rejects_invalid_atoms_and_unresolvable_bond_endpoints() {
    for atom in [
        "1 C 0 0",
        "x C 0 0 0 0",
        "1 C 0 0 NaN 0",
        "1 C 0 0 wrong 0",
        "1 Qq 0 0 0 0",
        "1 R0 0 0 0 0",
        "1 R# 0 0 0 0 RGROUPS=(2 1 2)",
        "1 R# 0 0 0 0 RGROUPS=(1 0)",
    ] {
        let error = parse_mol_v3000(&v3000(atom, "1 1 1 2")).err().unwrap();
        assert!(!error.to_string().is_empty(), "{atom}");
    }
    for bond in ["1 1", "1 x 1 2", "1 1 x 2", "1 1 1 x", "1 1 9 2", "1 1 1 9"] {
        let error = parse_mol_v3000(&v3000("1 C 0 0 0 0", bond))
            .err()
            .unwrap_or_else(|| panic!("accepted bond {bond}"));
        assert!(
            matches!(error, MolParseError::InvalidBondLine { .. }),
            "{bond}: {error}"
        );
    }
}

#[test]
fn chemical_json_reports_invalid_values_at_the_corresponding_field() {
    let original = json!({"atoms":{"elements":{"number":[6,8]},"formalCharges":[0,0],"isotopes":[0,0],"coords":{"3d":[0,0,0,1,0,0]}},"bonds":{"connections":{"index":[0,1]},"order":[1]}});
    for (pointer, value, diagnostic) in [
        ("/atoms/elements/number/0", json!(-1), "number"),
        ("/atoms/elements/number/0", json!(999), "atomic number"),
        ("/atoms/formalCharges/0", json!(128), "formalCharges"),
        ("/atoms/formalCharges/0", json!("x"), "formalCharges"),
        ("/atoms/isotopes/0", json!(65536), "isotopes"),
        ("/atoms/isotopes/0", json!(-1), "isotopes"),
        ("/atoms/coords/3d/0", json!("x"), "coords"),
        ("/atoms/coords/3d", json!([0]), "coords"),
        ("/bonds/connections/index/0", json!(-1), "index"),
        ("/bonds/connections/index", json!([0]), "even"),
        ("/bonds/connections/index/0", json!(2), "references"),
        ("/bonds/connections/index/1", json!(2), "references"),
        ("/bonds/order/0", json!(false), "order"),
    ] {
        let mut input = original.clone();
        *input.pointer_mut(pointer).unwrap() = value;
        let error = parse_cjson(&input.to_string()).err().unwrap();
        assert!(error.to_string().contains(diagnostic), "{pointer}: {error}");
    }
    for input in [
        "{",
        "{}",
        "{\"atoms\":{}}",
        "{\"atoms\":{\"elements\":{\"number\":[6]}},\"bonds\":{}}",
    ] {
        assert!(parse_cjson(input).is_err(), "{input}");
    }
    let text = original.to_string();
    for limits in [
        CjsonParseLimits {
            max_input_bytes: 1,
            ..Default::default()
        },
        CjsonParseLimits {
            max_json_depth: 0,
            ..Default::default()
        },
        CjsonParseLimits {
            max_array_items: 1,
            ..Default::default()
        },
        CjsonParseLimits {
            max_atoms: 1,
            ..Default::default()
        },
        CjsonParseLimits {
            max_bonds: 0,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            parse_cjson_with_limits(&text, &limits),
            Err(CjsonError::ResourceLimit { .. })
        ));
    }
    let strings = json!({"atoms":{"elements":{"number":[6]}},"metadata":"long"}).to_string();
    assert!(matches!(
        parse_cjson_with_limits(
            &strings,
            &CjsonParseLimits {
                max_string_bytes: 1,
                ..Default::default()
            }
        ),
        Err(CjsonError::ResourceLimit { .. })
    ));
}
