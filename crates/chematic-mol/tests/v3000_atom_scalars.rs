use chematic_core::AtomIdx;
use chematic_mol::{
    MolParseError, parse_mol_v3000_with_coords, roundtrip_mol_v3000_source, write_mol_v3000,
};
use serde_json::Value;

fn packet() -> Value {
    serde_json::from_str(include_str!(
        "../../../validation/v3000_atom_scalars/cases.json"
    ))
    .unwrap()
}

#[test]
fn explicit_scalar_meaning_and_quoted_text_survive_native_routes() {
    for c in packet()["cases"].as_array().unwrap() {
        let input = c["input"].as_str().unwrap();
        let (mol, meta, coords) = parse_mol_v3000_with_coords(input).unwrap();
        let atom = mol.atom(AtomIdx(0));
        let e = &c["expected"];
        assert_eq!(
            u64::from(atom.element.atomic_number()),
            e["number"].as_u64().unwrap(),
            "{}",
            c["id"]
        );
        assert_eq!(
            i64::from(atom.charge),
            e["charge"].as_i64().unwrap(),
            "{}",
            c["id"]
        );
        assert_eq!(
            u64::from(atom.isotope.unwrap_or(0)),
            e["isotope"].as_u64().unwrap(),
            "{}",
            c["id"]
        );
        assert_eq!(
            u64::from(atom.atom_map.unwrap_or(0)),
            e["atom_map"].as_u64().unwrap(),
            "{}",
            c["id"]
        );
        assert_eq!(roundtrip_mol_v3000_source(input).unwrap(), input);
        let output = write_mol_v3000(&mol, &meta, &coords);
        let (parsed, _, _) = parse_mol_v3000_with_coords(&output).unwrap();
        let after = parsed.atom(AtomIdx(0));
        assert_eq!(
            (after.charge, after.isotope, after.atom_map),
            (atom.charge, atom.isotope, atom.atom_map)
        );
        if c["scope"] == "quoted_compatibility" {
            assert!(output.contains("LABEL=\"memo CHG=2 MASS=13 end\""));
        }
    }
}

#[test]
fn malformed_unsupported_and_duplicate_scalars_fail_on_source_line() {
    for c in packet()["invalid"].as_array().unwrap() {
        let input = c["input"].as_str().unwrap();
        let field = c["field"].as_str().unwrap();
        for text in [
            input.to_string(),
            input.replacen("M  V30 1 ", "M  V30 1 -\nM  V30 ", 1),
        ] {
            for result in [
                parse_mol_v3000_with_coords(&text).map(|_| ()),
                roundtrip_mol_v3000_source(&text).map(|_| ()),
            ] {
                assert!(
                    matches!(result,Err(MolParseError::InvalidAtomLine{line:8,ref detail}) if detail.contains(field)),
                    "{}: {:?}",
                    c["id"],
                    result.err()
                );
            }
        }
    }
}

#[test]
fn field_storage_boundaries_and_grouped_or_escaped_text_are_explicit() {
    let source = packet()["cases"][0]["input"].as_str().unwrap().to_owned();
    for (charge, mass, mapping) in [(-128, 1, 0), (127, 65535, 65535)] {
        let text = source.replace(
            "C 0 0 0 0",
            &format!("C 0 0 0 {mapping} CHG={charge} MASS={mass}"),
        );
        let (mol, _, _) = parse_mol_v3000_with_coords(&text).unwrap();
        let atom = mol.atom(AtomIdx(0));
        assert_eq!(
            (atom.charge, atom.isotope, atom.atom_map),
            (charge, Some(mass), (mapping != 0).then_some(mapping))
        );
    }
    for property in [
        "FOO=(4 0 CHG=2 MASS=13 0)",
        "LABEL=\"memo \"\" CHG=2 MASS=13 end\"",
    ] {
        let text = source.replace("C 0 0 0 0", &format!("C 0 0 0 0 {property}"));
        let (mol, meta, coords) = parse_mol_v3000_with_coords(&text).unwrap();
        assert_eq!(
            (mol.atom(AtomIdx(0)).charge, mol.atom(AtomIdx(0)).isotope),
            (0, None)
        );
        assert!(write_mol_v3000(&mol, &meta, &coords).contains(property));
    }
}
