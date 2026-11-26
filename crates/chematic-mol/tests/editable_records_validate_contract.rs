use chematic_mol::{
    FractRecord, GprRecord, Gromos96Document, LpmdDocument, UnixyzRecord, ViewmolRecord,
};
use serde_json::Value;

fn first_fixture(source: &str, field: &str) -> Value {
    let packet: Value = serde_json::from_str(source).unwrap();
    packet["valid"][0][field].clone()
}

#[test]
fn edited_records_are_revalidated_before_use() {
    let mut fract: FractRecord = serde_json::from_value(first_fixture(
        include_str!("../../../validation/fract_preservation/cases.json"),
        "expected",
    ))
    .unwrap();
    assert!(fract.validate().is_ok());
    fract.atoms.clear();
    assert!(fract.validate().is_err());

    let mut gpr: GprRecord = serde_json::from_value(first_fixture(
        include_str!("../../../validation/gpr_preservation/cases.json"),
        "expected",
    ))
    .unwrap();
    assert!(gpr.validate().is_ok());
    gpr.atoms.clear();
    assert!(gpr.validate().is_err());

    let mut gromos: Gromos96Document = serde_json::from_value(first_fixture(
        include_str!("../../../validation/gromos96_preservation/cases.json"),
        "record",
    ))
    .unwrap();
    assert!(gromos.validate().is_ok());
    gromos.blocks.clear();
    assert!(gromos.validate().is_err());

    let mut lpmd: LpmdDocument = serde_json::from_value(first_fixture(
        include_str!("../../../validation/lpmd_preservation/cases.json"),
        "record",
    ))
    .unwrap();
    assert!(lpmd.validate().is_ok());
    lpmd.frames.clear();
    assert!(lpmd.validate().is_err());

    let mut unixyz: UnixyzRecord = serde_json::from_value(first_fixture(
        include_str!("../../../validation/unixyz_preservation/cases.json"),
        "expected",
    ))
    .unwrap();
    assert!(unixyz.validate().is_ok());
    unixyz.atoms[0].atomic_number = 0;
    assert!(unixyz.validate().is_err());

    let mut viewmol: ViewmolRecord = serde_json::from_value(first_fixture(
        include_str!("../../../validation/viewmol_preservation/cases.json"),
        "expected",
    ))
    .unwrap();
    assert!(viewmol.validate().is_ok());
    viewmol.factor_to_angstrom = "0".into();
    assert!(viewmol.validate().is_err());
}
