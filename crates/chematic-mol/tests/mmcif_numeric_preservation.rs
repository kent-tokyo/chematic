use chematic_mol::{MmcifAtomRecord, UnitCell, parse_mmcif, write_mmcif};
use serde_json::Value;

fn atom_bits(atom: &MmcifAtomRecord) -> Vec<u64> {
    [atom.x, atom.y, atom.z, atom.occupancy, atom.b_iso]
        .map(f64::to_bits)
        .to_vec()
}

fn cell_bits(cell: &UnitCell) -> Vec<u64> {
    [cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma]
        .map(f64::to_bits)
        .to_vec()
}

#[test]
fn authored_numeric_packet_round_trips_without_decimal_rounding() {
    let packet: Value = serde_json::from_str(include_str!(
        "../../../validation/mmcif_numeric_preservation/cases.json"
    ))
    .unwrap();
    for case in packet["cases"].as_array().unwrap() {
        let before = parse_mmcif(case["input"].as_str().unwrap()).unwrap();
        let text = write_mmcif(
            &before.atoms,
            before.cell.as_ref(),
            before.space_group.as_deref(),
            "test",
        );
        let after = parse_mmcif(&text).unwrap();
        assert_eq!(before.atoms, after.atoms, "{}", case["id"]);
        assert_eq!(before.space_group, after.space_group);
        for (i, (a, b)) in before.atoms.iter().zip(&after.atoms).enumerate() {
            let expected: Vec<_> = ["x", "y", "z", "occupancy", "b_iso"]
                .map(|key| {
                    case["expected"]["atoms"][i][key]
                        .as_f64()
                        .unwrap()
                        .to_bits()
                })
                .to_vec();
            assert_eq!(atom_bits(a), expected, "{} source", case["id"]);
            assert_eq!(atom_bits(b), expected, "{} output", case["id"]);
        }
        assert_eq!(
            before.cell.as_ref().map(cell_bits),
            after.cell.as_ref().map(cell_bits)
        );
        if let Some(cell) = after.cell {
            let expected: Vec<_> = ["a", "b", "c", "alpha", "beta", "gamma"]
                .map(|key| case["expected"]["cell"][key].as_f64().unwrap().to_bits())
                .to_vec();
            assert_eq!(cell_bits(&cell), expected);
        }
    }
}

#[test]
fn finite_binary64_storage_boundaries_and_signed_zero_are_retained() {
    let mut parsed = parse_mmcif("data_x\nloop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nN 0 0 0\n").unwrap();
    for value in [
        0.0,
        -0.0,
        f64::from_bits(1),
        -f64::from_bits(1),
        f64::MIN_POSITIVE,
        f64::MAX,
        f64::MIN,
        f64::from_bits(1.0_f64.to_bits() + 1),
        1e-200,
        -1e200,
    ] {
        let atom = &mut parsed.atoms[0];
        atom.x = value;
        atom.y = value;
        atom.z = value;
        atom.occupancy = value;
        atom.b_iso = value;
        let output = write_mmcif(&parsed.atoms, None, None, "storage");
        let row = output
            .lines()
            .find(|line| line.starts_with("ATOM "))
            .unwrap();
        assert!(row.len() <= 256, "numeric magnitude must not inflate a row");
        let after = parse_mmcif(&output).unwrap();
        assert_eq!(atom_bits(&after.atoms[0]), vec![value.to_bits(); 5]);
    }
}
