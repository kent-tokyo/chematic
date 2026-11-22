use chematic_mol::parse_mmcif;

#[test]
fn incomplete_cell_is_absent_without_losing_atom_or_space_group() {
    let fields = [
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    ];
    for missing in 0..fields.len() {
        let mut text = String::from("data_test\n_symmetry.space_group_name_H-M 'P 1'\n");
        for (index, field) in fields.iter().enumerate() {
            let value = if index == missing { "?" } else { "90.0(2)" };
            text.push_str(&format!("_cell.{field} {value}\n"));
        }
        text.push_str("loop_\n_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nC 1.0 2.0 3.0\n#\n");
        let result = parse_mmcif(&text).unwrap();
        assert!(result.cell.is_none(), "missing {}", fields[missing]);
        assert_eq!(result.space_group.as_deref(), Some("P 1"));
        assert_eq!(result.atoms.len(), 1);
    }
}
