// ===========================================================================
// Tests
//
// Only happy paths are exercised here: constructing a `JsValue` (as every
// `Err` branch above does) outside a real wasm/JS runtime aborts the native
// test process rather than returning an `Err` -- the same constraint
// documented on `crate::mol_io::extxyz_frame_from_json_args` and the
// `run_reactants` error-path note in `crate::tests`. Error-path coverage of
// each hand-written JSON<->struct mapping is implicit: every `Ok` path
// below round-trips through the exact same field-extraction code an error
// would come from.
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const CUBE_2X2X2: &str = "Water density\n\
Generated for chematic tests\n\
1    0.000000    0.000000    0.000000\n\
2    1.000000    0.000000    0.000000\n\
2    0.000000    1.000000    0.000000\n\
2    0.000000    0.000000    1.000000\n\
8    8.000000    0.500000    0.500000    0.500000\n\
0.0 1.0 2.0 3.0\n\
4.0 5.0 6.0 7.0\n";

    const OPENDX_2X2X2: &str = "object 1 class gridpositions counts 2 2 2\n\
origin -1.0 -1.0 -1.0\n\
delta 0.5 0.0 0.0\n\
delta 0.0 0.5 0.0\n\
delta 0.0 0.0 0.5\n\
object 2 class gridconnections counts 2 2 2\n\
object 3 class array type double rank 0 items 8 data follows\n\
0.0 1.0 2.0\n\
3.0 4.0 5.0\n\
6.0 7.0\n\
attribute \"dep\" string \"positions\"\n\
object \"regular positions regular connections\" class field\n\
component \"positions\" value 1\n\
component \"connections\" value 2\n\
component \"data\" value 3\n";

    #[test]
    fn cube_round_trip_via_json() {
        assert_eq!(mol_from_cube(CUBE_2X2X2).unwrap().atom_count(), 1);
        let json1 = cube_grid_json(CUBE_2X2X2).unwrap();
        let text2 = write_cube_json(&json1).unwrap();
        let json2 = cube_grid_json(&text2).unwrap();
        assert_eq!(json1, json2);
    }

    #[test]
    fn opendx_round_trip_via_json() {
        let json1 = opendx_grid_json(OPENDX_2X2X2).unwrap();
        // parse_opendx always yields Angstrom (no unit tag in the format),
        // so the fail-closed default writer must accept it directly.
        let text2 = write_opendx_json(&json1).unwrap();
        let json2 = opendx_grid_json(&text2).unwrap();
        assert_eq!(json1, json2);
        // The lossy writer must also accept an already-Angstrom grid.
        let text3 = write_opendx_lossy_json(&json1).unwrap();
        assert_eq!(opendx_grid_json(&text3).unwrap(), json1);
    }

    const MMCIF_ATOMS_JSON: &str = r#"[
        {"group_pdb":"ATOM","serial":1,"element":"C","atom_name":"CA","alt_loc":null,
         "res_name":"ALA","chain_id":"A","res_seq":1,"label_seq_id":1,"icode":null,
         "x":1.0,"y":2.0,"z":3.0,"occupancy":1.0,"b_iso":20.0,"formal_charge":null,
         "entity_id":"1","model_num":1},
        {"group_pdb":"HETATM","serial":2,"element":"Zn","atom_name":"ZN","alt_loc":null,
         "res_name":"ZN","chain_id":"A","res_seq":2,"label_seq_id":null,"icode":"A",
         "x":4.0,"y":5.0,"z":6.0,"occupancy":0.5,"b_iso":15.0,"formal_charge":2,
         "entity_id":null,"model_num":1}
    ]"#;

    #[test]
    fn mmcif_round_trip_via_json() {
        let text = write_mmcif_json(
            MMCIF_ATOMS_JSON,
            r#"{"a":10.0,"b":10.0,"c":10.0,"alpha":90.0,"beta":90.0,"gamma":90.0}"#,
            "P1",
            "TEST",
        )
        .unwrap();
        assert_eq!(mol_from_mmcif(&text).unwrap().atom_count(), 2);
        let coords: Vec<[f64; 3]> =
            serde_json::from_str(&mmcif_coords_json(&text).unwrap()).unwrap();
        assert_eq!(coords, vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let out = mmcif_to_json(&text).unwrap();
        assert!(out.contains("\"formal_charge\":2"));
        assert!(out.contains("\"space_group\":\"P1\""));
    }

    const PQR_ATOMS_JSON: &str = r#"[
        {"group_pdb":"ATOM","serial":1,"atom_name":"N","res_name":"ALA","chain_id":null,
         "res_seq":1,"icode":null,"x":1.0,"y":2.0,"z":3.0,"charge":-0.4,"radius":1.5,"element":"N"},
        {"group_pdb":"HETATM","serial":2,"atom_name":"ZN","res_name":"ZN","chain_id":"A",
         "res_seq":2,"icode":null,"x":4.0,"y":5.0,"z":6.0,"charge":2.0,"radius":1.09,"element":"Zn"}
    ]"#;

    #[test]
    fn pqr_round_trip_via_json() {
        let text = write_pqr_json(PQR_ATOMS_JSON).unwrap();
        assert_eq!(mol_from_pqr(&text).unwrap().atom_count(), 2);
        let out = pqr_to_json(&text).unwrap();
        assert!(out.contains("\"radius\":1.09"));
        assert!(out.contains("\"chain_id\":\"A\""));
    }

    #[test]
    fn pqr_infer_element_matches_rust_api() {
        assert_eq!(pqr_infer_element("ATOM", "ALA", "CA").as_deref(), Some("C"));
        assert_eq!(
            pqr_infer_element("HETATM", "ZN", "ZN").as_deref(),
            Some("Zn")
        );
    }

    const QCSCHEMA_WATER_JSON: &str = r#"{
        "schema_name": "qcschema_molecule",
        "schema_version": 1,
        "symbols": ["O", "H", "H"],
        "geometry": [0.0, 0.0, 0.0, 0.0, 0.0, 1.8, 0.0, 1.7, -0.5],
        "molecular_charge": 0.0,
        "molecular_multiplicity": 1
    }"#;

    #[test]
    fn qcschema_molecule_round_trip() {
        let handle = mol_from_qcschema_molecule(QCSCHEMA_WATER_JSON).unwrap();
        assert_eq!(handle.atom_count(), 3);
        let coords_json = qcschema_molecule_coords_json(QCSCHEMA_WATER_JSON).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&coords_json).unwrap();
        let coords_only = serde_json::to_string(&parsed["coords"]).unwrap();
        let written = to_qcschema_molecule_json(&handle, &coords_only, 0.0, 1).unwrap();
        let handle2 = mol_from_qcschema_molecule(&written).unwrap();
        assert_eq!(handle2.atom_count(), 3);
    }

    #[test]
    fn qcschema_validate_atomic_input_round_trip() {
        let json = format!(
            r#"{{"schema_name":"qcschema_input","schema_version":1,"molecule":{QCSCHEMA_WATER_JSON},
                "driver":"energy","model":{{"method":"b3lyp","basis":"def2-svp"}},"keywords":{{}}}}"#
        );
        let out = qcschema_validate_atomic_input(&json).unwrap();
        assert!(out.contains("b3lyp"));
    }

    #[test]
    fn qcschema_validate_atomic_input_preserves_unknown_fields() {
        // `extras` (a spec-defined open bag) and a wholly unrecognized
        // top-level key must both still be present after the
        // parse_atomic_input -> write_atomic_input round trip --
        // chematic_mol::AtomicInput::{extras,unknown_fields} already
        // implement this; this test pins that guarantee at the WASM
        // binding level specifically (not just "a known field survives",
        // which qcschema_validate_atomic_input_round_trip above already
        // covers).
        let json = format!(
            r#"{{"schema_name":"qcschema_input","schema_version":1,"molecule":{QCSCHEMA_WATER_JSON},
                "driver":"energy","model":{{"method":"b3lyp","basis":"def2-svp"}},
                "extras":{{"vendor_tag":"keep-me"}},
                "x_vendor_extension":{{"nested":42}}}}"#
        );
        let out = qcschema_validate_atomic_input(&json).unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["extras"]["vendor_tag"], "keep-me");
        assert_eq!(v["x_vendor_extension"]["nested"], 42);
    }

    #[test]
    fn qcschema_validate_atomic_result_preserves_unknown_fields() {
        // chematic_mol::AtomicResult has the same extras/unknown_fields
        // mechanism as AtomicInput (both are JsonObject fields collected
        // via collect_unknown) -- same guarantee, pinned independently
        // since AtomicResult has its own parse/write pair.
        let json = format!(
            r#"{{"schema_name":"qcschema_output","schema_version":1,"molecule":{QCSCHEMA_WATER_JSON},
                "driver":"energy","model":{{"method":"b3lyp","basis":"def2-svp"}},
                "provenance":{{"creator":"test"}},"success":true,"return_result":-76.3,
                "extras":{{"vendor_tag":"keep-me"}},
                "x_vendor_extension":{{"nested":42}}}}"#
        );
        let out = qcschema_validate_atomic_result(&json).unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["extras"]["vendor_tag"], "keep-me");
        assert_eq!(v["x_vendor_extension"]["nested"], 42);
    }

    const ORCA_INPUT_JSON: &str = r##"{
        "comments": ["# a comment line"],
        "keywords": ["B3LYP", "def2-SVP"],
        "blocks": [{"name":"pal","raw":"nprocs 4","has_end":true}],
        "coords": {"type":"xyz","charge":0,"multiplicity":1,"atoms":[
            {"element":"O","x":0.0,"y":0.0,"z":0.0,"frozen":[false,false,false],"extra":null},
            {"element":"H","x":0.0,"y":0.0,"z":0.96,"frozen":[false,false,false],"extra":null}
        ]}
    }"##;

    #[test]
    fn orca_input_round_trip_via_json() {
        let text = write_orca_input_json(ORCA_INPUT_JSON).unwrap();
        assert_eq!(mol_from_orca_input(&text).unwrap().atom_count(), 2);
        let out = orca_input_to_json(&text).unwrap();
        assert!(out.contains("\"B3LYP\""));
        assert!(out.contains("\"type\":\"xyz\""));
    }

    #[test]
    fn orca_output_to_json_extracts_energy_and_trajectory() {
        let text = "\
Total Charge           Charge          ....    0
Multiplicity           Mult            ....    1
---------------------------------
CARTESIAN COORDINATES (ANGSTROEM)
---------------------------------
  O     0.000000    0.000000    0.000000
  H     0.000000    0.000000    0.960000
-------------------------   --------------------
FINAL SINGLE POINT ENERGY       -76.320145981234
-------------------------   --------------------
****ORCA TERMINATED NORMALLY****
";
        let out = orca_output_to_json(text).unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["final_energy_hartree"].as_f64(), Some(-76.320145981234));
        assert_eq!(v["trajectory"].as_array().unwrap().len(), 1);
        assert_eq!(
            v["trajectory"][0]["elements"],
            serde_json::json!(["O", "H"])
        );
        assert_eq!(v["termination"]["kind"], "normal");
    }

    const LAMMPS_DATA_JSON: &str = r#"{
        "counts": [["atoms", 2], ["atom types", 1]],
        "atom_style": "atomic",
        "simulation_box": {"lo":[0.0,0.0,0.0],"hi":[10.0,10.0,10.0],"tilt":null},
        "masses": [{"atom_type":1,"mass":12.011}],
        "atoms": [
            {"id":1,"molecule_id":null,"atom_type":1,"charge":null,"x":1.0,"y":1.0,"z":1.0,"image":null},
            {"id":2,"molecule_id":null,"atom_type":1,"charge":null,"x":2.0,"y":2.0,"z":2.0,"image":[1,0,0]}
        ],
        "velocities": [],
        "bonds": [],
        "unparsed_sections": []
    }"#;

    #[test]
    fn lammps_data_round_trip_via_json() {
        let text = write_lammps_data_json(LAMMPS_DATA_JSON).unwrap();
        let out = lammps_data_to_json(&text, "atomic").unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["atoms"].as_array().unwrap().len(), 2);
        assert_eq!(v["atoms"][1]["image"], serde_json::json!([1, 0, 0]));
        assert_eq!(v["atom_style"], "atomic");
    }

    // Note: an "unsupported atom_style" error-path test is omitted here --
    // like `run_reactants`'s error-path tests in `crate::tests`,
    // constructing the resulting `JsValue` outside a real wasm/JS runtime
    // aborts the native test process. `chematic_mol::LammpsDataError::
    // UnsupportedAtomStyle` coverage lives in chematic-mol's own unit tests;
    // `lammps_atom_style_from_str` (this file's only new logic in that
    // path) is a direct 1:1 string match with no branch of its own to miss.

    const LAMMPS_DUMP_FRAME_JSON: &str = r#"{
        "timestep": 100,
        "num_atoms": 2,
        "box_bounds": {"lo":[0.0,0.0,0.0],"hi":[10.0,10.0,10.0],"tilt":null},
        "boundary_flags": ["pp","pp","pp"],
        "column_names": ["id","x","y","z"],
        "rows": [[1.0,1.0,1.0,1.0],[2.0,2.0,2.0,2.0]]
    }"#;

    #[test]
    fn lammps_dump_frame_round_trip_via_json() {
        let text = write_lammps_dump_frame_json(LAMMPS_DUMP_FRAME_JSON).unwrap();
        let out = lammps_dump_frame_to_json_str(&text).unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v["timestep"], 100);
        assert_eq!(v["rows"].as_array().unwrap().len(), 2);
    }

    #[test]
    fn lammps_trajectory_round_trip_multi_frame() {
        let frames_json = format!("[{LAMMPS_DUMP_FRAME_JSON},{LAMMPS_DUMP_FRAME_JSON}]");
        let text = write_lammps_trajectory_json(&frames_json).unwrap();
        let out = lammps_trajectory_to_json(&text).unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v.as_array().unwrap().len(), 2);
    }

    // Fixtures below reuse the exact box/tilt/xs values from
    // chematic_mol::lammps_dump's own `orthogonal_frame`/`triclinic_frame`
    // test fixtures and their `cartesian_positions_passes_through_x_y_z`/
    // `scaled_coordinates_triclinic_hand_computed`/
    // `cartesian_positions_none_for_unwrapped_only_columns` tests, so the
    // hand-computed expected values are shared with (not re-derived from)
    // that module's own oracle.

    const LAMMPS_DUMP_XYZ_FRAME_JSON: &str = r#"{
        "timestep": 1000,
        "num_atoms": 2,
        "box_bounds": {"lo":[0.0,0.0,0.0],"hi":[10.0,20.0,30.0],"tilt":null},
        "boundary_flags": ["pp","pp","pp"],
        "column_names": ["id","type","x","y","z"],
        "rows": [[1.0,1.0,1.0,2.0,3.0],[2.0,1.0,4.0,5.0,6.0]]
    }"#;

    #[test]
    fn lammps_dump_cartesian_positions_passthrough_xyz() {
        let out = lammps_dump_cartesian_positions_json(LAMMPS_DUMP_XYZ_FRAME_JSON).unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(v, serde_json::json!([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]));
    }

    // box lo=[0,0,0] hi=[10,10,10], tilt xy=2 xz=1 yz=0.5 (a genuinely
    // triclinic box -- all 3 tilt components nonzero), xs=ys=zs=0.5.
    const LAMMPS_DUMP_TRICLINIC_XS_FRAME_JSON: &str = r#"{
        "timestep": 2000,
        "num_atoms": 1,
        "box_bounds": {"lo":[0.0,0.0,0.0],"hi":[10.0,10.0,10.0],"tilt":[2.0,1.0,0.5]},
        "boundary_flags": ["pp","ff","ss"],
        "column_names": ["id","xs","ys","zs"],
        "rows": [[1.0,0.5,0.5,0.5]]
    }"#;

    #[test]
    fn lammps_dump_cartesian_positions_triclinic_hand_computed() {
        // x = xlo + xs*(xhi-xlo) + ys*xy + zs*xz = 0 + 0.5*10 + 0.5*2 + 0.5*1 = 5 + 1 + 0.5 = 6.5
        // y = ylo + ys*(yhi-ylo) + zs*yz         = 0 + 0.5*10 + 0.5*0.5       = 5 + 0.25    = 5.25
        // z = zlo + zs*(zhi-zlo)                 = 0 + 0.5*10                 = 5.0
        let out =
            lammps_dump_cartesian_positions_json(LAMMPS_DUMP_TRICLINIC_XS_FRAME_JSON).unwrap();
        let v: serde_json::Value = serde_json::from_str(&out).unwrap();
        let positions = v.as_array().unwrap();
        assert_eq!(positions.len(), 1);
        let p = positions[0].as_array().unwrap();
        assert!((p[0].as_f64().unwrap() - 6.5).abs() < 1e-9);
        assert!((p[1].as_f64().unwrap() - 5.25).abs() < 1e-9);
        assert!((p[2].as_f64().unwrap() - 5.0).abs() < 1e-9);
    }

    const LAMMPS_DUMP_UNWRAPPED_ONLY_FRAME_JSON: &str = r#"{
        "timestep": 1000,
        "num_atoms": 2,
        "box_bounds": {"lo":[0.0,0.0,0.0],"hi":[10.0,20.0,30.0],"tilt":null},
        "boundary_flags": ["pp","pp","pp"],
        "column_names": ["id","xu","yu","zu"],
        "rows": [[1.0,100.0,200.0,300.0],[2.0,1.0,2.0,3.0]]
    }"#;

    #[test]
    fn lammps_dump_cartesian_positions_none_for_unwrapped_only_is_json_null() {
        // "unwrapped" coordinates are a different physical quantity from a
        // scaled coordinate and must never be silently treated as one --
        // cartesian_positions() returns None here, which this binding must
        // surface as JSON `null`, not an error and not `[]`.
        let out =
            lammps_dump_cartesian_positions_json(LAMMPS_DUMP_UNWRAPPED_ONLY_FRAME_JSON).unwrap();
        assert_eq!(out, "null");
    }

    // -----------------------------------------------------------------
    // Typed-array accessor helpers -- exercised natively against their
    // pure-Rust helper functions (not the `#[wasm_bindgen]`-exposed
    // Float64Array/Uint32Array functions themselves, which construct real
    // `js_sys` values and would abort a native test process; see this
    // file's "Typed-array accessors" section header comment). Each
    // assertion cross-checks the helper's output against the JSON
    // sibling's own already-tested output, the strongest available
    // same-input cross-check.
    // -----------------------------------------------------------------

    #[test]
    fn cube_values_vec_matches_cube_grid_json() {
        let values = cube_values_vec(CUBE_2X2X2).unwrap();
        let json = cube_grid_json(CUBE_2X2X2).unwrap();
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let expected: Vec<f64> = serde_json::from_value(v["values"].clone()).unwrap();
        assert_eq!(values, expected);
        assert_eq!(values, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);
    }

    #[test]
    fn cube_shape_vec_matches_cube_grid_json() {
        let shape = cube_shape_vec(CUBE_2X2X2).unwrap();
        let json = cube_grid_json(CUBE_2X2X2).unwrap();
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let expected: [usize; 3] = serde_json::from_value(v["shape"].clone()).unwrap();
        assert_eq!(shape, [2u32, 2, 2]);
        assert_eq!(shape.map(|d| d as usize), expected);
    }

    #[test]
    fn opendx_values_vec_matches_opendx_grid_json() {
        let values = opendx_values_vec(OPENDX_2X2X2).unwrap();
        let json = opendx_grid_json(OPENDX_2X2X2).unwrap();
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let expected: Vec<f64> = serde_json::from_value(v["values"].clone()).unwrap();
        assert_eq!(values, expected);
        assert_eq!(values, vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]);
    }

    #[test]
    fn opendx_shape_vec_matches_opendx_grid_json() {
        let shape = opendx_shape_vec(OPENDX_2X2X2).unwrap();
        assert_eq!(shape, [2u32, 2, 2]);
    }

    #[test]
    fn dump_rows_flat_matches_frame_json_rows() {
        let flat = dump_rows_flat(LAMMPS_DUMP_FRAME_JSON).unwrap();
        // LAMMPS_DUMP_FRAME_JSON: column_names has 4 entries, 2 atoms.
        assert_eq!(flat, vec![1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0]);
    }

    #[test]
    fn dump_cartesian_flat_matches_cartesian_positions_json_xyz_passthrough() {
        let flat = dump_cartesian_flat(LAMMPS_DUMP_XYZ_FRAME_JSON).unwrap();
        assert_eq!(flat, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }

    #[test]
    fn dump_cartesian_flat_matches_cartesian_positions_json_triclinic_hand_computed() {
        // Same hand-computed expected values as
        // `lammps_dump_cartesian_positions_triclinic_hand_computed` above.
        let flat = dump_cartesian_flat(LAMMPS_DUMP_TRICLINIC_XS_FRAME_JSON).unwrap();
        assert_eq!(flat.len(), 3);
        assert!((flat[0] - 6.5).abs() < 1e-9);
        assert!((flat[1] - 5.25).abs() < 1e-9);
        assert!((flat[2] - 5.0).abs() < 1e-9);
    }

    #[test]
    fn dump_cartesian_flat_errs_for_unwrapped_only_columns() {
        // The JSON sibling (lammps_dump_cartesian_positions_json) returns
        // JSON `null` for this same fixture -- a `Float64Array` has no
        // `null`, so this helper (and the Float64Array-returning
        // #[wasm_bindgen] function built on it) returns `Err` instead. See
        // lammps_dump_cartesian_positions_f64's doc comment.
        let err = dump_cartesian_flat(LAMMPS_DUMP_UNWRAPPED_ONLY_FRAME_JSON).unwrap_err();
        assert!(err.contains("no recognized coordinate columns"));
    }
}
