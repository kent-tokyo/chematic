use chematic_core::{AtomIdx, BondOrder, Element};
use chematic_mol::{
    CcdError, CcdParseLimits, parse_ccd_component as parse, parse_ccd_component_with_limits,
};

const WATER: &str = include_str!("../../../validation/openbabel_issues/ccd/water.cif");

#[test]
fn explicit_graph_coordinates_identifiers_and_extra_fields_survive() {
    let c = parse(WATER).unwrap();
    assert_eq!(c.component_id, "SYN");
    assert_eq!((c.atoms.len(), c.bonds.len()), (3, 2));
    assert_eq!(c.atoms[0].atom_id, "O#1");
    assert_eq!(c.atoms[0].model_coordinates, Some([1., 2., 3.]));
    assert_eq!(c.atoms[0].ideal_coordinates, Some([0., 0., 0.]));
    assert_eq!(c.atoms[0].extra_fields["_chem_comp_atom.alt_atom_id"], "O1");
    assert_eq!(c.unhandled_tags, ["_chem_comp.name"]);
    let mol = c.to_molecule(false).unwrap();
    assert_eq!((mol.atom_count(), mol.bond_count()), (3, 2));
    assert_eq!(mol.total_formula(), "H2O");
    assert!(mol.atoms().all(|(_, a)| a.hydrogen_count == Some(0)));
}

#[test]
fn singletons_charges_and_isotopes_do_not_gain_inferred_hydrogens() {
    let c = parse("data_ion\n_chem_comp.id ION\n_chem_comp_atom.comp_id ION\n_chem_comp_atom.atom_id CL\n_chem_comp_atom.type_symbol Cl\n_chem_comp_atom.charge -1\n").unwrap();
    let mol = c.to_molecule(false).unwrap();
    assert_eq!(mol.atom(AtomIdx(0)).element, Element::CL);
    assert_eq!(mol.atom(AtomIdx(0)).charge, -1);
    assert_eq!(mol.total_formula(), "Cl");
    let d = parse(&WATER.replace("H1 H 0", "H1 D 0")).unwrap();
    assert_eq!(d.atoms[1].isotope, Some(2));
    assert_eq!(
        d.to_molecule(false).unwrap().atom(AtomIdx(1)).isotope,
        Some(2)
    );
}

#[test]
fn stereo_labels_are_retained_and_graph_loss_requires_consent() {
    let text = "data_s\n_chem_comp.id S\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_stereo_config\nS C C 0 R\nS F F 0 N\nS CL Cl 0 N\nS BR Br 0 N\nS H H 0 N\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\nS C F SING\nS C CL SING\nS C BR SING\nS C H SING\n";
    let c = parse(text).unwrap();
    assert_eq!(c.atoms[0].stereo_config.as_deref(), Some("R"));
    assert_eq!(c.to_molecule(false).err().unwrap(), CcdError::StereoLoss);
    let graph = c.to_molecule(true).unwrap();
    assert_eq!(graph.bond_count(), 4);
    assert!(
        graph
            .atoms()
            .all(|(_, a)| a.chirality == Default::default())
    );
    let c = parse(&template(
        "T C1 C N N\nT C2 C N N\nT F1 F N N\nT F2 F N N\nT H1 H N N\nT H2 H N N",
        "T C1 C2 DOUB N E\nT C1 F1 SING N N\nT C2 F2 SING N N\nT C1 H1 SING N N\nT C2 H2 SING N N",
    ))
    .unwrap();
    assert_eq!(c.bonds[0].stereo_config.as_deref(), Some("E"));
    assert_eq!(c.to_molecule(false).err().unwrap(), CcdError::StereoLoss);
}

#[test]
fn aromatic_annotation_projects_aromatic_bonds_but_retains_source_orders() {
    let mut atoms = String::new();
    let mut bonds = String::new();
    for i in 0..6 {
        atoms.push_str(&format!("T C{i} C Y N\nT H{i} H N N\n"));
        let next = (i + 1) % 6;
        let order = if i % 2 == 0 { "DOUB" } else { "SING" };
        bonds.push_str(&format!(
            "T C{i} C{next} {order} Y N\nT C{i} H{i} SING N N\n"
        ));
    }
    let c = parse(&template(&atoms, &bonds)).unwrap();
    assert_eq!(
        c.bonds.iter().filter(|b| b.value_order == "DOUB").count(),
        3
    );
    let mol = c.to_molecule(false).unwrap();
    assert_eq!(mol.total_formula(), "C6H6");
    assert_eq!(
        mol.bonds()
            .filter(|(_, b)| b.order == BondOrder::Aromatic)
            .count(),
        6
    );
    assert_eq!(mol.atoms().filter(|(_, a)| a.aromatic).count(), 6);
}

#[test]
fn validates_loop_shape_and_graph_identity_without_partial_results() {
    for text in [
        WATER.replace("H2 H 0", "H1 H 0"),
        WATER.replace("'O#1' H2 SING", "UNKNOWN H2 SING"),
        WATER.replace("'O#1' H2 SING", "'O#1' 'O#1' SING"),
        WATER.replace("'O#1' H2 SING", "H1 'O#1' SING"),
        WATER.replace("SYN H1 H 0", "OTHER H1 H 0"),
        WATER.replace("SING N N", "DELO N N"),
        WATER.replace("_chem_comp_bond.atom_id_2", "_chem_comp_bond.atom_id_1"),
        format!("{WATER}\nSYN"),
        format!("{WATER}\ndata_other\n_chem_comp.id OTHER"),
        WATER.replace("H1 H 0 N N 2", "H1 Xx 0 N N 2"),
        WATER.replace("H1 H 0 N N 2", "H1 H 128 N N 2"),
    ] {
        assert!(parse(&text).is_err(), "{text}");
    }
}

#[test]
fn coordinates_refuse_nonfinite_or_incomplete_triples() {
    for value in ["NaN", "inf", "1e309", "?", "."] {
        let text = WATER.replace("N N 1 2 3", &format!("N N {value} 2 3"));
        assert!(parse(&text).is_err(), "{value}");
    }
    let c = parse(&WATER.replace("N N 1 2 3", "N N ? ? ?")).unwrap();
    assert_eq!(c.atoms[0].model_coordinates, None);
    assert_eq!(c.atoms[0].ideal_coordinates, Some([0., 0., 0.]));
}

#[test]
fn limits_cover_input_tokens_atoms_and_bonds_including_singletons() {
    for limits in [
        CcdParseLimits {
            max_input_bytes: 1,
            ..Default::default()
        },
        CcdParseLimits {
            max_tokens: 1,
            ..Default::default()
        },
        CcdParseLimits {
            max_atoms: 2,
            ..Default::default()
        },
        CcdParseLimits {
            max_bonds: 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            parse_ccd_component_with_limits(WATER, &limits),
            Err(CcdError::ResourceLimit { .. })
        ));
    }
    assert!(
        parse_ccd_component_with_limits(
            WATER,
            &CcdParseLimits {
                max_input_bytes: WATER.len(),
                max_atoms: 3,
                max_bonds: 2,
                ..Default::default()
            }
        )
        .is_ok()
    );
}

#[test]
fn lexical_boundaries_quotes_comments_case_and_crlf() {
    let c = parse(&WATER.replace('\n', "\r\n")).unwrap();
    assert_eq!(c, parse(WATER).unwrap());
    let c = parse(
        &WATER
            .replace("'O#1'", "'loop_'")
            .replace("loop_\n_chem", "LOOP_\n_chem"),
    )
    .unwrap();
    assert_eq!(c.atoms[0].atom_id, "loop_");
    assert!(parse(&WATER.replace("'O#1' O", "'unterminated O")).is_err());
    assert!(parse(&WATER.replace("\n;\nloop_", "\nnot_closed\nloop_")).is_err());
    assert!(parse("").is_err());
    assert!(parse("data_empty\n_chem_comp.id E").is_err());
}

fn template(atoms: &str, bonds: &str) -> String {
    let atoms = atoms
        .lines()
        .map(|line| {
            let mut values: Vec<_> = line.split_whitespace().collect();
            values.insert(3, "0");
            values.join(" ")
        })
        .collect::<Vec<_>>()
        .join("\n");
    format!(
        "data_T\n_chem_comp.id T\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n_chem_comp_atom.pdbx_stereo_config\n{atoms}\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n_chem_comp_bond.pdbx_stereo_config\n{bonds}\n"
    )
}

#[test]
fn uppercase_two_letter_symbols_are_elements_not_atom_name_guesses() {
    for (symbol, canonical) in [
        ("NA", "Na"),
        ("CL", "Cl"),
        ("FE", "Fe"),
        ("BE", "Be"),
        ("CE", "Ce"),
        ("HE", "He"),
    ] {
        let text = format!(
            "data_I\n_chem_comp.id I\n_chem_comp_atom.comp_id I\n_chem_comp_atom.atom_id X\n_chem_comp_atom.type_symbol {symbol}\n_chem_comp_atom.charge 0\n"
        );
        let c = parse(&text).unwrap();
        assert_eq!(c.atoms[0].element, canonical);
        assert_eq!(
            c.to_molecule(false)
                .unwrap()
                .atom(AtomIdx(0))
                .element
                .symbol(),
            canonical
        );
    }
}

#[test]
fn unknown_charge_survives_in_record_but_is_not_imputed_in_graph() {
    for marker in ["?", "."] {
        let c = parse(&WATER.replace("'O#1' O 0", &format!("'O#1' O {marker}"))).unwrap();
        assert_eq!(c.atoms[0].charge, None);
        assert!(matches!(
            c.to_molecule(true),
            Err(CcdError::UnknownCharge { .. })
        ));
    }
}

#[test]
fn public_diagnostics_explain_refusals_and_cif_table_ambiguities() {
    for (text, reason) in [
        ("not_a_data_block", "expected data block"),
        ("data_X\n_chem_comp.id", "missing singleton value"),
        (
            "data_X\n_chem_comp.id X\n_chem_comp.id X",
            "duplicate component ID",
        ),
        ("data_X\n_chem_comp.id X\nloop_", "empty loop header"),
        (
            "data_X\n_chem_comp.id X\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_bond.comp_id\nX X",
            "mixed-category loop",
        ),
        (
            "data_X\n_chem_comp.id X\n_chem_comp_atom.atom_id A\n_chem_comp_atom.atom_id B",
            "duplicate singleton tag",
        ),
        ("data_X\n_chem_comp.id 'X", "unterminated quoted value"),
        (
            "data_X\n_chem_comp.name\n;value\n;suffix",
            "text delimiter suffix",
        ),
    ] {
        assert!(parse(text).unwrap_err().to_string().contains(reason));
    }
    let limit = parse_ccd_component_with_limits(
        WATER,
        &CcdParseLimits {
            max_atoms: 2,
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(limit.to_string().contains("atoms exceeds limit 2"));
    let c = parse(&WATER.replace("'O#1' O 0", "'O#1' O ?")).unwrap();
    assert!(
        c.to_molecule(false)
            .err()
            .unwrap()
            .to_string()
            .contains("O#1")
    );
    assert!(
        CcdError::StereoLoss
            .to_string()
            .contains("allow_stereo_loss=true")
    );
    let ion = "data_I\n_chem_comp.id I\n_chem_comp_atom.comp_id I\n_chem_comp_atom.atom_id X\n_chem_comp_atom.type_symbol CL\n";
    assert!(matches!(
        parse_ccd_component_with_limits(
            ion,
            &CcdParseLimits {
                max_atoms: 0,
                ..Default::default()
            }
        ),
        Err(CcdError::ResourceLimit { .. })
    ));
    for replacement in ["H1 H 0 Z N 2", "H1 H 0 N E 2"] {
        assert!(parse(&WATER.replace("H1 H 0 N N 2", replacement)).is_err());
    }
    assert!(
        parse(&format!("{WATER}\n_chem_comp_atom.atom_id OTHER"))
            .unwrap_err()
            .to_string()
            .contains("mixed singleton/loop")
    );
}

#[test]
fn long_repeated_column_names_refuse_before_row_expansion() {
    let tag = format!("_chem_comp_atom.extra_{}", "x".repeat(20_000));
    let mut text = format!(
        "data_A\n_chem_comp.id A\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n{tag}\n"
    );
    for i in 0..2000 {
        text.push_str(&format!("A C{i} C v\n"));
    }
    assert!(text.len() < 100_000);
    let error = parse(&text).unwrap_err();
    assert!(matches!(
        error,
        CcdError::ResourceLimit {
            resource: "expanded loop bytes",
            ..
        }
    ));
}
