use chematic_core::{AtomIdx, BondIdx, StereoGroup, StereoGroupKind};
use chematic_rxn::{
    Reaction, ReactionParseLimits, parse_reaction_cxsmiles, parse_reaction_cxsmiles_with_limits,
    write_reaction_cxsmiles,
};

fn fingerprint(reaction: &Reaction) -> Vec<String> {
    let mut result = Vec::new();
    for (role, molecules) in [&reaction.reactants, &reaction.agents, &reaction.products]
        .iter()
        .enumerate()
    {
        for mol in *molecules {
            let cip = chematic_chem::assign_cip(mol);
            for group in mol.stereo_groups() {
                for atom in &group.atom_indices {
                    assert!(mol.atom(*atom).chirality.is_tetrahedral());
                    result.push(format!(
                        "{role}:{:?}:{:?}:{:?}",
                        group.kind,
                        mol.atom(*atom).atom_map,
                        cip.get(*atom)
                    ));
                }
            }
        }
    }
    result.sort();
    result
}

#[test]
fn issue9732_group_follows_stereocentre_after_component_reordering() {
    let reaction = parse_reaction_cxsmiles("O.C[C@H](O)CC>>CCC(C)=O |&1:2|").unwrap();
    let output = write_reaction_cxsmiles(&reaction, true).unwrap();
    let reread = parse_reaction_cxsmiles(&output).unwrap();
    assert_eq!(fingerprint(&reaction), fingerprint(&reread), "{output}");
    assert!(reread.reactants[0].atom_count() > 1, "{output}");
    assert_eq!(reread.reactants[0].stereo_groups().len(), 1);
    assert!(reread.reactants[1].stereo_groups().is_empty());
}

#[test]
fn authored_packet_preserves_roles_map_identity_cip_and_group_membership() {
    let packet: serde_json::Value = serde_json::from_str(include_str!(
        "../../../validation/rdkit-9732-2026-10-24/cases.json"
    ))
    .unwrap();
    for case in packet["positive"].as_array().unwrap() {
        let input = case["input"].as_str().unwrap();
        let reaction = parse_reaction_cxsmiles(input).unwrap();
        for canonical in [false, true] {
            let output = write_reaction_cxsmiles(&reaction, canonical).unwrap();
            let reread = parse_reaction_cxsmiles(&output).unwrap();
            assert_eq!(
                fingerprint(&reaction),
                fingerprint(&reread),
                "{input} -> {output}"
            );
            assert_eq!(
                write_reaction_cxsmiles(&reread, canonical).unwrap(),
                output,
                "{}",
                case["id"]
            );
        }
    }
    for case in packet["negative"].as_array().unwrap() {
        assert!(
            parse_reaction_cxsmiles(case["input"].as_str().unwrap()).is_err(),
            "{}",
            case["id"]
        );
    }
}

#[test]
fn reordered_equivalent_input_and_atom_renumbering_keep_group_identity() {
    let a = parse_reaction_cxsmiles("O.C[C@H](O)CC>>CCC(C)=O |&1:2|").unwrap();
    let b = parse_reaction_cxsmiles("C[C@H](O)CC.O>>CCC(C)=O |&1:1|").unwrap();
    assert_eq!(
        write_reaction_cxsmiles(&a, true).unwrap(),
        write_reaction_cxsmiles(&b, true).unwrap()
    );
    let mut reaction =
        parse_reaction_cxsmiles("O.[CH3:10][C@H:11](F)[C@@H:12](O)Cl>>C |&3:2,4|").unwrap();
    let original = fingerprint(&reaction);
    let mol = &reaction.reactants[1];
    let order = (0..mol.atom_count())
        .rev()
        .map(|i| AtomIdx(i as u32))
        .collect::<Vec<_>>();
    reaction.reactants[1] = mol.renumber_atoms(&order).unwrap();
    let output = write_reaction_cxsmiles(&reaction, true).unwrap();
    assert_eq!(
        fingerprint(&parse_reaction_cxsmiles(&output).unwrap()),
        original
    );
}

#[test]
fn limits_include_extension_and_programmatic_invalid_groups_are_refused() {
    let input = "O.C[C@H](O)CC>>C |&1:2|";
    let limits = ReactionParseLimits {
        max_input_bytes: input.len() - 1,
        ..Default::default()
    };
    assert!(parse_reaction_cxsmiles_with_limits(input, &limits).is_err());
    for index in [AtomIdx(0), AtomIdx(999)] {
        let mut reaction = parse_reaction_cxsmiles("C[C@H](F)Cl>>C").unwrap();
        reaction.reactants[0]
            .set_stereo_groups(vec![StereoGroup::new(StereoGroupKind::And(1), vec![index])]);
        assert!(write_reaction_cxsmiles(&reaction, true).is_err());
    }
    let mut reaction = parse_reaction_cxsmiles("C[C@H](F)Cl>>C").unwrap();
    reaction.reactants[0].set_stereo_groups_with_bonds(vec![(
        StereoGroup::new(StereoGroupKind::And(1), vec![AtomIdx(1)]),
        vec![BondIdx(0)],
    )]);
    assert!(write_reaction_cxsmiles(&reaction, true).is_err());
}
