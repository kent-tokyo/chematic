use super::super::*;

// Complete baseline function; compare local stereo and the molecule graph.
fn frozen_perceive_tetrahedral_from_3d(mol: &mut Molecule, conformer: &Coords3D) {
    let mut ordinary_centers: std::collections::HashSet<_> =
        chematic_perception::stereo_centers(mol)
            .into_iter()
            .map(|(idx, _)| idx)
            .collect();
    // The CIP classifier keeps ring-path provenance and remote stereo.
    ordinary_centers.extend(chematic_chem::potential_stereocenter_indices(mol));
    ordinary_centers.extend(additional_3d_ring_centers(mol, false));
    let candidates: Vec<_> = mol
        .atoms()
        .filter(|(_, atom)| atom.chirality == Chirality::None && !atom.aromatic)
        .filter_map(|(idx, atom)| {
            let mut order: Vec<u32> = mol.neighbors(idx).map(|(neighbor, _)| neighbor.0).collect();
            order.sort_unstable();
            let neighbors: Vec<AtomIdx> = order.iter().map(|&neighbor| AtomIdx(neighbor)).collect();
            let sulfur_lone_pair = matches!(atom.element, Element::S | Element::SE)
                && order.len() == 3
                && chematic_core::implicit_hcount(mol, idx) == 0
                && chematic_perception::substituent_signatures_distinct(mol, idx, &neighbors)
                && {
                    let valence: u8 = mol
                        .neighbors(idx)
                        .map(|(_, bond)| mol.bond(bond).order.order_int())
                        .sum();
                    (atom.charge == 1 && valence == 3)
                        || (atom.charge == 0
                            && valence == 4
                            && mol.neighbors(idx).any(|(neighbor, bond)| {
                                mol.atom(neighbor).element == Element::O
                                    && mol.bond(bond).order == BondOrder::Double
                            }))
                };
            if !ordinary_centers.contains(&idx) && !sulfur_lone_pair {
                return None;
            }
            if ordinary_centers.contains(&idx)
                && !sulfur_lone_pair
                && mol.neighbors(idx).any(|(neighbor, bond)| {
                    !matches!(
                        mol.bond(bond).order,
                        BondOrder::Single | BondOrder::Up | BondOrder::Down | BondOrder::Dative
                    ) && !(atom.element == Element::P
                        && mol.bond(bond).order == BondOrder::Double
                        && mol.atom(neighbor).element == Element::O)
                })
            {
                return None;
            }
            if order.len() == 3 {
                order.push(STEREO_H_SENTINEL);
            } else if order.len() != 4 {
                return None;
            }
            Some((idx, order))
        })
        .collect();
    for (idx, order) in candidates {
        if let Some(chirality) = tetrahedral_chirality_from_3d(&order, idx, conformer) {
            mol.set_stereo_neighbor_order(idx, order);
            mol.set_chirality(idx, chirality);
        }
    }
}

fn snapshot(mol: &Molecule) -> serde_json::Value {
    serde_json::json!({
        "atoms": mol.atoms().map(|(i,a)| format!("{i:?}:{a:?}:{:?}:{:?}:{:?}",mol.stereo_neighbor_order(i),mol.atom_tag(i),mol.r_group_label(i))).collect::<Vec<_>>(),
        "bonds": mol.bonds().map(|(i,b)| format!("{i:?}:{b:?}:{:?}:{:?}",mol.bond_direction(i),mol.bond_direction_anchor(i))).collect::<Vec<_>>(),
        "groups": format!("{:?}",mol.stereo_groups()),
        "group_bonds": (0..mol.stereo_groups().len()).map(|i| format!("{:?}",mol.stereo_group_bond_indices(i))).collect::<Vec<_>>(),
    })
}

fn coordinates(count: usize, reflected: bool) -> Coords3D {
    Coords3D {
        points: (0..count)
            .map(|i| {
                let x = i as f64 * 0.1251;
                let y = ((i * 37) % 101) as f64 * 0.1234;
                let z = ((i * 73) % 97) as f64 * 0.2345;
                Point3::new(x, y, if reflected { -z } else { z })
            })
            .collect(),
    }
}

fn compare(mol: &Molecule, coords: &Coords3D) {
    let mut expected = mol.clone();
    let mut actual = mol.clone();
    frozen_perceive_tetrahedral_from_3d(&mut expected, coords);
    perceive_tetrahedral_from_3d(&mut actual, coords);
    // The intentional CTAB extension permits P=N as one tetrahedral slot.
    // Remove only newly recovered imino-P tags before comparing with the
    // frozen baseline; every other atom and graph field must still match.
    let added_imino_phosphorus: Vec<_> = actual
        .atoms()
        .filter(|(idx, atom)| {
            atom.element == Element::P
                && atom.chirality.is_tetrahedral()
                && expected.atom(*idx).chirality == Chirality::None
                && actual.neighbors(*idx).any(|(neighbor, bond)| {
                    actual.atom(neighbor).element == Element::N
                        && actual.bond(bond).order == BondOrder::Double
                })
        })
        .map(|(idx, _)| idx)
        .collect();
    for idx in added_imino_phosphorus {
        actual.set_chirality(idx, Chirality::None);
        actual.clear_stereo_neighbor_order(idx);
    }
    assert_eq!(snapshot(&actual), snapshot(&expected));
}

#[test]
fn only_imino_phosphorus_additions_change_pinned_14072_graphs() {
    let corpora = [
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
        )),
        include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../validation/cosmolkit_comparison/corpora/stress-v1.smi"
        )),
    ];
    let mut count = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let mol =
                chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            compare(&mol, &coordinates(mol.atom_count(), false));
            count += 1;
        }
    }
    assert_eq!(count, 14072);
}

#[test]
fn all_five_vertex_graphs_and_reflections_keep_original_assignments() {
    let edges: Vec<_> = (0..5)
        .flat_map(|a| ((a + 1)..5).map(move |b| (a, b)))
        .collect();
    for mask in 0..(1u32 << edges.len()) {
        let mut b = MoleculeBuilder::new();
        for _ in 0..5 {
            b.add_atom(Atom::new(Element::C));
        }
        for (i, &(a, c)) in edges.iter().enumerate() {
            if mask >> i & 1 == 1 {
                b.add_bond(AtomIdx(a), AtomIdx(c), BondOrder::Single)
                    .unwrap();
            }
        }
        let mol = b.build();
        for reflected in [false, true] {
            compare(&mol, &coordinates(5, reflected));
        }
    }
}

#[test]
fn candidate_boundaries_lone_pairs_rings_and_existing_stereo_are_preserved() {
    for smiles in [
        "C",
        "CCCCCC",
        "c1ccccc1",
        "[H][C@](F)(Cl)Br",
        "FC(Cl)(Br)I",
        "F[C@](Cl)(Br)I",
        "F[C@@](Cl)(Br)I",
        "C[N+](F)(Cl)Br",
        "N1CC2CCC1C2",
        "C1CC2CCC1C2",
        "C1C2CC3CC1CC(C2)C3",
        "ClC1CC(F)(Br)CC1I",
        "C[S+](F)Cl",
        "CS(=O)CC",
        "C[Se+](F)Cl",
        "CP(=O)(O)CC",
        "[Pt@SP1](Cl)(Br)(I)F",
        "[P](F)(Cl)(Br)(I)N",
        "[Fe](C)(N)(O)(F)(Cl)Br",
        "[13CH3][C@H](F)Cl",
        "[C@H](F)(Cl)Br",
        "C.CC.F",
        "[Na+].[Cl-]",
        "[CH2][CH2]",
        "*C(*)*",
    ] {
        let mol = chematic_smiles::parse_template(smiles).unwrap();
        for reflected in [false, true] {
            compare(&mol, &coordinates(mol.atom_count(), reflected));
        }
        compare(&mol, &Coords3D { points: Vec::new() });
        compare(
            &mol,
            &Coords3D {
                points: vec![Point3::zero(); mol.atom_count()],
            },
        );
        let mut existing = mol.clone();
        for (i, _) in mol.atoms() {
            existing.set_chirality(i, Chirality::Clockwise);
        }
        compare(&existing, &coordinates(mol.atom_count(), false));
    }
    compare(
        &MoleculeBuilder::new().build(),
        &Coords3D { points: Vec::new() },
    );
}
