use super::*;
use chematic_core::{Atom, Element, MoleculeBuilder};

fn frozen_stereo_centers(mol: &Molecule) -> Vec<(AtomIdx, bool)> {
    let ring_set = find_sssr(mol);
    let mut centers = Vec::new();

    for (idx, atom) in mol.atoms() {
        // Skip aromatics and obvious non-centers.
        if atom.aromatic {
            continue;
        }

        let heavy_nbs: Vec<AtomIdx> = mol
            .neighbors(idx)
            .filter(|(nb, _)| mol.atom(*nb).element.atomic_number() != 1)
            .map(|(nb, _)| nb)
            .collect();
        let implicit_h = chematic_core::implicit_hcount(mol, idx) as usize;
        let is_ring_nitrogen = atom.element.atomic_number() == 7
            && heavy_nbs.len() == 3
            && implicit_h == 0
            && ring_set.contains_atom(idx);
        let groups = heavy_nbs.len() + implicit_h + usize::from(is_ring_nitrogen);

        if groups != 4 {
            continue;
        } // only tetrahedral candidates

        // Compare heavy-atom branches from this center. Implicit H and the
        // lone-pair group of a constrained ring nitrogen are fourth groups by
        // construction and do not need a graph signature.
        if !substituent_signatures_distinct(mol, idx, &heavy_nbs) {
            continue;
        }
        if implicit_h > 0 && heavy_nbs.len() == 3 {
            // An implicit H is a fourth group and is distinct from every
            // heavy-atom branch. A ring nitrogen's fourth group is its lone
            // pair, handled by the same branch-distinctness condition.
        } else if heavy_nbs.len() != 4 && !is_ring_nitrogen {
            continue;
        }

        centers.push((idx, atom.chirality != Chirality::None));
    }

    centers
}

fn compare(mol: &Molecule) {
    let before: Vec<_> = mol
        .atoms()
        .map(|(i, a)| format!("{i:?}:{a:?}:{:?}", mol.stereo_neighbor_order(i)))
        .collect();
    assert_eq!(stereo_centers(mol), frozen_stereo_centers(mol));
    assert_eq!(
        before,
        mol.atoms()
            .map(|(i, a)| format!("{i:?}:{a:?}:{:?}", mol.stereo_neighbor_order(i)))
            .collect::<Vec<_>>()
    );
}

#[test]
fn all_14072_pinned_graphs_match_frozen_classifier() {
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
            compare(
                &chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap(),
            );
            count += 1;
        }
    }
    assert_eq!(count, 14072);
}

#[test]
fn every_five_atom_graph_and_hydrogen_override_matches() {
    let edges: Vec<_> = (0..5)
        .flat_map(|a| ((a + 1)..5).map(move |b| (a, b)))
        .collect();
    for mask in 0..(1u32 << edges.len()) {
        for hydrogen_count in [None, Some(0), Some(1), Some(2), Some(3), Some(4)] {
            let mut b = MoleculeBuilder::new();
            for _ in 0..5 {
                let mut atom = Atom::new(Element::C);
                atom.hydrogen_count = hydrogen_count;
                b.add_atom(atom);
            }
            for (i, &(a, c)) in edges.iter().enumerate() {
                if mask >> i & 1 == 1 {
                    b.add_bond(AtomIdx(a), AtomIdx(c), BondOrder::Single)
                        .unwrap();
                }
            }
            compare(&b.build());
        }
    }
}

#[test]
fn stereocenter_boundary_groups_and_ring_nitrogen_match() {
    for smiles in [
        "C",
        "CCCCCC",
        "CC(F)CCCCCCCCCCCCCCCC",
        "c1ccccc1",
        "N[C@@H](C)C(=O)O",
        "[H][C@](F)(Cl)Br",
        "FC(Cl)(Br)I",
        "F[C@](Cl)(Br)I",
        "F[C@@](Cl)(Br)I",
        "C[N+](F)(Cl)Br",
        "N1CC2CCC1C2",
        "C1CC2CCC1C2",
        "C1C2CC3CC1CC(C2)C3",
        "C[S+](CC)CCC",
        "CS(=O)CC",
        "C[Se+](CC)CCC",
        "C[Se](=O)CC",
        "[13CH3]C(F)C",
        "[CH2]",
        "[C]",
        "[NH4+]",
        "[H]N1CC2CCC1C2",
        "[Fe]",
        "ClC1CC(F)(Br)CC1I",
    ] {
        compare(&chematic_smiles::parse_template(smiles).unwrap());
    }
}
