//! Frozen component facts; existing pinned 10k and stress inputs only.
use super::*;

struct ReferenceInfo {
    cyclic: std::sync::Arc<Vec<bool>>,
    /// Component label of every atom with a cyclic bond (`u32::MAX` otherwise).
    root: Vec<u32>,
    /// Per label: `COMP_*` bits.
    kind: Vec<u8>,
    /// Per label: cyclic bonds and atoms.
    edges: Vec<u32>,
    atoms: Vec<u32>,
    has_aromatic_bond: bool,
    /// Every aromatic bond joins two flagged atoms (and one exists), and no
    /// flagged neutral group-16 atom is hypercoordinate.
    explicit: bool,
    non_ring_aromatic_bond: bool,
    any_flagged: bool,
}

fn reference_info(mol: &Molecule) -> ReferenceInfo {
    let n = mol.atom_count();
    let cyclic = crate::sssr::ring_bond_flags_shared(mol);
    let components = crate::sssr::ring_components_shared(mol);
    let c = components.count as usize;
    let mut root = vec![u32::MAX; n];
    let mut kind = vec![0u8; c];
    let mut edges = vec![0u32; c];
    let mut atoms = vec![0u32; c];
    let mut has_aromatic_bond = false;
    let mut consistent = true;
    let mut non_ring_aromatic_bond = false;
    for (bond_idx, bond) in mol.bonds() {
        if bond.order == BondOrder::Aromatic {
            has_aromatic_bond = true;
            if !mol.atom(bond.atom1).aromatic || !mol.atom(bond.atom2).aromatic {
                consistent = false;
            }
        }
        if cyclic[bond_idx.0 as usize] {
            let r = components.label[bond.atom1.0 as usize];
            root[bond.atom1.0 as usize] = r;
            root[bond.atom2.0 as usize] = r;
            let r = r as usize;
            edges[r] += 1;
            kind[r] |= match bond.order {
                BondOrder::Aromatic => COMP_AROMATIC,
                BondOrder::Single | BondOrder::Up | BondOrder::Down => COMP_NON_AROMATIC,
                _ => COMP_NON_AROMATIC | COMP_NON_AROMATIC_MULTIPLE,
            };
        } else if bond.order == BondOrder::Aromatic {
            non_ring_aromatic_bond = true;
        }
    }
    let mut any_flagged = false;
    for (idx, atom) in mol.atoms() {
        any_flagged |= atom.aromatic;
        if atom.aromatic
            && atom.charge == 0
            && matches!(atom.element.atomic_number(), 8 | 16 | 34 | 52)
            && mol.degree(idx) + chematic_core::implicit_hcount(mol, idx) as usize > 2
        {
            consistent = false;
        }
        let r = root[idx.0 as usize];
        if r != u32::MAX {
            atoms[r as usize] += 1;
            if atom.aromatic {
                kind[r as usize] |= COMP_FLAGGED;
            }
        }
    }
    for (_, bond) in mol.bonds() {
        if matches!(
            bond.order,
            BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
        ) {
            for end in [bond.atom1, bond.atom2] {
                let r = root[end.0 as usize];
                if r != u32::MAX {
                    kind[r as usize] |= COMP_INCIDENT_MULTIPLE;
                }
            }
        }
    }
    ReferenceInfo {
        cyclic,
        root,
        kind,
        edges,
        atoms,
        has_aromatic_bond,
        explicit: has_aromatic_bond && consistent,
        non_ring_aromatic_bond,
        any_flagged,
    }
}

#[test]
fn contiguous_component_facts_match_original() {
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
    let mut rows = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let input =
                chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            let expected = reference_info(&input);
            for _ in 0..2 {
                let actual = ComponentInfo::new(&input);
                let n = input.atom_count();
                let c = expected.kind.len();
                assert_eq!(actual.atom_count, n);
                assert_eq!(actual.facts.len(), n + 2 * c);
                assert_eq!(actual.cyclic, expected.cyclic);
                assert_eq!(&actual.facts[..n], expected.root.as_slice());
                assert_eq!(actual.kind, expected.kind);
                assert_eq!(&actual.facts[n..n + c], expected.edges.as_slice());
                assert_eq!(&actual.facts[n + c..], expected.atoms.as_slice());
                assert_eq!(
                    actual.roots().collect::<Vec<_>>(),
                    (0..c)
                        .filter(|&i| expected.edges[i] > 0)
                        .collect::<Vec<_>>()
                );
                assert_eq!(actual.has_aromatic_bond, expected.has_aromatic_bond);
                assert_eq!(actual.explicit, expected.explicit);
                assert_eq!(
                    actual.non_ring_aromatic_bond,
                    expected.non_ring_aromatic_bond
                );
                assert_eq!(actual.any_flagged, expected.any_flagged);
            }
            rows += 1;
        }
    }
    assert_eq!(rows, 14072);
}
