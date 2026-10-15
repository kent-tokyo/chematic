//! Differential oracle for reuse of the cyclic-component partition.
//! Uses only the existing licensed 10k corpus, with its pinned manifest.
use super::*;

fn reference_typing(
    mol: &Molecule,
    typing: impl Fn(AtomIdx, &[bool]) -> (ElectronDonorType, bool),
) -> ParityTyping {
    let atom_count = mol.atom_count();
    let cyclic = crate::sssr::ring_bond_flags_shared(mol);
    let mut ring_atom = vec![false; atom_count];
    let mut components = UnionFind::new(atom_count);
    for (index, bond) in mol.bonds() {
        if cyclic[index.0 as usize] {
            ring_atom[bond.atom1.0 as usize] = true;
            ring_atom[bond.atom2.0 as usize] = true;
            components.union(bond.atom1.0, bond.atom2.0);
        }
    }
    let mut donor = vec![ElectronDonorType::None; atom_count];
    let mut candidate = vec![false; atom_count];
    for atom in 0..atom_count {
        if ring_atom[atom] {
            (donor[atom], candidate[atom]) = typing(AtomIdx(atom as u32), cyclic.as_slice());
        }
    }
    let keep = reference_component_mask(mol, cyclic.as_slice(), &ring_atom, &candidate, components);
    ParityTyping {
        donor,
        candidate,
        keep,
    }
}

fn reference_component_mask(
    mol: &Molecule,
    cyclic: &[bool],
    ring_atom: &[bool],
    candidate: &[bool],
    mut components: UnionFind,
) -> Vec<bool> {
    let mut candidate_components = UnionFind::new(mol.atom_count());
    let mut keep_root = vec![false; mol.atom_count()];
    for (index, bond) in mol.bonds() {
        let (left, right) = (bond.atom1.0 as usize, bond.atom2.0 as usize);
        if cyclic[index.0 as usize]
            && candidate[left]
            && candidate[right]
            && !candidate_components.union(left as u32, right as u32)
        {
            keep_root[components.find(left as u32) as usize] = true;
        }
    }
    (0..mol.atom_count())
        .map(|atom| ring_atom[atom] && keep_root[components.find(atom as u32) as usize])
        .collect()
}

fn compare(input: &Molecule) {
    let type_atom = |idx, cyclic: &[bool]| {
        let donor = get_atom_electron_donor_type(input, idx, cyclic);
        (donor, is_atom_candidate_for_aromaticity(input, idx, donor))
    };
    let expected = reference_typing(input, type_atom);
    for _ in 0..2 {
        let actual = ParityTyping::new(input, type_atom);
        assert_eq!(actual.donor, expected.donor);
        assert_eq!(actual.candidate, expected.candidate);
        assert_eq!(actual.keep, expected.keep);
    }
}

#[test]
fn fixed_ten_thousand_typings_match_original_partition() {
    let corpus = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
    ));
    let mut rows = 0;
    for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
        let input =
            chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
        compare(&input);
        rows += 1;
    }
    assert_eq!(rows, 10000);
}

#[test]
fn individually_seeded_ring_data_preserves_original_partition() {
    let corpus = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi"
    ));
    for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
        let input =
            chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
        let flags = crate::sssr::ring_bond_flags_shared(&input);
        let components = crate::sssr::ring_components_shared(&input);
        let flags_only = input.clone();
        flags_only.seed_derived(chematic_core::DerivedSlot::RingBondFlags, flags);
        compare(&flags_only);
        let components_only = input.clone();
        components_only.seed_derived(chematic_core::DerivedSlot::RingComponents, components);
        compare(&components_only);
    }
}
