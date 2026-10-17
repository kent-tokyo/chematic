//! Frozen atom-invariant construction; share only unchanged old helpers.
use super::{connectivity_invariant, hash_vec, rdkit_isotope_delta};
use crate::ecfp::{rdkit_total_degree, rdkit_total_h_count};
use chematic_core::{AtomIdx, Molecule};
use smallvec::SmallVec;

fn reference_connectivity_invariant(mol: &Molecule, idx: AtomIdx, ring_atoms: &[bool]) -> u32 {
    let atom = mol.atom(idx);
    let mut components: SmallVec<[u32; 6]> = SmallVec::new();
    // Dummy atoms have atomic number 0 (chematic stores them as C).
    components.push(if atom.wildcard {
        0
    } else {
        atom.element.atomic_number() as u32
    });
    components.push(rdkit_total_degree(mol, idx) as u32);
    components.push(rdkit_total_h_count(mol, idx) as u32);
    components.push((atom.charge as i32) as u32);
    components.push(rdkit_isotope_delta(mol, idx) as u32);
    if ring_atoms.get(idx.0 as usize).copied().unwrap_or(false) {
        components.push(1);
    }
    hash_vec(&components)
}

fn compare(mol: &Molecule, label: &str) {
    let rings = chematic_perception::ring_atom_flags(mol);
    let alternating: Vec<bool> = (0..mol.atom_count()).map(|i| i % 2 == 0).collect();
    for flags in [&rings[..], &[][..], &alternating[..]] {
        for (idx, _) in mol.atoms() {
            assert_eq!(
                connectivity_invariant(mol, idx, flags),
                reference_connectivity_invariant(mol, idx, flags),
                "{label} atom {}",
                idx.0
            );
        }
    }
}

#[test]
fn connectivity_components_match_frozen_all_corpus_atoms() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut checked = 0;
    for path in [
        "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi",
        "validation/cosmolkit_comparison/corpora/stress-v1.smi",
    ] {
        for line in std::fs::read_to_string(root.join(path))
            .unwrap()
            .lines()
            .filter(|s| !s.trim().is_empty())
        {
            let smiles = line.split_whitespace().next().unwrap();
            let mol = chematic_smiles::parse(smiles).unwrap();
            compare(&mol, smiles);
            if let Ok(view) = chematic_perception::apply_aromaticity_rdkit_parity_experimental(&mol)
            {
                compare(&view, smiles);
            }
            checked += 1;
        }
    }
    assert_eq!(checked, 14072);
}

#[test]
fn bracket_wildcard_isotopes_charges_and_h_graphs_match_frozen() {
    for smiles in [
        "*",
        "[*H]",
        "[*H2]",
        "[H]",
        "[2H]C([H])([H])[H]",
        "[H]C([H])([H])[H]",
        "[13CH4]",
        "[NH4+]",
        "[O-]",
        "[Na+].Cl",
        "c1cc[nH]c1",
        "[C:7]([H])([H])([H])[H]",
        "[C@H](F)(Cl)Br",
    ] {
        let mol = chematic_smiles::parse(smiles).unwrap();
        compare(&mol, smiles);
    }
    // The missing ring-flag entry is false, including an isolated atom.
    let mol = chematic_smiles::parse("CC").unwrap();
    assert_eq!(
        connectivity_invariant(&mol, AtomIdx(1), &[true]),
        reference_connectivity_invariant(&mol, AtomIdx(1), &[true])
    );
}

#[test]
fn numeric_h_count_casts_and_saturation_match_frozen_builder_graphs() {
    use chematic_core::{Atom, BondOrder, Element, MoleculeBuilder};
    for wildcard in [false, true] {
        for count in [0, 1, 254, 255] {
            for explicit in [0, 1, 255, 256, 260] {
                let mut builder = MoleculeBuilder::new();
                let mut center = Atom::new(Element::C);
                center.wildcard = wildcard;
                center.hydrogen_count = Some(count);
                let center = builder.add_atom(center);
                for i in 0..explicit {
                    let mut h = Atom::new(Element::H);
                    h.wildcard = i == 1;
                    let h = builder.add_atom(h);
                    builder.add_bond(center, h, BondOrder::Single).unwrap();
                }
                compare(&builder.build(), "builder H cast/saturation");
            }
        }
    }
}
