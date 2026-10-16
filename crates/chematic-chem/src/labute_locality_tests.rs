//! Frozen source-main Labute oracle, including nonfinite output bits.
use super::*;

fn frozen_labute_asa_core(mol: &Molecule, collect_per_atom: bool) -> (Vec<f64>, f64, f64) {
    // RDKit's getLabuteAtomContribs (MolSurf.cpp) with includeHs, operation
    // for operation; bond aromaticity as RDKit perceives it.
    let n = mol.atom_count();
    if n == 0 {
        return (Vec::new(), 0.0, 0.0);
    }
    let view = crate::descriptors::descriptor_aromaticity(mol);
    let view: &Molecule = &view;
    // Most drug-like molecules fit inline.  Keeping radius and overlap in one
    // buffer avoids two heap allocations and improves locality without
    // changing RDKit's accumulation order.
    let mut atom_data: smallvec::SmallVec<[(f64, f64); 64]> =
        mol.atoms().map(|(_, a)| (rb0(a), 0.0)).collect();
    let mut accumulate_bond = |bidx| {
        let bond = view.bond(bidx);
        let i = bond.atom1.0 as usize;
        let j = bond.atom2.0 as usize;
        let (ri, rj) = (atom_data[i].0, atom_data[j].0);
        let mut bij = ri + rj;
        bij -= bond_scale(view, bidx);
        let dij = (ri - rj).abs().max(bij).min(ri + rj);
        atom_data[i].1 += rj * rj - (ri - dij) * (ri - dij) / dij;
        atom_data[j].1 += ri * ri - (rj - dij) * (rj - dij) / dij;
    };
    // Bonds in RDKit's numbering (ring closures last), without materializing
    // the complete order vector for every molecule.
    mol.for_each_bond_in_rdkit_order(&mut accumulate_bond);
    // One hydrogen per atom, pooled into a single term (as RDKit does).
    let rh = crate::descriptors::RDKIT_RB0[0];
    let mut h_contrib = 0.0f64;
    for i in 0..n {
        let ri = atom_data[i].0;
        let bij = ri + rh;
        let dij = (ri - rh).abs().max(bij).min(ri + rh);
        atom_data[i].1 += rh * rh - (ri - dij) * (ri - dij) / dij;
        h_contrib += ri * ri - (rh - dij) * (rh - dij) / dij;
    }
    let mut per_atom = if collect_per_atom {
        Vec::with_capacity(n)
    } else {
        Vec::new()
    };
    let mut total = 0.0;
    for i in 0..n {
        let (ri, overlap) = atom_data[i];
        let contribution = PI * ri * (4.0 * ri - overlap);
        total += contribution;
        if collect_per_atom {
            per_atom.push(contribution);
        }
    }
    let h_pool_area = if h_contrib.abs() > 1e-4 {
        PI * rh * (4.0 * rh - h_contrib)
    } else {
        0.0
    };
    (per_atom, h_pool_area, total + h_pool_area)
}

fn compare(input: &Molecule) {
    for collect in [false, true] {
        let reference = input.clone();
        let candidate = input.clone();
        for pass in 0..2 {
            let expected = frozen_labute_asa_core(&reference, collect);
            let actual = labute_asa_core(&candidate, collect);
            assert_eq!(
                actual.0.len(),
                expected.0.len(),
                "collect={collect}, pass={pass}"
            );
            for (index, (actual, expected)) in actual.0.iter().zip(&expected.0).enumerate() {
                assert_eq!(
                    actual.to_bits(),
                    expected.to_bits(),
                    "atom={index}, collect={collect}, pass={pass}"
                );
            }
            assert_eq!(
                actual.1.to_bits(),
                expected.1.to_bits(),
                "hydrogen pool, collect={collect}, pass={pass}"
            );
            assert_eq!(
                actual.2.to_bits(),
                expected.2.to_bits(),
                "total, collect={collect}, pass={pass}"
            );
        }
    }
}

#[test]
fn fused_labute_loop_keeps_all_pinned_per_atom_and_total_bits() {
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
    let mut inline = 0;
    let mut spilled = 0;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let input =
                chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            inline += usize::from(input.atom_count() <= 64);
            spilled += usize::from(input.atom_count() > 64);
            compare(&input);
            rows += 1;
        }
    }
    assert_eq!(rows, 14072);
    assert!(inline > 0 && spilled > 0);
}

#[test]
fn fused_labute_loop_keeps_empty_dummy_nonfinite_and_capacity_boundaries() {
    compare(&chematic_core::MoleculeBuilder::new().build());
    for size in [1, 2, 63, 64, 65, 127, 128, 129, 256] {
        compare(&chematic_smiles::parse_template(&"C".repeat(size)).unwrap());
    }
    for smiles in [
        "*",
        "**",
        "[H]",
        "[H][H]",
        "[2H]O[2H]",
        "[Fe]",
        "[Na+].[Cl-]",
        "N->[Cu+2]<-N",
        "c1ccccc1",
        "C1=CC=CC=C1",
        "C1CC1",
        "C%10CCCCC%10",
        "c1ccccc1.c1ccccc1",
    ] {
        compare(&chematic_smiles::parse_template(smiles).unwrap());
    }
}
