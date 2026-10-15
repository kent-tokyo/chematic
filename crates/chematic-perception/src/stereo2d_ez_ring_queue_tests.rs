//! Baseline breadth-first traversal with return-only visited-count instrumentation.
use super::*;

fn reference_ring_search(mol: &Molecule, bond: BondIdx, size: usize) -> (bool, usize) {
    let (start, goal) = (mol.bond(bond).atom1, mol.bond(bond).atom2);
    let mut seen = HashSet::from([start]);
    let mut frontier = vec![start];
    for _ in 0..size.saturating_sub(2) {
        let mut next = Vec::new();
        for &a in &frontier {
            for (nb, bi) in mol.neighbors(a) {
                if bi == bond {
                    continue;
                }
                if nb == goal {
                    return (true, seen.len());
                }
                if seen.insert(nb) {
                    next.push(nb);
                }
            }
        }
        frontier = next;
    }
    (false, seen.len())
}

#[test]
fn ring_queue_keeps_original_depth_and_membership_on_all_pinned_bonds() {
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
    let sizes = (0..=12).chain([33, 65]).collect::<Vec<_>>();
    let mut rows = 0usize;
    let mut queries = 0usize;
    let mut inline_cases = 0usize;
    let mut hashed_cases = 0usize;
    let mut found = 0usize;
    let mut not_found = 0usize;
    for corpus in corpora {
        for line in corpus.lines().filter(|line| !line.trim().is_empty()) {
            let mol =
                chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            for (bi, _) in mol.bonds() {
                for &size in &sizes {
                    let (expected, visited) = reference_ring_search(&mol, bi, size);
                    assert_eq!(
                        in_ring_smaller_than(&mol, bi, size),
                        expected,
                        "row {rows}, bond {bi:?}, size {size}"
                    );
                    queries += 1;
                    if visited <= 32 {
                        inline_cases += 1;
                    } else {
                        hashed_cases += 1;
                    }
                    if expected {
                        found += 1;
                    } else {
                        not_found += 1;
                    }
                }
            }
            rows += 1;
        }
    }
    assert_eq!(rows, 14072);
    assert!(queries > 1_000_000);
    assert!(inline_cases > 0 && hashed_cases > 0 && found > 0 && not_found > 0);
    eprintln!(
        "ring_queue_oracle {}",
        serde_json::json!({"rows":rows,"queries":queries,"inline_cases":inline_cases,"hashed_cases":hashed_cases,"found":found,"not_found":not_found})
    );
}
