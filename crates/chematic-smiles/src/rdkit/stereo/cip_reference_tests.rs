//! Frozen ranking code from a07a1294 (identical CIP bodies to c84c68f2); only the function and pool names change.

use super::super::{
    mol::{BondType, Mol},
    periodic,
};
use super::iterate_cip_ranks;

/// `buildCIPInvariants`.
fn build_cip_invariants(mol: &Mol) -> Vec<i64> {
    const N_MASS_BITS: i64 = 10;
    const MAX_MASS: i64 = 1 << N_MASS_BITS;
    mol.atoms
        .iter()
        .map(|atom| {
            let num = i64::from(atom.anum % 128);
            let mut mass: i64 = 0;
            if atom.isotope != 0 {
                mass =
                    i64::from(atom.isotope) - i64::from(periodic::most_common_isotope(atom.anum));
                if mass >= 0 {
                    mass += 1;
                }
            }
            mass += MAX_MASS / 2;
            if mass < 0 {
                mass = 0;
            } else {
                mass %= MAX_MASS;
            }
            let mut invariant = num;
            invariant = (invariant << N_MASS_BITS) | mass;
            let mapnum = atom.map.map_or(-1i64, i64::from);
            let mapnum = (mapnum + 1) % 1024;
            invariant = (invariant << 10) | mapnum;
            invariant
        })
        .collect()
}

/// `findSegmentsToResort`: assigns ranks along `sorted`, returns tied
/// sections and the number of distinct ranks.
fn find_segments_to_resort(
    sorted: &[usize],
    cip: &[Vec<i32>],
    curr_rank: &mut [u32],
) -> (Vec<(usize, usize)>, usize) {
    let mut res: Vec<(usize, usize)> = Vec::new();
    let mut num_independent = sorted.len();
    let mut current = 0usize;
    let mut running = 0u32;
    curr_rank[sorted[0]] = 0;
    let mut in_equal = false;
    for i in 1..sorted.len() {
        if cip[sorted[current]] == cip[sorted[i]] {
            curr_rank[sorted[i]] = running;
            num_independent -= 1;
            if !in_equal {
                in_equal = true;
                res.push((i - 1, 0));
            }
        } else {
            running += 1;
            curr_rank[sorted[i]] = running;
            current = i;
            if in_equal {
                res.last_mut().expect("open section").1 = i;
                in_equal = false;
            }
        }
    }
    if in_equal {
        res.last_mut().expect("open section").1 = sorted.len() - 1;
    }
    (res, num_independent)
}

/// `iterateCIPRanks`.
fn reference_iterate_cip_ranks(
    mol: &Mol,
    invars: &[i64],
    ranks: &mut Vec<u32>,
    seed_with_invars: bool,
) {
    let n = mol.atoms.len();
    ranks.resize(n, 0);
    if n == 0 {
        return;
    }
    // Entry buffers are reused across calls (RDKit reserves 16 per entry).
    let mut pool: Vec<Vec<i32>> =
        REFERENCE_ENTRY_POOL.with(|pool| std::mem::take(&mut *pool.borrow_mut()));
    while pool.len() < n {
        pool.push(Vec::with_capacity(16));
    }
    let cip = &mut pool[..n];
    for (e, &v) in cip.iter_mut().zip(invars) {
        e.clear();
        e.push(v as i32);
    }
    let mut sorted: Vec<usize> = (0..n).collect();
    sorted.sort_by(|&x, &y| cip[x].cmp(&cip[y]));
    let mut curr_rank = vec![0u32; n];
    let (mut needs_sorting, mut num_ranks) = find_segments_to_resort(&sorted, cip, &mut curr_rank);
    ranks.copy_from_slice(&curr_rank);
    for i in 0..n {
        if seed_with_invars {
            cip[i][0] = invars[i] as i32;
        } else {
            cip[i][0] = mol.atoms[i].anum as i32;
            cip[i].push(ranks[i] as i32);
        }
    }
    let cip_rank_index = if seed_with_invars { 1 } else { 2 };
    let max_its = n / 2 + 1;
    let mut num_its = 0;
    let mut last_num_ranks: i64 = -1;
    // Bond features: per atom, (count, neighbour) in bond order, stored
    // back to back (atom `a`'s are `features[feature_start[a]..feature_start[a + 1]]`).
    let mut feature_start: Vec<usize> = Vec::with_capacity(n + 1);
    let mut features: Vec<(u32, usize)> = Vec::with_capacity(2 * mol.bonds.len());
    feature_start.push(0);
    for a in 0..n {
        for &b in &mol.atom_bonds[a] {
            let bond = &mol.bonds[b];
            let nbr = bond.other(a);
            let special = bond.bt == BondType::Double
                && mol.atoms[nbr].anum == 15
                && matches!(mol.degree(nbr), 3 | 4);
            let count = if special { 1 } else { bond.bt.twice() };
            features.push((count, nbr));
        }
        feature_start.push(features.len());
    }
    let total_hs: Vec<usize> = (0..n).map(|a| mol.total_num_hs(a) as usize).collect();
    let mut vals: Vec<i32> = Vec::with_capacity(16);
    while !needs_sorting.is_empty()
        && num_its < max_its
        && (last_num_ranks < 0 || (last_num_ranks as usize) < num_ranks)
    {
        for index in 0..n {
            vals.clear();
            for &(count, nbr) in &features[feature_start[index]..feature_start[index + 1]] {
                for _ in 0..count {
                    vals.push(ranks[nbr] as i32 + 1);
                }
            }
            vals.sort_unstable_by(|x, y| y.cmp(x));
            let entry = &mut cip[index];
            entry.reserve(vals.len() + total_hs[index]);
            entry.extend_from_slice(&vals);
            entry.resize(entry.len() + total_hs[index], 0);
        }
        last_num_ranks = num_ranks as i64;
        for &(first, last) in &needs_sorting {
            sorted[first..=last].sort_by(|&x, &y| cip[x].cmp(&cip[y]));
        }
        let (ns, nr) = find_segments_to_resort(&sorted, cip, &mut curr_rank);
        needs_sorting = ns;
        num_ranks = nr;
        ranks.copy_from_slice(&curr_rank);
        if last_num_ranks as usize != num_ranks {
            for i in 0..n {
                cip[i].resize(cip_rank_index + 1, 0);
                cip[i][cip_rank_index] = ranks[i] as i32;
            }
        }
        num_its += 1;
    }
    REFERENCE_ENTRY_POOL.with(|p| *p.borrow_mut() = pool);
}

thread_local! {
 static REFERENCE_ENTRY_POOL: std::cell::RefCell<Vec<Vec<i32>>> = const { std::cell::RefCell::new(Vec::new()) };
}

fn compare_modes(mol: &Mol, label: &str) {
    let natural = build_cip_invariants(mol);
    let tied = vec![0; mol.atoms.len()];
    let unique: Vec<i64> = (0..mol.atoms.len()).map(|i| (i as i64) << 5).collect();
    // Exercise the i64 -> i32 truncation in both old and current seeding.
    let mixed: Vec<i64> = (0..mol.atoms.len())
        .map(|i| ((i.wrapping_mul(7919) % 7) as i64 - 3) * 0x1_0000_0001)
        .collect();
    for invars in [&natural, &tied, &unique, &mixed] {
        for seed in [false, true] {
            let mut reference = vec![99; mol.atoms.len() + 3];
            let mut actual = reference.clone();
            reference_iterate_cip_ranks(mol, invars, &mut reference, seed);
            iterate_cip_ranks(mol, invars, &mut actual, seed);
            assert_eq!(actual, reference, "{label}, seed={seed}, invars={invars:?}");
        }
    }
}

#[test]
fn frozen_cip_all_corpus_rows_and_seed_modes() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut checked = 0;
    let mut unsupported = 0;
    for path in [
        "validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi",
        "validation/cosmolkit_comparison/corpora/stress-v1.smi",
    ] {
        let corpus = std::fs::read_to_string(root.join(path)).expect("checked-in corpus");
        for (index, line) in corpus
            .lines()
            .filter(|line| !line.trim().is_empty())
            .enumerate()
        {
            let smiles = line.split_whitespace().next().unwrap();
            let graph = crate::parse(smiles).expect("corpus parses");
            let Ok(mut mol) = super::super::parse::from_chematic(&graph) else {
                unsupported += 1;
                continue;
            };
            let label = format!("{path}:{index}:{smiles}");
            compare_modes(&mol, &label);
            if super::super::sanitize::remove_hs_and_sanitize(&mut mol).is_ok() {
                compare_modes(&mol, &label);
            }
            checked += 1;
        }
    }
    assert_eq!(checked + unsupported, 14072);
    assert!(
        checked >= 14000,
        "checked={checked}, unsupported={unsupported}"
    );
    eprintln!("frozen CIP: {checked} rows compared, {unsupported} unsupported model conversions");
}

#[test]
fn refinement_survives_empty_shrinking_and_thread_local_calls() {
    let cases = [
        "C[C@]12CC[C@H]3[C@@H](CC=C4C[C@@H](O)CC[C@@]43C)[C@@H]1CC[C@@H]2O",
        "[Na+]",
        "",
        "C",
        "[13CH3][C@H](F)Cl",
        "F/C=C/F",
        "O=P(O)(O)O",
        "C1CCCCC1",
        "[C:7]([H])([H])([H])[H]",
        "",
    ];
    let workers: Vec<_> = (0..4)
        .map(|_| {
            std::thread::spawn(move || {
                for _ in 0..5 {
                    for smiles in cases {
                        let mol = if smiles.is_empty() {
                            Mol::default()
                        } else {
                            let graph = crate::parse(smiles).expect("boundary parses");
                            super::super::parse::from_chematic(&graph).unwrap()
                        };
                        compare_modes(&mol, smiles);
                    }
                }
            })
        })
        .collect();
    for worker in workers {
        worker.join().unwrap();
    }
}

#[test]
fn adjacent_tied_sections_share_boundaries_without_double_refinement() {
    let keys: Vec<Vec<i32>> = [0, 0, 1, 1, 2, 2, 3, 4]
        .into_iter()
        .map(|k| vec![k])
        .collect();
    let order: Vec<_> = (0..8).collect();
    let mut ranks = vec![0; 8];
    assert_eq!(
        find_segments_to_resort(&order, &keys, &mut ranks).0,
        vec![(0, 2), (2, 4), (4, 6)]
    );
    let graph = crate::parse("CCCCCCCC").unwrap();
    let mol = super::super::parse::from_chematic(&graph).unwrap();
    for seed in [false, true] {
        let mut old = Vec::new();
        let mut actual = Vec::new();
        reference_iterate_cip_ranks(&mol, &[0, 0, 1, 1, 2, 2, 3, 4], &mut old, seed);
        iterate_cip_ranks(&mol, &[0, 0, 1, 1, 2, 2, 3, 4], &mut actual, seed);
        assert_eq!(actual, old);
    }
}

#[test]
fn generated_bond_orders_maps_isotopes_and_graph_sizes_match_frozen_ranks() {
    use super::super::mol::{Atom, Bond};
    for n in 0_usize..65 {
        for seed in 0_usize..4 {
            let mut mol = Mol::default();
            for a in 0..n {
                let mut atom = Atom::new(if a.is_multiple_of(9) { 15 } else { 6 });
                atom.isotope = if a.is_multiple_of(5) { 13 } else { 0 };
                atom.map = (a.is_multiple_of(3)).then_some((a % 7) as u32);
                atom.num_explicit_hs = (a % 3) as u32;
                atom.implicit_valence = (a % 4) as i32;
                atom.no_implicit = a.is_multiple_of(2);
                mol.add_atom(atom);
            }
            let types = [
                BondType::Single,
                BondType::Double,
                BondType::Triple,
                BondType::Aromatic,
                BondType::Dative,
                BondType::Quadruple,
            ];
            for a in 0..n {
                for b in a..n {
                    if b != a + 1 && !(a * 17 + b * 7 + seed).is_multiple_of(19) {
                        continue;
                    }
                    mol.add_bond(Bond::new(a, b, types[(a + b + seed) % types.len()]));
                }
            }
            compare_modes(&mol, &format!("generated n={n},seed={seed}"));
        }
    }
}
