//! Frozen complete model-gate functions from main 3e7f04b2.
use super::*;

fn frozen_rdkit_model_may_disagree_uncached(mol: &Molecule) -> bool {
    let n = mol.atom_count();
    let mut has_aromatic_atom = false;
    for i in 0..n {
        let idx = AtomIdx(i as u32);
        let atom = mol.atom(idx);
        if atom.wildcard {
            return true;
        }
        let el = atom.element;
        let charge = i32::from(atom.charge);
        // Twice the bond-order sum (aromatic bonds count 1.5).
        let mut twice = 0u32;
        let mut aromatic_bonds = 0;
        for &(nb, b) in mol.neighbor_slice(idx) {
            let order = mol.bond(b).order;
            let Some(t) = twice_order(order) else {
                return true;
            };
            twice += t;
            if order == BondOrder::Aromatic {
                aromatic_bonds += 1;
            }
            if atom.aromatic
                && mol.atom(nb).aromatic
                && matches!(order, BondOrder::Double | BondOrder::Triple)
            {
                return true;
            }
        }
        let bracket = atom.hydrogen_count.is_some();
        let hs = u32::from(atom.hydrogen_count.unwrap_or(0));
        if atom.aromatic {
            has_aromatic_atom = true;
            if aromatic_bonds < 2 {
                return true;
            }
            if bracket {
                let conn = mol.degree(idx) + hs as usize;
                let ok = match (el, charge) {
                    (Element::C, 0) | (Element::B, 0) => conn == 3,
                    (Element::N | Element::P, 0) => conn == 2 || conn == 3,
                    (Element::N | Element::P, 1) => conn == 3 || conn == 2,
                    (Element::N | Element::P, -1) => conn == 2,
                    (Element::O | Element::S | Element::SE, 0 | 1) => conn == 2,
                    _ => false,
                };
                if !ok {
                    return true;
                }
            }
            if !matches!(
                el,
                Element::C
                    | Element::N
                    | Element::O
                    | Element::S
                    | Element::P
                    | Element::SE
                    | Element::B
            ) {
                return true;
            }
            continue;
        }
        if aromatic_bonds > 0 {
            return true;
        }
        let val = (twice / 2 + hs) as i32;
        match el {
            Element::C => {
                if charge != 0 || (bracket && val < 4) || val > 4 {
                    return true;
                }
            }
            Element::B => {
                if charge != 0 || (bracket && val < 3) || val > 3 {
                    return true;
                }
            }
            Element::N | Element::P => {
                let base = 3 + charge;
                if (bracket && val < base) || (el == Element::N && val > base) {
                    return true;
                }
            }
            Element::O | Element::S | Element::SE => {
                let base = 2 + charge;
                if (bracket && val < base) || (el == Element::O && val > base) {
                    return true;
                }
            }
            Element::F | Element::CL | Element::BR | Element::I => {
                if (bracket && val < 1 + charge) || val > 1 + charge.max(0) {
                    return true;
                }
            }
            Element::H | Element::SI => {
                if charge != 0 {
                    return true;
                }
            }
            _ => return true,
        }
    }
    // A ring triple bond (`C1=CSC#C1`): RDKit may flag it aromatic while
    // keeping its type, which the parity view does not model.
    if has_cyclic_triple_bond(mol) {
        return true;
    }
    // With no aromatic atoms there is no aromatic component to check.
    // The atom/valence and cyclic-triple checks above still run unchanged.
    has_aromatic_atom && !frozen_aromatic_systems_are_huckel(mol)
}

fn frozen_aromatic_systems_are_huckel(mol: &Molecule) -> bool {
    let n = mol.atom_count();
    let mut seen = vec![false; n];
    let mut electrons_of = vec![0u32; n];
    let mut stack = Vec::new();
    for start in 0..n {
        if seen[start] || !mol.atom(AtomIdx(start as u32)).aromatic {
            continue;
        }
        seen[start] = true;
        stack.push(start);
        let mut comp = Vec::new();
        let mut electrons = 0u32;
        let mut twice_bonds = 0usize;
        while let Some(a) = stack.pop() {
            let idx = AtomIdx(a as u32);
            match pi_electrons(mol, idx) {
                Some(e) => {
                    electrons += e;
                    electrons_of[a] = e;
                }
                None => return false,
            }
            comp.push(a);
            for &(nb, b) in mol.neighbor_slice(idx) {
                if mol.bond(b).order != BondOrder::Aromatic {
                    continue;
                }
                twice_bonds += 1;
                let j = nb.0 as usize;
                if !seen[j] {
                    seen[j] = true;
                    stack.push(j);
                }
            }
        }
        if electrons < 2 || !(electrons - 2).is_multiple_of(4) {
            return false;
        }
        let bonds = twice_bonds / 2;
        if bonds < comp.len() {
            // a tree: aromatic atoms outside any aromatic ring
            return false;
        }
        if bonds > comp.len() && !frozen_fused_rings_are_huckel(mol, &comp, &electrons_of) {
            return false;
        }
    }
    true
}

fn frozen_fused_rings_are_huckel(mol: &Molecule, comp: &[usize], electrons_of: &[u32]) -> bool {
    let mut local = vec![u32::MAX; mol.atom_count()];
    for (k, &a) in comp.iter().enumerate() {
        local[a] = k as u32;
    }
    // Flat aromatic adjacency of the component in local numbering.
    let mut start = Vec::with_capacity(comp.len() + 1);
    let mut flat: Vec<usize> = Vec::with_capacity(comp.len() * 3);
    for &a in comp {
        start.push(flat.len());
        flat.extend(
            mol.neighbor_slice(AtomIdx(a as u32))
                .iter()
                .filter(|&&(_, b)| mol.bond(b).order == BondOrder::Aromatic)
                .map(|&(nb, _)| local[nb.0 as usize] as usize),
        );
    }
    start.push(flat.len());
    let adj = |x: usize| &flat[start[x]..start[x + 1]];
    let m = comp.len();
    let mut parent = vec![usize::MAX; m];
    let mut queue = std::collections::VecDeque::with_capacity(m);
    for u in 0..m {
        for &v in adj(u) {
            if v < u {
                continue;
            }
            // BFS from u to v without the edge u-v.
            parent.iter_mut().for_each(|p| *p = usize::MAX);
            parent[u] = u;
            queue.clear();
            queue.push_back(u);
            'bfs: while let Some(x) = queue.pop_front() {
                for &y in adj(x) {
                    if (x == u && y == v) || parent[y] != usize::MAX {
                        continue;
                    }
                    parent[y] = x;
                    if y == v {
                        break 'bfs;
                    }
                    queue.push_back(y);
                }
            }
            if parent[v] == usize::MAX {
                return false;
            }
            let mut electrons = 0u32;
            let mut x = v;
            loop {
                electrons += electrons_of[comp[x]];
                if x == u {
                    break;
                }
                x = parent[x];
            }
            if electrons < 2 || !(electrons - 2).is_multiple_of(4) {
                return false;
            }
        }
    }
    true
}

fn compare(mol: &Molecule) {
    let baseline = mol.clone();
    let candidate = mol.clone();
    for _ in 0..2 {
        assert_eq!(
            rdkit_model_may_disagree_uncached(&candidate),
            frozen_rdkit_model_may_disagree_uncached(&baseline)
        );
        assert_eq!(
            aromatic_systems_are_huckel(&candidate),
            frozen_aromatic_systems_are_huckel(&baseline)
        );
    }
}

#[test]
fn all_pinned_model_gate_outcomes_are_unchanged() {
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
            let mol =
                chematic_smiles::parse_template(line.split_whitespace().next().unwrap()).unwrap();
            inline += usize::from(mol.atom_count() <= 64);
            spilled += usize::from(mol.atom_count() > 64);
            compare(&mol);
            rows += 1;
        }
    }
    assert_eq!(rows, 14072);
    assert!(inline > 0 && spilled > 0);
}

#[test]
fn unusual_and_disconnected_graphs_keep_gate_decisions() {
    compare(&chematic_core::MoleculeBuilder::new().build());
    for smiles in [
        "*",
        "[H]",
        "[C]",
        "[nH]1cccc1",
        "c1cccc1",
        "c1ccccccc1",
        "c1ccccc#1",
        "C1=CSC#C1",
        "c1ccccc1.c1ccccc1",
        "[Na+].[Cl-]",
        "[O-][Cl](=O)(=O)=O",
        "N->[Cu]<-N",
        "[13CH3][NH3+]",
        "o1cccccccc1",
        "s1cccccccc1",
        "O=C1NNC(=O)N1",
    ] {
        compare(&chematic_smiles::parse_template(smiles).unwrap());
    }
}

fn fused_hexagons(rings: usize, spectators: usize) -> Molecule {
    use chematic_core::{Atom, MoleculeBuilder};
    let mut b = MoleculeBuilder::new();
    let add = |b: &mut MoleculeBuilder| {
        let mut atom = Atom::new(Element::C);
        atom.aromatic = true;
        b.add_atom(atom)
    };
    let first: Vec<_> = (0..6).map(|_| add(&mut b)).collect();
    for i in 0..6 {
        b.add_bond(first[i], first[(i + 1) % 6], BondOrder::Aromatic)
            .unwrap();
    }
    let mut edge = (first[0], first[1]);
    for _ in 1..rings {
        let fresh: Vec<_> = (0..4).map(|_| add(&mut b)).collect();
        let path = [edge.0, fresh[0], fresh[1], fresh[2], fresh[3], edge.1];
        for pair in path.windows(2) {
            b.add_bond(pair[0], pair[1], BondOrder::Aromatic).unwrap();
        }
        edge = (fresh[1], fresh[2]);
    }
    for _ in 0..spectators {
        b.add_atom(Atom::new(Element::C));
    }
    b.build()
}

#[test]
fn fused_bfs_keeps_paths_across_inline_and_spill_boundaries() {
    for rings in [2, 8, 15, 16, 17, 31, 32] {
        for spectators in [0, 64] {
            let mol = fused_hexagons(rings, spectators);
            let n = 6 + 4 * (rings - 1);
            let comp: Vec<_> = (0..n).collect();
            for style in 0..3 {
                let electrons: Vec<_> = (0..mol.atom_count())
                    .map(|i| match style {
                        0 => 1,
                        1 => (i % 3) as u32,
                        _ => 2,
                    })
                    .collect();
                assert_eq!(
                    fused_rings_are_huckel(&mol, &comp, &electrons),
                    frozen_fused_rings_are_huckel(&mol, &comp, &electrons),
                    "rings={rings}, spectators={spectators}, style={style}"
                );
            }
            compare(&mol);
        }
    }
}
