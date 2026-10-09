//! A hook through which an exact port of RDKit's sanitization (registered
//! by chematic-smiles, which depends on this crate) corrects the
//! RDKit-parity aromatic view ([`crate::apply_aromaticity_rdkit_parity_shared`])
//! on molecules where that view may disagree with RDKit.
//!
//! The parity view agrees with RDKit on ordinary molecules (every row of
//! the ChEMBL and rdkit-js corpora), so the hook only runs for molecules
//! [`rdkit_model_may_disagree`] flags: unusual elements or dummy atoms,
//! radicals and valence states RDKit's cleanup rewrites, charged carbon or
//! boron, aromatic atoms outside an aromatic ring system or with exocyclic
//! multiple bonds RDKit's electron counting treats specially, and aromatic
//! systems whose Hückel electron count is not 4n+2 (antiaromatic rings,
//! annulenes). Everything else keeps the parity view's cost.

use std::sync::OnceLock;

use chematic_core::{AtomIdx, BondOrder, Element, Molecule};

/// `hook(mol, view)`: `view` corrected with RDKit's model of `mol`, or
/// `None` when they agree (or the model cannot be computed).
pub type RdkitModelHook = fn(&Molecule, &Molecule) -> Option<Molecule>;

static HOOK: OnceLock<RdkitModelHook> = OnceLock::new();

/// Registers the RDKit-model hook (first registration wins).
pub fn set_rdkit_model_hook(hook: RdkitModelHook) {
    let _ = HOOK.set(hook);
}

/// The registered hook, when `mol` needs it.
pub(crate) fn hook_for(mol: &Molecule) -> Option<RdkitModelHook> {
    let hook = *HOOK.get()?;
    rdkit_model_may_disagree(mol).then_some(hook)
}

fn twice_order(order: BondOrder) -> Option<u32> {
    Some(match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => 2,
        BondOrder::Double => 4,
        BondOrder::Triple => 6,
        BondOrder::Aromatic => 3,
        _ => return None,
    })
}

/// Whether the RDKit-parity view of `mol` may differ from RDKit's
/// sanitized molecule (see the module documentation). Linear time.
pub fn rdkit_model_may_disagree(mol: &Molecule) -> bool {
    let n = mol.atom_count();
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
    if mol.bonds().any(|(_, b)| b.order == BondOrder::Triple) {
        let ring = crate::sssr::ring_bond_flags_shared(mol);
        if mol.bonds().any(|(bi, b)| {
            b.order == BondOrder::Triple && ring.get(bi.0 as usize).copied().unwrap_or(false)
        }) {
            return true;
        }
    }
    !aromatic_systems_are_huckel(mol)
}

/// Pi electrons an aromatic atom contributes to its ring system, `None`
/// where the count is not obvious.
fn pi_electrons(mol: &Molecule, idx: AtomIdx) -> Option<u32> {
    let atom = mol.atom(idx);
    let degree = mol.degree(idx);
    let hs = atom.hydrogen_count.unwrap_or(0);
    let mut exo_double = false;
    for &(nb, b) in mol.neighbor_slice(idx) {
        let order = mol.bond(b).order;
        if order == BondOrder::Double && !mol.atom(nb).aromatic {
            let nel = mol.atom(nb).element;
            if !matches!(nel, Element::O | Element::S | Element::N) || atom.element != Element::C {
                return None;
            }
            exo_double = true;
        }
    }
    Some(match (atom.element, atom.charge) {
        (Element::C, 0) => {
            if exo_double {
                0
            } else {
                1
            }
        }
        (Element::N | Element::P, 0) => {
            if degree + usize::from(hs) >= 3 {
                2
            } else {
                1
            }
        }
        (Element::N | Element::P, 1) => 1,
        (Element::N | Element::P, -1) => 2,
        (Element::O | Element::S | Element::SE, 0) => {
            if degree > 2 {
                return None;
            }
            2
        }
        (Element::O | Element::S | Element::SE, 1) => 1,
        (Element::B, 0) => 0,
        _ => return None,
    })
}

/// Whether every aromatic component (atoms joined by aromatic bonds) has a
/// Hückel count of 4n+2 pi electrons, is not a tree, and, when fused, has
/// 4n+2-electron smallest rings ([`fused_rings_are_huckel`]).
fn aromatic_systems_are_huckel(mol: &Molecule) -> bool {
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
        if electrons < 2 || (electrons - 2) % 4 != 0 {
            return false;
        }
        let bonds = twice_bonds / 2;
        if bonds < comp.len() {
            // a tree: aromatic atoms outside any aromatic ring
            return false;
        }
        if bonds > comp.len() && !fused_rings_are_huckel(mol, &comp, &electrons_of) {
            return false;
        }
    }
    true
}

/// For a fused aromatic component: the smallest aromatic ring through each
/// aromatic bond (RDKit judges fused systems ring by ring) has 4n+2 pi
/// electrons, and every aromatic bond is in such a ring.
fn fused_rings_are_huckel(mol: &Molecule, comp: &[usize], electrons_of: &[u32]) -> bool {
    let mut local = std::collections::HashMap::with_capacity(comp.len());
    for (k, &a) in comp.iter().enumerate() {
        local.insert(a, k);
    }
    let adj: Vec<Vec<usize>> = comp
        .iter()
        .map(|&a| {
            mol.neighbor_slice(AtomIdx(a as u32))
                .iter()
                .filter(|&&(_, b)| mol.bond(b).order == BondOrder::Aromatic)
                .map(|&(nb, _)| local[&(nb.0 as usize)])
                .collect()
        })
        .collect();
    let m = comp.len();
    let mut parent = vec![usize::MAX; m];
    let mut queue = std::collections::VecDeque::with_capacity(m);
    for u in 0..m {
        for &v in &adj[u] {
            if v < u {
                continue;
            }
            // BFS from u to v without the edge u-v.
            parent.iter_mut().for_each(|p| *p = usize::MAX);
            parent[u] = u;
            queue.clear();
            queue.push_back(u);
            'bfs: while let Some(x) = queue.pop_front() {
                for &y in &adj[x] {
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
            if electrons < 2 || (electrons - 2) % 4 != 0 {
                return false;
            }
        }
    }
    true
}
