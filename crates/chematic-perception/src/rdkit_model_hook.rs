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
    // Memoized: several consumers ask for one molecule.
    *mol.derived(chematic_core::DerivedSlot::RdkitModelGate, || {
        rdkit_model_may_disagree_uncached(mol)
    })
}

fn rdkit_model_may_disagree_uncached(mol: &Molecule) -> bool {
    let mut has_aromatic_atom = false;
    for raw in 0..mol.atom_count() {
        let (aromatic, may_disagree) = atom_may_disagree(mol, AtomIdx(raw as u32));
        if may_disagree {
            return true;
        }
        has_aromatic_atom |= aromatic;
    }
    // A ring triple bond (`C1=CSC#C1`): RDKit may flag it aromatic while
    // keeping its type, which the parity view does not model.
    if has_cyclic_triple_bond(mol) {
        return true;
    }
    // With no aromatic atoms there is no aromatic component to check.
    // The atom/valence and cyclic-triple checks above still run unchanged.
    has_aromatic_atom && !aromatic_systems_are_huckel(mol)
}

fn atom_may_disagree(mol: &Molecule, idx: AtomIdx) -> (bool, bool) {
    let atom = mol.atom(idx);
    if atom.wildcard {
        return (atom.aromatic, true);
    }
    let mut twice_bond_order = 0u32;
    let mut aromatic_bonds = 0;
    for &(neighbor, bond) in mol.neighbor_slice(idx) {
        let order = mol.bond(bond).order;
        let Some(twice) = twice_order(order) else {
            return (atom.aromatic, true);
        };
        twice_bond_order += twice;
        aromatic_bonds += usize::from(order == BondOrder::Aromatic);
        if atom.aromatic
            && mol.atom(neighbor).aromatic
            && matches!(order, BondOrder::Double | BondOrder::Triple)
        {
            return (true, true);
        }
    }
    let explicit_hydrogens = atom.hydrogen_count.is_some();
    let hydrogens = u32::from(atom.hydrogen_count.unwrap_or(0));
    let disagrees = if atom.aromatic {
        aromatic_atom_may_disagree(
            atom.element,
            i32::from(atom.charge),
            mol.degree(idx) + hydrogens as usize,
            explicit_hydrogens,
            aromatic_bonds,
        )
    } else {
        aliphatic_atom_may_disagree(
            atom.element,
            i32::from(atom.charge),
            (twice_bond_order / 2 + hydrogens) as i32,
            explicit_hydrogens,
            aromatic_bonds,
        )
    };
    (atom.aromatic, disagrees)
}

fn aromatic_atom_may_disagree(
    element: Element,
    charge: i32,
    connectivity: usize,
    bracketed: bool,
    aromatic_bonds: usize,
) -> bool {
    if aromatic_bonds < 2 {
        return true;
    }
    let supported = matches!(
        element,
        Element::C | Element::N | Element::O | Element::S | Element::P | Element::SE | Element::B
    );
    let valid_bracket_valence = !bracketed
        || match (element, charge) {
            (Element::C, 0) | (Element::B, 0) => connectivity == 3,
            (Element::N | Element::P, 0 | 1) => connectivity == 2 || connectivity == 3,
            (Element::N | Element::P, -1) => connectivity == 2,
            (Element::O | Element::S | Element::SE, 0 | 1) => connectivity == 2,
            _ => false,
        };
    !supported || !valid_bracket_valence
}

fn aliphatic_atom_may_disagree(
    element: Element,
    charge: i32,
    valence: i32,
    bracketed: bool,
    aromatic_bonds: usize,
) -> bool {
    if aromatic_bonds > 0 {
        return true;
    }
    match element {
        Element::C => charge != 0 || (bracketed && valence < 4) || valence > 4,
        Element::B => charge != 0 || (bracketed && valence < 3) || valence > 3,
        Element::N | Element::P => {
            let base = 3 + charge;
            (bracketed && valence < base) || (element == Element::N && valence > base)
        }
        Element::O | Element::S | Element::SE => {
            let base = 2 + charge;
            (bracketed && valence < base) || (element == Element::O && valence > base)
        }
        Element::F | Element::CL | Element::BR | Element::I => {
            (bracketed && valence < 1 + charge) || valence > 1 + charge.max(0)
        }
        Element::H | Element::SI => charge != 0,
        _ => true,
    }
}

fn has_cyclic_triple_bond(mol: &Molecule) -> bool {
    if !mol.bonds().any(|(_, bond)| bond.order == BondOrder::Triple) {
        return false;
    }
    let ring = crate::sssr::ring_bond_flags_shared(mol);
    mol.bonds().any(|(idx, bond)| {
        bond.order == BondOrder::Triple && ring.get(idx.0 as usize).copied().unwrap_or(false)
    })
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
        if electrons < 2 || !(electrons - 2).is_multiple_of(4) {
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
