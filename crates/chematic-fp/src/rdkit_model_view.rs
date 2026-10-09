//! The molecule RDKit fingerprints: chematic's memoized RDKit-parity view
//! ([`chematic_perception::with_rdkit_parity_view`]) corrected, where it
//! disagrees, with the aromaticity, bond orders, charges and hydrogen counts
//! of chematic-smiles' port of RDKit's `MolFromSmiles` sanitization (the
//! model the RDKit-compatible SMILES writer uses). The parity view agrees
//! with the port on ordinary molecules, so the common path only pays for
//! the port's sanitization; inputs such as `c1ccc1`, which RDKit
//! dearomatizes, get RDKit's graph.

use chematic_core::{AtomIdx, BondIdx, BondOrder, Molecule};
use chematic_perception::AromaticityError;

/// Runs `f` on the RDKit-model view of `mol` (see the module docs); `f`
/// receives the parity view's error where chematic cannot perceive `mol`.
pub(crate) fn with_rdkit_model_view<R>(
    mol: &Molecule,
    f: impl FnOnce(Result<&Molecule, &AromaticityError>) -> R,
) -> R {
    chematic_perception::with_rdkit_parity_view(mol, |view| match view {
        Ok(v) => match patched(mol, v) {
            Some(p) => f(Ok(&p)),
            None => f(Ok(v)),
        },
        Err(e) => f(Err(e)),
    })
}

/// The bond-order class RDKit's model distinguishes.
fn class(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single | BondOrder::Up | BondOrder::Down => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 12,
        BondOrder::Dative => 17,
        _ => 0,
    }
}

/// `view` with the port's model where they disagree; `None` when they
/// agree (or the port cannot model `mol`, or the view is not on `mol`'s
/// graph).
pub(crate) fn patched(mol: &Molecule, view: &Molecule) -> Option<Molecule> {
    if view.atom_count() != mol.atom_count() || view.bond_count() != mol.bond_count() {
        return None;
    }
    let model = chematic_smiles::rdkit_sanitized_model(mol).ok()?;
    let atom_differs = |i: usize| -> bool {
        match model.atoms[i] {
            Some((arom, charge, _)) => {
                let a = view.atom(AtomIdx(i as u32));
                a.aromatic != arom || a.charge != charge
            }
            None => false,
        }
    };
    let bond_differs = |b: usize| -> bool {
        match model.bonds[b] {
            Some(order) => class(view.bond(BondIdx(b as u32)).order) != class(order),
            None => false,
        }
    };
    if !(0..mol.atom_count()).any(atom_differs) && !(0..mol.bond_count()).any(bond_differs) {
        return None;
    }
    let mut out = view.clone();
    let mut touched = vec![false; mol.atom_count()];
    for b in 0..mol.bond_count() {
        if bond_differs(b) {
            let order = model.bonds[b].expect("modelled bond");
            let bond = view.bond(BondIdx(b as u32));
            touched[bond.atom1.0 as usize] = true;
            touched[bond.atom2.0 as usize] = true;
            out.set_bond_order(BondIdx(b as u32), order);
        }
    }
    for i in 0..mol.atom_count() {
        let Some((arom, charge, hs)) = model.atoms[i] else {
            continue;
        };
        if atom_differs(i) || touched[i] {
            let idx = AtomIdx(i as u32);
            out.set_atom_aromatic(idx, arom);
            out.set_charge(idx, charge);
            if !view.atom(idx).wildcard {
                out.set_hydrogen_count(idx, Some(hs.saturating_sub(model.removed_h_neighbors[i])));
            }
        }
    }
    Some(out)
}

/// Whether the RDKit model corrects chematic's parity view of `mol`
/// (diagnostics).
#[doc(hidden)]
pub fn rdkit_model_view_differs(mol: &Molecule) -> bool {
    chematic_perception::with_rdkit_parity_view(mol, |view| match view {
        Ok(v) => patched(mol, v).is_some(),
        Err(_) => false,
    })
}
