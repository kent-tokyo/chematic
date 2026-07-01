//! Atom/bond mutation bindings (add/remove atom, add/remove bond, charge/element edits).

use crate::MolHandle;
use wasm_bindgen::prelude::*;

/// Return a new `MolHandle` with one atom appended.
///
/// The second return value is the new atom's index (as a JS number).
/// Use `with_atom_added_idx` to retrieve the index.
#[wasm_bindgen]
pub fn mol_with_atom_added(mol: &MolHandle, element_symbol: &str) -> Result<MolHandle, JsValue> {
    let element = chematic_core::Element::from_symbol(element_symbol)
        .ok_or_else(|| JsValue::from_str(&format!("unknown element: {element_symbol}")))?;
    let atom = chematic_core::Atom::new(element);
    let (new_mol, _) = mol.inner.with_atom_added(atom);
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

/// Return the index that would be assigned to an atom appended to `mol`.
#[wasm_bindgen]
pub fn mol_next_atom_idx(mol: &MolHandle) -> u32 {
    mol.inner.atom_count() as u32
}

/// Return a new `MolHandle` with one bond added between `a` and `b`.
///
/// `order` — 1 = single, 2 = double, 3 = triple.
/// Returns a JS error if the bond already exists or `a == b`.
#[wasm_bindgen]
pub fn mol_with_bond_added(
    mol: &MolHandle,
    a: u32,
    b: u32,
    order: u32,
) -> Result<MolHandle, JsValue> {
    let bond_order = match order {
        2 => chematic_core::BondOrder::Double,
        3 => chematic_core::BondOrder::Triple,
        _ => chematic_core::BondOrder::Single,
    };
    let (new_mol, _bond_idx) = mol
        .inner
        .with_bond_added(
            chematic_core::AtomIdx(a),
            chematic_core::AtomIdx(b),
            bond_order,
        )
        .map_err(|e| JsValue::from_str(&e.to_string()))?;
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

/// Return a new `MolHandle` with the formal charge of atom `idx` changed.
///
/// Returns a JS error if `idx` is out of range.
#[wasm_bindgen]
pub fn mol_with_atom_charge(mol: &MolHandle, idx: u32, charge: i32) -> Result<MolHandle, JsValue> {
    if idx as usize >= mol.inner.atom_count() {
        return Err(JsValue::from_str(&format!(
            "atom index {idx} out of range (molecule has {} atoms)",
            mol.inner.atom_count()
        )));
    }
    let new_mol = mol
        .inner
        .with_atom_charge(chematic_core::AtomIdx(idx), charge as i8);
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

/// Return a new `MolHandle` with the element of atom `idx` changed.
///
/// `element_symbol` — periodic-table symbol, e.g. `"N"`, `"O"`, `"Cl"`.
/// Returns a JS error if `idx` is out of range or the symbol is unknown.
#[wasm_bindgen]
pub fn mol_with_atom_element(
    mol: &MolHandle,
    idx: u32,
    element_symbol: &str,
) -> Result<MolHandle, JsValue> {
    if idx as usize >= mol.inner.atom_count() {
        return Err(JsValue::from_str(&format!(
            "atom index {idx} out of range (molecule has {} atoms)",
            mol.inner.atom_count()
        )));
    }
    let el = chematic_core::Element::from_symbol(element_symbol)
        .ok_or_else(|| JsValue::from_str(&format!("unknown element: {element_symbol}")))?;
    let new_mol = mol.inner.with_atom_element(chematic_core::AtomIdx(idx), el);
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

/// Return a new `MolHandle` with atom `idx` and all its bonds removed.
///
/// Atom indices above `idx` shift down by 1.  Returns a JS error if `idx`
/// is out of range.
#[wasm_bindgen]
pub fn mol_with_atom_removed(mol: &MolHandle, idx: u32) -> Result<MolHandle, JsValue> {
    if idx as usize >= mol.inner.atom_count() {
        return Err(JsValue::from_str(&format!(
            "atom index {idx} out of range (molecule has {} atoms)",
            mol.inner.atom_count()
        )));
    }
    let (new_mol, _remap) = mol.inner.with_atom_removed(chematic_core::AtomIdx(idx));
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

/// Return a new `MolHandle` with bond `idx` removed.
///
/// Atom indices are unchanged; bond indices above `idx` shift down.
/// Returns a JS error if `idx` is out of range.
#[wasm_bindgen]
pub fn mol_with_bond_removed(mol: &MolHandle, idx: u32) -> Result<MolHandle, JsValue> {
    if idx as usize >= mol.inner.bond_count() {
        return Err(JsValue::from_str(&format!(
            "bond index {idx} out of range (molecule has {} bonds)",
            mol.inner.bond_count()
        )));
    }
    let new_mol = mol.inner.with_bond_removed(chematic_core::BondIdx(idx));
    Ok(MolHandle {
        inner: std::rc::Rc::new(new_mol),
    })
}

// ---------------------------------------------------------------------------
// SDF / V3000 write
// ---------------------------------------------------------------------------
