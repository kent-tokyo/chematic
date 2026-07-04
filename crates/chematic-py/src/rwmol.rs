//! `RWMol` — RDKit-compatible editable molecule (subset: Add/Remove atom & bond).
//!
//! Wraps an owned `chematic_core::Molecule` (not the `Arc`-shared one `Mol`
//! uses) so it can be mutated in place. `GetMol()` snapshots the current
//! state into an independent, read-only `Mol`.

use std::collections::HashMap;
use std::sync::Arc;

use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;

use chematic_core::{Atom, AtomIdx, BondOrder, Element, MoleculeBuilder};

use crate::Mol;

fn parse_bond_order(order: &str) -> PyResult<BondOrder> {
    match order {
        "SINGLE" => Ok(BondOrder::Single),
        "DOUBLE" => Ok(BondOrder::Double),
        "TRIPLE" => Ok(BondOrder::Triple),
        "AROMATIC" => Ok(BondOrder::Aromatic),
        other => Err(PyValueError::new_err(format!(
            "unsupported bond order: {other:?} (expected SINGLE/DOUBLE/TRIPLE/AROMATIC)"
        ))),
    }
}

/// RDKit-compatible editable molecule. See module docs for supported subset.
#[pyclass(name = "RWMol")]
pub struct RWMol {
    inner: chematic_core::Molecule,
    props: HashMap<String, String>,
}

#[pymethods]
#[allow(non_snake_case)] // method names mirror RDKit's RWMol API (AddAtom, GetMol, ...)
impl RWMol {
    /// Create an editable molecule, optionally copy-constructed from `mol`
    /// (which is left unmodified — `RWMol` always edits an independent copy).
    #[new]
    #[pyo3(signature = (mol=None))]
    fn new(mol: Option<&Mol>) -> Self {
        match mol {
            Some(m) => RWMol {
                inner: MoleculeBuilder::from_molecule(&m.inner).build(),
                props: m.props.clone(),
            },
            None => RWMol {
                inner: MoleculeBuilder::new().build(),
                props: HashMap::new(),
            },
        }
    }

    /// Add an atom by atomic number, returning its (0-based) atom index.
    fn AddAtom(&mut self, atomic_num: u8) -> PyResult<u32> {
        let element = Element::from_atomic_number(atomic_num)
            .ok_or_else(|| PyValueError::new_err(format!("invalid atomic number: {atomic_num}")))?;
        Ok(self.inner.add_atom(Atom::new(element)).0)
    }

    /// Add a bond between two atom indices with the given order
    /// (``"SINGLE"``/``"DOUBLE"``/``"TRIPLE"``/``"AROMATIC"``).
    ///
    /// Returns the molecule's bond count *after* adding, matching RDKit's
    /// ``RWMol.AddBond`` return convention (not the new bond's own index).
    fn AddBond(&mut self, begin_idx: u32, end_idx: u32, order: &str) -> PyResult<usize> {
        let order = parse_bond_order(order)?;
        self.inner
            .add_bond(AtomIdx(begin_idx), AtomIdx(end_idx), order)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(self.inner.bond_count())
    }

    /// Remove atom `idx` and all bonds involving it. Remaining atom indices
    /// above `idx` shift down by one, matching RDKit's `RemoveAtom`.
    fn RemoveAtom(&mut self, idx: u32) -> PyResult<()> {
        if idx as usize >= self.inner.atom_count() {
            return Err(PyIndexError::new_err(format!(
                "atom index out of range: {idx}"
            )));
        }
        self.inner.remove_atom(AtomIdx(idx));
        Ok(())
    }

    /// Remove the bond between two atom indices (a no-op if none exists),
    /// matching RDKit's atom-index-pair `RemoveBond` signature.
    fn RemoveBond(&mut self, begin_idx: u32, end_idx: u32) -> PyResult<()> {
        if begin_idx as usize >= self.inner.atom_count()
            || end_idx as usize >= self.inner.atom_count()
        {
            return Err(PyIndexError::new_err("atom index out of range"));
        }
        if let Some((bond_idx, _)) = self
            .inner
            .bond_between(AtomIdx(begin_idx), AtomIdx(end_idx))
        {
            self.inner.remove_bond(bond_idx);
        }
        Ok(())
    }

    fn GetNumAtoms(&self) -> usize {
        self.inner.atom_count()
    }

    fn GetNumBonds(&self) -> usize {
        self.inner.bond_count()
    }

    /// Snapshot the current state into an independent, read-only `Mol`.
    /// Further edits to this `RWMol` do not affect the returned `Mol`.
    fn GetMol(&self) -> Mol {
        Mol {
            inner: Arc::new(MoleculeBuilder::from_molecule(&self.inner).build()),
            props: self.props.clone(),
        }
    }
}
