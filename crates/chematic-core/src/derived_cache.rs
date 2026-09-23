//! Per-molecule memoization of derived perception data.
//!
//! Many descriptors, fingerprints and SMARTS matches need the same derived
//! views of a molecule (the SSSR, an aromaticity-perceived copy, …). Without
//! memoization each public call recomputed them from scratch, which made
//! repeated calls on the same [`Molecule`](crate::Molecule) pay the full
//! perception cost every time.
//!
//! The cache is invisible to callers:
//!
//! * it is cleared by every `&mut self` method on `Molecule`, so a mutated
//!   graph can never observe stale derived data;
//! * cloning a `Molecule` yields an *empty* cache, so code that clones and
//!   then edits (including crate-internal `with_*` helpers) stays correct;
//! * values are computed outside the cell and published with `set`, so a
//!   computation may itself consult other slots of the same molecule.

use std::any::Any;
use std::panic::{RefUnwindSafe, UnwindSafe};
use std::sync::{Arc, OnceLock};

/// Well-known derived-data slots. Each slot stores exactly one Rust type,
/// chosen by the perception crate that owns it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
#[doc(hidden)]
pub enum DerivedSlot {
    /// `chematic_perception::RingSet` from `find_sssr`.
    Sssr = 0,
    /// Cyclic-bond flags from `ring_bond_flags`.
    RingBondFlags = 1,
    /// Aromaticity-perceived copy used by descriptor calculations.
    DescriptorAromatic = 2,
    /// Wildman–Crippen atom-type index per atom (on the descriptor aromatic view).
    CrippenTypes = 3,
    /// Winning canonical ranks and canonical SMILES string.
    CanonicalRanks = 4,
    /// Result of the RDKit-parity aromaticity perception.
    RdkitParityAromatic = 5,
}

const SLOT_COUNT: usize = 6;

/// Stored values keep `Molecule`'s auto traits: `Send + Sync` for Rayon and
/// the unwind-safety traits so `Molecule` stays usable across
/// `std::panic::catch_unwind` (as it was before the cache existed).
type Erased = dyn Any + Send + Sync + RefUnwindSafe + UnwindSafe;

type Slot = OnceLock<Arc<Erased>>;

pub(crate) struct DerivedCache {
    slots: [Slot; SLOT_COUNT],
}

impl Default for DerivedCache {
    fn default() -> Self {
        Self {
            slots: [const { OnceLock::new() }; SLOT_COUNT],
        }
    }
}

impl Clone for DerivedCache {
    /// Clones start empty; see the module documentation.
    fn clone(&self) -> Self {
        Self::default()
    }
}

impl DerivedCache {
    pub(crate) fn get_or_compute<T, F>(&self, slot: DerivedSlot, compute: F) -> Arc<T>
    where
        T: Any + Send + Sync + RefUnwindSafe + UnwindSafe,
        F: FnOnce() -> T,
    {
        let cell = &self.slots[slot as usize];
        if let Some(existing) = cell.get() {
            let any: Arc<dyn Any + Send + Sync> = Arc::<Erased>::clone(existing);
            if let Ok(value) = any.downcast::<T>() {
                return value;
            }
            // A slot is owned by exactly one type; a mismatch is a
            // programming error, but never return wrong data — just recompute.
            debug_assert!(false, "derived-cache slot {slot:?} type mismatch");
            return Arc::new(compute());
        }
        let value = Arc::new(compute());
        let erased: Arc<Erased> = value.clone();
        // A concurrent writer may have won; either value is identical.
        let _ = cell.set(erased);
        value
    }
}
