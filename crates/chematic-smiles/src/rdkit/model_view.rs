//! Sanitized core views shared by RDKit-compatible model callbacks.

use super::{Molecule, RdkitSmilesError, has_added_hydrogens, parse, pdb_read, sanitize, stereo};

thread_local! {
    static IN_RDKIT_MODEL: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Memoize a scalar RDKit-model calculation on its original molecule.
/// Nested model calls deliberately bypass this cache: they run on the
/// supplied graph, while an outer call may first sanitize or remove H.
#[doc(hidden)]
#[inline]
pub fn with_rdkit_model_scalar(
    mol: &Molecule,
    slot: chematic_core::derived_cache::DerivedScalar,
    f: impl FnOnce(&Molecule) -> f64,
) -> f64 {
    if IN_RDKIT_MODEL.with(|c| c.get()) {
        return f(mol);
    }
    mol.derived_scalar(slot, || with_rdkit_model_molecule_impl(mol, f, true))
}

/// Runs `f` on RDKit's own molecule for `mol` (`Chem.MolFromSmiles(s)`:
/// hydrogen graph atoms removed, RDKit's aromaticity, bond orders, charges
/// and per-atom hydrogen counts, as [`rdkit_hydrogen_suppressed`] builds it)
/// when chematic's reading of `mol` may differ from RDKit's
/// ([`chematic_perception::rdkit_model_may_disagree`], or hydrogen graph
/// atoms), and on `mol` itself otherwise. Nested calls run on the molecule
/// they are given.
pub fn with_rdkit_model_molecule<R>(mol: &Molecule, f: impl FnOnce(&Molecule) -> R) -> R {
    with_rdkit_model_molecule_impl(mol, f, false)
}

fn with_rdkit_model_molecule_impl<R>(
    mol: &Molecule,
    f: impl FnOnce(&Molecule) -> R,
    reuse_aromatic_view: bool,
) -> R {
    if IN_RDKIT_MODEL.with(|c| c.get()) {
        return f(mol);
    }
    let needed = chematic_perception::rdkit_model_may_disagree(mol)
        || mol
            .atoms()
            .any(|(_, a)| !a.wildcard && a.element == chematic_core::Element::H);
    // Scalar descriptors on the same immutable graph share sanitization.
    // Clone and every graph edit clear this slot together with other derived
    // data. Cache failures too, retaining the original graph fallback.
    let cached = (needed && reuse_aromatic_view).then(|| {
        mol.derived(chematic_core::DerivedSlot::RdkitScalarModel, || {
            rdkit_model_molecule(mol, true).ok()
        })
    });
    let owned = if needed && !reuse_aromatic_view {
        rdkit_model_molecule(mol, false).ok()
    } else {
        None
    };
    let model = cached
        .as_ref()
        .and_then(|cached| cached.as_ref().as_ref())
        .or(owned.as_ref());
    match model {
        Some(m) => {
            struct Reset;
            impl Drop for Reset {
                fn drop(&mut self) {
                    IN_RDKIT_MODEL.with(|c| c.set(false));
                }
            }
            IN_RDKIT_MODEL.with(|c| c.set(true));
            let _reset = Reset;
            f(m)
        }
        None => f(mol),
    }
}

#[cfg(test)]
#[path = "model_view_cache_tests.rs"]
mod cache_tests;

/// Descriptor/fingerprint callbacks consume a core molecule, which cannot
/// carry CIP ranks or the port's potential-stereocenter flags. Keep stereo
/// cleanup for specified tags/directions, but skip flag-only ranking when
/// there is no specified stereo. Sanitization and hydrogen handling are the
/// same as in [`super::rdkit_mol_from_smiles`].
pub(super) fn rdkit_model_molecule(
    mol: &Molecule,
    reuse_aromatic_view: bool,
) -> Result<Molecule, RdkitSmilesError> {
    let mut m = parse::from_chematic(mol)?;
    if has_added_hydrogens(mol) {
        sanitize::sanitize_keeping_hs(&mut m)?;
    } else {
        sanitize::remove_hs_and_sanitize(&mut m)?;
    }
    stereo::legacy_stereo_perception_unflagged(&mut m);
    let view = pdb_read::to_chematic(&m, true)?;
    if reuse_aromatic_view && !chematic_perception::rdkit_parity_view_is_identity(&view) {
        // This graph has just been sanitized by the same RDKit model used
        // by the parity hook. Descriptor callers may reuse it instead of
        // perceiving and sanitizing it again. Keep the usual descriptor
        // normalization (including aromatic ring triples) in its consumer.
        // Identity views need no copy; generic callbacks need no such cache.
        view.seed_derived(
            chematic_core::DerivedSlot::RdkitParityAromatic,
            std::sync::Arc::new(Ok::<_, chematic_perception::AromaticityError>(view.clone())),
        );
    }
    Ok(view)
}
