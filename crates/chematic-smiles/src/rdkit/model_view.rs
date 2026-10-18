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

/// Memoize a scalar with a proof-based constant shortcut on the original
/// graph, before model preparation. `shortcut` may return a value only when
/// it equals `f` on the sanitized/H-normalized model as well. Nested model
/// callbacks still evaluate `f` directly, bypassing outer caches/shortcuts.
#[doc(hidden)]
#[inline]
pub fn with_rdkit_model_scalar_shortcut(
    mol: &Molecule,
    slot: chematic_core::derived_cache::DerivedScalar,
    shortcut: impl FnOnce(&Molecule) -> Option<f64>,
    f: impl FnOnce(&Molecule) -> f64,
) -> f64 {
    if IN_RDKIT_MODEL.with(|c| c.get()) {
        return f(mol);
    }
    mol.derived_scalar(slot, || {
        shortcut(mol).unwrap_or_else(|| with_rdkit_model_molecule_impl(mol, f, true))
    })
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

/// Prepare the sanitized model shared by scalar descriptor calls when the
/// input graph needs one.
///
/// Returns `true` when model preparation was required (including a cached
/// failure) and `false` when scalar descriptors can consume `mol` directly.
/// This is a binding-level preparation hook: the low-level Rust parser stays
/// lazy, while bindings whose parser contract includes sanitization can put
/// that work on the parse side of a parse/descriptor boundary.
#[doc(hidden)]
pub fn prepare_rdkit_scalar_model(mol: &Molecule) -> bool {
    let needed = chematic_perception::rdkit_model_may_disagree(mol)
        || mol
            .atoms()
            .any(|(_, atom)| !atom.wildcard && atom.element == chematic_core::Element::H);
    if needed {
        let _ = mol.derived(chematic_core::DerivedSlot::RdkitScalarModel, || {
            rdkit_model_molecule(mol, true).ok()
        });
    }
    needed
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

#[cfg(test)]
mod shortcut_tests {
    use super::*;
    use chematic_core::derived_cache::DerivedScalar;
    use std::cell::Cell;

    fn no_polar_atoms(mol: &Molecule) -> bool {
        !mol.atoms()
            .any(|(_, atom)| matches!(atom.element.atomic_number(), 7 | 8))
    }

    #[test]
    fn scalar_shortcut_preserves_cache_clone_and_edit_rules() {
        let mut mol = crate::parse_template("c1ccccc1").unwrap();
        let shortcuts = Cell::new(0);
        let callbacks = Cell::new(0);
        for _ in 0..2 {
            let value = with_rdkit_model_scalar_shortcut(
                &mol,
                DerivedScalar::RdkitTpsa,
                |input| {
                    shortcuts.set(shortcuts.get() + 1);
                    no_polar_atoms(input).then_some(0.0)
                },
                |input| {
                    callbacks.set(callbacks.get() + 1);
                    input
                        .atoms()
                        .filter(|(_, a)| matches!(a.element.atomic_number(), 7 | 8))
                        .count() as f64
                },
            );
            assert_eq!(value.to_bits(), 0.0f64.to_bits());
        }
        assert_eq!(shortcuts.get(), 1);
        assert_eq!(callbacks.get(), 0);
        assert!(
            mol.derived_if_computed::<bool>(chematic_core::DerivedSlot::RdkitModelGate)
                .is_none()
        );
        let cloned = mol.clone();
        let value = with_rdkit_model_scalar_shortcut(
            &cloned,
            DerivedScalar::RdkitTpsa,
            |_| {
                shortcuts.set(shortcuts.get() + 1);
                Some(0.0)
            },
            |_| panic!("zero callback"),
        );
        assert_eq!(value, 0.0);
        assert_eq!(shortcuts.get(), 2);
        let _ = mol.add_atom(chematic_core::Atom::new(chematic_core::Element::O));
        let value = with_rdkit_model_scalar_shortcut(
            &mol,
            DerivedScalar::RdkitTpsa,
            |input| {
                shortcuts.set(shortcuts.get() + 1);
                no_polar_atoms(input).then_some(0.0)
            },
            |input| {
                callbacks.set(callbacks.get() + 1);
                input
                    .atoms()
                    .filter(|(_, a)| matches!(a.element.atomic_number(), 7 | 8))
                    .count() as f64
            },
        );
        assert_eq!(value, 1.0);
        assert_eq!(shortcuts.get(), 3);
        assert_eq!(callbacks.get(), 1);
    }

    #[test]
    fn nested_model_callback_bypasses_shortcut_and_outer_scalar_cache() {
        let mol = crate::parse("C").unwrap();
        assert_eq!(
            with_rdkit_model_scalar_shortcut(
                &mol,
                DerivedScalar::RdkitTpsa,
                |_| Some(0.0),
                |_| panic!("shortcut")
            ),
            0.0
        );
        let previous = IN_RDKIT_MODEL.with(|cell| cell.replace(true));
        let value = with_rdkit_model_scalar_shortcut(
            &mol,
            DerivedScalar::RdkitTpsa,
            |_| panic!("nested shortcut"),
            |_| 42.0,
        );
        IN_RDKIT_MODEL.with(|cell| cell.set(previous));
        assert_eq!(value, 42.0);
        assert_eq!(
            mol.derived_scalar(DerivedScalar::RdkitTpsa, || panic!("outer cache")),
            0.0
        );
    }
}
