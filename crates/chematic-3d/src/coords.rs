//! 3D coordinate types for molecular structures.
//!
//! `Point3`/`Coords3D` are now defined once in `chematic_core::coords3d` and
//! re-exported here as a compatibility surface, per the 3D Breakthrough
//! Program's Coords3D unification (`docs/rfcs/3d_breakthrough_master_plan.md`,
//! decision 1a and the Wave 1->2 integration note in
//! `chematic_core::coords3d`'s own doc comment). Every existing
//! `crate::coords::{Point3, Coords3D}` / `chematic_3d::{Point3, Coords3D}`
//! call site in this workspace keeps compiling unchanged: the two types were
//! kept field-for-field and method-for-method identical (same `x`/`y`/`z`
//! layout, same `points: Vec<Point3>` layout, same `new_zeroed`/`get`/`set`/
//! `atom_count` names and panics-on-out-of-range indexing convention) exactly
//! so this bridge could be a plain re-export with zero call-site changes.
//!
//! The `chematic-core` type is a strict superset (adds `is_finite()` on both
//! types, `Default`/`PartialEq` on `Coords3D`) so nothing that worked before
//! this bridge can stop working after it.

pub use chematic_core::{Coords3D, Point3};
