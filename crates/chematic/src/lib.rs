#![forbid(unsafe_code)]
//! `chematic` — umbrella crate for the chematic ecosystem.
//!
//! Enable sub-crates via feature flags:
//!
//! ```toml
//! [dependencies]
//! chematic = { version = "0.1.0", features = ["full"] }
//! ```
//!
//! Available features: `smiles`, `perception`, `mol`, `depict`, `fp`,
//! `chem`, `smarts`, `rxn`, `threed`, `full`.

#[cfg(feature = "smiles")]
pub use chematic_core as core;
#[cfg(feature = "smiles")]
pub use chematic_smiles as smiles;
#[cfg(feature = "perception")]
pub use chematic_perception as perception;
#[cfg(feature = "mol")]
pub use chematic_mol as mol;
#[cfg(feature = "depict")]
pub use chematic_depict as depict;
#[cfg(feature = "fp")]
pub use chematic_fp as fp;
#[cfg(feature = "chem")]
pub use chematic_chem as chem;
#[cfg(feature = "smarts")]
pub use chematic_smarts as smarts;
#[cfg(feature = "rxn")]
pub use chematic_rxn as rxn;
#[cfg(feature = "threed")]
pub use chematic_3d as threed;

#[cfg(test)]
mod tests {
    #[test]
    fn test_crate_compiles() {
        // Umbrella crate compiles with no features enabled.
        assert!(true);
    }
}
